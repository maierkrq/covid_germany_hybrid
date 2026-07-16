/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2017-2020 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <numeric>

#include "dune/grid/config.h"

#include "fem/fixdune.hh"
#include "fem/boundaryInterpolation.hh"
#include "linalg/dynamicMatrix.hh"
#include "utilities/combinatorics.hh"

namespace Kaskade
{
  namespace 
  {                                       
    // TODO: use utilities/bezier.hh
    // Bezier basis on [0,1]
    template <class Real>
    Real B(int n, int i, Real t)
    {
      return binomial(n,i)*power(1-t,n-i)*power(t,i);
    }
    
    // Derivative of Bezier basis on [0,1]
    template <class Real>
    Real Bdt(int n, int i, Real t)
    {
      Real f1dt = n==i? 0: -(n-i)*power(1-t,n-i-1);
      Real f2dt = i==0? 0: i*power(t,i-1);
      return binomial(n,i)*(f1dt*power(t,i) + power(1-t,n-i)*f2dt);
    }

    
    // This evaluates a parametric curve between x and y with tangent vectors tx and ty, respectively,
    // at t in [0,1], and returns the point on the curve as well as the curve tangent. The tangent vectors
    // provided need to be normalized.
    //
    // The curve is a rational cubic Bezier curve. In the symmetric convex setting, i.e.
    // if both tangent vectors point to the same side of the xy-line, are coplanar, and form the same angle with the 
    // xy-line, the resulting curve is a circular arc.
    //
    // For straight lines, i.e. tx = (y-x)/|y-x| and ty=-tx, the arc is not only a straight line, but also uniformly
    // parametrized, i.e. p(t) = (1-t)x + t*y.
    //
    // The construction follows the idea of Hamann/Farin/Nielson: A parametric triangular patch based on
    // generalized conics, 1991. Note, however, that the formulas given there are not completely correct.
    // Assuming for the moment a symmetric convex setting, a circular arc is constructed as a rational
    // quadratic Bezier curve. This is then degree-elevated to a cubic rational Bezier curve. One 
    // notices that this curve is also a reasonable interpolation curve in the nonconvex and nonsymmetric
    // case.
    template <class Real, int d>
    std::pair<Dune::FieldVector<Real,d>,Dune::FieldVector<Real,d>>
    arc(Dune::FieldVector<Real,d> const& x, Dune::FieldVector<Real,d> const& y, 
        Dune::FieldVector<Real,d> tx, Dune::FieldVector<Real,d> ty,
        Real t)
    {
      using Vector = Dune::FieldVector<Real,d>;
      
      if (tx*(y-x) < 0)                               // make sure the tangent vectors
        tx = -tx;                                     // point inwards
      if (ty*(x-y) < 0)
        ty = -ty;

      assert(std::abs(tx.two_norm()-1)<1e-13);
      assert(std::abs(ty.two_norm()-1)<1e-13);

      // Compute angles alpha and beta between the base line xy and the tangent vectors.
      Real D = (x-y).two_norm();                    // length of base line
      Real cosAlpha = std::min(tx*(y-x)/D,1.0);     // make sure that rounding error does not lead to cosAlpha>1
      Real alpha = std::acos(cosAlpha);             // angle between tangent at x and baseline
      Real cosBeta = std::min(ty*(x-y)/D,1.0);      // make sure that rounding error does not lead to cosBeta>1
      Real beta = std::acos(cosBeta);               // angle between tangent at y and baseline
      assert(cosAlpha>0 && cosBeta>0);
      
      // In the coplanar setting, the angle gamma between the tangent vectors is that of the 
      // triangle's corner opposite to the base line.
      Real gamma = M_PI - alpha - beta;             // angle between tangents
      Real om = sin(gamma/2);                       // weight for circular arc by rational quadratic Bezier
      
      // Perform degree elevation, representing a quadratic rational Bezier curve by
      // a cubic one.
      Real d1 = om*D/(1+2*om)/cosAlpha;             // offset of b1 along tx
      Real d2 = om*D/(1+2*om)/cosBeta;              // offset of b2 along ty
      om = (2*om+1)/3;                              // new weight shared by b1 and b2
      std::array<Vector,4> b{x, x+d1*tx, y+d2*ty, y};
      std::array<Real,4> omega{1.0, om, om, 1.0};
      
      Vector nom(0.0), nom_dt(0.0);
      Real denom = 0.0, denom_dt = 0.0;
      for (int i=0; i<=3; ++i)
      {
        Real tmp = omega[i]*B(3,i,t);
        nom += tmp*b[i];
        denom += tmp;
        
        Real tmp_dt = omega[i]*Bdt(3,i,t);
        nom_dt += tmp_dt*b[i];
        denom_dt += tmp_dt;
      }
      
      return std::make_pair(nom/denom,(nom_dt*denom-nom*denom_dt)/square(denom));
    }
  }
  
  
  // --------------------------------------------------------------------------------------------
  // --------------------------------------------------------------------------------------------

  // 2D boundary interpolation. That's easy, we only need to consider edges.
  
  
  template <class ctype>
  Dune::FieldVector<ctype,2> boundaryInterpolation(std::array<Dune::FieldVector<ctype,2>,2> const& vertices,
                                                   std::array<Dune::FieldVector<ctype,2>,2> const& normals,
                                                   std::array<Dune::FieldVector<ctype,2>,2> const& cornerNormals,
                                                   Dune::FieldVector<ctype,2> const& b)
  {
      using Vector = Dune::FieldVector<ctype,2>;
      ctype t = b[1];
      
      if (t==0)               // that's easy: we're sitting on the starting point. In that case
        return vertices[0];   // the second point may actually be invalid (because it's unclear
                              // which edge was meant), so this corner case (pun intended) is covered here as well.

      // Compute the curve tangent vectors.
      Vector const& nx = cornerNormals[0];
      Vector const& ny = cornerNormals[1];
      Vector tx; tx[0] = nx[1]; tx[1] = -nx[0];   
      Vector ty; ty[0] = -ny[1]; ty[1] = ny[0];   
      
      return arc(vertices[0],vertices[1],normalize(tx),normalize(ty),t).first;
  }
  
  // --------------------------------------------------------------------------------------------
  // --------------------------------------------------------------------------------------------  
  
  // Now 3D boundary interpolation. For filling boundary faces by superimposing boundary arcs, we need 
  // more infrastructure, and in particular have to interpolate also the corner normals along the edges.
  
  namespace
  {
    
    // For 3D, computes an edge arc and interpolates the given normal vectors along the edge.
    // A pair of position on the arc and normal is returned. 
    //
    // x,y: the edge end points
    // tx,ty: the end point boundary curve tangents, pointing inwards
    // cnx,cny: the corner normals, i.e. the surface normals at the edge endpoints
    template <class Real>
    std::pair<Dune::FieldVector<Real,3>,Dune::FieldVector<Real,3>>
    edgeByTangent(Dune::FieldVector<Real,3> const& x,  Dune::FieldVector<Real,3> const& y, 
                  Dune::FieldVector<Real,3> const& tx, Dune::FieldVector<Real,3> const& ty,
                  Dune::FieldVector<Real,3> const& cnx, Dune::FieldVector<Real,3> const& cny,
                  Real t)
    {
      using Vector = Dune::FieldVector<Real,3>;
      
      if (t==0)                       // that's easy: we're sitting on the starting point. In that case
        return std::make_pair(x,cnx); // the second point may actually be invalid (because it's unclear
      if (t==1.0)                     // which edge was meant), so this corner case (pun intended) is covered here as well.
        return std::make_pair(y,cny);
                          
      // Compute point on arc and the tangent vector
      Vector z, tz;                  
      std::tie(z,tz) = arc(x,y,tx,ty,t);
      
      // Finally compute the surface normal at the evaluated point. The surface normal vector should be approximately
      // (1-t)*cnx+t*cny, but that vector is not orthogonal to the curve tangent vector. We therefore rotate it in
      // the plane spanned by these two vectors, such that it is perpendicular also to the tangent vector.
      // We take care that the result points into the same direction as nx and ny.
      Vector nz = normalize(vectorProduct(vectorProduct(tz,(1-t)*cnx+t*cny),tz));
      assert(nz.two_norm2()>0);
    
      return std::make_pair(z,nz);
    }
    
    // For 3D, computes an edge arc and interpolates the given normal vectors along the edge.
    // A pair of position on the arc and normal is returned. 
    // preconditoin: nx+ny \ne 0
    // Based on Hamann/Farin/Nielson 1991
    template <class Real>
    std::pair<Dune::FieldVector<Real,3>,Dune::FieldVector<Real,3>>
    edge(Dune::FieldVector<Real,3> const& x,  Dune::FieldVector<Real,3> const& y, 
         Dune::FieldVector<Real,3> const& nx, Dune::FieldVector<Real,3> const& ny,
         Dune::FieldVector<Real,3> const& cnx, Dune::FieldVector<Real,3> const& cny,
         Real t)
    {
      using Vector = Dune::FieldVector<Real,3>;
      
      if (t==0)                       // that's easy: we're sitting on the starting point. In that case
        return std::make_pair(x,cnx); // the second point may actually be invalid (because it's unclear
      if (t==1.0)                     // which edge was meant), so this corner case (pun intended) is covered here as well.
        return std::make_pair(y,cny);
                    
      // First obtain the plane in which the curve lies. nx+ny is an in-plane vector, as well as y-x.
      Vector npx = vectorProduct(nx,y-x);  // normal from cross product
      Vector npy = vectorProduct(ny,x-y);  // normal from cross product
      
      // Compute the curve tangent vectors.
      Vector tx = normalize(vectorProduct(cnx,npx));    
      Vector ty = normalize(vectorProduct(cny,npy));    
      
      return edgeByTangent(x,y,tx,ty,cnx,cny,t);
    }
    
    // Weight function for blending (corner-edge) surface patches
    // i: number of vertex opposite the considered edge
    // b: barycentric coordinates of point inside patch
    template <class Real>
    Real w(int i, Dune::FieldVector<Real,3> const& b)
    {
      std::array<Real,3> B{b[1]*b[2],b[0]*b[2],b[0]*b[1]};
      
      return B[i] / (B[0]+B[1]+B[2]);
    }
    
  }
    
  // --------------------------------------------------------------------------------------------
  // --------------------------------------------------------------------------------------------
      
  
  template <class ctype>
  Dune::FieldVector<ctype,2> boundaryInterpolationByEdge(std::array<Dune::FieldVector<ctype,2>,2> const& vertices,
                                                         std::array<BoundaryEdge<int,Dune::FieldVector<ctype,2>>,1> const& edges,
                                                         bool /* isBoundaryFace */,
                                                         Dune::FieldVector<ctype,2> const& b)
  {
    constexpr int dim = 2;
    ctype eps = std::numeric_limits<ctype>::epsilon();

    // Check trivial case: we're on one corner of the edge
    for (int i=0; i<dim; ++i)
      if (b[i]>1-100*eps)
        return vertices[i];

    // We're inside the edge (there's only one). Compute the cubic rational Bezier curve.
    auto const& edge = edges[0];
    return arc(vertices[edge.vertex(0)],vertices[edge.vertex(1)],edge.tangent(0),edge.tangent(1),b[edge.vertex(1)]).first;
  }


  template <class ctype>
  Dune::FieldMatrix<ctype,2,2> boundaryInterpolationByEdgeDerivative(std::array<Dune::FieldVector<ctype,2>,2> const& vertices,
                                                                     std::array<BoundaryEdge<int,Dune::FieldVector<ctype,2>>,1> const& edges,
                                                                     bool /* isBoundaryFace */,
                                                                     Dune::FieldVector<ctype,2> const& b)
  {
    // We're on one edge. Compute the cubic rational Bezier curve.
//     return arc(vertices[edges[0].vertex(0)],vertices[edges[0].vertex(1)],edges[0].tangent(0),edges[0].tangent(1),b[edges[0].vertex(1)]).first;
    auto const& edge = edges[0];
    Dune::FieldMatrix<ctype,2,2> dx_db(0.0);
    auto dx_dt = arc(vertices[edge.vertex(0)],vertices[edge.vertex(1)],edge.tangent(0),edge.tangent(1),b[edge.vertex(1)]).second;
    for (int i=0; i<2; ++i)
      dx_db[i][ edge.vertex(1)] = dx_dt[i];
    return dx_db;
  }



  
  template <class ctype>
  Dune::FieldVector<ctype,3> boundaryInterpolationByEdge(std::array<Dune::FieldVector<ctype,3>,3> const& vertices,
                                                         std::array<BoundaryEdge<int,Dune::FieldVector<ctype,3>>,3> const& edges,
                                                         bool isBoundaryFace,
                                                         Dune::FieldVector<ctype,3> const& b)
  {
    constexpr int dim = 3;
    using Vector = Dune::FieldVector<ctype,dim>;
    using Edge = BoundaryEdge<int,Vector>;
    ctype eps = std::numeric_limits<ctype>::epsilon();
    
    // Check trivial case: we're on one corner of the triangle
    for (int i=0; i<dim; ++i)
      if (b[i]>1-100*eps)
        return vertices[i];
      
    // Check simple case: we're on one edge. Compute the cubic rational Bezier curve.
    for (int i=0; i<dim; ++i)
      if (b[i]<100*eps)
      {
        int k = (i+1)%3;
        int l = (i+2)%3;
        assert(std::find(begin(edges),end(edges),Edge(k,l)) != end(edges));
        Edge const& edge = *std::find(begin(edges),end(edges),Edge(k,l));
        return arc(vertices[edge.vertex(0)],vertices[edge.vertex(1)],edge.tangent(0),edge.tangent(1),b[edge.vertex(1)]).first;
      }
    
    // Now the real thing: we're inside the patch. For each vertex, compute a cubic rational Bezier curve
    // to the opposite edge, through the desired point, and blend them together.
    
    // First get the corner normals as the vector product of the tangent vectors...
    Vector n[3];
    if (isBoundaryFace)                                             // ... at least for boundary faces
    {
      Vector ts[3][3];
      for (int i=0; i<3; ++i)
      {
        ts[edges[i].vertex(0)][edges[i].vertex(1)] = edges[i].tangent(0);
        ts[edges[i].vertex(1)][edges[i].vertex(0)] = edges[i].tangent(1);
      }
        
      for (int i=0; i<3; ++i)
      {
        Vector ta = ts[i][(i+1)%3];
        Vector tb = ts[i][(i+2)%3];
        Vector nab = vectorProduct(ta,tb);                          // corner surface normal is orthogonal to both tangents
        if (nab.two_norm()<1e-3)                                    // ... unless these are (virtually) parallel, which is the 
        {                                                           //     case if, e.g., two edges form a differentiable feature curve
          Vector nf = vectorProduct(vertices[(i+1)%3]-vertices[i],
                                    vertices[(i+2)%3]-vertices[i]); // In that case we take the face normal and rotate it to be
          Vector nn = vectorProduct(ta,nf);                         // orthogonal to the edge tangents.
          nab = vectorProduct(nn,ta);
        }
        n[i] = normalize(nab);          
      }
    }    
    else                                                            // ... but for interior faces, we just take the face normal.
    {                                                               // Unless there is a boundary edge, this leads to 
      n[0] = normalize(vectorProduct(vertices[1]-vertices[0],       // simple linear interpolation.
                                     vertices[2]-vertices[0]));
      n[1] = n[0];
      n[2] = n[0];
    }
    
    
    Vector r;
    for (int i=0; i<3; ++i)
    {
      Vector const& x = vertices[i];                                // the corner point x
      Edge const& e = *std::find_if(begin(edges),end(edges),        // and the opposite edge
                                    [&](auto const& e) { return !e.incident(i); });
      Vector const& nx = n[i];                                      // with its normal
      Vector y, ny;                                                 // Now compute the corresponding point y
      int k = e.vertex(0);
      int l = e.vertex(1);
      std::tie(y,ny) = edgeByTangent(vertices[k],vertices[l],       // on the opposite edge and its
                                     e.tangent(0),e.tangent(1),     // surface normal vector.
                                     n[k],n[l],
                                     b[l]/(1-b[i]));
      Vector z = edge(x,y,nx,ny,nx,ny,1-b[i]).first;                // the point on the curve from x to y TODO: map this to edge by tangent
      r += w(i,b)*z;
    }
    return r;
  }
  
  

  
  // --------------------------------------------------------------------------------------------
  // --------------------------------------------------------------------------------------------
    
  namespace Simplex
  {
    // Extending the boundary interpolation into the domain interior. This is for 2D and 3D.

    // TODO: Idea: do not distinguish between boundary faces and interior faces. Just define interior edges such that all
    //       tangents are parallel to the edges themselves.

    template <class ctype, int dim>
    Vector<ctype,dim> cellInterpolation(VertexPositions<ctype,dim> const& vertices,
                                        EdgeDirections<ctype,dim> const& boundaryEdges,
                                        FaceFlags<dim> const& isBoundaryFace,
                                        Barycentric<ctype,dim> b)
    {
      using Vector = Vector<ctype,dim>;

      constexpr ctype eps = std::numeric_limits<ctype>::epsilon();


      // Check trivial case: we're on one corner of the simplex
      for (int i=0; i<=dim; ++i)
        if (b[i]>1-200*eps)
          return Vector(0.0);

      // Make sure that even if we're on some face, we don't get
      // division by zero below, and a proper weighting of boundary faces.
      // For this we move the point slightly into the interior.
      for (ctype& bi: b)
        bi = std::max(100*eps,bi);                            // make sure we're a *little* bit off the face
      b /= std::accumulate(b.begin(),b.end(),0.0);            // and afterwards normalize the barycentric coordinates


      // We will extend the boundary face displacements linearly into the cell interior, and blend them together
      // such that on every boundary face only its own contribution remains. This reproduces the G1 boundary interpolation.
      // On non-boundary faces, a (somewhat arbitrary but consistent) blending is obtained.
      ctype sumOfWeights = 0;                             // sum of blending weights, for later normalization
      Vector u(0.0);                                      // the displacement

      for (int i=0; i<dim+1; ++i)                         // step through all faces
      {

        // The vertices of face i are all except for vertex i itself.
        std::array<int,dim> vidx;                         // Note that vidx is monotone
        for (int k=0; k<i; ++k)      vidx[k] = k;
        for (int k=i+1; k<=dim; ++k) vidx[k-1] = k;

        // Extract the vertex positions
        std::array<Vector,dim> vs;                        // the vertices of the boundary face
        Vector bs;                                        // the barycentric coordinates of the corresponding boundary point
        ctype total = eps;                                // sum of barycentric face coordinates, for barycentric normalization
        ctype w = 1.0;                                    // the blending weight is just the barycentric face bubble
        for (int k=0; k<dim; ++k)
        {
          vs[k] = vertices[vidx[k]];
          bs[k] = b[vidx[k]];
          total += bs[k];                                 // sum up the barycentric coordinates for later normalization
          w *= bs[k]+eps;                                     // blending weight is just the product of barycentric coordinates
        }
        assert(total>0);
        bs /= total;                                      // renormalize to barycentric coordinates on the face
        sumOfWeights += w;                                // add up the blending weights for later normalization

        if (w>1000*eps*sumOfWeights)                      // only compute point if it's weight is nonzero
        {
          // Extract the boundary edges incident to the face
          std::array<BoundaryEdge<int,Vector>,dim*(dim-1)/2> edges;
          int idx = 0;
          for (int k=0; k<dim-1; ++k)
            for (int l=k+1; l<dim; ++l)
            {
              auto const& e = *std::find(begin(boundaryEdges),end(boundaryEdges),BoundaryEdge<int,Vector>(vidx[k],vidx[l]));
              edges[idx++] = BoundaryEdge<int,Vector>(k,l,e.tangent(0),e.tangent(1)); // vertices and tangents match since vidx is monotone
            }


          Vector z = boundaryInterpolationByEdge(vs,edges,isBoundaryFace[i],bs);   // get the displaced boundary point

          for (int k=0; k<dim; ++k)
            z -= bs[k]*vs[k];                             // and the displacement by subtracting the boundary point itself
          z *= total;                                     // linear scale down in the interior

          u += w*z;                                       // blend the displacements together
        }
      }
      u /= sumOfWeights;                                  // normalize to a convex combination of the displacements

      return u;
    }

    // explicit instantiation
    template
    Vector<double,3> cellInterpolation(VertexPositions<double,3> const& vertices,
                                       EdgeDirections<double,3> const& edges,
                                       FaceFlags<3> const& isBoundaryFace,
                                       Barycentric<double,3> b);
    template
    Vector<double,2> cellInterpolation(VertexPositions<double,2> const& vertices,
                                       EdgeDirections<double,2> const& edges,
                                       FaceFlags<2> const& isBoundaryFace,
                                       Barycentric<double,2> b);
  

    template <class ctype, int dim>
    Dune::FieldMatrix<ctype,dim,dim+1> cellInterpolationDerivative(VertexPositions<ctype,dim> const& vertices,
                                                                    EdgeDirections<ctype,dim> const& boundaryEdges,
                                                                    FaceFlags<dim> const& isBoundaryFace,
                                                                    Barycentric<ctype,dim> b)
    {
      using Vector = Vector<ctype,dim>;

      constexpr ctype eps = std::numeric_limits<ctype>::epsilon();

      for (ctype& bi: b)
        bi = std::max(100*eps,bi);                            // make sure we're a *little* bit off the face
      b /= std::accumulate(b.begin(),b.end(),0.0);            // and afterwards normalize the barycentric coordinates


      ctype sumOfWeights = 0;                             // sum of blending weights, for later normalization
      Dune::FieldVector<ctype,dim+1> dsumOfWeights_db(0.0);

      Vector u(0.0);                                      // the displacement
      Dune::FieldMatrix<ctype,dim,dim+1> du_db(0.0);     // and its derivative

      for (int i=0; i<dim+1; ++i)                         // step through all faces
      {

        // The vertices of face i are all except for vertex i itself.
        std::array<int,dim> vidx;                         // Note that vidx is monotone
        for (int k=0; k<i; ++k)      vidx[k] = k;
        for (int k=i+1; k<=dim; ++k) vidx[k-1] = k;

        // Extract the vertex positions
        std::array<Vector,dim> vs;                        // the vertices of the boundary face
        Vector bs;                                        // the barycentric coordinates of the corresponding boundary point
        Dune::FieldMatrix<ctype,dim,dim+1> dbs_db(0.0);

        ctype total = eps;                                // sum of barycentric face coordinates, for barycentric normalization
        Dune::FieldVector<ctype,dim+1> dtotal_db(0.0);

        ctype w = 1.0;                                    // the blending weight is just the barycentric face bubble
        Dune::FieldVector<ctype,dim+1> dw_db(0.0);

        for (int k=0; k<dim; ++k)
        {
          vs[k] = vertices[vidx[k]];

          bs[k] = b[vidx[k]];
          dbs_db[k][vidx[k]] = 1;

          total += bs[k];                                 // sum up the barycentric coordinates for later normalization
          dtotal_db[vidx[k]] += 1;

          dw_db = (bs[k]+eps)*dw_db + w*dbs_db[k];
          w *= bs[k]+eps;                                     // blending weight is just the product of barycentric coordinates

        }
        assert(total>0);

        dbs_db = (total*dbs_db - Dune::outerProduct(bs,dtotal_db)) / power(total,2);
        bs /= total;                                      // renormalize to barycentric coordinates on the face

        sumOfWeights += w;                                // add up the blending weights for later normalization
        dsumOfWeights_db += dw_db;

        if (w>1000*eps*sumOfWeights)                      // only compute point if it's weight is nonzero
        {
          // Extract the boundary edges incident to the face
          std::array<BoundaryEdge<int,Vector>,dim*(dim-1)/2> edges;
          int idx = 0;
          for (int k=0; k<dim-1; ++k)
            for (int l=k+1; l<dim; ++l)
            {
              auto const& e = *std::find(begin(boundaryEdges),end(boundaryEdges),BoundaryEdge<int,Vector>(vidx[k],vidx[l]));
              edges[idx++] = BoundaryEdge<int,Vector>(k,l,e.tangent(0),e.tangent(1)); // vertices and tangents match since vidx is monotone
            }


          Vector z = boundaryInterpolationByEdge(vs,edges,isBoundaryFace[i],bs);   // get the displaced boundary point
          auto dz_dbs = boundaryInterpolationByEdgeDerivative(vs,edges,isBoundaryFace[i],bs);
//           {
//             auto bs2 = bs;
//             bs2[0] += 1e-6; bs2[1] -= 1e-6;
//             auto z2 = boundaryInterpolationByEdge(vs,edges,isBoundaryFace[i],bs2);
//             auto dznum = (z2-z)/1e-6;
//             auto dz = dz_dbs*Vector{1,-1};
//             std::cout << "dznum = " << dznum << "\ndz   = " << dz << "\n\n";
//           }

          for (int k=0; k<dim; ++k)
          {
            z -= bs[k]*vs[k];                             // and the displacement by subtracting the boundary point itself
            for (int i=0; i<dim; ++i)
              dz_dbs[i][k] -= vs[k][i];
          }

          auto dz_db = dz_dbs * dbs_db;

          dz_db = total*dz_db + outerProduct(z,dtotal_db);
          z *= total;                                     // linear scale down in the interior

          u += w*z;                                       // blend the displacements together
          du_db += w*dz_db + outerProduct(z,dw_db);
       }
      }
      du_db = (sumOfWeights*du_db - outerProduct(u,dsumOfWeights_db)) / power(sumOfWeights,2);
      u /= sumOfWeights;                                  // normalize to a convex combination of the displacements

      return du_db;
    }


    // explicit instantiation
    template
    Dune::FieldMatrix<double,2,3> cellInterpolationDerivative(VertexPositions<double,2> const& vertices,
                                                              EdgeDirections<double,2> const& edges,
                                                              FaceFlags<2> const& isBoundaryFace,
                                                              Barycentric<double,2> b);
  }

  // ----------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------------------------------------------
  
  
  template <class Index, class ctype, int dim>
  AngleFeatureDetector<Index,ctype,dim>::AngleFeatureDetector(ctype featureCurveAngle, ctype featureEdgeAngle)
  {
    // Compute thresholds such that
    // (i)  t1+t2<featureCurveThreshold <=> angle(t1,t2)<featureCurveAngle
    // (ii) n1*n2<featureEdgeThreshold <=> angle(n1,n2)>featureEdgeAngle
    ctype pi = boost::math::constants::pi<ctype>();
    featureEdgeAngle = std::max(0.0,std::min(featureEdgeAngle,pi));
    featureCurveAngle = std::max(0.0,std::min(featureCurveAngle,pi));
    featureEdgeThreshold = std::cos(featureEdgeAngle);
    featureCurveThreshold = 2*std::sin(featureCurveAngle/2);
  }

  template <class Index, class ctype, int dim>
  bool AngleFeatureDetector<Index,ctype,dim>::isFeatureVertex(Index /* vidx */, Vector const& t1, Vector const& t2,
                                                              int /* bsid1 */, int /* bsid2 */) const
  {
    return (t1+t2).two_norm() >= featureCurveThreshold;
  }

  template <class Index, class ctype, int dim>
  bool AngleFeatureDetector<Index,ctype,dim>::isFeatureEdge(Index /* vidx */, Vector const& n1, Vector const& n2,
                                                            int /* bsid1 */, int /* bsid2 */) const
  {
    return n1*n2 < featureEdgeThreshold;
  }

  template <class Index, class ctype, int dim>
  bool AngleFeatureDetector<Index,ctype,dim>::isFeatureCurve(Index /* vidx */, Vector const& t1, Vector const& t2) const
  {
    return (t1+t2).two_norm() >= featureCurveThreshold;
  }

  // explicit instantiation
  template class AngleFeatureDetector<unsigned int,double,2>;
  template class AngleFeatureDetector<unsigned int,double,3>;

  namespace
  {
    /**
     * Returns a vector that represents the "center line" of the convex polyhedral cone
     * \f$ K = \{ x \mid \forall i: n_i^T x \ge 0 \} \f$. This is defined as the solution of
     * the following QP:
     * \f[ \min b + x^T x \quad \text{s.t.} \forall i: n_i^T x \le b \f]
     * If the cone is trivial, i.e. \f$ K = \{ 0 \} \f$ due to incompatible constraints, the result somehow
     * represents the "least violation" solution.
     */
    template <class Vector>
    Vector centerOfPolyhedralCone(std::vector<Vector> const& ns)
    {
      using ctype = typename Vector::field_type;
      
      // TODO: quick hack - we should use a QP solver
      Vector x = std::accumulate(begin(ns),end(ns),Vector(0.0));
      ctype spmin, spmax;
      
      for (int k=0; k<1000; ++k)
      {
        spmin = std::numeric_limits<ctype>::max();
        spmax = std::numeric_limits<ctype>::min();
        int imin = 0;
        for (int i=0; i<ns.size(); ++i)
        {
          ctype tmp = ns[i]*x;
          if (tmp < spmin)
          {
            spmin = tmp;
            imin = i;
          }
          spmax = std::max(spmax,tmp);
        }
        if (spmax < 2*spmin)
          break;
        x += ns[imin];
      }
      
      return normalize(x);
    }
    
  }
  
  // 3D algorithm
  template <class Index, class ctype>
  void adjustEdgeTangents(BoundaryStar<Index,ctype,3>& star,
                          std::map<std::array<Index,2>,Dune::FieldVector<ctype,3>>& orientedEdges,
                          GeometricFeatureDetector<Index,ctype,3> const& features)
  {
    using Vector = Dune::FieldVector<ctype,3>;
    
    int const ncorners = star.corners.size();
    
    // Sort the facets such that they are arranged in a circular fashion around the central vertex.
    for (int i=0; i<ncorners-1; ++i)
    {
      Index const commonVertex = star.corners[i].vidx[1];   // the far vertex that shall be shared by the corners
      auto& next = star.corners[i+1];                       // the next corner that shall share the far vertex

      // Among the remaining corners, find the one which shares the far vertex with the current corner.
      for (int j=i+2; j<ncorners; ++j)
        if (star.corners[j].vidx[0]==commonVertex || star.corners[j].vidx[1]==commonVertex)
        {
          std::swap(next,star.corners[j]);
          break;
        }

      // Orient the corner in such a way that current.vidx[1] == next.vidx[0].
      if (next.vidx[1]==commonVertex)
      {
        std::swap(next.vidx[0],next.vidx[1]);
        std::swap(next.vx[0],next.vx[1]);
      }
    }
    
    // Group the corners in one or several contiguous groups, such that different groups are separated
    // by a feature edge and an edge incident to corners of the same group is not a feature edge.
    // The remarkable exception here is the presence of just one corner group and one feature edge.
    std::vector<int> groups(ncorners,-1);
    groups[0] = 0;
    for (int i=0; i<ncorners-1; ++i)
      if (features.isFeatureEdge(star.vidx,star.corners[i].fn,star.corners[i+1].fn,
                                 star.corners[i].bsi, star.corners[i+1].bsi))
        groups[i+1] = groups[i]+1;
      else
        groups[i+1] = groups[i];

    // Now we've looked forward - but looking backwards may reveal that there was no feature edge where we started.
    if (!features.isFeatureEdge(star.vidx,star.corners[0].fn,star.corners.back().fn,
                                star.corners[0].bsi, star.corners.back().bsi))    // if that is the case
      for (int& g: groups)                                              // the final group (at the end)
        if (g==groups.back())                                           // is the same as the start group
          g = 0;

    int const nGroups = *std::max_element(begin(groups),end(groups)) + 1;

    // If there is only one group, there is no feature edge and all corners should share the same normal. That's 
    // the simple case, and we're done.
    if (nGroups==1)
    {
      // Get the (unique) normal at the vertex. We provide an average of the surrounding
      // face normals as proposal, but leave the final value to the feature detector.
      std::vector<Vector> fns(star.corners.size());
      for (int i=0; i<fns.size(); ++i)
        fns[i] = star.corners[i].fn;
      Vector vn = features.vertexNormal(star.vidx,centerOfPolyhedralCone(fns));

      
      for (auto& corner: star.corners)
      {
        Vector np = vectorProduct(vn,corner.vx[0]-star.vx);
        Vector t = normalize(vectorProduct(np,vn));
        orientedEdges[std::array<Index,2>{star.vidx,corner.vidx[0]}] = t;
        
        np = vectorProduct(vn,corner.vx[1]-star.vx);
        t = normalize(vectorProduct(np,vn));
        orientedEdges[std::array<Index,2>{star.vidx,corner.vidx[1]}] = t;
      }
    }
    else if (nGroups==2)                                  // exactly two groups indicates there might be a
    {                                                     // feature curve passing through this vertex
      std::array<Vector,2> ts;
      int off = 0;
      for (int i=0; i<ncorners; ++i)
        if (groups[i] != groups[(i+1)%ncorners])          // this is one of the two feature edges
        {
          ts[off] = normalize(star.corners[i].vx[1]-star.vx);
          ++off;
        }
          
      for (int i=0; i<ncorners; ++i)                                  
        if (groups[i] != groups[(i+1)%ncorners])
        {
          if (!features.isFeatureCurve(star.vidx,ts[0],ts[1]))        // the feature edges are nearly collinear but opposite
            orientedEdges[std::array<Index,2>{star.vidx,star.corners[i].vidx[1]}] 
                = normalize(normalize(star.corners[i].vx[1]-star.vx)-(ts[0]+ts[1])/2.0);    // make them collinear by averaging such that a G1-continuous feature curve results
          else                                                          // but if not almost collinear
            orientedEdges[std::array<Index,2>{star.vidx,star.corners[i].vidx[1]}] 
                = normalize(star.corners[i].vx[1]-star.vx);                     // take them as they are given
        }
      std::vector<Vector> ns(nGroups,Vector(0.0));                      // face normals within each group are averaged to a common surface normal
      for (int i=0; i<ncorners; ++i)
        ns[groups[i]] += star.corners[i].fn;
      
      for (int i=0; i<ncorners; ++i)
        if (groups[i]==groups[(i+1)%ncorners])
        {
          Vector np = vectorProduct(ns[groups[i]],star.corners[i].vx[1]-star.vx);
          orientedEdges[std::array<Index,2>{star.vidx,star.corners[i].vidx[1]}]
              = normalize(vectorProduct(np,ns[groups[i]]));
        }
    }
    else
    {
      // If there are more than two facet groups, at least three feature edges meet here.
      // There's no need to meddle with the tangents at all.
      for (auto& corner: star.corners)
      {
        orientedEdges[std::array<Index,2>{star.vidx,corner.vidx[0]}] = normalize(corner.vx[0]-star.vx);
        orientedEdges[std::array<Index,2>{star.vidx,corner.vidx[1]}] = normalize(corner.vx[1]-star.vx);
      }
    }
  }
  
  // explicit instantiation
  template
  void adjustEdgeTangents(BoundaryStar<unsigned int,double,3>& star, 
                          std::map<std::array<unsigned int,2>,Dune::FieldVector<double,3>>& orientedEdges,
                          GeometricFeatureDetector<unsigned int,double,3> const& features);


  
  // 2D algorithm
  template <class Index, class ctype>
  void adjustEdgeTangents(BoundaryStar<Index,ctype,2>& star, 
                          std::map<std::array<Index,2>,Dune::FieldVector<ctype,2>>& orientedEdges,
                          GeometricFeatureDetector<Index,ctype,2> const& features)
  {
    using Vector = Dune::FieldVector<ctype,2>;

    assert(star.corners.size()==2);                               // in 2D, a boundary vertex has exactly two incident boundary edges (corners)
    Vector e0 = normalize(star.corners[0].vx[0] - star.vx);       // and each corner only one edge, with opposite vertex 0
    Vector e1 = normalize(star.corners[1].vx[0] - star.vx);       // Now we have the two edge tangents e0 and e1 pointing in opposite direction

    bool isFeatureVertex = features.isFeatureVertex(star.vidx,e0,e1,star.corners[0].bsi,star.corners[1].bsi);

    if (isFeatureVertex)                                          // There's a large angle between the edges - a feature vertex
    {
      orientedEdges[std::array<Index,2>{star.vidx,star.corners[0].vidx[0]}] = e0;
      orientedEdges[std::array<Index,2>{star.vidx,star.corners[1].vidx[0]}] = e1;
    }
    else                                                          // Almost collinear edges - their difference is an average tangent.
    {
      // Compute the average normal.
      Vector n = star.corners[0].fn + star.corners[1].fn;
      Vector nn = normalize(features.vertexNormal(star.vidx,n));
      assert(n*nn>0);
      // Now we have to rotate the joint normal back as a tangent vector (that points along e0).
      Vector t; t[0] = nn[1]; t[1] = -nn[0];
      if (t*e0 < 0)
        t = -t;
      orientedEdges[std::array<Index,2>{star.vidx,star.corners[0].vidx[0]}] = t;  // make sure they point into  the
      orientedEdges[std::array<Index,2>{star.vidx,star.corners[1].vidx[0]}] = -t; // right direction
    }
  }
  
  // explicit instantiation
  template
  void adjustEdgeTangents(BoundaryStar<unsigned int,double,2>& star, 
                          std::map<std::array<unsigned int,2>,Dune::FieldVector<double,2>>& orientedEdges,
                          GeometricFeatureDetector<unsigned int,double,2> const& features);

}
