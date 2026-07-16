/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2020-2020 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef BOUNDARYFACE_HPP
#define BOUNDARYFACE_HPP

#include <iterator>
#include <mutex>
#include <utility>
#include <vector>

#include "dune/common/fvector.hh"
#include "dune/grid/config.h"

#include "fem/barycentric.hh"
#include "fem/boundaryFace.hh"
#include "fem/fixdune.hh"
#include "fem/gridBasics.hh"
#include "utilities/detailed_exception.hh"

namespace Kaskade
{
  template <class Grid, class Face, class Displacement,int dimw>
  BoundaryFace<Grid,Face,Displacement,dimw>::BoundaryFace(Face const& f, Displacement const* d)
    : disp(d)
  {
    // Some grids do not tolerate concurrent read access.
    // On construction, we access the grid face often, so we
    // simply lock the whole constructor.
    std::lock_guard lock(GridLocking<Grid>::mutex());

    face = std::make_shared<Face>(f);

    auto const& geo = face->geometry();    // face geometry
    gt = geo.type();
    auto const& refElem = Dune::referenceElement<ctype,facedimension>(gt);

    faceAffine = geo.affine();
    displacementAffine = !disp || disp->order(cell())<=1;

    // Precompute some quantities in case of simple (i.e. affine) geometry of the
    // deformed face.
    if (faceAffine && displacementAffine)
    {
      LocalPosition xi = LocalPosition(0.0);
      std::tie(deformationJacobian,gridJacobian) = jacobiansUncached(xi);
      origin = globalUncached(xi);
      normal = unitOuterNormalUncached(xi);
    }

    // Compute the area of the deformed face.
    if (disp)
    {
      // Compute the Jacobian of the mapping from reference element to deformed face.
      // This is a d x (d-1) matrix J. We evaluate it at the center.
      // TODO: for planar meshes, geo.jacobianTransposed() must be extended with zeros,
      //       as it is only (d-1) x (d-1)
      LocalPosition xi = refElem.position(0,0);
      auto J = jacobians(xi).first;

      // Now the area is vol(refElem) * det(J^T J)^{1/2}
      vol = refElem.volume()*std::sqrt( (transpose(J)*J).determinant() );
    }
    else
      vol = geo.volume();
  }

  // ---------------------------------------------------------------------------------------------------

  template <class Grid, class Face, class Displacement,int dimw>
  BoundaryFace<Grid,Face,Displacement,dimw>::~BoundaryFace()
  {
    // The destructor of grid intersections need not be thread-safe
    // (this includes DeformingGrid as of Dune 2.7.0), hence we need
    // to lock the face destruction.
    std::lock_guard lock(GridLocking<Grid>::mutex());
    face.reset();
  }

  // ---------------------------------------------------------------------------------------------------

  template <class Grid, class Face, class Displacement,int dimw>
  std::pair<typename BoundaryFace<Grid,Face,Displacement,dimw>::GlobalPosition,
            typename BoundaryFace<Grid,Face,Displacement,dimw>::GlobalPosition>
  BoundaryFace<Grid,Face,Displacement,dimw>::boundingBox() const
  {
    // Compute bounding box. If the deformed face is affine (i.e. planar),
    // it is convex and we may just compute the bounding box of its corners.
    GlobalPosition minCorner(std::numeric_limits<ctype>::max());
    GlobalPosition maxCorner(std::numeric_limits<ctype>::lowest());

    auto refElem = Dune::referenceElement<ctype,facedimension>(gt);
    int const nCorners = refElem.size(facedimension);
    for (int i=0; i<nCorners; ++i)
    {
      LocalPosition xi = refElem.position(i,facedimension); // get corner i (of codimension=facedimension)
      GlobalPosition x = global(xi);
      minCorner = min(minCorner,x);
      maxCorner = max(maxCorner,x);
    }

    // If it is not affine, computing an exact bounding box is highly nontrivial.
    // We assume (or hope) it's not too curved, and just add the center point.
    // TODO: consider using more points (quadrature rule) or subdividing the face.
    if (! (faceAffine && displacementAffine))
    {
      LocalPosition xi = refElem.position(0,0);             // get face center (of codimension=0)
      GlobalPosition x = global(xi);
      minCorner = min(minCorner,x);
      maxCorner = max(maxCorner,x);
    }

    return std::pair(minCorner,maxCorner);
  }

  // ---------------------------------------------------------------------------------------------------

  template <class Grid, class Face, class Displacement,int dimw>
  typename BoundaryFace<Grid,Face,Displacement,dimw>::GlobalPosition
  BoundaryFace<Grid,Face,Displacement,dimw>::global(LocalPosition const& xi) const
  {
    if (faceAffine && displacementAffine)
      return origin + deformationJacobian*xi;
    else
      return globalUncached(xi);
  }

  template <class Grid, class Face, class Displacement,int dimw>
  typename BoundaryFace<Grid,Face,Displacement,dimw>::GlobalPosition
  BoundaryFace<Grid,Face,Displacement,dimw>::globalUncached(LocalPosition const& xi) const
  {
    std::lock_guard lock(GridLocking<Grid>::mutex());

    // First get the global position of the undeformed grid face.
    GlobalPosition x = extend(face->geometry().global(xi),PositionExtension(0.0));

    // ... and then displace it.
    if (disp)                         // displacement given, now evaluate
    {
      if constexpr(kind==VOLUMEMESH)  // we need to evaluate the displacement on the cell
        x += disp->value(face->inside(),
                         face->geometryInInside().global(xi));
      else                            // in surface meshes, the displacement is given on the face
        x += disp->value(*face,xi);   // (which actually is the cell)
    }

    return x;
  }

  // ---------------------------------------------------------------------------------------------------

  template <class Grid, class Face, class Displacement,int dimw>
  typename BoundaryFace<Grid,Face,Displacement,dimw>::GlobalPosition
  BoundaryFace<Grid,Face,Displacement,dimw>::unitOuterNormal(LocalPosition const& xi) const
  {
    if (faceAffine && displacementAffine)
      return normal;
    else
      return unitOuterNormalUncached(xi);
  }

  template <class Grid, class Face, class Displacement,int dimw>
  typename BoundaryFace<Grid,Face,Displacement,dimw>::GlobalPosition
  BoundaryFace<Grid,Face,Displacement,dimw>::unitOuterNormalUncached(LocalPosition const& xi) const
  {
    static_assert(2<=worlddimension && worlddimension<=3,"boundary face only defined for 2D and 3D");
    static_assert(kind==VOLUMEMESH,"not yet implemented for other cases");

    std::lock_guard lock(GridLocking<Grid>::mutex());

    // Compute the undeformed face position and normal.
    GlobalPosition n;

    auto [J,gJ] = jacobians(xi);

    // Tangents of the deformed face are obtained by applying the Jacobian of the
    // mapping from the reference domain to the deformed face. (This does not hold
    // for the normals, though, such that we must go via the tangents.)
    // Finally, the normal is obtained from the tangent vectors (by turning a tangent by
    // 90° or by taking the vector product of two tangents). For volume meshes, we need to make sure
    // that the orientation is correct, i.e. the normal actually points outwards.
    // Here we use that the deformation preserves orientation (though not orthogonality).
    if constexpr (kind==VOLUMEMESH)
    {
      if constexpr (worlddimension==2 && facedimension==1)
      {
        GlobalPosition t = J * LocalPosition{1.0};     // obtain single tangent vector
        n = GlobalPosition{t[1],-t[0]};                // turn by 90 degrees

        CellGlobalPosition gt = gJ * LocalPosition{1.0};
        CellGlobalPosition gN = CellGlobalPosition{gt[1],-gt[0]};
        if (face->outerNormal(xi)*gN < 0)              // turning by 90° points inwards
          n = -n;
      }
      else if (worlddimension==3 && facedimension==2)
      {
        GlobalPosition t0 = J*LocalPosition{1.0,0.0}, t1 = J*LocalPosition{0.0,1.0};
        n = vectorProduct(t0,t1);           // take the cross product

        CellGlobalPosition gt0 = gJ*LocalPosition{1.0,0.0}, gt1 = gJ*LocalPosition{0.0,1.0};
        CellGlobalPosition gN = vectorProduct(gt0,gt1);

        if (face->outerNormal(xi)*gN < 0)
          n = -n;
      }
    }


    return n/n.two_norm();                // return the *unit* outer normal vector
  }

  // ---------------------------------------------------------------------------------------------------

  template <class Grid, class Face, class Displacement,int dimw>
  std::optional<std::pair<typename BoundaryFace<Grid,Face,Displacement,dimw>::LocalPosition,
                          typename BoundaryFace<Grid,Face,Displacement,dimw>::ctype>>
  BoundaryFace<Grid,Face,Displacement,dimw>::intersection(GlobalPosition const& a,
                                                          GlobalPosition const& b) const
  {
    using ReturnType = std::optional<std::pair<LocalPosition,ctype>>;

    // The line segment direction.
    GlobalPosition const d = b-a;

    // Check for the easy case.
    if (faceAffine && displacementAffine)
    {
      // We look for global(xi) = o+J*xi = a+t*d, i.e. [J -d] [xi t]' = a-o
      GlobalPosition xit;
      horzcat(deformationJacobian,-d).solve(xit,a-origin);

      // Extract solution.
      LocalPosition xi; std::copy(xit.begin(),xit.begin()+facedimension,xi.begin());
      ctype t = xit[facedimension];

      if (0<=t && t<=1 && checkInside(gt,xi)<=0)
        return std::pair(xi,t);
      else
        return ReturnType();
    }

    // Now cover the complex case: the deformed face is not affine. We need to
    // solve a nonlinear equation global(xi) = a+t*d to find the intersection.

    ctype dNorm2 = d.two_norm2();

    // Compute initial guess for the intersection point: face center ...
    auto refElem = Dune::referenceElement<ctype,facedimension>(gt);
    LocalPosition xi = refElem.position(0,0);

    // Now do some Newton steps for solving global(xi) = a+td. Use F=global(xi)-a-td and
    // compute the Jacobian dF = [J -d]. The Newton step is [dxi,dt] = - dF \ F;
    for (int i=0; i<5; ++i)
    {
      GlobalPosition x = global(xi);
      ctype t = d*(x-a) / dNorm2;   // best approximation in [a,b], i.e. min |a+td - x| => t = d^T(x-a)/|d|^2.
      GlobalPosition F = x-a-t*d;

      auto dF = horzcat(jacobians(xi).first,-d);
      GlobalPosition dxit; dF.solve(dxit,-F);
      LocalPosition dxi; std::copy(dxit.begin(),dxit.begin()+facedimension,dxi.begin());

      // Apply Newton correction.
      xi += dxi;
      t += dxit[facedimension];

      // Project onto reference simplex
      xi = projectToUnitSimplex(xi);

      // Check for termination.
      if (dxit.two_norm()<1e-6)     // apparently converged, both xi and t are of order 1, hence 1e-6 "makes sense"
      {
        if (0<=t && t<=1)           // and intersection in the right position on the line
          return ReturnType(std::pair(xi,t));
        else
          return ReturnType();
      }
    }

    return ReturnType();
  }

  // ---------------------------------------------------------------------------------------------------

  template <class Grid, class Face, class Displacement,int dimw>
  typename BoundaryFace<Grid,Face,Displacement,dimw>::GlobalPosition
  BoundaryFace<Grid,Face,Displacement,dimw>::extend(CellGlobalPosition const& x,
                                                    PositionExtension const& z)
  {
    if constexpr(kind==VOLUMEMESH || kind==SURFACEMESH) // grid-global coordinates are already world coordinates
      return x;
    else
      return vertcat(x,z);
  }

  // ---------------------------------------------------------------------------------------------------

  template <class Grid, class Face, class Displacement,int dimw>
  typename std::pair<typename BoundaryFace<Grid,Face,Displacement,dimw>::DeformationJacobian,
                     typename BoundaryFace<Grid,Face,Displacement,dimw>::GridJacobian>
  BoundaryFace<Grid,Face,Displacement,dimw>::jacobians(LocalPosition const& xi) const
  {
    if (faceAffine && displacementAffine)
      return std::pair(deformationJacobian,gridJacobian);
    else
      return jacobiansUncached(xi);
  }

  template <class Grid, class Face, class Displacement,int dimw>
  typename std::pair<typename BoundaryFace<Grid,Face,Displacement,dimw>::DeformationJacobian,
                     typename BoundaryFace<Grid,Face,Displacement,dimw>::GridJacobian>
  BoundaryFace<Grid,Face,Displacement,dimw>::jacobiansUncached(LocalPosition const& xi) const
  {
    // TODO: cover the case of surface and planar meshes
    static_assert(kind == VOLUMEMESH,"not yet implemented");

    GridJacobian dF = transpose(face->geometry().jacobianTransposed(xi));

    if (disp)
    {
      // deformation is I + displacement
      auto dJ = unitMatrix<ctype,worlddimension>() + disp->derivative(cell(),face->geometryInInside().global(xi));
      return std::pair(dJ*dF,dF);
    }
    else
      return std::pair(dF,dF);
  }

  // ---------------------------------------------------------------------------------------------------

  template <class Grid, class Face, class Displacement,int dimw>
  typename BoundaryFace<Grid,Face,Displacement,dimw>::Cell
  BoundaryFace<Grid,Face,Displacement,dimw>::cell() const
  {
    if constexpr (kind==VOLUMEMESH)    // We are a boundary face in the grid ...
      return face->inside();           // ... and return the corresponding cell.
    else   // surface or planar mesh
      return *face;                    // We are already a cell.
  }
}

#endif
