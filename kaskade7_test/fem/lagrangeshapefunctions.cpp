/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2014-2018 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "dune/grid/config.h"

#include "fem/lagrangeshapefunctions.hh"
#include "utilities/gaussNodes.hh"


namespace Kaskade
{
  
  // --------------------------------------------------------------------------
  // --------------------------------------------------------------------------
  
  template <int dim>
  LagrangeBase<dim>::LagrangeBase(int order_)
  : myOrder(order_)
  , ls(SimplexLagrangeDetail::size(dim,myOrder))
  {
    // Create barycentric indices for all shape functions
    for (int i=0; i<ls.size(); ++i)
    {
      auto ti = SimplexLagrangeDetail::tupleIndex<dim>(myOrder,i);
      std::copy(begin(ti),end(ti),begin(ls[i]));
      ls[i][dim] = myOrder - std::accumulate(begin(ti),end(ti),0);
    }
  }

  // --------------------------------------------------------------------------
  
  template <int dim, class Real>
  EquidistantLagrange<dim,Real>::EquidistantLagrange(int order)
  : LagrangeBase<dim>(order)
  , normalization(this->size(),1.0)
  {
    // Compute normalization factors
    for (int i=0; i<this->size(); ++i)
      normalization[i] = 1 / value(i,nodalPosition(i));
  }

  template <int dim, class Real>
  Dune::FieldVector<Real,dim> EquidistantLagrange<dim,Real>::nodalPosition(int k) const
  {
    // construct equidistant interpolation nodes. max  in denominator is just to prevent GCC from
    // choking on division by 0.
    Dune::FieldVector<Real,dim> xi;
    for (int i=0; i<dim; ++i)
      xi[i] =  this->order()==0? 1.0/(dim+1): static_cast<Real>(this->ls[k][i])/this->order();
    return xi;
  }

  template <int dim, class Real>
  Real EquidistantLagrange<dim,Real>::value(int k, Vector const& xi) const
  {
    Dune::FieldVector<Real,dim+1> zeta = barycentric(xi);
    
    // Remember that sum(xi)=order. Thus, phi below is the product of
    // order linear factors, i.e. the value of a polynomial of degree
    // order. Moreover, due to the form of the linear factors, it is
    // zero on each affine facet that lies below the associated node
    // in barycentric coordinate direction i. Since the barycentric
    // coordinate directions "enclose" the node, it is the only one
    // that is not covered by any such facet. 
    Real phi = 1.0;
    for (int i=0; i<=dim; ++i)
      for (int j=0; j<this->ls[k][i]; ++j)            // the inner loop generates l[i] factors
      {
        Real pos = static_cast<Real>(j)/this->order();  // Equidistant points. This will not be reached for order==0.
        phi *= zeta[i]-pos;    
      }
    return normalization[k] * phi;
  }
  
  template <int dim, class Real>
  Real EquidistantLagrange<dim,Real>::derivative(int k, int dir, Vector const& xi) const
  {
    Dune::FieldVector<Real,dim+1> zeta = barycentric(xi);
    
    // The derivative is simply based on the product rule and the
    // chain rule for zeta (algorithmic differentiation). Note that only two components of zeta have
    // a nonvanishing derivative, these are zeta[dir] (derivative 1)
    // and zeta[dim] (derivative -1).
    Dune::FieldVector<Real,dim+1> dzeta(0); dzeta[dir] = 1; dzeta[dim] = -1;
    
    Real phi = 1, dphi = 0;                         // constant and its derivative (which is of course zero)
    for (int i=0; i<=dim; ++i)
      for (int j=0; j<this->ls[k][i]; ++j)            // the inner loop generates l[i] factors
      {
        Real pos = static_cast<Real>(j)/this->order();  // Equidistant points. This will not be reached for order==0.
                                                    
        dphi = dphi*(zeta[i]-pos) + phi*dzeta[i];     // compute phi * (zeta[i]-pos) and its derivative 
        phi *= zeta[i]-pos;                           // dphi*(zeta[i]-pos) + phi*dzeta[i]
      }
    
    return normalization[k] * dphi;
  }

  template <int dim, class Real>
  Real EquidistantLagrange<dim,Real>::derivative2(int k, int dir1, int dir2, Vector const& xi) const
  {
    Dune::FieldVector<Real,dim+1> zeta = barycentric(xi);
    Dune::FieldVector<Real,dim+1> dzeta1(0); dzeta1[dir1] = 1; dzeta1[dim] = -1; // direction for derivative 
    Dune::FieldVector<Real,dim+1> dzeta2(0); dzeta2[dir2] = 1; dzeta2[dim] = -1; // direction for derivative
    
    Real phi = 1, dphi1 = 0, dphi2 = 0, dphi12 = 0;   // constant and its derivatives (which are of course zero)
    for (int i=0; i<=dim; ++i)
      for (int j=0; j<this->ls[k][i]; ++j)            // the inner loop generates l[i] factors
      {
        Real pos = static_cast<Real>(j)/this->order();  // Equidistant points. This will not be reached for order==0.
        // compute phi * (zeta[i]-pos) and its derivatives dphi1*(zeta[i]-pos) + phi*dzeta1[i],  
        // dphi2*(zeta[i]-pos) + phi*dzeta2[i], and second derivative 
        // dphi12*(zeta[i]-pos) + dphi1*dzeta2[i] + dphi2*dzeta1[i] + phi*dzeta12[i]. 
        // Due to linearity of zeta, the last term vanishes.
        dphi12 = dphi12*(zeta[i]-pos) + dphi1*dzeta2[i] + dphi2*dzeta1[i];
        dphi1 = dphi1*(zeta[i]-pos) + phi*dzeta1[i];
        dphi2 = dphi2*(zeta[i]-pos) + phi*dzeta2[i];
        phi *= zeta[i]-pos;
      }
    
    return normalization[k] * dphi12;
  }
  
  // explicit instantiation
  template class EquidistantLagrange<1,double>;
  template class EquidistantLagrange<2,double>;
  template class EquidistantLagrange<3,double>;
  
  // --------------------------------------------------------------------------
  // --------------------------------------------------------------------------

  
  template <int dim, class Real>
  LobattoGridLagrange<dim,Real>::LobattoGridLagrange(int order)
  : LagrangeBase<dim>(order)
  , lobatto(lobattoNodes(order+1))
  , basis(order)
  , vInv(this->size(),this->size())
  {
    for (double& t: lobatto)            // transform Lobatto nodes from [-1,1] 
      t = (t+1)/2;                      // to [0,1]
      
    // Create interpolation matrix: the inverse of the generalized Vandermonde matrix
    // of basis values at our interpolation nodes.
    //
    // We use the equidistant Lagrange basis, which is suboptimal, but simple and 
    // easily available. In 2D, the Vandermonde condition number grows to 
    // about 4000 at order 14, which is a factor 40 larger than the (much more complex)
    // Proriol basis as reported in Warburton 2006. That's bad, but we opt for simplicity
    // here.
    for (int i=0; i<this->size(); ++i)
    {
      Vector xi = nodalPosition(i);
      for (int j=0; j<this->size(); ++j)
        vInv[i][j] = basis.value(j,xi);
    }

    invert(vInv);
  }

  template <int dim, class Real>
  Dune::FieldVector<Real,dim> LobattoGridLagrange<dim,Real>::nodalPosition(int k) const
  {
    // Construct interpolation nodes from a Lobatto grid similar to MG Blyth, C Pozrikidis. 
    // A Lobatto interpolation grid over the triangle. IMA J. Appl. Math., 1-17, 2005.
    // This is done by selecting the *barycentric* coordinates according to the Lobatto nodes 
    // (which leads to non-normalized "barycentric" coordinates. 
    Dune::FieldVector<Real,dim+1> lambda;
    for (int i=0; i<=dim; ++i)
      lambda[i] = lobatto[this->ls[k][i]];
    lambda /= lambda.one_norm();
    
    Vector xi;
    std::copy_n(lambda.begin(),dim,xi.begin()); // extract the cartesian coordinates
    return xi;
  }

  template <int dim, class Real>
  Real LobattoGridLagrange<dim,Real>::value(int k, Vector const& xi) const
  {
    Real phi = 0;
    for (int i=0; i<this->size(); ++i)
      phi += vInv[i][k] * basis.value(i,xi);
    return phi;
  }
  
  template <int dim, class Real>
  Real LobattoGridLagrange<dim,Real>::derivative(int k, int dir, Vector const& xi) const
  {
    Real dphi = 0;
    for (int i=0; i<this->size(); ++i)
      dphi += vInv[i][k] * basis.derivative(i,dir,xi);
    return dphi;
  }

  template <int dim, class Real>
  Real LobattoGridLagrange<dim,Real>::derivative2(int k, int dir1, int dir2, Vector const& xi) const
  {
    Real ddphi = 0;
    for (int i=0; i<this->size(); ++i)
      ddphi += vInv[i][k] * basis.derivative2(i,dir1,dir2,xi);
    return ddphi;
  }
  
  // explicit instantiation
  template class LobattoGridLagrange<1,double>;
  template class LobattoGridLagrange<2,double>;
  template class LobattoGridLagrange<3,double>;
  
  // --------------------------------------------------------------------------
  // --------------------------------------------------------------------------

  
  template <class ctype, int dimension, class Scalar, int O>
  LagrangeCubeShapeFunction<ctype,dimension,Scalar,O>::LagrangeCubeShapeFunction(Dune::FieldVector<int,dimension> const& xi_)
  : xi(xi_)
  {
    for (int i=0; i<dim; ++i)
      assert(xi[i]>=0 && xi[i]<=order);

    // construct interpolation nodes (here: equidistant nodes). ?: in
    // denominator is just to prevent division by 0 for order 0.
    for (int i=0; i<=order; ++i)
      if (order==0)
        node[i] = 0.5;
      else
        node[i] = static_cast<CoordType>(i)/(order? order: 1);

    // compute normalization factor such that the shape functions have
    // value 1 at their node
    Dune::FieldVector<CoordType,dim> x;
    for (int i=0; i<dim; ++i)
      x[i] = node[xi[i]];
    normalization = 1.0 / evaluateNonNormalized(x);

    // compute index local in element
    local_ = 0;
    for (int i=0; i<dim; ++i)
      local_ = (order+1)*local_ + xi[i];

    // compute codimension of subentity in which the node resides. If
    // an index is 0 or order, the node lies on a face in that
    // coordinate direction, hence the codimension has to be increased
    // by 1.
    codim_ = 0;
    for (int i=0; i<dim; ++i)
      if (xi[i]==0 || xi[i]==order)
        ++codim_;

    // compute subentity index. See the Dune reference element page for subentity numbering.
    entity_ = 0;
    if (codim_ == 0) // the cell itself
      entity_ = 0;
    else if (codim_ == dim) {
      // vertices are numbered according to the binary number of their
      // reversed coordinates. All of xi's entries are either 0 or
      // order.
      entity_ = 0;
      for (int i=dim-1; i>=0; --i)
        entity_ = 2*entity_ + (xi[i]==0? 0: 1);
    } else if (codim_ == 1) {
      // faces are numbered according to which coordinate direction they delimit.
      entity_ = 0;
      for (int i=0; i<dim; ++i)
        if (xi[i]==0 || xi[i]==order) // only one i passes this test!
          entity_ = 2*i + (xi[i]==0? 0: 1);
    } else if (codim_ == dim-1) {
      entity_ = 0;
      int off = 0;
      for (int i=0; i<dim; ++i) {
        if (xi[i]>0 && xi[i]<order) // exactly one i passes this test!
          off = (1<<(dim-1))*i;
        else
          entity_ = 2*entity_ + (xi[i]==0? 0: 1);
      }
      entity_ += off;
    } else
      assert("Unknown codimension!"==0);

    // Precompute completely local subentity indices. Computation of
    // those subentity indices which depend on global orientation of
    // the subentity are postponed to actual call.
    if (codim()==0) {
      // We are in the interior of the cell. Interpete the interior
      // nodes as indices to the reduced order order-1
      // by subtracting 1 from all coordinates.
      entityIndex_ = 0;
      for (int i=0; i<dim; ++i)
        if (xi[i]>0 && xi[i]<order)
          entityIndex_ = (order-1)*entityIndex_ + xi[i]-1;
    } else if (order<=2) {
      // Ok, since there is only one degree of freedom on this
      // subentity, a globally unique numbering is trivial...
      entityIndex_ = 0;
    } else
      // Oops, this is the hard case. Here the numbering depends on
      // the actual element. Schedule the computation to be performed
      // lateron.
      entityIndex_ = -1;
  }
  // explicit instantiation
  template LagrangeCubeShapeFunction<double,2,double,0>::LagrangeCubeShapeFunction(Dune::FieldVector<int,2> const& xi_);
  template LagrangeCubeShapeFunction<double,2,double,1>::LagrangeCubeShapeFunction(Dune::FieldVector<int,2> const& xi_);
  template LagrangeCubeShapeFunction<double,2,double,2>::LagrangeCubeShapeFunction(Dune::FieldVector<int,2> const& xi_);
  
  template LagrangeCubeShapeFunction<double,3,double,0>::LagrangeCubeShapeFunction(Dune::FieldVector<int,3> const& xi_);
  template LagrangeCubeShapeFunction<double,3,double,1>::LagrangeCubeShapeFunction(Dune::FieldVector<int,3> const& xi_);
  template LagrangeCubeShapeFunction<double,3,double,2>::LagrangeCubeShapeFunction(Dune::FieldVector<int,3> const& xi_);

  template LagrangeCubeShapeFunction<double,2,float,0>::LagrangeCubeShapeFunction(Dune::FieldVector<int,2> const& xi_);
  template LagrangeCubeShapeFunction<double,2,float,1>::LagrangeCubeShapeFunction(Dune::FieldVector<int,2> const& xi_);
  template LagrangeCubeShapeFunction<double,2,float,2>::LagrangeCubeShapeFunction(Dune::FieldVector<int,2> const& xi_);

  template LagrangeCubeShapeFunction<double,3,float,0>::LagrangeCubeShapeFunction(Dune::FieldVector<int,3> const& xi_);
  template LagrangeCubeShapeFunction<double,3,float,1>::LagrangeCubeShapeFunction(Dune::FieldVector<int,3> const& xi_);
  template LagrangeCubeShapeFunction<double,3,float,2>::LagrangeCubeShapeFunction(Dune::FieldVector<int,3> const& xi_);
  
  
  template <class ctype, int dimension, class Scalar, int O>
  Scalar LagrangeCubeShapeFunction<ctype,dimension,Scalar,O>::evaluateDerivativeNonNormalized(int dir, Dune::FieldVector<ctype,dimension> const& x) const 
  {
    ResultType independentFactor = 1;
    for (int i=0; i<dim; ++i)
      if (i!=dir)
        for (int k=0; k<=order; ++k)
          if (k!=xi[i])
            independentFactor *= x[i]-node[k];

    ResultType result = 0;
    for (int j=0; j<=order; ++j) {
      if (j!=xi[dir]) {
        ResultType dphi = 1;
        for (int  k=0; k<=order; ++k)
          if (k!=xi[dir] && k!=j)
            dphi *= x[dir]-node[k];
        result += dphi;
      }
    }

    return independentFactor * result;
  }
  // explicit instantiation
  template double LagrangeCubeShapeFunction<double,2,double,0>::evaluateDerivativeNonNormalized(int dir, Dune::FieldVector<double,2> const& x) const;
  template double LagrangeCubeShapeFunction<double,2,double,1>::evaluateDerivativeNonNormalized(int dir, Dune::FieldVector<double,2> const& x) const;
  template double LagrangeCubeShapeFunction<double,2,double,2>::evaluateDerivativeNonNormalized(int dir, Dune::FieldVector<double,2> const& x) const;
  
  template double LagrangeCubeShapeFunction<double,3,double,0>::evaluateDerivativeNonNormalized(int dir, Dune::FieldVector<double,3> const& x) const;
  template double LagrangeCubeShapeFunction<double,3,double,1>::evaluateDerivativeNonNormalized(int dir, Dune::FieldVector<double,3> const& x) const;
  template double LagrangeCubeShapeFunction<double,3,double,2>::evaluateDerivativeNonNormalized(int dir, Dune::FieldVector<double,3> const& x) const;

  template float LagrangeCubeShapeFunction<double,2,float,0>::evaluateDerivativeNonNormalized(int dir, Dune::FieldVector<double,2> const& x) const;
  template float LagrangeCubeShapeFunction<double,2,float,1>::evaluateDerivativeNonNormalized(int dir, Dune::FieldVector<double,2> const& x) const;
  template float LagrangeCubeShapeFunction<double,2,float,2>::evaluateDerivativeNonNormalized(int dir, Dune::FieldVector<double,2> const& x) const;

  template float LagrangeCubeShapeFunction<double,3,float,0>::evaluateDerivativeNonNormalized(int dir, Dune::FieldVector<double,3> const& x) const;
  template float LagrangeCubeShapeFunction<double,3,float,1>::evaluateDerivativeNonNormalized(int dir, Dune::FieldVector<double,3> const& x) const;
  template float LagrangeCubeShapeFunction<double,3,float,2>::evaluateDerivativeNonNormalized(int dir, Dune::FieldVector<double,3> const& x) const;
  
  
  template <class ctype, int dimension, class Scalar, int Ord>
  LagrangeCubeShapeFunctionSet<ctype,dimension,Scalar,Ord>::LagrangeCubeShapeFunctionSet()
  : ShapeFunctionSet<ctype,dimension,Scalar>(Dune::GeometryType(Dune::GeometryType::cube,dimension))
  {
    Dune::FieldVector<int,dimension> xi(0);
    for (int i=0; i<power(Ord+1,dimension); ++i) 
    {
      sf.push_back(value_type(xi));
      this->iNodes.push_back(sf.back().position());
      
      // shift to next node index
      for (int j=0; j<dimension; ++j)
        if (xi[j]<Ord) {
          ++xi[j];
          break;
        } else
          xi[j] = 0;
    }
    
    this->size_ = sf.size();
    this->order_ = dimension * Ord;
  }
  // explicit instantiation
  template LagrangeCubeShapeFunctionSet<double,1,double,0>::LagrangeCubeShapeFunctionSet();
  template LagrangeCubeShapeFunctionSet<double,1,double,1>::LagrangeCubeShapeFunctionSet();
  template LagrangeCubeShapeFunctionSet<double,1,double,2>::LagrangeCubeShapeFunctionSet();
  
  template LagrangeCubeShapeFunctionSet<double,2,double,0>::LagrangeCubeShapeFunctionSet();
  template LagrangeCubeShapeFunctionSet<double,2,double,1>::LagrangeCubeShapeFunctionSet();
  template LagrangeCubeShapeFunctionSet<double,2,double,2>::LagrangeCubeShapeFunctionSet();
  
  template LagrangeCubeShapeFunctionSet<double,3,double,0>::LagrangeCubeShapeFunctionSet();
  template LagrangeCubeShapeFunctionSet<double,3,double,1>::LagrangeCubeShapeFunctionSet();
  template LagrangeCubeShapeFunctionSet<double,3,double,2>::LagrangeCubeShapeFunctionSet();

  template LagrangeCubeShapeFunctionSet<double,1,float,0>::LagrangeCubeShapeFunctionSet();
  template LagrangeCubeShapeFunctionSet<double,1,float,1>::LagrangeCubeShapeFunctionSet();
  template LagrangeCubeShapeFunctionSet<double,1,float,2>::LagrangeCubeShapeFunctionSet();

  template LagrangeCubeShapeFunctionSet<double,2,float,0>::LagrangeCubeShapeFunctionSet();
  template LagrangeCubeShapeFunctionSet<double,2,float,1>::LagrangeCubeShapeFunctionSet();
  template LagrangeCubeShapeFunctionSet<double,2,float,2>::LagrangeCubeShapeFunctionSet();

  template LagrangeCubeShapeFunctionSet<double,3,float,0>::LagrangeCubeShapeFunctionSet();
  template LagrangeCubeShapeFunctionSet<double,3,float,1>::LagrangeCubeShapeFunctionSet();
  template LagrangeCubeShapeFunctionSet<double,3,float,2>::LagrangeCubeShapeFunctionSet();
  
}
