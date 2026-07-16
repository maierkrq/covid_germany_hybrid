/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2019 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef HEATTRANSFER_HH
#define HEATTRANSFER_HH

#include <cmath>

#include "fem/fixdune.hh"
#include "fem/variables.hh"

// Example simple heat transfer equation.
// This defines the problem   -\kappa\Delta u + q u = f,
// where f is constructed such that the solution is u = x*(x-1)*exp(-(x-0.5)^2).
// In addition, homogeneous Dirichlet boundary conditions are prescribed at x=0.

template <class RType, class VarSet>
class StatHeatTransferFunctional: public FunctionalBase<VariationalFunctional>
{
  using Self = StatHeatTransferFunctional<RType,VarSet>;

public:
  using Scalar = RType;
  using OriginVars = VarSet;
  using AnsatzVars = VarSet;
  using TestVars = VarSet;

private:
  static int const uSIdx = spaceIndex<AnsatzVars,0>;

public:

  class DomainCache 
  {
  public:
    DomainCache(Self const& F_,
                typename AnsatzVars::VariableSet const& vars_,
                int flags=7):
      F(F_), data(vars_) 
    {}

    template <class Entity>
    void moveTo(Entity const &entity) { e = &entity; }

    template <class Position, class Evaluators>
    void evaluateAt(Position const& xi, Evaluators const& evaluators)
    {
      using namespace boost::fusion;
      auto x = e->geometry().global(xi);

      u   = component<0>(data).value(at_c<uSIdx>(evaluators));
      du  = component<0>(data).derivative(at_c<uSIdx>(evaluators))[0];
      
      // manufacture the right hand side f = -\kappa\Delta u + q u,
      // where u = v*w with v=x*(x-1) and w=exp(-(x-0.5)^2).
      double v, w, vX, vXX, wX, wXX, uXX;
      v   = x[0]*(x[0] - 1);
      w   = exp(-(x[0] - 0.5)*(x[0] - 0.5));
      
      vX  =  2*x[0] - 1;
      vXX =  2;
      wX  = -2*(x[0] - 0.5)*w;
      wXX = -2*w - 2*(x[0] - 0.5)*wX;
      
      uXX = vXX*w + 2*vX*wX + v*wXX;
      f   = -F.kappa*uXX + F.q*v*w;
    }

    Scalar
    d0() const 
    {
      return (F.kappa*du*du + F.q*u*u)/2 - f*u;
    }
    
    template<int row, int dim> 
    Dune::FieldVector<Scalar, TestVars::template Components<row>::m>
    d1 (VariationalArg<Scalar,dim> const& arg) const 
    {
      return du*arg.derivative[0] + F.q*u*arg.value - f*arg.value; 
    }

    template<int row, int col, int dim> 
    Dune::FieldMatrix<Scalar, TestVars::template Components<row>::m, AnsatzVars::template Components<col>::m>
    d2 (VariationalArg<Scalar,dim> const &argT, VariationalArg<Scalar,dim> const &argA) const
    {
      return argT.derivative[0]*argA.derivative[0] + F.q*argT.value*argA.value;
    }

  private:
    Self const& F;
    typename AnsatzVars::VariableSet const& data;
    typename AnsatzVars::Grid::template Codim<0>::Entity const* e;
    Scalar f;
    Scalar u;
    Dune::FieldVector<Scalar,AnsatzVars::Grid::dimension> du;
  };

/// \class BoundaryCache
///

  class BoundaryCache 
  {
  public:
    static const bool hasInteriorFaces = false;
    using FaceIterator = typename AnsatzVars::Grid::LeafIntersectionIterator;

    BoundaryCache(Self const&,
                  typename AnsatzVars::VariableSet const& vars_,
                  int flags=7):
      data(vars_)
    {}

    void moveTo(FaceIterator const& entity)
    {
      Scalar x = entity->geometry().center()[0];
      dirichlet = (x<=1e-12) || (x>=(1-1e-12));
    }

    template <class Evaluators>
    void evaluateAt(Dune::FieldVector<typename AnsatzVars::Grid::ctype,
                                      AnsatzVars::Grid::dimension-1>
                   const& x, Evaluators const& evaluators) 
    {
      u = component<0>(data).value(boost::fusion::at_c<uSIdx>(evaluators));
    }

    Scalar
    d0() const 
    {
      if (dirichlet)
        return penalty*u*u/2;
      else return 0.0;
    }
    
    template<int row, int dim> 
    Dune::FieldVector<Scalar, TestVars::template Components<row>::m>
    d1 (VariationalArg<Scalar,dim> const& arg) const 
    {
      if (dirichlet)
        return penalty*u*arg.value[0];
      else return 0.0;
    }

    template<int row, int col, int dim> 
    Dune::FieldMatrix<Scalar, TestVars::template Components<row>::m,
                      AnsatzVars::template Components<col>::m>
    d2 (VariationalArg<Scalar,dim> const &argT,
                       VariationalArg<Scalar,dim> const &argA) const
    {
      if ( dirichlet )
        return penalty*argT.value*argA.value;
      else return 0.0;
    }

  private:
    typename AnsatzVars::VariableSet const& data;
    bool dirichlet = false;
    Scalar penalty = 1e12;
    Scalar u = 0;
  };


/// \struct StatHeatTransferFunctional constructor
///

  StatHeatTransferFunctional(Scalar kappa_, Scalar q_): kappa(kappa_), q(q_)
  {
  }
  
  
  
/// \struct D2
///

  template <int row>
  struct D1: public FunctionalBase<VariationalFunctional>::D1<row> 
  {
    static bool const present   = true;
    static bool const constant  = false;
  };
  
public:
  template <int row, int col>
  struct D2: public FunctionalBase<VariationalFunctional>::D2<row,col>
  {
    static bool const present   = true;
    static bool const symmetric = true;
    static bool const lumped    = false;
  };

/// \fn integrationOrder
///

  template <class Cell>
  int integrationOrder(Cell const& /* cell */,
                       int shapeFunctionOrder, bool boundary) const 
  {
    if (boundary) 
      return 2*shapeFunctionOrder;
    else
      return 2*shapeFunctionOrder;
  }
  
private:
  Scalar kappa, q;
};

/// \example statHeatTransfer.cpp
/// show the usage of StatHeatTransferFunctional describing a stationary heat transfer problem,
/// no adaptive grid refinement.
///
#endif
