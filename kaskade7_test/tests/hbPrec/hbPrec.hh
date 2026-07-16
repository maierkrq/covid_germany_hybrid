/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2020 Zuse Institute Berlin                            */
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

double const PI = 3.141592653589793238462;

/// Example simple heat transfer equation
///

template <class RType, class VarSet>
class HBPrecFunctional: public FunctionalBase<VariationalFunctional>
{
  using Self = HBPrecFunctional<RType,VarSet>;

public:
  using Scalar = RType;
  using OriginVars = VarSet;
  using AnsatzVars = VarSet;
  using TestVars = VarSet;
  using GridView = typename AnsatzVars::GridView;
  //static ProblemType const type = VariationalFunctional;

/// \class DomainCache
///

  class DomainCache 
  {
  public:
    DomainCache(Self const& F_,
                typename AnsatzVars::VariableSet const& vars_,
                int flags=7):
      F(F_), data(vars_) 
    {}

    void moveTo(Cell<GridView> const& c) { cell = c; }

    template <class Evaluators>
    void evaluateAt(LocalPosition<GridView> const& xi, Evaluators const& evaluators)
    {
      auto x = cell.geometry().global(xi);

      u   = data.template value<0>(evaluators);
      du  = data.template derivative<0>(evaluators)[0];
      
      //right hand side f = -laplace(u)+8u , where u is known solution
      f = (4*PI*PI+8)*cos(2*PI*x[0])-(8*x[1]-4)*(8*x[1]-4)*exp(4*(x[1]-0.5)*(x[1]-0.5));
    }

    Scalar d0() const 
    {
      return (du*du + 8*u*u)/2 - f*u;
    }
    
    template<int row, int dim> 
    Dune::FieldVector<Scalar, TestVars::template Components<row>::m>
    d1 (VariationalArg<Scalar,dim> const& arg) const 
    {
      return du*arg.derivative[0] + 8*u*arg.value - f*arg.value; 
    }

    template<int row, int col, int dim> 
    Dune::FieldMatrix<Scalar, TestVars::template Components<row>::m,
                      AnsatzVars::template Components<col>::m>
    d2 (VariationalArg<Scalar,dim> const &arg1, VariationalArg<Scalar,dim> const &arg2) const 
    {
      return arg1.derivative[0]*arg2.derivative[0] + 8*arg1.value*arg2.value;
    }

  private:
    Self const& F;
    typename AnsatzVars::VariableSet const& data;
    Cell<GridView> cell;
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

    void moveTo(FaceIterator const& f)
    {
      face = f;
    }

    template <class Position, class Evaluators>
    void evaluateAt(Position const& xi, Evaluators const& evaluators)
    {
      x = face->geometry().global(xi);

      u = data.template value<0>(evaluators);
      u0 = 0;   
    }

    // inhomogeneous Neumann condition on one part of the boundary and homogeneous Neumann condition on the rest
    Scalar
    d0() const 
    {
      return 0.0;
    }
    
    template<int row, int dim> 
    Dune::FieldVector<Scalar, TestVars::template Components<row>::m>
    d1 (VariationalArg<Scalar,dim> const& arg) const 
    {
      if ( (x[1]<=1e-12) || (x[1]>=(1-1e-12)) )
        return -4*exp(1)*arg.value[0];
      else 
        return 0.0;
    }

    template<int row, int col, int dim> 
    Dune::FieldMatrix<Scalar, TestVars::template Components<row>::m,
                      AnsatzVars::template Components<col>::m>
    d2 (VariationalArg<Scalar,dim> const &arg1, 
                       VariationalArg<Scalar,dim> const &arg2) const 
    {
      return 0.0;
    }

  private:
    typename AnsatzVars::VariableSet const& data;
    FaceIterator face;
    Dune::FieldVector<typename AnsatzVars::Grid::ctype,AnsatzVars::Grid::dimension> x;
    Scalar u, u0;
  };


  
  
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
      return 2*shapeFunctionOrder-1;
  }
  
};

/// \example hbPrec.cpp
/// show the usage of HBPrecFunctional describing a stationary heat transfer problem,
/// no adaptive grid refinement.
///
#endif
