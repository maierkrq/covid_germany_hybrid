/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2019-2020 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef ELASTOMECHANICS_HH
#define ELASTOMECHANICS_HH

#include <type_traits>

#include <boost/math/constants/constants.hpp>

#include "fem/functional_aux.hh"
#include "fem/diffops/elastoVariationalFunctionals.hh"


// Deriving from FunctionalBase introduces default D1 and D2 structures.
template <class VarSet>
class ElasticityFunctional: public Kaskade::FunctionalBase<VariationalFunctional>
{
public:
  using Scalar = double;
  using AnsatzVars = VarSet;
  using TestVars = VarSet;
  using OriginVars = VarSet;
  using Grid = typename AnsatzVars::Grid;
  using GridView = typename AnsatzVars::GridView;

  static int const dim = Grid::dimension;
  using ElasticEnergy = Kaskade::Elastomechanics::LameNavier<dim,Scalar>;
  using Vector = Dune::FieldVector<Scalar, dim>;
  using Matrix = Dune::FieldMatrix<Scalar, dim, dim>;
  static int constexpr u_Idx = 0;
  static int constexpr u_Space_Idx = spaceIndex<AnsatzVars,u_Idx>;


  class DomainCache 
  {
  public:
    DomainCache(ElasticityFunctional const& functional, typename AnsatzVars::VariableSet const& vars_, int flags=7)
    : vars(vars_), energy(functional.moduli), length(functional.length)
    {}

    void moveTo(Cell<GridView> const& c)
    {
      cell = c;
    }

    template <class Position, class Evaluators>
    void evaluateAt(Position const& xi, Evaluators const& evaluators)
    {
      energy.setLinearizationPoint( component<u_Idx>(vars).derivative(boost::fusion::at_c<u_Space_Idx>(evaluators)) );
      auto x = cell.geometry().global(xi)[0];
      u  = vars.template value<u_Idx>(evaluators);
      double pi = boost::math::constants::pi<double>();
      ddf2 = 4*pi*pi/length/length*sin(2*pi*x/length);
      ddf3 = 2*length/(x+0.5*length)/(x+0.5*length)/(x+0.5*length);
    }

    Scalar d0() const
    {
      return energy.d0()-ddf2*u[1]-ddf3*u[2];
    }

    template<int row>
    Vector d1 (VariationalArg<Scalar,dim> const& arg) const
    {
      Vector v{0,ddf2*arg.value,ddf3*arg.value}; // (0); v[1] = ddf2*arg.value; v[2] = ddf3*arg.value;
      return energy.d1(arg)-v;
    }

    template<int row, int col>
    Matrix d2 (VariationalArg<Scalar,dim> const& argTest, VariationalArg<Scalar,dim> const& argAnsatz) const
    {
      return energy.d2(argTest,argAnsatz);
    }

  private:
    typename AnsatzVars::VariableSet const& vars;
    ElasticEnergy energy;
    Cell<GridView> cell;
    Vector u;
    Scalar ddf2, ddf3, length;
  };

  class BoundaryCache : public CacheBase<ElasticityFunctional,BoundaryCache>
  {
  public:
    using FaceIterator = typename GridView::IntersectionIterator;

    BoundaryCache(ElasticityFunctional const& f_, typename AnsatzVars::VariableSet const& vars_, int flags=7)
    : vars(vars_), length(f_.length), radius(f_.radius)
    {}

    void moveTo(FaceIterator const& face)
    {
	  
      Vector inX(0); inX[0] = 1;                        // unit vector in direction x
      auto n = face->centerUnitOuterNormal();         // unit outer normal
      e= &face;
      if (n*inX > 0.5)          // top face: boundary conditions
      {
        alpha = 1e8;             // penalty factor for not fulfilling Dirichlet u=u0, requires large penalty for hard materials such as steel with a Young's modulus around 7e9
        beta  = 0;             // Force vector for condition
      }
      else if (n*inX < -0.5)    // bottom face: boundary condition 
      {
        alpha = 1e8;
        beta  = 0;
      }
      else                   // sides: boundary condition
      {
        alpha = 0; 
        beta = 0;         
        beta[1]  = n[1]; beta[2] = n[2];
      
      }
    }

    template <class Evaluators>
    void evaluateAt(Dune::FieldVector<typename Grid::ctype,dim-1> const& x, Evaluators const& evaluators)
    {
      using namespace boost::fusion;

      u = at_c<u_Idx>(vars.data).value(at_c<u_Space_Idx>(evaluators));
      xglob = (*e)->geometry().global(x);
      double pi = boost::math::constants::pi<double>();
      df2 = 2*pi/length*cos(2*pi*xglob[0]/length);
      df3 = length/(xglob[0]+0.5*length)/(xglob[0]+0.5*length)-4.0/3.0/length;
    }

    Scalar
    d0() const
    { 
      Vector T(0); T[0]=beta[1]*df2+beta[2]*df3;
      return alpha*(u)*(u) - (T*u);
    }

    template<int row>
    Scalar d1_impl (VariationalArg<Scalar,dim,dim> const& arg) const
    {
	  Vector T(0); T[0]=beta[1]*df2+beta[2]*df3;
      return 2*alpha*(u*arg.value) - T*arg.value;
    }

    template<int row, int col>
    Scalar d2_impl (VariationalArg<Scalar,dim,dim> const &arg1, VariationalArg<Scalar,dim,dim> const &arg2) const
    {
      return 2*alpha*(arg1.value*arg2.value);
    }

  private:
    typename AnsatzVars::VariableSet const& vars;
    FaceIterator const* e;
    Vector u, u0, beta,  xglob;
    Scalar alpha, strainfactor, length, radius, df2, df3;
  };

  ElasticityFunctional(Kaskade::ElasticModulus const& moduli_,Scalar length_=4.0, Scalar radius_=1.0)
  : moduli(moduli_), length(length_), radius(radius_)
  {
  }

  template <class Cell>
  int integrationOrder(Cell const& /* cell */, int shapeFunctionOrder, bool boundary) const
  {
    if (boundary)
      return 2*shapeFunctionOrder;      // mass term u*u on boundary
    else
      return 2*(shapeFunctionOrder-1);  // energy term "u_x * u_x" in interior
  }

private:
  Kaskade::ElasticModulus moduli;
  Scalar length;
  Scalar radius;
};

#endif /* ELASTOMECHANICS_HH_ */
