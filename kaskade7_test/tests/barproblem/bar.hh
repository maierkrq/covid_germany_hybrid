/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2013 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef ELASTOMECHANICS_HH
#define ELASTOMECHANICS_HH

#include <type_traits>
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
  static int const dim = AnsatzVars::Grid::dimension;
  using ElasticEnergy = Kaskade::Elastomechanics::LameNavier<dim,Scalar>;
  using Vector = Dune::FieldVector<Scalar, dim>;
  using Matrix = Dune::FieldMatrix<Scalar, dim, dim>;
  static int constexpr u_Idx = 0;
  static int constexpr u_Space_Idx = boost::fusion::result_of::value_at_c<typename AnsatzVars::Variables, u_Idx>::type::spaceIndex;


  class DomainCache 
  {
  public:
    DomainCache(ElasticityFunctional const& functional, typename AnsatzVars::VariableSet const& vars_, int flags=7)
    : vars(vars_), energy(functional.moduli), length(functional.length)
    {}

    template <class Entity>
    void moveTo(Entity const& entity) {e= &entity;}

    template <class Position, class Evaluators>
    void evaluateAt(Position const& x, Evaluators const& evaluators)
    {
      energy.setLinearizationPoint( boost::fusion::at_c<u_Idx>(vars.data).derivative(boost::fusion::at_c<u_Space_Idx>(evaluators)) );
      Vector xglob = e->geometry().global(x);
      u  = boost::fusion::at_c<u_Idx>(vars.data).value(boost::fusion::at_c<u_Space_Idx>(evaluators));
      double pi=3.141592653589793238462;
      ddf2=4*pi*pi/length/length*sin(2*pi*xglob[0]/length);
      
    }

    Scalar d0() const
    {
      return energy.d0()-ddf2*u[1];
    }

    template<int row>
    Vector d1 (VariationalArg<Scalar,dim> const& arg) const
    {
	  Vector v(0); v[1]=ddf2*arg.value;
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
    typename AnsatzVars::Grid::template Codim<0>::Entity const* e;
    Vector u;
    Vector xglob;
    Scalar ddf2, length;
  };

  class BoundaryCache : public CacheBase<ElasticityFunctional,BoundaryCache>
  {
  public:
    using FaceIterator = typename AnsatzVars::Grid::LeafIntersectionIterator;

    BoundaryCache(ElasticityFunctional const& f_, typename AnsatzVars::VariableSet const& vars_, int flags=7)
    : vars(vars_), length(f_.length)
    {}

    void moveTo(FaceIterator const& face)
    {
	  
      Vector inX(0); inX[0] = 1;                        // unit vector in direction x
      auto n = face->centerUnitOuterNormal();         // unit outer normal
      e= &face;
      if (n*inX > 0.5)          // top face: boundary conditions
      {
        alpha = 1e12;             // penalty factor for not fulfilling Dirichlet u=u0, requires large penalty for hard materials such as steel with a Young's modulus around 7e9
        beta  = 0;             // Force vector for condition
      }
      else if (n*inX < -0.5)    // bottom face: boundary condition 
      {
        alpha = 1e12;
        beta  = 0;
      }
      else                   // sides: boundary condition
      {
        alpha = 0; 
        beta = 0;         
        beta[1]  = n[1];
      
      }

    }

    template <class Evaluators>
    void evaluateAt(Dune::FieldVector<typename AnsatzVars::Grid::ctype,dim-1> const& x, Evaluators const& evaluators)
    {
      using namespace boost::fusion;

      u = at_c<u_Idx>(vars.data).value(at_c<u_Space_Idx>(evaluators));
      xglob = (*e)->geometry().global(x);
      double pi=3.141592653589793238462;
      df2=2*pi/length*cos(2*pi*xglob[0]/length);
      
    }

    Scalar
    d0() const
    { 
      Vector T(0); T[0]=beta[1]*df2;
      return alpha*(u)*(u) - (T*u);
    }

    template<int row>
    Scalar d1_impl (VariationalArg<Scalar,dim,dim> const& arg) const
    {
	  Vector T(0); T[0]=beta[1]*df2;
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
    Scalar alpha, strainfactor, length,  df2;
  };

  ElasticityFunctional(Kaskade::ElasticModulus const& moduli_,Scalar length_=4.0, Scalar radius_=1.0): moduli(moduli_), length(length_)
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
};

#endif /* ELASTOMECHANICS_HH_ */
