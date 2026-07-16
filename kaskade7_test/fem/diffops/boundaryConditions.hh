/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2019-2022 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef BOUNDARYCONDITIONS_HH
#define BOUNDARYCONDITIONS_HH

#include "fem/shapefunctioncache.hh"

namespace Kaskade
{
  /**
   * \ingroup boundaryconditions
   * \brief A simple boundary cache implementing homogeneous Neumann boundary conditions.
   *
   * This realizes \f$ u_x \kappa n = 0 \f$ for all variables, i.e. the boundary integral term vanishes.
   *
   * Within the variational functional (or weak formulation), specify global homogeneous Neumann
   * boundary conditions by
   * \code
   * using Self = ...; // type of variational functional
   * using BoundaryCache = Kaskade::HomogeneousNeumannBoundary<Self>;
   * \endcode
   */
  template <class Functional>
  class HomogeneousNeumannBoundary
  {
    using Scalar = typename Functional::Scalar;
    static constexpr int dim = Functional::AnsatzVars::Grid::dimension;

  public:
    /**
     * \brief Constructor.
     *
     * This is intended for nonlinear variational functionals.
     */
    HomogeneousNeumannBoundary(Functional const&, typename Functional::OriginVars::VariableSet const&, int) {}

    /**
     * \brief Constructor.
     *
     * This is intended for linear variational functionals.
     */
    HomogeneousNeumannBoundary(int) {}

    template <class FaceIterator>
    void moveTo(FaceIterator const&) {}

    template <class Evaluators>
    void evaluateAt(Dune::FieldVector<typename Functional::AnsatzVars::Grid::ctype, dim-1> const&, Evaluators const&) {}

    Scalar d0() const { return 0; }

    template <int row>
    auto d1(VariationalArg<Scalar,dim> const& ) const
    {
      return Dune::FieldVector<Scalar,Functional::TestVars::template Components<row>::m>(0);
    }

    template <int row, int col>
    auto d2(VariationalArg<Scalar,dim> const& , VariationalArg<Scalar,dim> const& ) const
    {
      return Dune::FieldMatrix<Scalar,Functional::TestVars::template Components<row>::m,Functional::AnsatzVars::template Components<row>::m>(0);
    }
  };

  /**
   * \ingroup boundaryconditions
   * \brief Dirichlet boundary conditions by the penalty approach.
   *
   * For the Dirichlet boundary condition \f$ u|_{\partial \Omega} = u_D \f$, this defines the term
   * \f$ \frac{\gamma}{2} h_F^{-2} |u-u_D|^2 \f$ (and its derivatives), which can be added to the variational
   * functional. The scaling by the face diameter \f$ h_F \f$ is supposed to provide nearly optimal convergence
   * rates for linear finite elements.
   *
   * The dependence on \f$ h_F^{-2} \f$ can lead to ill-conditioned stiffness matrices on mesh refinement.
   * Consider using \ref DirichletNitscheBoundary instead.
   */
  template <class GridView, int components, class ScalarType=typename GridView::ctype>
  class DirichletPenaltyBoundary
  {
    using Scalar = ScalarType;
    static constexpr int dim = GridView::dimension;
    static constexpr int dimworld = GridView::dimensionworld;

  public:

    using Vector = Dune::FieldVector<Scalar,components>;

    /**
     * \brief Moves the boundary condition to a new face.
     *
     * This is called before calling d0,d1,d2 on a new face, and here the
     * weight \f$ h_F^{-2} \f$ is computed once.
     */
    void moveTo(typename GridView::IntersectionIterator const& fi)
    {
      // Compute h^{-2}. h is volume^(1/(dim-1)). Take care of special case for dim=1.
      if (dim > 1)
        weight = std::pow(fi->geometry().volume(),-2/(dim-1));
      else
        weight = 1;
    }

    /**
     * \brief Defines the data for the boundary condition.
     * 
     * Call this from the `evaluateAt` method in the derived class.
     * \param gamma the nominal penalty factor (which is internally multiplied by h^{-2} to ensure convergence).
     * \param u the displacement around which is linearized
     * \param ud the given Dirichlet boundary value
     */
    void setBoundaryData(Scalar gamma_, Vector const& u_, Vector const& ud)
    {
      gamma = gamma_;
      assert(gamma>=0);
      u = u_;
      deviation = u-ud;
    }

    Scalar d0() const { return 0.5*gamma*weight * (deviation*deviation); }

    Vector d1(VariationalArg<Scalar,dim> const& argT) const
    {
      return gamma*weight*argT.value[0]*deviation;
    }

    Dune::FieldMatrix<Scalar,components,components>
    d2(VariationalArg<Scalar,dim> const& argTest, VariationalArg<Scalar,dim> const& argAnsatz) const
    {
      return gamma*weight*argTest.value[0]*argAnsatz.value[0] * unitMatrix<Scalar,components>();
    }

  private:
    Scalar gamma;
    Scalar weight; // this is h^{-2}
    Vector u, deviation;
  };

  /**
   * \ingroup boundaryconditions
   * \brief Dirichlet boundary conditions by Nitsche's method.
   *
   * For the Dirichlet boundary condition \f$ u|_{\partial \Omega} = u_D \f$, this defines the term
   * \f[  \frac{\gamma}{2} h_F^{-1} |u-u_D|^2 - u^T u_x \kappa n\f] (and its derivatives), which can be added to the variational
   * functional. The scaling by the face diameter \f$ h_F \f$ is supposed to provide nearly optimal convergence
   * rates for linear finite elements.
   *
   * This can lead to indefinite stiffness matrices if \f$ \gamma \f$ is not chosen sufficiently large,
   * but does not cause ill-conditioning on mesh refinement, in contrast to \ref DirichletPenaltyBoundary.
   *
   * This follows R. Stenberg: On some techniques for approximating boundary conditions in the finite element method,
   * J Comp. Appl. Math. 63:139-148, 1995.
   *
   * \tparam GridView the grid view on which the problem is defined (used for ctype and face type)
   * \tparam components the number of vectorial components of the FE function
   * \tparam Scalar the type of scalar entries of the FE function
   *
   * \todo Consider extending to M. Juntunen, R. Stenberg: NITSCHE’S METHOD FOR GENERAL BOUNDARY CONDITIONS,
   * Math. Comp. 78(267):1353-1374, 2009.

   */
  template <class GridView, int components, class ScalarType=typename GridView::ctype>
  class DirichletNitscheBoundary
  {
    using Scalar = ScalarType;
    static constexpr int dim = GridView::dimension;
    static constexpr int dimworld = GridView::dimensionworld;

  public:

    using Vector = Dune::FieldVector<Scalar,components>;

    /**
     * \brief Moves the boundary condition to a new face.
     *
     * This is called before calling d0,d1,d2 on a new face, and here the
     * weight \f$ h_F^{-1} \f$ is computed once.
     */
    void moveTo(typename GridView::IntersectionIterator const& fi)
    {
      face = *fi;

      // Compute h^{-1}. h is volume^(1/(dim-1)). Take care of special case for dim=1.
      if (dim > 1)
        weight = std::pow(face.geometry().volume(),-1/(dim-1));
      else
        weight = 1;
    }

    /**
     * \brief Defines the data for the boundary condition.
     *
     * \param xi The local position \f$ \xi \f$ within the face.
     * \param gamma The penalty factor.
     * \param u The current boundary value around which to linearize.
     * \param ud The Dirichlet value.
     * \param ux The current solution's derivative at the boundary.
     * \param kappa The diffusivity tensor.
     */
    void setBoundaryData(Dune::FieldVector<typename GridView::ctype, dim-1> const& xi,
                         Scalar gamma_, Vector const& u, Vector const& ud, Dune::FieldMatrix<Scalar,components,dim> ux_,
                         Dune::FieldMatrix<Scalar,dim,dim> const& kappa)
    {
      gamma = gamma_;
      assert(gamma>=0);
      deviation = u-ud;
      ux = ux_;
      kn = kappa * face.unitOuterNormal(xi);
    }

    Scalar d0() const { return 0.5*weight*gamma*(deviation*deviation) - deviation*(ux*kn); }

    Vector d1(VariationalArg<Scalar,dim> const& arg) const
    {
      return weight*gamma*arg.value[0]*deviation - arg.value[0]*(ux*kn) - deviation*(arg.derivative[0]*kn);
    }

    Dune::FieldMatrix<Scalar,components,components>
    d2(VariationalArg<Scalar,dim> const& argTest, VariationalArg<Scalar,dim> const& argAnsatz) const
    {
      Dune::FieldMatrix<Scalar,components,components> ones(1.0);

      Scalar t1 = argTest.value[0]*(argAnsatz.derivative[0]*kn);
      Scalar t2 = argAnsatz.value[0]*(argTest.derivative[0]*kn);

      return weight*gamma*argTest.value[0]*argAnsatz.value[0] * unitMatrix<Scalar,components>()
             - (t1+t2)*ones;
    }

  private:
    Scalar gamma;
    Scalar weight; // this is h^{-1}
    typename GridView::Intersection face;
    Dune::FieldVector<Scalar,dim> kn; // this is kappa*n
    Vector deviation;
    Dune::FieldMatrix<Scalar,components,dim> ux;
  };


}

#endif


