/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2019-2019 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef SHELLKINEMATICS
#define SHELLKINEMATICS

#include "dune/common/fmatrix.hh"
#include "dune/common/fvector.hh"

#include "fem/diffops/materialLaws.hh"
#include "fem/fixdune.hh"
#include "fem/shapefunctioncache.hh"


namespace Kaskade::Elastomechanics
{
  /**
   * \ingroup diffopsElasto
   * \brief A lower-dimensional kinematics ansatz useful for shell models
   *
   * This represents the spatial derivative \f$ \Phi' = d\Phi / d[x,z] \f$ of the parametrized deformation
   * \f$ \Phi: \Omega\subset \R^{d+1} \to \R^{d+1} \f$ on \f$\Omega = \omega \times \R \f$
   * defined as \f[ \Phi(x,z;p) = \phi(x;p) + z t(x;p) n(x;p) \f]
   * with \f$ \phi, n: \omega\to\R^{d+1} \f$, \f$ t:\omega\to\R \f$, and \f$ n^T\phi' = 0\f$.
   * For the shell case dim=2, \f$ n = \phi_{x_1} \times \phi_{x_2} \f$ holds.
   *
   * The spatial derivative is of the form \f$ \Phi'(x,z;p) = A(x;p) + z B(x;p) \f$ with
   * \f$ A,B\in \R^{(d+1)\times(d+1)} \f$.
   *
   * \tparam dim the dimension of the parametrization domain, usually 2 (or 1 for beams in 2D, currently not implemented)
   * \tparam Real a real field type, usually double
   *
   * \relates ShellKinematics
   */
  template <int dim=2, class Real=double>
  class ShellKinematicsDerivative
  {
  public:
    /**
     * \brief The type of deformation vectors in \f$ \R^{d+1} \f$
     */
    using Vector = Dune::FieldVector<Real,dim+1>;

    using Matrix = Dune::FieldMatrix<Real,dim+1,dim+1>;

    /**
     * \brief Constructor defining a (pretty useless) default linearization point.
     */
    ShellKinematicsDerivative();

    /**
     * \brief Constructor defining a linearization point.
     */
    ShellKinematicsDerivative(Dune::FieldMatrix<Real,dim+1,dim> const& dphi_dx_,
                              Tensor<Real,dim+1,dim,dim> const& ddphi_ddx_,
                              Real t, Dune::FieldVector<Real,dim> const& dt_dx_);

    /**
     * \brief Defines a new evaluation/linearization point.
     */
    void setLinearizationPoint(Dune::FieldMatrix<Real,dim+1,dim> const& dphi_dx_,
                               Tensor<Real,dim+1,dim,dim> const& ddphi_ddx_,
                               Real t, Dune::FieldVector<Real,dim> const& dt_dx_);

    /**
     * \brief Evaluates \f$ (A,B) \f$ making up \f$ \Phi' = A+zB \f$.
     */
    std::pair<Matrix,Matrix> d0() const;

    /**
     * \brief Computes the directional parametrization derivative of \f$ \Phi'(x,z;p) \f$.
     *
     * \f[ \frac{\partial (A,B)}{\partial p} v \f]
     *
     * \param v the parametrization variation direction in which to take the derivative
     * \param i for i=0,...,dim-1, the derivative wrt to the i-th component of \f$ \phi \f$
     *          is taken. For i=dim, the derivative wrt to \f$ t \f$ is computed.
     */
    std::pair<Matrix,Matrix> d1(VariationalArg<Real,dim,1> const& v, int i) const;

    /**
     * \brief Computes the second directional parametrization derivative of \f$ \Phi' \f$.
     *
     * \f[ \frac{\partial^2 (A,B)}{\partial p_i \partial p_j} [v,w], \quad p_i, p_j\in\{\phi,t\} \f]
     *
     * \param v the parametrization variation direction in which to take the first derivative
     * \param i for i=0,...,dim-1, the derivative wrt to the i-th component of \f$ \phi \f$
     *          is taken. For i=dim, the derivative wrt to \f$ t \f$ is computed.
     * \param w the parametrization variation direction in which to take the second derivative
     * \param j for j=0,...,dim-1, the derivative wrt to the j-th component of \f$ \phi \f$
     *          is taken. For j=dim, the derivative wrt to \f$ t \f$ is computed.
     */
    std::pair<Matrix,Matrix> d2(VariationalArg<Real,dim,1> const& v, int i,
                                VariationalArg<Real,dim,1> const& w, int j) const;

  private:
    Dune::FieldVector<Real,dim+1>     phi, n;
    Dune::FieldMatrix<Real,dim+1,dim> dphi_dx, dn_dx;
    Tensor<Real,dim+1,dim,dim>        ddphi_ddx;
    Real                              t;
    Dune::FieldVector<Real,dim>       dt_dx;
  };

  /**
   * \ingroup diffopsElasto
   * \brief A lower-dimensional kinematics ansatz useful for shell models
   *
   * This represents the deformation \f$ \Phi: \Omega\subset \R^{d+1} \to \R^{d+1} \f$
   * on \f$\Omega = \omega \times \R \f$
   * as \f[ \Phi(x,z) = \phi(x) + z t(x) n(x) \f]
   * with \f$ \phi, n: \omega\to\R^{d+1} \f$, \f$ t:\omega\to\R \f$, and \f$ n^T\phi' = 0\f$.
   *
   * \tparam dim the dimension of the parametrization domain, usually 2 (or 1 for beams in 2D)
   * \tparam Real a real field type, usually double
   */
  template <int dim=2, class Real=double>
  class ShellKinematics
  {
  public:
    /**
     * \brief The type of deformation vectors in \f$ \R^{d+1} \f$
     */
    using Vector = Dune::FieldVector<Real,dim+1>;

    ShellKinematics(Vector const& phi, Dune::FieldMatrix<Real,dim+1,dim> const& dphi_dx, Real t);

    /**
     * \brief Returns an object for evaluating \f$ \Phi' \f$ and its parametrization derivatives.
     */
    ShellKinematicsDerivative<dim,Real> derivative() const;

    Vector d0() const
    {
    }

    /**
     * \brief Computes the directional parametrization derivative.
     *
     * \f[ \frac{\partial\Phi}{\partial p} v, \quad p\in\{\phi,t\} \f]
     *
     * \tparam row for `row`=0, the derivative wrt \f$ \phi \f$ is evaluated, for `row`=1 wrt \f$ t \f$
     */
    template <int row>
    auto d1(VariationalArg<Real,dim,1> const& v) const;

  private:
    Dune::FieldVector<Real,dim+1>     phi, n;
    Dune::FieldMatrix<Real,dim+1,dim> dphi_dx;
    Real                              t;
  };

  // ----------------------------------------------------------------------------------------------

  /**
   * \ingroup diffopsElasto
   * \brief A stored energy density class for Taylor-based shell models.
   *
   * Assume on the reference configuration \f$\Omega = \omega\times\mathopen]-\tau,\tau\mathclose[ \f$ with
   * \f$ \omega\subset\R^d \f$ is a deformation \f$ \Phi:\Omega\to\R^{d+1}\f$ defined and there is a
   * stored energy density \f$ W: \R^{(d+1)\times(d+1)}_{\rm sym} \to \R \f$
   * formulated in terms of the Euler-Lagrange strain tensor \f$ E = \frac{1}{2}(C-I) \f$ with
   * \f$ C = \Phi'^T \Phi' \f$. Integrating over \f$ t\in \mathopen]-\tau,\tau\mathclose[ \f$, we define the
   * marginal stored energy density
   * \f[ w:\omega\to\R, \quad w(x) = \int_{-\tau}^\tau W(E(x,t)) \, dt. \f]
   *
   * For Kirchhoff shells (\f$ d=2\f$) with \f$ \Phi(x,t) = \phi(x) + t(x)n(x) \f$ and
   * \f$ n(x) = \phi_{x_0}(x) \times \phi_{x_1}(x) \f$, the derivative is of the form \f$ \Phi'(x,t) = A(x)+tB(x) \f$.
   * Then, the marginal stored energy is
   * \f[ w(x) = 2\tau W(E(x,0)) + \frac{\tau^3}{3} W'(E(x,0))[B^TB] +
   *     \frac{\tau^3}{12} W''(E(x,0))[B^TA+A^TB,B^TA+A^TB] + \mathcal{O}(\tau^5), \f]
   * where \f$ E(x,0) = \frac{1}{2} (A(x)^TA(x)-I) \f$. The first order term represents in-plane strain energy contributions,
   * and the third order term bending energy. If \f$ \tau \f$ is small, the remainder term of
   * order \f$ \tau^5 \f$ can be neglected.
   *
   * This class implements this marginal stored energy density and its derivatives wrt \f$ \phi \f$ and \f$ t \f$
   * (and their spatial derivatives).
   *
   * \tparam Energy the type of stored energy density \f$ W \f$ (currently, only 3D material laws can be used)
   * \tparam Real the scalar field type of real numbers
   *
   * Explicit instantiations are provided for dim=2, Real=double, and Energy=StVenantKirchhoff<3,double>. If you
   * want to instantiate the ShellEnergy for different material laws, include fem/diffops/shellKinematics.hpp,
   * where the template definition is provided.
   */
  template <class Energy=MaterialLaws::StVenantKirchhoff<3,double>,
            class Real=double>
  class ShellEnergy
  {
  public:
    /**
     * \brief The dimension of the parameter domain.
     */
    static int const dim = Energy::dim-1;

    ShellEnergy(Energy const& W, Real tau);

    /**
     * \brief Constructor defining a linearization point.
     * \param W the volume stored energy density
     * \param dphi_dx the spatial derivative of deformation \f$ \phi \f$
     * \param ddphi_ddx the second spatial derivative of deformation \f$ \phi \f$
     * \param t the local shell thickness \f$ t \f$
     * \param dt_dx the spatial derivative of the thickness \f$ t \f$
     */
    ShellEnergy(Energy const& W, Real tau,
                Dune::FieldMatrix<Real,dim+1,dim> const& dphi_dx,
                Tensor<Real,dim+1,dim,dim> const& ddphi_ddx,
                Real t, Dune::FieldVector<Real,dim> const& dt_dx)
    : ShellEnergy(W,tau)
    {
      setLinearizationPoint(dphi_dx,ddphi_ddx,t,dt_dx);
    }

    /**
     * \brief Defines a new evaluation/linearization point.
     */
    void setLinearizationPoint(Dune::FieldMatrix<Real,dim+1,dim> const& dphi_dx,
                               Tensor<Real,dim+1,dim,dim> const& ddphi_ddx,
                               Real t, Dune::FieldVector<Real,dim> const& dt_dx);

    /**
     * \brief Evaluates the marginal stored energy density \f$ w(x) \f$.
     */
    Real d0() const;

    /**
     * \brief Evaulates the parametric directional derivative of the marginal stored energy density \f$ w \f$ wrt changes.
     * in \f$ \phi\f$ and \f$ t \f$.
     *
     * \tparam row either 0 or 1. For row=0, the derivative wrt \f$ \phi \f$ is computed, for row=1 the
     *             derivative wrt \f$ t \f$.
     * \param v the direction in which to take the derivative
     *
     * \return for row=0, a dim+1 vector \f$ r_i \f$ with the derivative values wrt \f$ \phi_i \f$,
     *         for row=1, a scalar, the derivative wrt \f$ t \f$
     */
    template <int row>
    Dune::FieldVector<Real,row==0?dim+1:1> d1(VariationalArg<Real,dim,1> const& v) const;

    /**
     * \brief Evaulates the second parametric directional derivative of the marginal stored energy density \f$ w \f$ wrt changes.
     * in \f$ \phi\f$ and \f$ t \f$.
     */
    template <int row, int col>
    Dune::FieldMatrix<Real,row==0?dim+1:1,col==0?dim+1:1> d2(VariationalArg<Real,dim,1> const& v,
                                                             VariationalArg<Real,dim,1> const& w) const;

  private:
    Energy                              W;
    Real                                tau;
    ShellKinematicsDerivative<dim,Real> dPhi;

    using Matrix = typename ShellKinematicsDerivative<dim,Real>::Matrix;

    Matrix A, B, E, BtB, BtA, BtA_AtB;

    // compute first directional derivative of energy, given the first directional derivatives of contributing matrices A and B
    Real d1(std::pair<Matrix,Matrix> const& dAdB) const;

    // compute second directional derivative of energy, given the first and second directional derivatives of contributing matrices A and B
    // first derivatives need to be given in both directions, whereas the second derivative needs only one value
    // (order of differentiation irrelevant)
    Real d2(std::pair<Matrix,Matrix> const& dAdB1, std::pair<Matrix,Matrix> const& dAdB2,
            std::pair<Matrix,Matrix> const& ddAddB) const;
  };
}

#endif
