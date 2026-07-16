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

#include "fem/diffops/materialLaws.hh"
#include "fem/diffops/shellKinematics.hpp"

namespace Kaskade::Elastomechanics
{
  template <int dim, class Real>
  using Result = std::pair<typename ShellKinematicsDerivative<dim,Real>::Matrix,
                           typename ShellKinematicsDerivative<dim,Real>::Matrix>;

  template <int dim, class Real>
  ShellKinematicsDerivative<dim,Real>::ShellKinematicsDerivative()
  : dphi_dx(0.0)
  , ddphi_ddx(0.0)
  , t(1.0)
  , dt_dx(0.0)
  {
    static_assert(dim==2); // currently only implemented for dim=2.
    dphi_dx[0][0] = 1;
    dphi_dx[1][1] = 1;
  }

  template <int dim, class Real>
  ShellKinematicsDerivative<dim,Real>::ShellKinematicsDerivative(Dune::FieldMatrix<Real,dim+1,dim> const& dphi_dx,
                                                                 Tensor<Real,dim+1,dim,dim> const& ddphi_ddx,
                                                                 Real t, Dune::FieldVector<Real,dim> const& dt_dx)
  {
    static_assert(dim==2); // currently only implemented for dim=2.
    setLinearizationPoint(dphi_dx,ddphi_ddx,t,dt_dx);
  }

  template <int dim, class Real>
  void ShellKinematicsDerivative<dim,Real>::setLinearizationPoint(Dune::FieldMatrix<Real,dim+1,dim> const& dphi_dx_,
                                                                  Tensor<Real,dim+1,dim,dim> const& ddphi_ddx_,
                                                                  Real t_, Dune::FieldVector<Real,dim> const& dt_dx_)
  {
    using TensorAccess::_;

    dphi_dx = dphi_dx_;
    ddphi_ddx = ddphi_ddx_;
    t = t_;
    dt_dx = dt_dx_;

    n = vectorProduct(column(dphi_dx,0),column(dphi_dx,1));
    Vector n0 = vectorProduct(ddphi_ddx(_,0,0),column(dphi_dx,1)) + vectorProduct(column(dphi_dx,0),ddphi_ddx(_,1,0));
    Vector n1 = vectorProduct(ddphi_ddx(_,0,1),column(dphi_dx,1)) + vectorProduct(column(dphi_dx,0),ddphi_ddx(_,1,1));
    dn_dx = horzcat(n0,n1);
  }


  template <int dim, class Real>
  Result<dim,Real> ShellKinematicsDerivative<dim,Real>::d0() const
  {
    Matrix A = horzcat(column(dphi_dx,0),           column(dphi_dx,1),            t*n);
    Matrix B = horzcat(t*column(dn_dx,0)+dt_dx[0]*n,t*column(dn_dx,1)+dt_dx[1]*n, 0*n);
    return std::make_pair(A,B);
  }

  template <int dim, class Real>
  Result<dim,Real> ShellKinematicsDerivative<dim,Real>::d1(VariationalArg<Real,dim,1> const& v, int i) const
  {
    using TensorAccess::_;

    if (i<dim+1) // results verified by numerical differentiation
    {
      // algorithmic differentiation by chain rule
      auto dphi_dx_dv = 0.0*dphi_dx; dphi_dx_dv[i] = v.derivative[0];
      Tensor<Real,dim+1,dim,dim> ddphi_ddx_dv = 0.0*ddphi_ddx; ddphi_ddx_dv[i] = v.hessian[0];

      // Compute the first derivative of the normal.
      auto dn_dv = vectorProduct(column(dphi_dx_dv,0),column(dphi_dx,1))
      + vectorProduct(column(dphi_dx,0),column(dphi_dx_dv,1));
      auto dn0_dv = vectorProduct(ddphi_ddx_dv(_,0,0),column(dphi_dx,1)) + vectorProduct(ddphi_ddx(_,0,0),column(dphi_dx_dv,1))
      + vectorProduct(column(dphi_dx_dv,0),ddphi_ddx(_,1,0)) + vectorProduct(column(dphi_dx,0),ddphi_ddx_dv(_,1,0));
      auto dn1_dv = vectorProduct(ddphi_ddx_dv(_,0,1),column(dphi_dx,1)) + vectorProduct(ddphi_ddx(_,0,1),column(dphi_dx_dv,1))
      + vectorProduct(column(dphi_dx_dv,0),ddphi_ddx(_,1,1)) + vectorProduct(column(dphi_dx,0),ddphi_ddx_dv(_,1,1));
      auto dn_dx_dv = horzcat(dn0_dv,dn1_dv);

      Matrix dA = horzcat(column(dphi_dx_dv,0), column(dphi_dx_dv,1), t*dn_dv);
      Matrix dB = horzcat(t*column(dn_dx_dv,0)+dt_dx[0]*dn_dv,
                          t*column(dn_dx_dv,1)+dt_dx[1]*dn_dv,
                          0.0*n);
      return std::make_pair(dA,dB);
    }
    else // results verified by numerical differentiation
    {
      Real dt_dv = v.value[0];
      auto dt_dx_dv = v.derivative[0];

      Matrix dA = horzcat(Dune::FieldMatrix<Real,dim+1,dim>(0.0), dt_dv*n);
      Matrix dB = horzcat(dt_dv*column(dn_dx,0)+dt_dx_dv[0]*n,dt_dv*column(dn_dx,1)+dt_dx_dv[1]*n, 0*n);
      return std::make_pair(dA,dB);
    }
  }

  template <int dim, class Real>
  Result<dim,Real> ShellKinematicsDerivative<dim,Real>::d2(VariationalArg<Real,dim,1> const& v, int i,
                                                           VariationalArg<Real,dim,1> const& w, int j) const
  {
    using TensorAccess::_;

    if (i<dim+1 && j<dim+1) // results verified by numerical differentiation
    {
      // Compute second derivative wrt phi. In phi, Phi is quadratic, i.e. the second derivative is constant
      // (but depends on t, of course).

      auto dphi_dx_dv = 0.0*dphi_dx; dphi_dx_dv[i] = v.derivative[0]; // obviously the second derivative
      auto dphi_dx_dw = 0.0*dphi_dx; dphi_dx_dw[j] = w.derivative[0]; // dphi_dx_dv_dw vanishes
      Tensor<Real,dim+1,dim,dim> ddphi_ddx_dv = 0.0*ddphi_ddx; ddphi_ddx_dv[i] = v.hessian[0]; // same here ... this simplifies
      Tensor<Real,dim+1,dim,dim> ddphi_ddx_dw = 0.0*ddphi_ddx; ddphi_ddx_dw[j] = w.hessian[0]; // the remaining expressions a lot

      // Compute the required derivatives of the normal. It turns out that only the second derivative
      // is needed. The first derivatives are stated (for reference), but are commented out as they are
      // not required. TODO: exploit the special sparsity of dphi_dx_dv and dphi_dx_dw to short-cut the v
      //                     vector products below.
      // auto dn_dv = vectorProduct(column(dphi_dx_dv,0),column(dphi_dx,1))
      //            + vectorProduct(column(dphi_dx,0),column(dphi_dx_dv,1));
      // auto dn_dw = vectorProduct(column(dphi_dx_dw,0),column(dphi_dx,1))
      //            + vectorProduct(column(dphi_dx,0),column(dphi_dx_dw,1));
      auto dn_dv_dw = vectorProduct(column(dphi_dx_dv,0),column(dphi_dx_dw,1))
      + vectorProduct(column(dphi_dx_dw,0),column(dphi_dx_dv,1));

      // auto dn0_dv = vectorProduct(ddphi_ddx_dv(_,0,0),column(dphi_dx,1)) + vectorProduct(ddphi_ddx(_,0,0),column(dphi_dx_dv,1))
      //             + vectorProduct(column(dphi_dx_dv,0),ddphi_ddx(_,1,0)) + vectorProduct(column(dphi_dx,0),ddphi_ddx_dv(_,1,0));
      // auto dn0_dw = vectorProduct(ddphi_ddx_dw(_,0,0),column(dphi_dx,1)) + vectorProduct(ddphi_ddx(_,0,0),column(dphi_dx_dw,1))
      //             + vectorProduct(column(dphi_dx_dw,0),ddphi_ddx(_,1,0)) + vectorProduct(column(dphi_dx,0),ddphi_ddx_dw(_,1,0));
      auto dn0_dv_dw = vectorProduct(ddphi_ddx_dv(_,0,0),column(dphi_dx_dw,1)) + vectorProduct(ddphi_ddx_dw(_,0,0),column(dphi_dx_dv,1))
      + vectorProduct(column(dphi_dx_dv,0),ddphi_ddx_dw(_,1,0)) + vectorProduct(column(dphi_dx_dw,0),ddphi_ddx_dv(_,1,0));

      // auto dn1_dv = vectorProduct(ddphi_ddx_dv(_,0,1),column(dphi_dx,1)) + vectorProduct(ddphi_ddx(_,0,1),column(dphi_dx_dv,1))
      //             + vectorProduct(column(dphi_dx_dv,0),ddphi_ddx(_,1,1)) + vectorProduct(column(dphi_dx,0),ddphi_ddx_dv(_,1,1));
      // auto dn1_dw = vectorProduct(ddphi_ddx_dw(_,0,1),column(dphi_dx,1)) + vectorProduct(ddphi_ddx(_,0,1),column(dphi_dx_dw,1))
      //             + vectorProduct(column(dphi_dx_dw,0),ddphi_ddx(_,1,1)) + vectorProduct(column(dphi_dx,0),ddphi_ddx_dw(_,1,1));
      auto dn1_dv_dw = vectorProduct(ddphi_ddx_dv(_,0,1),column(dphi_dx_dw,1)) + vectorProduct(ddphi_ddx_dw(_,0,1),column(dphi_dx_dv,1))
      + vectorProduct(column(dphi_dx_dv,0),ddphi_ddx_dw(_,1,1)) + vectorProduct(column(dphi_dx_dw,0),ddphi_ddx_dv(_,1,1));

      // auto dn_dx_dv = horzcat(dn0_dv,dn1_dv);
      // auto dn_dx_dw = horzcat(dn0_dw,dn1_dw);
      auto dn_dx_dv_dw = horzcat(dn0_dv_dw,dn1_dv_dw);

      Matrix ddA = horzcat(Dune::FieldMatrix<Real,dim+1,2>(0.0), t*dn_dv_dw);
      Matrix ddB = horzcat(t*column(dn_dx_dv_dw,0)+dt_dx[0]*dn_dv_dw,
                           t*column(dn_dx_dv_dw,1)+dt_dx[1]*dn_dv_dw,
                           0.0*n);
      return std::make_pair(ddA,ddB);
    }
    else if (i<dim+1 && j==dim+1) // results verified by numerical differentiation
    {
      // Compute mixed second derivative, first wrt phi, then wrt t
      auto dphi_dx_dv = 0.0*dphi_dx; dphi_dx_dv[i] = v.derivative[0]; // obviously the second derivative
      Tensor<Real,dim+1,dim,dim> ddphi_ddx_dv = 0.0*ddphi_ddx; ddphi_ddx_dv[i] = v.hessian[0]; // same here ... this simplifies

      Real dt_dw = w.value[0];
      auto dt_dx_dw = w.derivative[0];

      // Compute the first derivative of the normal.
      auto dn_dv = vectorProduct(column(dphi_dx_dv,0),column(dphi_dx,1))
      + vectorProduct(column(dphi_dx,0),column(dphi_dx_dv,1));
      auto dn0_dv = vectorProduct(ddphi_ddx_dv(_,0,0),column(dphi_dx,1)) + vectorProduct(ddphi_ddx(_,0,0),column(dphi_dx_dv,1))
      + vectorProduct(column(dphi_dx_dv,0),ddphi_ddx(_,1,0)) + vectorProduct(column(dphi_dx,0),ddphi_ddx_dv(_,1,0));
      auto dn1_dv = vectorProduct(ddphi_ddx_dv(_,0,1),column(dphi_dx,1)) + vectorProduct(ddphi_ddx(_,0,1),column(dphi_dx_dv,1))
      + vectorProduct(column(dphi_dx_dv,0),ddphi_ddx(_,1,1)) + vectorProduct(column(dphi_dx,0),ddphi_ddx_dv(_,1,1));
      auto dn_dx_dv = horzcat(dn0_dv,dn1_dv);

      Matrix dA = horzcat(Dune::FieldMatrix<Real,dim+1,2>(0.0), dt_dw*dn_dv);
      Matrix dB = horzcat(dt_dw*column(dn_dx_dv,0)+dt_dx_dw[0]*dn_dv,
                          dt_dw*column(dn_dx_dv,1)+dt_dx_dw[1]*dn_dv,
                          0.0*n);
      return std::make_pair(dA,dB);


    }
    else if (i==dim+1 && j<dim+1)
    {
      // Compute mixed second derivative, first wrt phi, then wrt t
      // second derivatives are symmetric - exploit this
      return d2(w,j,v,i);
    }
    else // results verified by numerical differentiation
    {
      // dPhi is linear in t, thus the second derivative wrt t vanishes.
      Matrix zero(0.0);
      return std::make_pair(zero,zero);
    }
  }

  // explicit instantiation
  template class ShellKinematicsDerivative<2,double>;

  // ----------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------

  using Energy = MaterialLaws::StVenantKirchhoff<3,double>;

  // explicit instantiation of ShellEnergy
  template class ShellEnergy<Energy,double>;
  template Dune::FieldVector<double,3> ShellEnergy<Energy,double>::d1<0>(VariationalArg<double,2,1> const&) const;
  template Dune::FieldVector<double,1> ShellEnergy<Energy,double>::d1<1>(VariationalArg<double,2,1> const&) const;
  template Dune::FieldMatrix<double,3,3> ShellEnergy<Energy,double>
            ::d2<0,0>(VariationalArg<double,2,1> const&,VariationalArg<double,2,1> const&) const;
  template Dune::FieldMatrix<double,3,1> ShellEnergy<Energy,double>
            ::d2<0,1>(VariationalArg<double,2,1> const&,VariationalArg<double,2,1> const&) const;
  template Dune::FieldMatrix<double,1,3> ShellEnergy<Energy,double>
            ::d2<1,0>(VariationalArg<double,2,1> const&,VariationalArg<double,2,1> const&) const;
  template Dune::FieldMatrix<double,1,1> ShellEnergy<Energy,double>
            ::d2<1,1>(VariationalArg<double,2,1> const&,VariationalArg<double,2,1> const&) const;

}

