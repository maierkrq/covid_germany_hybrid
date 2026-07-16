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

#ifndef SHELLKINEMATICS_HPP
#define SHELLKINEMATICS_HPP

#include "fem/diffops/shellKinematics.hh"

namespace Kaskade::Elastomechanics
{

  namespace ShellKinematics_Detail
  {
    template <class Matrix>
    Matrix symmetricSum(Matrix const& A)
    {
      return A + transpose(A);
    }
  }

  template <class Energy, class Real>
  ShellEnergy<Energy,Real>::ShellEnergy(Energy const& W_, Real tau_)
  : W(W_)
  , tau(tau_)
  {
    static_assert(dim==2);
  }


  template <class Energy, class Real>
  void ShellEnergy<Energy,Real>::setLinearizationPoint(Dune::FieldMatrix<Real,dim+1,dim> const& dphi_dx,
                                                       Tensor<Real,dim+1,dim,dim> const& ddphi_ddx,
                                                       Real t, Dune::FieldVector<Real,dim> const& dt_dx)
  {
    using namespace ShellKinematics_Detail;

    dPhi.setLinearizationPoint(dphi_dx,ddphi_ddx,t,dt_dx);

    std::tie(A,B) = dPhi.d0();
    E = 0.5*(transpose(A)*A-unitMatrix<Real,dim+1>());
    BtB = transpose(B)*B;
    BtA = transpose(B)*A;
    BtA_AtB = symmetricSum(BtA);

    W.setLinearizationPoint(E);
  }

  template <class Energy, class Real>
  Real ShellEnergy<Energy,Real>::d0() const
  {
    Real const tau3 = tau*tau*tau;
    return 2*tau*W.d0() + /*0.5*/1.0/3.0*tau3*W.d1(BtB) + tau3/12*W.d2(BtA_AtB,BtA_AtB);
  }

  template <class Energy, class Real>
  Real ShellEnergy<Energy,Real>::d1(std::pair<Matrix,Matrix> const& dAdB) const
  {
    using namespace ShellKinematics_Detail;
    auto const& [dA,dB] = dAdB;

    Real const tau3 = tau*tau*tau;
    auto dE = 0.5*symmetricSum(transpose(dA)*A);
    auto dBtB = symmetricSum(transpose(dB)*B);
    auto dBtA = transpose(dB)*A + transpose(B)*dA;
    auto dBtA_AtB = symmetricSum(dBtA);

    return 2*tau*W.d1(dE) + /*0.5*/1.0/3.0*tau3*(W.d2(dE,BtB)+W.d1(dBtB)) + tau3/12*(W.d3(dE,BtA_AtB,BtA_AtB)+2*W.d2(BtA_AtB,dBtA_AtB));
  }

  template <class Energy, class Real>
  Real ShellEnergy<Energy,Real>::d2(std::pair<Matrix,Matrix> const& dAdB1, std::pair<Matrix,Matrix> const& dAdB2,
                                    std::pair<Matrix,Matrix> const& ddAddB) const
  {
    using namespace ShellKinematics_Detail;
    auto const& [dA1,dB1] = dAdB1;
    auto const& [dA2,dB2] = dAdB2;
    auto const& [ddA,ddB] = ddAddB;

    Real const tau3 = tau*tau*tau;
    auto dE1 = 0.5*symmetricSum(transpose(dA1)*A);
    auto dE2 = 0.5*symmetricSum(transpose(dA2)*A);
    auto ddE = 0.5*symmetricSum(transpose(ddA)*A + transpose(dA1)*dA2);

    auto dBtB1 = symmetricSum(transpose(dB1)*B);
    auto dBtB2 = symmetricSum(transpose(dB2)*B);
    auto ddBtB = symmetricSum(transpose(ddB)*B + transpose(dB1)*dB2);

    auto dBtA1 = transpose(dB1)*A + transpose(B)*dA1;
    auto dBtA2 = transpose(dB2)*A + transpose(B)*dA2;
    auto ddBtA = transpose(ddB)*A + transpose(dB1)*dA2 + transpose(dB2)*dA1 + transpose(B)*ddA;

    auto dBtA_AtB1 = dBtA1 + transpose(dBtA1);
    auto dBtA_AtB2 = dBtA2 + transpose(dBtA2);
    auto ddBtA_AtB = ddBtA + transpose(ddBtA);

    return 2*tau*(W.d2(dE1,dE2) + W.d1(ddE))
           + /*0.5*/1.0/3.0*tau3*( W.d3(dE1,BtB,dE2)+W.d2(ddE,BtB)+W.d2(dE1,dBtB2) + W.d2(dBtB1,dE2)+W.d1(ddBtB) )
           + tau3/12*( W.d4(dE1,BtA_AtB,BtA_AtB,dE2)+W.d3(ddE,BtA_AtB,BtA_AtB)+W.d3(dE1,dBtA_AtB2,BtA_AtB)+W.d3(dE1,BtA_AtB,dBtA_AtB2)
                    + 2*(W.d3(BtA_AtB,dBtA_AtB1,dE2)+W.d2(dBtA_AtB2,dBtA_AtB1)+W.d2(BtA_AtB,ddBtA_AtB))  );
  }

  template <class Energy, class Real>
  template <int row>
  Dune::FieldVector<Real,row==0?ShellEnergy<Energy,Real>::dim+1:1> ShellEnergy<Energy,Real>::d1(VariationalArg<Real,dim,1> const& v) const
  {
    // implementation verified by numerical differentiation
    if constexpr (row==0)
    {
      Dune::FieldVector<Real,dim+1> r;
      for (int i=0; i<dim+1; ++i)
        r[i] = d1(dPhi.d1(v,i));
      return r;
    }
    else
      return d1(dPhi.d1(v,dim+1));
  }



  template <class Energy, class Real>
  template <int row, int col>
  Dune::FieldMatrix<Real,row==0?ShellEnergy<Energy,Real>::dim+1:1,col==0?ShellEnergy<Energy,Real>::dim+1:1>
  ShellEnergy<Energy,Real>::d2(VariationalArg<Real,dim,1> const& v, VariationalArg<Real,dim,1> const& w) const
  {
    // results have been verified by numerical differentiation (but only for quadratic stored energy)
    if constexpr (row==0)
      if constexpr (col==0)
      {
        Dune::FieldMatrix<Real,dim+1,dim+1> r;
        for (int i=0; i<dim+1; ++i)
          for (int j=0; j<dim+1; ++j)
            r[i][j] = d2(dPhi.d1(v,i),dPhi.d1(w,j),dPhi.d2(v,i,w,j));
        return r;
      }
      else
      {
        Dune::FieldMatrix<Real,dim+1,1> r;
        for (int i=0; i<dim+1; ++i)
          r[i][0] = d2(dPhi.d1(v,i),dPhi.d1(w,dim+1),dPhi.d2(v,i,w,dim+1));
        return r;
      }
    else
      if constexpr (col==0)
      {
        Dune::FieldMatrix<Real,1,dim+1> r;
        for (int j=0; j<dim+1; ++j)
          r[0][j] = d2(dPhi.d1(v,dim+1),dPhi.d1(w,j),dPhi.d2(v,dim+1,w,j));
        return r;
      }
      else
        return d2(dPhi.d1(v,dim+1),dPhi.d1(w,dim+1),dPhi.d2(v,dim+1,w,dim+1));
  }

}
#endif
