/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2017-2020 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <algorithm>
#include <cassert>
// #include <chrono>
#include <iostream>
#include <iterator>

#include "dune/grid/config.h"

#include "fem/fixdune.hh"
#include "io/matlab.hh"
#include "linalg/conjugation.hh"
#include "linalg/dynamicMatrixOps.hh"
#include "linalg/matrixProduct.hh"
#include "linalg/qp.hh"
#include "linalg/symmetricOperators.hh"
#include "linalg/direct.hh"
#include "linalg/qpLinesearch.hh"
#include "linalg/threadedMatrix.hh"
#include "linalg/triplet.hh"
#include "linalg/crsutil.hh"
#include "utilities/enums.hh"
#include "utilities/detailed_exception.hh"
#include "utilities/power.hh"
#include "utilities/threading.hh"
#include "utilities/timing.hh"

namespace Kaskade
{
  template <int d, QPStructure sparsity, class Real>
  QPPenalizedSolver<d,sparsity,Real>::QPPenalizedSolver(MatrixA const& A_, MatrixB const& B_)
  : A(A_), B(B_), allCols(B.M())
  {
    std::iota(begin(allCols),end(allCols),0);
  }
  template <int d, QPStructure sparsity, class Real>
  QPPenalizedSolver<d,sparsity,Real>::QPPenalizedSolver(MatrixA const& A_, MatrixB const& B_, Real gamma_)
  : QPPenalizedSolver(A_,B_)
  {
    setPenalty(gamma_);
  }

  template <int d, QPStructure sparsity, class Real>
  QPPenalizedSolver<d,sparsity,Real>& QPPenalizedSolver<d,sparsity,Real>::setMinSteps(int i)
  {
    assert(i>=0);
    minSteps = i;
    return *this;
  }


  template <int d, QPStructure sparsity, class Real>
  QPPenalizedSolver<d,sparsity,Real>& QPPenalizedSolver<d,sparsity,Real>::setMaxSteps(int i)
  {
    assert(i>=0);
    maxSteps = i;
    return *this;
  }

  template <int d, QPStructure sparsity, class Real>
  QPPenalizedSolver<d,sparsity,Real>& QPPenalizedSolver<d,sparsity,Real>::setPenalty(Real g)
  {
    assert(g>=0);
    gamma = g;
    return *this;
  }

  template <int d, QPStructure sparsity, class Real>
  std::tuple<typename QPPenalizedSolver<d,sparsity,Real>::VectorX,typename QPPenalizedSolver<d,sparsity,Real>::VectorB,int>
  QPPenalizedSolver<d,sparsity,Real>::solve(VectorX const& c, VectorB const& b, VectorB const& lambda, double tol) const
  {
    VectorX x(c.N());
    x = 0;
    return solve(c,b,lambda,tol,x);
  }


  // *************************************************************************


  template <int d, QPStructure sparsity, class Real>
  std::tuple<typename QPPenalizedSolver<d,sparsity,Real>::VectorX,int>
  QPPenalizedSolver<d,sparsity,Real>::solve(VectorX const& c, VectorB const& b, double tol) const
  {
    // Solve with initial iterate 0.
    VectorX x(c.N());
    x = 0;
    return solve(c,b,tol,x);
  }

  template <int d, QPStructure sparsity, class Real>
  std::tuple<typename QPPenalizedSolver<d,sparsity,Real>::VectorX,typename QPPenalizedSolver<d,sparsity,Real>::VectorB,int>
  QPPenalizedSolver<d,sparsity,Real>::solve(VectorX const& c, VectorB const& b, VectorB lambda, double tol, VectorX const& x) const
  {
    // Substitute bl <- b - lambda/gamma
    auto bl = b;
    bl.axpy(-1.0/gamma,lambda);

    auto [xNew,iter] = solve(c,bl,tol,x);

    // Compute frist order multiplier update lambda <- (lambda + gamma*(B*x-b))_+, 
    // i.e. lambda <- lambda + max(-lambda, gamma*(B*x-b))
    // but take care that the resulting updated multiplier lambda remains nonegative.
    // This makes only sense, however, if we have actually minimized - otherwise
    // there is no new information in the constraint violation.
    if (iter>0)
    {
      B.mv(xNew,bl);                              // B*x
      bl.axpy(-1,b);                              // B*x-b (don't take B*x-b+lambda/gamma here)
      for (int i=0; i<bl.N(); ++i)
        lambda[i] = std::max(-1.*lambda[i][0],gamma*bl[i][0]);   // compute first order update (stored in lambda)
    }

    return std::make_tuple(xNew,lambda,iter);
  }

  //-----------------------------------------------------------------------------------------------

  template <int d, QPStructure sparsity, class Real>
  std::tuple<typename QPPenalizedSolver<d,sparsity,Real>::VectorX,int>
  QPPenalizedSolver<d,sparsity,Real>::solve(VectorX const& c, VectorB const& b, double tol, VectorX const& x) const
  {
    // The iteration vector.
    auto xa = x;

    std::vector<int> active, oldActive{-1};


    int iter = 0;
    for (; iter<maxSteps; ++iter)
    {
      // check which constraints are active, i.e. Bx-b > 0
      auto bx = b;
      B.mv(xa,bx);
      for (int i=0; i<bx.N(); ++i)
        if (bx[i]>b[i])
          active.push_back(i);

      if (iter>=minSteps && active==oldActive) // we did a Newton step on that active set before - and found the solution
        break;

      // form Hessian H = A + gamma*Ba'*Ba, where Ba is the active matrix,
      auto Ba = B(active,allCols);
      MatrixA H;

      if constexpr (sparsity==QPStructure::DENSE)
        H = A + gamma*transpose(Ba)*Ba;
      else
      {
        auto Ia = sparseUnitMatrix<Real,1>(Ba.N());
        Ia *= gamma;
        H = A + conjugation(Ba,Ia);
      }

      // form derivative A*x+c + gamma*Ba'*(Ba*x-b)
      VectorB tmp(active.size());
      for (int i=0; i<active.size(); ++i)
        tmp[i] = bx[active[i]]-b[active[i]];
      VectorX df = c;
      A.umv(xa,df);
      Ba.usmtv(gamma,tmp,df);

      // solve for Newton step
      VectorX dx(x.N());
      if constexpr (sparsity==QPStructure::DENSE)
        dx = posv(H,df);
      else
      {
        dx = 0;
        MatrixAsTriplet<Real> H_trip(H);                                  // setup for direct solve
        MatrixRepresentedOperator<MatrixAsTriplet<Real>, VectorX, VectorX> Hp_(H_trip);
        auto Hinv = directInverseOperator(Hp_, DirectType::UMFPACK, MatrixProperties::SYMMETRIC);
        Hinv.applyscaleadd(1.0, df, dx);

        H.usmv(-1,dx,df);             // update residual
        Hinv.applyscaleadd(1,df,dx);  // improve solution by iterative refinement
      }

      dx *= -1;

      auto t = qpLinesearch(A,B,gamma,c,b,xa,dx);

      if (t==0)                       // Step length zero (probably bad correction due to rounding errors).
        break;                        // As we don't change anything at all, further iterations are useless.
                                      // (We'd break early next iteration anyways, but can do it already here.)
      // Take step
      xa.axpy(t,dx);

      if (dx.two_norm() < tol*xa.two_norm())  // If the correction drops below the requested tolerance, we stop.
        break;                                // TODO: Remember this is dangerous if convergence is slow.

      std::swap(active,oldActive);
      active.clear();
    }

    return std::make_tuple(xa,maxSteps);
  }

  template <int d, QPStructure sparsity, class Real>
  typename QPPenalizedSolver<d,sparsity,Real>::VectorB
  QPPenalizedSolver<d,sparsity,Real>::multiplierEstimate(VectorX const& x, VectorX c) const
  {
    // Solve B B' \lambda = B (Ax+c).

//     // First compute right hand side r = B(Ax+c)
//     A.umv(x,c);
//     VectorB r(B.N());
//     B.mv(c,r);
//
//     // Now compute B B' and solve the system
//     VectorB lambda(B.N());
//     Real regularize = 100*std::numeric_limits<Real>::epsilon()*B.frobenius_norm2();
//     if constexpr (sparsity==QPStructure::DENSE)
//     {
//       auto BBt = B*transpose(B) + dynamicUnitMatrix(regularize*unitMatrix<Real,1>(),B.N());
//     }
//     else
//     {
//       auto Ib = sparseUnitMatrix<Real,1>(B.N());
//       Ib *= regularize*B.frobenius_norm2();
//       auto BBt = B*transpose(B) + Ib;
//
//       MatrixAsTriplet<Real> BBtTrip(BBt);                                  // setup for direct solve
//       MatrixRepresentedOperator<MatrixAsTriplet<Real>, VectorB, VectorB> BBtOp(BBtTrip);
//       auto invBBt = directInverseOperator(BBtOp, DirectType::UMFPACK, MatrixProperties::SYMMETRIC);
//       invBBt.apply(r,lambda);
//     }



    std::cerr << "NOT YET IMPLEMENTED: QPPenalizedSolver::multiplierEstimate\n";
    abort();
  }
  

  template class QPPenalizedSolver<1,QPStructure::DENSE,double>;
  template class QPPenalizedSolver<2,QPStructure::DENSE,double>;
  template class QPPenalizedSolver<3,QPStructure::DENSE,double>;
  template class QPPenalizedSolver<1,QPStructure::SPARSE,double>;
  template class QPPenalizedSolver<2,QPStructure::SPARSE,double>;
  template class QPPenalizedSolver<3,QPStructure::SPARSE,double>;

  //--------------------------------------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------------------------------------

  template <int d, QPStructure sparsity, class Real>
  QPALSolver<d,sparsity,Real>::QPALSolver(MatrixA const& A, MatrixB const& B)
  : qp(A,B)
  {
    // Estimate a reasonable penalty factor gamma. A reasonable choice is based on the following
    // idea: AL is in essence a gradient method for the Lagrange dual problem. The contraction
    // rate is bounded by r = (Ec-ec)/(Ec+ec), where Ec and ec are maximal and minimal
    // eigenvalues of B(A+gamma B^T B)B^T, see Bertsekas, "Constrained optimization and Lagrange
    // multiplier methods", p 127.
    // We assume B is well-conditioned and roughly estimate Ec = 1/gamma, ec = b^2/(a+gamma*b^2)
    // with a = |A|_F and b = |B|_F (a very coarse estimate). Aiming at r=0.01 then yields
    // gamma = 50 a / b^2. According to first tests, this estimate of the contraction rate is
    // rather optimistic (by a factor of 10), so we take a 10 times larger value.
    Real a = A.frobenius_norm();
    Real b2 = B.frobenius_norm2();
    qp.setPenalty(500*a/b2);
  }
  
  template <int d, QPStructure sparsity, class Real>
  QPALSolver<d,sparsity,Real>& QPALSolver<d,sparsity,Real>::setMinSteps(int i)
  {
    assert(i>=0);
    minSteps = std::max(0,i);
    maxSteps = std::max(minSteps,maxSteps);
    return *this;
  }
  
  template <int d, QPStructure sparsity, class Real>
  QPALSolver<d,sparsity,Real>& QPALSolver<d,sparsity,Real>::setMaxSteps(int i)
  {
    assert(i>=0);
    maxSteps = std::max(0,i);
    minSteps = std::min(minSteps,maxSteps);
    return *this;
  }
  
  template <int d, QPStructure sparsity, class Real>
  QPALSolver<d,sparsity,Real>& QPALSolver<d,sparsity,Real>::setPenalty(Real gamma)
  {
    qp.setPenalty(gamma);
    return *this;
  }
    
  template <int d, QPStructure sparsity, class Real>
  std::tuple<typename QPALSolver<d,sparsity,Real>::VectorX,
             typename QPALSolver<d,sparsity,Real>::VectorB,int> 
  QPALSolver<d,sparsity,Real>::solve(VectorX const& c, VectorB const& b, double tol) const
  {
    VectorX x(c.N()); x = 0;
    return solve(c,b,tol,x);
  }
  
  template <int d, QPStructure sparsity, class Real>
  std::tuple<typename QPALSolver<d,sparsity,Real>::VectorX,
             typename QPALSolver<d,sparsity,Real>::VectorB,int> 
  QPALSolver<d,sparsity,Real>::solve(VectorX const& c, VectorB const& b, double tol, VectorX x) const
  {
    VectorB lambda(b.N()); lambda = 0;
    
    int i = 0;
    for ( ; i<maxSteps; ++i)
    {
      auto [xNew,dLambda,iter] = qp.solve(c,b,lambda,0.0,x);
      x = xNew;
      lambda += dLambda;
      
      if (dLambda.two_norm()<=tol*lambda.two_norm() && i>=minSteps)
        break;
    }
    
    return std::tuple(x,lambda,i);
  }
  
  
  template class QPALSolver<1,QPStructure::DENSE,double>;
  template class QPALSolver<2,QPStructure::DENSE,double>;
  template class QPALSolver<3,QPStructure::DENSE,double>;
  template class QPALSolver<1,QPStructure::SPARSE,double>;
  template class QPALSolver<2,QPStructure::SPARSE,double>;
  template class QPALSolver<3,QPStructure::SPARSE,double>;

  //--------------------------------------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------------------------------------
}
