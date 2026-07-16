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
#include <chrono>
#include <iostream>
#include <iterator>
#include <stdlib.h>     // for rand()

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

  template <class MatrixA, class MatrixB, class VectorX, class VectorB>
  std::array<typename ScalarTraits<typename MatrixA::field_type>::Real,4>
  checkKKT(MatrixA const& A, MatrixB const& B,VectorX const& c, VectorB const& b,
           VectorX const& x, VectorB const& lambda)
  {
    using Real = typename ScalarTraits<typename MatrixA::field_type>::Real;
    std::array<Real,4> result{0,0,0,0};

    int const m = b.N();

    // compute adjoint residual Ax+c+B^T*lambda
    VectorX r = c;
    A.umv(x,r);
    B.usmtv(1.0,lambda,r);
    result[0] = r.two_norm();

    // compute the constraint residual (b-Bx)_-
    // equivalent to (Bx-b)_+
    VectorB g = b;
    B.usmv(-1.0,x,g);
    for (int i=0; i<m; ++i)
      if (g[i][0] < 0)
        result[1] += g[i].two_norm2();
    result[1] = std::sqrt(result[1]);

    // compute the multiplier "residual" lambda_-
    for (int i=0; i<m; ++i)
      if (lambda[i][0] < 0)
        result[2] += lambda[i].two_norm2();
    result[2] = std::sqrt(result[2]);

    // compute complementarity (Bx-b)^T*lambda
    result[3] = -(g*lambda);

    return result;
  }

  // explicit instantiation
  template <int d, int s>
  using DenseMatrix = DynamicMatrix<Dune::FieldMatrix<double,d,s>>;

  template <int d, int s>
  using SparseMatrix = NumaBCRSMatrix<Dune::FieldMatrix<double,d,s>>;

  template <int d>
  using DenseVector = Dune::BlockVector<Dune::FieldVector<double,d>>;

  template
  std::array<double,4> checkKKT(DenseMatrix<1,1> const& A, DenseMatrix<1,1> const& B,
                                DenseVector<1> const& c, DenseVector<1> const& b,
                                DenseVector<1> const& x, DenseVector<1> const& lambda);
  template
  std::array<double,4> checkKKT(DenseMatrix<2,2> const& A, DenseMatrix<1,2> const& B,
                                DenseVector<2> const& c, DenseVector<1> const& b,
                                DenseVector<2> const& x, DenseVector<1> const& lambda);
  template
  std::array<double,4> checkKKT(DenseMatrix<3,3> const& A, DenseMatrix<1,3> const& B,
                                DenseVector<3> const& c, DenseVector<1> const& b,
                                DenseVector<3> const& x, DenseVector<1> const& lambda);
  template
  std::array<double,4> checkKKT(SparseMatrix<1,1> const& A, SparseMatrix<1,1> const& B,
                                DenseVector<1> const& c, DenseVector<1> const& b,
                                DenseVector<1> const& x, DenseVector<1> const& lambda);
  template
  std::array<double,4> checkKKT(SparseMatrix<2,2> const& A, SparseMatrix<1,2> const& B,
                                DenseVector<2> const& c, DenseVector<1> const& b,
                                DenseVector<2> const& x, DenseVector<1> const& lambda);
  template
  std::array<double,4> checkKKT(SparseMatrix<3,3> const& A, SparseMatrix<1,3> const& B,
                                DenseVector<3> const& c, DenseVector<1> const& b,
                                DenseVector<3> const& x, DenseVector<1> const& lambda);


  //-----------------------------------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------------
  template <class Real, class Implementation>
  Implementation& QPBoundSolverBase<Real,Implementation>::setMinSteps(int count)
  {
    assert(count >= 0);
    minSteps = count;
    return static_cast<Implementation&>(*this);
  }

  template <class Real, class Implementation>
  Implementation& QPBoundSolverBase<Real,Implementation>::setMaxSteps(int count)
  {
    assert(count >= 0);
    maxSteps = count;
    return static_cast<Implementation&>(*this);
  }


  template <class Real, class Implementation>
  std::tuple<typename QPBoundSolverBase<Real,Implementation>::Vector,int>
  QPBoundSolverBase<Real,Implementation>::solve(Vector const& c, Vector const& b, double tol) const
  {
    // Create initial guess. Tempting values are 0 (but this may be infeasible if b[i]<0 for some i)
    // and b (which is feasible, but can be arbitrarily bad if unconstrained variables are denoted
    // by b[i] = +inf). Thus we start with min(0,b).
    Vector x(dimension);
    x = 0;
    return solve(x,c,b,tol);
  }

  // TODO: Cf. Nocedal/Wright 2006 (Chap 16.6) or Axehill/Hansson 2008
  template <class Real, class Implementation>
  std::tuple<typename QPBoundSolverBase<Real,Implementation>::Vector,int>
  QPBoundSolverBase<Real,Implementation>::solve(Vector const& x_, Vector const& c, Vector const& b, double tol) const
  {
    assert(x_.N()==dimension);
    assert(c.N()==dimension);
    assert(b.N()==dimension);
    
    checkNaN(c);
    checkNaN(b);
    
    //  solve  QP min 0.5*x'*A*x + c'*x   s.t. x<=b
    Vector g(dimension), s(dimension), As(dimension);
    
    Vector x(dimension);
    for (int i=0; i<dimension; ++i)
      x[i] = std::min(x_[i][0],b[i][0]);
    
    Real normA = static_cast<Implementation const&>(*this).frobenius_norm();
    
    // Compute value of objective
    static_cast<Implementation const&>(*this).mv(x,As);
    As *= 0.5;
    As += c;
    Real objective = x*As;
    
    int iter = 0;
    while (true)
    {
      ++iter;

      double errorlevel = std::numeric_limits<Real>::epsilon() * (normA*x.two_norm() + c.two_norm());

      // Start with a Newton step on the currently inactive variables and perform a line search
      // (because the full Newton step could hit a constraint).
      getProjectedNewton(x,c,b,g,s);
      bool newtonEncounteredConstraint = projectedLineSearch(x,b,g,s,&c);
      static_cast<Implementation const&>(*this).mv(s,As);
      Real stepEnergy = s*As;

      // Next compute the projected gradient. If all active constraints remain active on
      // following the negative projected gradient direction, and the Newton step was taken full,
      // we have found the solution:
      // We came here by a Newton step on the inactive variables, such that their gradients are zero.
      // All the active variables have vanishing projected gradients as well.
      // TODO: This assumes Newton steps are exact. Provide a means to allow inexact solves and then skip this exit path.
      //
      // As we have found the exact solution, waiting for minSteps to be reached makes no sense.
      bool pendingInactivity = getProjectedGradient(x,c,b,g,s);
      if (!newtonEncounteredConstraint && !pendingInactivity)
      {
// std::cout << "terminating QP due to successful Newton step after " << iter << " iterations\n";
        break;
      }
      
      static_cast<Implementation const&>(*this).mv(x,As);
      As *= 0.5;
      As += c;
      Real newObjective = x*As;
      if (newObjective >= objective) // no progress - probably due to roundoff error 
      {
        std::cout << "terminating QP due to nonmonotonicity after " << iter << " iterations: Jold=" << objective 
                  << ", Jnew=" << newObjective << "\n";
        getProjectedGradient(x,c,b,g,s);
        std::cout << "(Gradient entries commented out)" << std::endl;
        //std::cout << "gradient: " << g << "\nprojected gradient: " << s << "\n";
        break;
      }
      objective = newObjective;

      bool gradientEncounteredConstraint = projectedLineSearch(x,b,g,s,&c);

      if( (iter>=minSteps && stepEnergy < tol && !gradientEncounteredConstraint)  || iter>=maxSteps )
      {
        if (iter>=maxSteps)
          std::cerr << "step energy: " << stepEnergy << " [limit: " << tol
                    << "], constraint encountered on gradient linesearch: " << gradientEncounteredConstraint 
                    << ", steps: " << iter << "\n";
if (stepEnergy < tol)
   std::cout << "terminating QP solver due to step energy " << stepEnergy << " < " << tol  << "\n";
        break;
      }
      
      if (iter >= 5000)
      {
        std::cerr.precision(15);
        std::cerr << "x = \n" << x << "\n";
        std::cerr << "s was \n" << s << " and ||s||_2 = " << s.two_norm() << "\n----------\n";
        if (iter > 5020) 
        {
          std::cerr << "c =[\n" << c << "]\nb = [\n" << b << "]\nerrorlevel = " << errorlevel << "\n";
          abort();
        }
      }
    }   

    return std::make_tuple(x,iter);
  }
  
  template <class Real, class Implementation>
  void QPBoundSolverBase<Real,Implementation>::gradientStep(Vector const& c, Vector const& b, Vector& x) const
  {
    assert(c.N()==dimension);
    assert(b.N()==dimension);
    assert(x.N()==dimension);
    
    checkNaN(c);
    checkNaN(b);
    checkNaN(x);
    
    for (int i=0; i<dimension; ++i)
      assert(x[i][0] <= b[i][0]);
    
    Vector g(dimension), s(dimension);
    getProjectedGradient(x,c,b,g,s);
    projectedLineSearch(x,b,g,s);
  }
  
  namespace 
  {
    template <class Real>
    bool active(Real a, Real b, Real epsFact)
    {
      using namespace std;
      constexpr Real eps = numeric_limits<Real>::epsilon();
      double inf = std::numeric_limits<Real>::infinity();
      if( b == inf)
        return false;
      else
        return b-a <= epsFact*eps*max(abs(a),abs(b));
    }
  }
  
  template <class Real, class Implementation>
  bool QPBoundSolverBase<Real,Implementation>::getProjectedGradient(Vector const& x, Vector const& c,
                                                                    Vector const& b, Vector& g, Vector& s) const
  {
    // gradient g = A*x+c
    static_cast<Implementation const&>(*this).mv(x,g);
    checkNaN(g);
    g += c;
    
    // project the negative gradient at active components not to point outside the negative quadrant
    // this gives the search direction
    bool pendingInactivity = false;
    for (int i=0; i<g.N(); ++i)
      if (active(x[i][0],b[i][0],Real(10)))
        if (g[i]<=0)                    // this active variable will remain active
          s[i] = 0;
        else                            // but this not. Flag the pending switch from active to inactive.
        {
          s[i] = -g[i];
          pendingInactivity = true;
        }
      else
        s[i] = -g[i];

    return pendingInactivity;
  }
  
  // Compute Newton step dx on the subspace of inactive components of x
  template <class Real, class Implementation>
  void QPBoundSolverBase<Real,Implementation>::getProjectedNewton(Vector const& x, Vector const& c, Vector const& b, Vector& g, Vector& dx) const
  {
    checkNaN(x);
    checkNaN(c);
    
    // gradient g = A*x+c
    static_cast<Implementation const&>(*this).mv(x,g);
    checkNaN(g);
    for (int i=0; i<g.N(); ++i)
      assert(!std::isnan(g[i]));
    g += c;
    
    // extract inactive components
    std::vector<int> inactive;
    for (int i=0; i<x.N(); ++i)
      if (!active(x[i][0],b[i][0],Real(10))) 
        inactive.push_back(i);          
                                        
    // check if there is anything to do at all
    if (inactive.empty())
      return;
      
    // solve min 1/2 (x+dx)^T A (x+dx) + c^T (x+dx)
    // first project the gradient
    Vector gp(inactive.size());
    for (int i=0; i<inactive.size(); ++i)
      gp[i] = g[inactive[i]];
        
    // solve for the Newton step F' dx = -F with F = Ax+c = g and F' = A
    Vector xp = static_cast<Implementation const&>(*this).solveSubmatrix(inactive,gp);
    checkNaN(xp);
    dx = 0;
    for (int i=0; i<inactive.size(); ++i)
      dx[inactive[i]] = -xp[i];
  }
  
  template <class Real, class Implementation>
  bool QPBoundSolverBase<Real,Implementation>::projectedLineSearch(Vector& x, Vector const& b, Vector& g, Vector& s,
    Vector const* c
  ) const
  {
    Vector As(b.N());
    static_cast<Implementation const&>(*this).mv(s,As);
    bool activeSetChanged = false;

    // Check that the provided search direction is actually a descent direction.
    // For positive semi-definite QPs, rounding errors can easily lead to
    // computing completely wrong Newton directions.
    if (g*s > 0)
      throw NonpositiveMatrixException("ascent search direction computed",__FILE__,__LINE__);

    while (true)
    {
      // nonvanishing search direction: compute step length t = -g^Ts / s^TAs (unconstrained)
      Real sAs = s*As;
      if (sAs < 0)
      {
        std::cerr << "ev: " << sAs / s.two_norm2() << "\n";
        sAs = 0;
      }
      
      Real t = sAs==0? std::numeric_limits<Real>::infinity()  // if no curvature - follow gradient to infinity
                     : -(g*s) / sAs;                          // otherwise exact linesearch

      // restrict step length to remain in admissible set
      int constraintEncountered = -1;
      for (int i=0; i<s.N(); ++i)
        if (s[i] > 0 && x[i]+t*s[i]>=b[i])  // If we hit a constraint 
        {                                   // this becomes active (we limit t to reach the constraint)
          t = (b[i]-x[i])/s[i];             // and since it has not been active before (s[i]>0),
          activeSetChanged = true;          // the active set changed - explore the new one with Newton
          constraintEncountered = i;        // note that we didn't find an unbound minimizer on this leg - one more iteration
        }
        
      if (std::isinf(t))                    // infinite step size means QP is undbounded 
      {
        std::cerr << "s = [" << s << "]\nx = [\n" << x << "]\nb = [\n" << b << "]\n";
        printMatrix(std::cerr);
        throw InfeasibleProblemException("QP unbounded.",__FILE__,__LINE__);
      } 
        
      // apply step x = x + ts
      x.axpy(t,s);


      if (constraintEncountered<0)          // linsearch found an unbound minimizer on a leg
        return activeSetChanged;            // and that's the end of it.
      else
      {
        // We encountered a constraint and continue the line search, projecting the search direction 
        // on the newly active constraint plane. For continuing the line search, we also need to 
        // update the gradient g(x) = Ax+c to g(x+ts) = g(x) + tAs.
        g.axpy(t,As);
        updateWithColumn(-s[constraintEncountered],constraintEncountered,As);
        s[constraintEncountered] = 0;
        // TODO: a complete matrix-vector multiplication is wasteful here, since s changes only in one 
        //       component and As can be updated by adding a multiple of a single column of A.
//         static_cast<Implementation const&>(*this).mv(s,As);
      }
    }
  }

  template <class Real, class Implementation>
  void QPBoundSolverBase<Real,Implementation>::updateWithColumn(Real a, int i, Vector& s) const
  {
    Vector aei(s.N());
    aei = 0;
    aei[i] = a;
    Vector aAei(s.N());
    static_cast<Implementation const&>(*this).mv(aei,aAei);
    s += aAei;
  }
  
  template <class Real, class Implementation>
  void QPBoundSolverBase<Real,Implementation>::checkNaN(Vector const& v) const
  {
    for (int i=0; i<v.N(); ++i)
      assert(!std::isnan(v[i]));
  }
  
//   template class QPBoundSolverBase<float, QPBoundSolver<float>>;
//   template class QPBoundSolverBase<float, QPDirectSparse<float>>;
  
  template class QPBoundSolverBase<double, QPBoundSolver<double>>;
  template class QPBoundSolverBase<double, QPDirectSparse<double>>;
  
  // ----------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------

  template <class Real>
  QPBoundSolver<Real>::QPBoundSolver(Matrix const& A_, QPConvexificationStrategy convex)
  : QPBoundSolverBase<Real,QPBoundSolver<Real>>(A_.N()), A(A_)
  {
  }
 
  template <class Real>
  typename QPBoundSolver<Real>::Vector QPBoundSolver<Real>::solveSubmatrix(std::vector<int> const& idx, Vector const& b) const
  {
    Matrix Ap = A(idx,idx);
    Real nap = Ap.frobenius_norm();
    for (int i=0; i<Ap.N(); ++i)
      Ap[i][i] += 1e-8*nap;   // regularization in case Ap is singular (semidefinite)
    return gesv(Ap,b);
  }

  template <class Real>
  void QPBoundSolver<Real>::updateWithColumn(Real a, int i, Vector& s) const
  {
    for (int j=0; j<s.N(); ++j)
      s[j] += a*A[j][i];
  }


  template class QPBoundSolver<float >;
  template class QPBoundSolver<double>;
  

  // ----------------------------------------------------------------------------------------------------------
  
  template <class Real>
  QPDirectSparse<Real>::QPDirectSparse(SparseMatrix const& A_, QPConvexificationStrategy convex)
  : QPBoundSolverBase<Real,QPDirectSparse<Real>>(A_.N()), A(A_)
  {
  }
  
  template <class Real>
  Real QPDirectSparse<Real>::frobenius_norm() const 
  {
    Real sum=0;
    // loop over all (non-zero) entires
    for (auto ri=A.begin(); ri!=A.end(); ++ri)
      for (auto ci=ri->begin(); ci!=ri->end(); ++ci)
        sum += (*ci)*(*ci);

    return std::sqrt(sum); 
  }
  
  template <class Real>
  typename QPDirectSparse<Real>::Vector QPDirectSparse<Real>::solveSubmatrix(std::vector<int> const& idx, Vector const& b) const
  {
    // get subindexed NumaBCRS-matrix A(idx, idx)
    SparseMatrix Ap = submatrix<SparseMatrix> (A, idx, idx);
    
    MatrixAsTriplet<Real> Ap_trip(Ap);                                  // setup for direct solve
    MatrixRepresentedOperator<MatrixAsTriplet<Real>, Vector, Vector> B(Ap_trip);

    Vector sol(idx.size()); sol = 0;
    directInverseOperator(B, DirectType::UMFPACK, MatrixProperties::SYMMETRIC).applyscaleadd(-1.0, b, sol);
    return sol;
  }
  
  template class QPDirectSparse<float >;
  template class QPDirectSparse<double>;

  
  // ----------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------

  template <QPStructure sparsity, class Real>
  QPSlackSolver<sparsity,Real>::QPSlackSolver(Matrix const& A_, Matrix const& B_, Real gamma_, QPConvexificationStrategy convex)
  : QPBoundSolverBase<Real,QPSlackSolver<sparsity,Real>>(A_.N()+B_.N()), A(A_), B(B_), gamma(gamma_)
  {
    assert(gamma>0);
    assert(A.M()==B.M());
    

#ifndef NDEBUG
    if constexpr (sparsity==QPStructure::DENSE)
    {
      Matrix I = dynamicUnitMatrix(gamma*unitMatrix<Real,1>(),B.N());
      Matrix Q = vertcat(horzcat(A+gamma*transpose(B)*B,transpose(B)),
                         horzcat(B,I));
      auto e = eigenvalues(Q);
      double lmin = *std::min_element(begin(e),end(e));
      if (lmin <= 0)
      {
        std::cout.precision(17);
        std::cout << "nonpositive QP!\nmax EV: " << *std::max_element(begin(e),end(e)) << ", min EV: "
                  << lmin << "\nQ = [\n" << Q << "];\n";
        std::cout.flush();

        if (false)
        {
          std::ofstream out("qp.m");
          out << "A=[\n" << std::setprecision(16) << A << "\n];\nB=[\n" << B << "\n];\ngamma=" << gamma
              << ";\nQ = [A B'; B gamma*eye(size(B,1))];\nQ1 = [" << Q << "];\n"
              << "lambdamin=" << lmin << ";\nlambda=[";
          for(auto l: e)
            out << l << " ";
          out << "];\n";
          out.close();
        }
      }
    }
#endif
 }
  
  template <QPStructure sparsity, class Real>
  void QPSlackSolver<sparsity,Real>::mv(Vector const& x, Vector& Ax) const
  { 
    int const n = A.N();
    int const m = B.N();
    assert(x.N() == n+m);
    
    // Compute block matrix vector product
    // [A+g*B'*B -g*B'][u]
    // [-g*B      g*I ][s]
    Vector u(n), s(m), ru(n), rs(m);
    for (int i=0; i<n; ++i)
      u[i] = x[i];
    for (int i=0; i<m; ++i)
      s[i] = x[i+n];
    
    B.mv(u,rs);
    rs *= -gamma;         // rs = -gamma*B*u
    B.mtv(rs,ru);
    ru *= -1;             // ru = gamma*B'*B*u
    A.umv(u,ru);          // ru = (A+gamma*B'*B)*u
    B.usmtv(-gamma,s,ru); // ru = (A+gamma*B'*B)*u - gamma*B'*s
    s *= gamma;
    rs += s;              // rs = -gamma*B*u + gamma*s


    for (int i=0; i<n; ++i)
      Ax[i] = ru[i];
    for (int i=0; i<m; ++i)
      Ax[i+n] = rs[i];
  }

//   TO BE IMPLEMENTED
//   template <QPStructure sparsity, class Real>
//   void QPSlackSolver<sparsity,Real>::updateWithColumn(Real a, int i, Vector& x) const
//   {
//     int n const = A.N();
//
//     if (i<n)      // first block column [A+gamma*B'*B; -gamma*B]
//     {
//     }
//     else          // second block column [-gamma*B'; gamma*I]
//     {
//     }
//   }

  template <QPStructure sparsity, class Real>
  Real QPSlackSolver<sparsity,Real>::frobenius_norm() const
  {
    Real normA2 = A.frobenius_norm2();
    Real normB2 = B.frobenius_norm2();

    return std::sqrt(normA2+gamma*gamma*normB2*normB2 + 2*gamma*gamma*normB2 + B.N()*gamma*gamma);
  }
  
  template <QPStructure sparsity, class Real>
  typename QPSlackSolver<sparsity,Real>::Vector
  QPSlackSolver<sparsity,Real>::solveSubmatrix(std::vector<int> const& idx, Vector const& r) const
  {
    assert(idx.size()==r.N());
      
    // The indices to be included are with respect to the block system Q=[A+g*B'*B -g*B'; -g*B g*I].
    // We know (by construction bu = +inf) that u components are inactive and hence included. Thus,
    // we only need to decide which s components are active or inactive.
    int const n = A.N();
    std::vector<int> isActive(B.N(),true), inactive;
    for (int i: idx)
      if (i>=n)
      {
        inactive.push_back(i-n);
        isActive[i-n] = false;
      }

    std::vector<int> idxU(n), active;
    std::iota(begin(idxU),end(idxU),0);
    for (int i=0; i<B.N(); ++i)
      if (isActive[i])
        active.push_back(i);

    int const m = inactive.size();
    assert(idx.size()==n+m);
    
    // We form the Schur complement wrt the lower right block gamma*I.
    // This yields A + B'*G*B with G_ii = gamma for active and 0 for inactive (included).
    // Hence we need to extract the active rows of B. In other terms, only active constraints
    // are penalized.
    Matrix Ba, Bi, S;
    if constexpr (sparsity==QPStructure::DENSE)
    {
      Ba = B(active,idxU);
      Bi = B(inactive,idxU);
      S = A + gamma*(transpose(Ba)*Ba);
    }
    else
    {
      Ba = submatrix<Matrix>(B,active,idxU);
      Bi = submatrix<Matrix>(B,inactive,idxU);
      Matrix Ia = sparseUnitMatrix<Real,1>(Ba.N());
      Ia *= gamma;
      S = A + conjugation(Ba,Ia);
    }

    Real nap = S.frobenius_norm();
    for (int i=0; i<n; ++i)
      S[i][i] += 1e-14*nap;                                // regularization in case S is singular (semidefinite)
      

    // Split the right hand side into u and s parts.
    Vector rU(n), rS(m);
    for (int i=0; i<n; ++i)
      rU[i] = r[i];
    for (int i=0; i<m; ++i)
      rS[i] = r[i+n];
    
    Bi.usmtv(1.0,rS,rU);                 // ru+B'*rs
    Vector u(n);
    if constexpr (sparsity==QPStructure::DENSE)
    {
      u = gesv(S,rU);
    }
    else
    {
      u = 0;
      MatrixAsTriplet<Real> S_trip(S);                                  // setup for direct solve
      MatrixRepresentedOperator<MatrixAsTriplet<Real>, Vector, Vector> Sp_(S_trip);
      auto Sinv = directInverseOperator(Sp_, DirectType::UMFPACK, MatrixProperties::SYMMETRIC);
      Sinv.applyscaleadd(1.0, rU, u);

      S.usmv(-1,u,rU); // residual
      Sinv.applyscaleadd(1,rU,u); // Nachiteration
    }
    Bi.usmv(gamma,u,rS);                  // rs+gamma*B*u
    rS /= gamma;

    
    // combine the result into a single vector
    Vector x(r.N()); x = 0;
    for (int i=0; i<n; ++i)
      x[i] = u[i];
    for (int i=0; i<m; ++i)
      x[i+n] = rS[i];

    return x;
  }

  template class QPSlackSolver<QPStructure::DENSE,double>;
  template class QPSlackSolver<QPStructure::SPARSE,double>;


  // ----------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------

  
  template <int d, class Real>
  QPSolver<d,Real>::QPSolver(MatrixA const& A, MatrixB const& B_, QPConvexificationStrategy convex)
  : Ainv(A), B(B_)
  {
    ScopedTimingSection tsec("QPSolver Constructor");

    if (convex == QPConvexificationStrategy::DONOTHING)
    {
      invertSpd(Ainv);
    }
    else if (convex == QPConvexificationStrategy::INCREMENTALADDITION)
    {
      Real Anorm = A.frobenius_norm();
      Real addDiag = std::numeric_limits<Real>::epsilon() * Anorm / (A.N()*A.N());
      bool positive = false;
      do
      {
        try 
        {
          invertSpd(Ainv);
          positive = true;
        }
        catch (NonpositiveMatrixException const& ex)
        {
          addDiag = 2*addDiag;
          Ainv = A;
          for (int i=0; i<Ainv.N(); ++i)
            Ainv[i][i] += addDiag*unitMatrix<Real,d>();
        }
      } while(!positive);
    }
    else
    {
      auto ev = eigenvalues(Ainv);
      Real addDiag = - std::min(1e-14 * *std::max_element(begin(ev),end(ev))/1e14,
                               1.1   * *std::min_element(begin(ev),end(ev)));
      for (int i=0; i<Ainv.N(); ++i)
           Ainv[i][i] += addDiag*unitMatrix<Real,d>();
      invertSpd(Ainv);      
    }

    auto M = B*Ainv*transpose(B);

    // If B.N() > B.M(), M is singular and may end up being indefinite due to rounding errors.
    // This can lead to the dual QP to be unbounded, such that the iteration does not terminate.
    // We prevent this here.
    if (B.N() > B.M())
    {
      std::cout << "overconstrained QP\n";
      auto e = eigenvalues(M);
      Real lmin = std::min(Real(0),*std::min_element(begin(e),end(e)));
      Real lmax = *std::max_element(begin(e),end(e));
      for (int i=0; i<M.N(); ++i)
        M[i][i] += -2*lmin + 1e-9*lmax;
    }

    // Create the negative Lagrangian dual problem (partially, specifying only the Hessian at this point)
    dual = std::make_unique<QPBoundSolver<Real>>(M);
  }

  template <int d, class Real>
  QPSolver<d,Real>::QPSolver(Self const& qp)
  : Ainv(qp.Ainv)
  , B(qp.B)
  , dual(std::make_unique<QPBoundSolver<Real>>(*qp.dual))
  {}

  template <int d, class Real>
  QPSolver<d,Real>& QPSolver<d,Real>::operator=(Self const& qp)
  {
    Ainv = qp.Ainv;
    B = qp.B;
    dual = std::make_unique<QPBoundSolver<Real>>(*qp.dual);
    return *this;
  }


  template <int d, class Real>
  std::tuple<typename QPSolver<d,Real>::VectorX,typename QPSolver<d,Real>::VectorB,int> 
  QPSolver<d,Real>::solve(VectorX const& c, VectorB const& b, double tol) const
  {
    // STEP I: build dual QP 
    // m = B*Ainv*c + b;
    VectorB m = b;
    VectorX tmpX = c;
    Ainv.mv(c,tmpX);
    B.umv(tmpX,m);

    // STEP II: solve dual QP max -0.5*y'*M*y + m'*y   s.t. y<=0
    //          equivalent to min 0.5*y'*M*y - m'*y s.t. y<=0
    m *= -1;
    VectorB zero = m; zero = 0;
    auto [y,iter] = dual->solve(m,zero,tol);



    // STEP III: reconstruct primal solution
    // x = Ainv*(B'*l-c);
    B.mtv(y,tmpX);
    tmpX -= c;
    VectorX x = c;
    Ainv.mv(tmpX,x);


    // Check that solution is ok. Feasibility: B*x-b <= 0.
    // But remember that solving the dual problem *approximately* can lead to infeasible primal solutions.
    // Thus, it's ok if the primal feasibility is violated within the tolerance (scaled by
    // |BA^{-1}B^T|, which is overestimated here).
    B.mv(x,m);
    m -= b;
    Real consMin = std::numeric_limits<Real>::max();
    Real consMax = -consMin;
    for (int i=0; i<m.N(); ++i)
    {
     consMin = std::min(consMin,m[i][0]);
     consMax = std::max(consMax,m[i][0]);
    }
std::cout << "QPSolver: max constraint violation = " << consMax << "\n";
    if (consMax > pow(B.frobenius_norm(),2)*Ainv.frobenius_norm()* tol*std::abs(consMin))
    {
     std::cerr << "Solution appears to be infeasible: constraints:\n" << m << "\n";
    //       throw InfeasibleProblemException("Problem appears to be infeasible.",__FILE__,__LINE__);
    }

    if (x.two_norm()>1e8)
    {
     std::cerr << "Ainv = [\n" << Ainv << "]\nB = [\n" << B << "]\n";
     std::cerr << "y = [\n" << y <<  "]\nBx-b = [\n" << m << "]\n";
     std::cerr << "x = [\n" << x << "];c = [\n" << c << "]\nb = [\n" << b << "]\n";
     throw UnboundedProblemException("QP appears to be unbounded.",__FILE__,__LINE__);
    }

    return std::make_tuple(x,y,iter);
  }

  template class QPSolver<1,double>;
  template class QPSolver<2,double>;
  template class QPSolver<3,double>;

  // ----------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------




}

#ifdef UNITTEST

#include <iostream>

int main(void)
{
  using namespace Kaskade;
  
  // work with dense matrices
  using QPDense = QPSolver<1,double>;
  using QPSlack = QPSlackSolver<QPStructure::DENSE,double>;
 
  int N = 2;
  int M = 3;
  QPDense::MatrixA A(N,N);
  QPDense::MatrixB B(M,N);
  QPDense::VectorX c(N);
  QPDense::VectorX b(M);
  
  A[0][0]= 1; A[0][1]=-1; A[1][0]=-1; A[1][1]= 2; 
  B[0][0]= 1; B[0][1]= 1; B[1][0]=-1; B[1][1]= 2; B[2][0]= 2; B[2][1]= 1;
  c[0]=-2; c[1]=-6;
  b[0]= 2; b[1]= 2; b[2]= 3;                      
  
  QPDense::VectorX uKnown(N);         // known solution: u = [0.6667 1.3333]
  uKnown[0] = 0.6667; uKnown[1] = 1.3333;
    
  
  // setup slack formulation
  double gamma = 100000;
  auto Bt = transpose(B);
  auto Aslack = A + gamma*Bt*B;
  QPSlack qpSlack(Aslack,-1.0*gamma*B,gamma);
  
  // rhs and bound for slack formulation  
  auto bS = qpSlack.bSlack(b);        
  auto cS = qpSlack.cSlack(c);
  
  // solve slack formulation
  auto sol = qpSlack.solve(cS, bS, 1e-5);
  auto u   = std::get<0>(sol);
  
  
  if(1)
  {
    double err = 0;
    for(int i=0; i<N; ++i)
    {
      err += fabs(u[i] - uKnown[i]);
    }
    std::cout << "|uKnown - u|_1 = " << err << "\n"; 
      
    for(int j=N; j<N+M; ++j)
      std::cout << "slack(" << j-N << ") = " << u[j] << "\n";
  }
  
  return 0;
};

#endif
