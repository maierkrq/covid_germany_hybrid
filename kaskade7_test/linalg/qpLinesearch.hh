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

#ifndef QPLINESEARCH_HH
#define QPLINESEARCH_HH

#include <memory>
#include <tuple>

#include "dune/istl/bvector.hh"

#include "linalg/dynamicMatrix.hh"
#include "linalg/threadedMatrix.hh"

namespace Kaskade
{
  /**
   * \ingroup qp
   * \brief Performs linesearch for penalized QPs.
   *
   * For a penalized convex QP of the form
   * \f[ \min_{\delta x} J(\delta x) = \frac{1}{2} \delta x^T A\delta x + c^T \delta x
   *     + \frac{\gamma}{2}\left\|B\delta x-b\right\|_+^2, \f]
   * and given a search direction \f$ \delta x \f$, this computes the (unique) step size \f$ t \ge 0 \f$
   * minimizing \f$ J(t\delta x) \f$.
   *
   * \param gamma the penalty factor \f$ \gamma \f$. If `gamma == std::numeric_limits<Real>::infinity()`
   *              holds, strict feasibility is enforced, i.e. \f$ tB\delta x \le b \f$. However,
   *              numerical roundoff errors may prevent to take a reasonable step in that case.
   *
   * This function is equivalent to \ref qpLinesearch(MatrixA const& A, MatrixB const& B, Real gamma,
   * VectorX const& c, VectorB const& b, VectorX const& x, VectorX const& dx) with \f$ x = 0 \f$.
   */
  template <class MatrixA, class MatrixB, class VectorX, class VectorB, class Real>
  Real qpLinesearch(MatrixA const& A, MatrixB const& B, Real gamma,
                    VectorX const& c, VectorB const& b, VectorX const& dx)
  {
    // Capture the trivial and useless cases:
    if (dx.two_norm2() == 0)
      return 1;

    if (c.two_norm2() == 0)
      return 0;


    // Select step size by minimizing J(t) = t^2/2 dx^T Adx + t c^T dx + gamma/2 |t Bdx - b|_+^2
    // over t. This is a piecewise quadratic function in t. Switching points can be induced only by
    // constraints i with (Bdx)_i>0 and b_i>0 (activation for growing t) or (Bdx)_i<0
    // and b_i<0 (inactivation for growing t).

    // J(t) = alpha t^2 + beta t + sum_j [eta_j (t-r_j)]_+^2.
    // Gather all quadratic contributions into alpha and all linear terms into beta. First the objective.
    VectorX Adx(dx.N());
    A.mv(dx,Adx);
    double alpha = dx*Adx / 2;
    double beta = dx*c;

    // TODO: RECONSIDER. EQUALITY MAY BE OK, IF UNIQUENESS AND CONVEXITY IS ENSURED BY ACTIVE CONSTRAINTS...
    //  if (alpha <= 0)
    //  {
      //  std::cout << "A = " << A << "\n";
      //  throw NonpositiveMatrixException("nonconvex direction encountered in QP linesearch: dx^T A dx = "
                                       //  +std::to_string(alpha),__FILE__,__LINE__);
    //  }
    if (alpha <= 0)
    {
      std::cout << "Nonconvex direction encountered in QP linesearch: dx^T A dx = " +std::to_string(alpha) << "\n";
      std::cout << "set steplength to 0.\n";
      return 0;
    }

    // Compute the unconstrained minimizer satisfying 2 alpha t + beta = 0.
    Real t = -beta/alpha/2;

    // If gamma==0, the constraints play no role at all. Then we're done here.
    // (Setting gamma to 0 can be used as a quick hack for improving performance when
    //  one knows for sure that the given direction dx does not interfere with the constraints.)
    if (gamma == (Real)0)
      return t;

    // Now consider constraints and distinguish between those that can lead to switching
    // and those that are always active or always inactive (for t>=0). The ones that remain
    // inactive can be ignored altogether.
    VectorB Bdx(b.N());
    B.mv(dx,Bdx);

    // For strict constraints (gamma = inf), we need not minimize after we hit the first
    // constraint - this simplifies things (and we avoid numerical difficulties) in a
    // special case implementation.
    if (gamma == std::numeric_limits<Real>::infinity())
    {
      // For each constraint, limit the stepsize to be feasible, i.e. t B dx <= b.
      for (int i=0; i<b.N(); ++i)
        if (Bdx[i] > 0)               // only consider constraints that might become (more) active
          t = std::min(t,b[i]/Bdx[i]);

      // If t<0, then either dx is no unconstrained descent direction, or we are already
      // infeasible and would violate a constraint even more when moving on (so becoming feasible
      // requires us to go backwards). This is definitely a bad situation.
      if (t < 0)
        throw InfeasibleProblemException("Infeasible initial starting point for linesearch (step "
                                         + std::to_string(t) + ").",
                                           __FILE__,__LINE__);
      return t;                       // And here, we're done.
    }

    // Now comes the tedious case: We have to look through all the segments
    // and find the minimizer.

    std::vector<std::pair<Real,Real>> rEtas;

    for (int i=0; i<b.N(); ++i)
      if (Bdx[i]>=0 && b[i]<=0)       // this constraint will remain active for all t>=0,
      {                               // thus we may simply include it in the quadratic objective
        alpha += gamma/2 * Bdx[i]*Bdx[i];
        beta -= gamma/2 * 2*Bdx[i]*b[i];
      }
      else if (Bdx[i]*b[i]>0)         // this constraint can switch
      {
        double r = b[i]/Bdx[i];       // intercept of t*dx where this constraint switches
        double eta = Bdx[i];          // change of constraint along dx

        rEtas.push_back(std::make_pair(r,eta));

        if (b[i]<0)              // This is already active at t=0 -> include its coefficients into the quadratic objective
        {                        // We may skip the constant terms, though, as we don't compare objective values
          alpha += gamma/2 * eta*eta;
          beta -= gamma/2 * 2*eta*eta*r;
        }
      }

    // Sort the constraints such that we can inspect one segment after the other
    // Sort it in reverse order, such that we can easily pop the encountered switches off the back
    std::sort(begin(rEtas),end(rEtas),FirstGreater());

    // Now look for a minimizer along t>=0.
    t = std::max(0.0,-beta / (2*alpha));         // minimizer of currently active quadratic objective on [0,oo[

    if (t==0)
    {
      // std::cout << "*** no descent in PQP linesearch -> t = 0 (|dx|=" << dx.two_norm() << ")\n";
      return 0;
    }

    while (!rEtas.empty())
    {
      if (t < rEtas.back().first)                     // The minimizer comes before the first switch and is therefore the
        break;                                        // global minimizer. We're done.

      // We ran into a constraint switch before reaching the minimizer. This means that for the remaining part,
      // the quadratic objective is changed.
      auto [r,Bdx] = rEtas.back();
      rEtas.pop_back();

      if (Bdx > 0)                                    // Constraint becomes active
      {                                               // -> add new quadratic contribution [Bdx*(t-r)]^2 to the objective
        alpha += gamma/2 * Bdx*Bdx;                   //    (as before, we can omit the constant term)
        beta -= gamma/2 * 2*Bdx*Bdx*r;
      }
      else                                            // Constraint becomes inactive
      {                                               // -> subtract quadratic contribution [Bdx*(t-r)]^2 that had been
        alpha -= gamma/2 * Bdx*Bdx;                   //    included in the objective up to here
        beta += gamma/2 * 2*Bdx*Bdx*r;
        if (alpha <= 0)
          throw NonpositiveMatrixException("nonconvex direction encountered in QP linesearch",__FILE__,__LINE__);
      }

      t = std::max(r,-beta / (2*alpha));              // Continue linesearch from the constraint switching point.
    }

    return t;
  }

  /**
   * \ingroup qp
   * \anchor qpLinesearchOffset
   * \brief Performs linesearch for penalized QPs.
   *
   * For a penalized convex QP of the form
   * \f[ \min_{x} J(x) = \frac{1}{2}  x^T A x + c^T  x
   *     + \frac{\gamma}{2}\left\|B x-b\right\|_+^2, \f]
   * and given a nonzero starting point \f$ x \f$ and a search direction \f$ \delta x \f$,
   * this computes the (unique) step size \f$ t \ge 0 \f$
   * minimizing \f$ J(x+t\delta x) \f$.
   */
  template <class MatrixA, class MatrixB, class VectorX, class VectorB, class Real>
  Real qpLinesearch(MatrixA const& A, MatrixB const& B, Real gamma,
                    VectorX const& c, VectorB const& b, VectorX const& x, VectorX const& dx)
  {
    // Shift problem such that it is centered around x.
    VectorX c0 = c;
    A.umv(x,c0);
    VectorB b0 = b;
    B.usmv(-1.0,x,b0);


    Real t = 0;
    try
    {
      t = qpLinesearch(A,B,gamma,c0,b0,dx);
    }
    catch (LinesearchException const& ex)
    {
      VectorB con(b.N());  // Bx-b
      B.mv(x,con);
      con -= b;
      for (int i=0; i<b.N(); ++i)
        con[i] = std::max(0.0,con[i][0]);

      B.usmtv(gamma,con,c0); // Ax+c+gamma*B'*(Bx-b)_+, the derivative at x

      if (c0.two_norm() > 1e3* (c.two_norm()+con.two_norm()) * std::numeric_limits<Real>::epsilon())
      {
        std::cout << "|grad| = " << c0.two_norm() << ", |c| = " << c.two_norm() << ", |Bx-b|_+ = " << con.two_norm() << "\n";
        throw LinesearchException(ex.what() + std::string("Does not appear to be due to roundoff."),__FILE__,__LINE__);
      }
    }

    return t;
  }
}

#endif
