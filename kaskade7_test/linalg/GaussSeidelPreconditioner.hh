/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2020-2020 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef GAUSS_SEIDEL_PRECONDITIONER_HH
#define GAUSS_SEIDEL_PRECONDITIONER_HH

#include "dune/istl/preconditioners.hh"
#include "dune/istl/solvercategory.hh"

#include "fem/fixdune.hh"
#include "fem/istlinterface.hh"

#include "linalg/matrixTraits.hh"
#include "linalg/symmetricOperators.hh"

#include "utilities/threading.hh"

namespace Kaskade
{

  
  // ----------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------

  enum class GaussSeidelPreconditionerMode { FORWARD=0x1, BACKWARD=0x2, SYMMETRIC=0x3 };

  /**
   * \ingroup preconditioners
   * \brief A simple single-level symmetric Gauss-Seidel preconditioner.
   *
   * For a square matrix \f$ A = L + D + R \f$, this solves
   * \f$ (wL+D) x = wy \f$ in forward mode, \f$ (D+wR) x = wy \f$ in backward mode, or
   * \f$  (wL+D)(D+wR) x = w^2y \f$ in symmetric mode for some stepsize (or damping factor)
   * \f$ w \f$, and thus realizes a SOR or SSOR method as a generalization of the Gauss-Seidel
   * method.
   *
   * \tparam Matrix
   */
  template <class Matrix,
            class Domain=typename MatrixTraits<Matrix>::NaturalDomain,
            class Range=typename MatrixTraits<Matrix>::NaturalRange>
  class GaussSeidelPreconditioner: public SymmetricPreconditioner<Domain,Range>
  {
    using Base = SymmetricPreconditioner<Domain,Range>;
    using field_type = typename Base::field_type;
    
  public:
    /**
     * \brief Constructor.
     * \param A the matrix to be preconditioned. This is copied and stored internally,
     *          such that the provided matrix need not exist after the constructor call.
     *          The matrix must be square (not necessarily symmetric, though) and must
     *          have invertible diagonal elements.
     * \param w the stepsize or damping factor
     */
    GaussSeidelPreconditioner(Matrix const& A_,
                              field_type w_ = 1.0,
                              GaussSeidelPreconditionerMode mode_=GaussSeidelPreconditionerMode::SYMMETRIC)
    : A(A_)
    , w(w_)
    , mode(mode_)
    {
      static_assert(std::is_same_v<Domain,Range>,
                    "Gauss Seidel works only for square matrices with coinciding domain and range types.");
      assert(MatrixTraits<Matrix>::isSquare(A));
    }

    void apply(Domain& x, Range const& y) override
    {
      Domain z = y;

      // forward sweep if requested: (L+D) z = y
      if ((int)mode & (int)GaussSeidelPreconditionerMode::FORWARD)
        for (int i=0; i<A.N(); ++i)
        {
          auto const& row = A[i];
          auto rhs = y[i];
          auto ci = row.begin();
          for (; ci!=row.end() && ci.index()<i; ++ci)                     // update right hand side entry
            ci->usmv(-1.0,z[ci.index()],rhs);                             // with L*z
          assert(ci!=row.end());                                          // make sure D exists
          ci->solve(z[i],rhs);                                            // solve with D
          z[i] *= w;
        }

      // backward sweep if requested: (D+R) x = z
      if ((int)mode & (int)GaussSeidelPreconditionerMode::BACKWARD)
        for (int i=A.N()-1; i>=0; --i)
        {
          auto const& row = A[i];
          auto rhs = z[i];
          auto ci = row.begin();                                          // move to D by skipping over L
          while (ci.index()<i && ci!=row.end())
            ++ci;
          assert(ci.index()==i);                                          // and make sure D exists
          auto a = *ci;
          for (; ci!=row.end(); ++ci)                                     // update right hand side entry
            ci->usmv(-1.0,x[ci.index()],rhs);                             // with R*x
          a.solve(x[i],rhs);                                              // solve with D
          x[i] *= w;
        }
    }

    field_type applyDp(Domain& x, Range const& y) override
    {
      apply(x,y);
      return x*y;
    }

    bool requiresInitializedInput() const override { return true; }

  private:
    Matrix A;
    field_type w;
    GaussSeidelPreconditionerMode mode;
  };
  


} // namespace Kaskade
#endif
