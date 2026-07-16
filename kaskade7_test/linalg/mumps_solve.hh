/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2022 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef MUMPS_SOLVE_HH
#define MUMPS_SOLVE_HH

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "dmumps_c.h"
#include "smumps_c.h"
#define ICNTL(I) icntl[(I)-1]
#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654

#include "linalg/factorization.hh"
#include "utilities/detailed_exception.hh"

namespace Kaskade
{

/**
 * \ingroup direct
 * \brief Factorization of sparse linear systems with mumps
 *
 */
template <class Scalar,class SparseIndexInt=int, class DIL=int>
class MUMPSFactorization: public Factorization<Scalar,SparseIndexInt>
{
private:
  using Base = Factorization<Scalar,SparseIndexInt>;

  void printMatrixInfo() {
    if (this->options.verbosity>=2)
    {
      const char *s = "UNKNOWN";
      if (this->options.matrixProperty==MatrixProperties::GENERAL) s = "GENERAL";
      else if (this->options.matrixProperty==MatrixProperties::SYMMETRIC) s = "MatrixProperties::SYMMETRIC";
      else if (this->options.matrixProperty==MatrixProperties::SYMMETRICSTRUCTURE) s = "SYMMETRICSTRUCTURE";
      else if (this->options.matrixProperty==MatrixProperties::POSITIVEDEFINITE) s = "MatrixProperties::POSITIVEDEFINITE";
      std::cout << "MUMPS" << " solver, n=" << n << ", nnz=" << N << ", matrix is " << s << std::endl;
    }
  }

public:

  /**
   * \brief Constructor keeping input data in triplet format (aka coordinate format) constant.
   *
   * Construction is factorization!
   * @param n size of the (square) matrix, i.e. the number of rows
   * @param ridx row indices, will be modified to 1-based indexing
   * @param cidx column indices, will be modified to 1-based indexing
   * @param values entry values (may be modified)
   */
  MUMPSFactorization(SparseIndexInt n_,
                     std::vector<SparseIndexInt>& ridx,
                     std::vector<SparseIndexInt>& cidx,
                     std::vector<Scalar>& values,
                     MatrixProperties property = MatrixProperties::GENERAL,
                     int verb = 0)
  : MUMPSFactorization(n_,ridx,cidx,values,typename Base::Options(property, verb))
  {}


  /**
   * \brief Constructor taking input data in triplet format (aka coordinate format).
   *
   * Construction is factorization!
   * @param n size of the (square) matrix, i.e. the number of rows and columns
   * @param ridx row indices
   * @param cidx column indices
   * @param values entry values
   */
  MUMPSFactorization(SparseIndexInt n_,
                     std::vector<SparseIndexInt>& ridx,
                     std::vector<SparseIndexInt>& cidx,
                     std::vector<Scalar>& values,
                     typename Base::Options options)
  : Base(options), N(ridx.size()), n(n_)
  {
    assert(cidx.size()==N && values.size()==N);

    printMatrixInfo();

    // MUMPS expects symmetric matrices to be given as only one triangular half-matrix. If both a_ij and a_ji are
    // given, they are summed up (which may not be what we want). Thus, if we tell MUMPS that the matrix is symmetric,
    // we set the upper right half to zero.
    if (this->options.matrixProperty==MatrixProperties::SYMMETRIC
        || this->options.matrixProperty==MatrixProperties::POSITIVEDEFINITE)
      for (SparseIndexInt i=0; i<N; ++i)
        if (cidx[i] > ridx[i])      // found a super-diagonal entry
          values[i] = 0;            // ... that we don't wont to fool MUMPS

    init();
    analyze(&ridx[0], &cidx[0]);
    factorize(&values[0]);
  }

  ~MUMPSFactorization()
  {
    mumps_par.job = JOB_END;
    mumps_par.ICNTL(3) = 0;
    callMumps(&mumps_par);
  }

  void init()
  {
    mumps_par.job = JOB_INIT;
    mumps_par.par = 1;

    // Tell MUMPS the structure of the matrix to be factorized. MUMPS distinguishes between
    // unsymmetric(0), general symmetric(2), and positive definite(1) matrices. For the latter,
    // no pivoting is performed.
    switch (this->options.matrixProperty)
    {
      case MatrixProperties::GENERAL:
        mumps_par.sym = 0;
        break;
      case MatrixProperties::SYMMETRICSTRUCTURE:
        mumps_par.sym = 0;
        break;
      case MatrixProperties::SYMMETRIC:
        mumps_par.sym = 2;
        break;
      case MatrixProperties::POSITIVEDEFINITE:
        mumps_par.sym = 1;
        break;
    }
    mumps_par.comm_fortran = USE_COMM_WORLD;
    callMumps(&mumps_par);
    if (this->options.verbosity < 1)        // If no verbose reporting is asked for
      mumps_par.ICNTL(3) = 0;               // reset the output stream
    mumps_par.ICNTL(14) += 5; // 5% additonal workspace
  }

  void analyze(SparseIndexInt* irn, SparseIndexInt* jcn)
  {
    mumps_par.job = 1;
    mumps_par.n = n;
    mumps_par.nz = N;

    // Convert to 1-based indexing
    for (int k=0; k<N; k++)
    {
      irn[k]++;
      jcn[k]++;
    }
    mumps_par.irn = irn;
    mumps_par.jcn = jcn;
    if(this->options.matrixProperty == MatrixProperties::SYMMETRIC
       && this->options.mumps.inertia)
    {
      mumps_par.ICNTL(24) = 1;
      mumps_par.cntl[2] = this->options.mumps.zeroPivotThreshold;
    }
    callMumps(&mumps_par);
    if (mumps_par.info[0] < 0)
      throw DirectSolverException("MUMPS analyze failed with error " + std::to_string(mumps_par.info[0]),__FILE__,__LINE__);
  }

  void factorize(Scalar* a)
  {
    mumps_par.job = 2;
    mumps_par.a = a;
    callMumps(&mumps_par);
    if (mumps_par.info[0] < 0)
    {
      std::string errmsg;
      switch(mumps_par.info[0])
      {
        case -13:
        {
          size_t entries = mumps_par.info[1]<0? -1000000*mumps_par.info[1]: mumps_par.info[1];
          errmsg = "\nAllocation of " + std::to_string(entries) + " scalar entries (" + std::to_string(entries*sizeof(Scalar)/1024/1024) + "MB) failed.";
          break;
        }
        default:
          break;
      }
      throw DirectSolverException("MUMPS factorization failed with error " + std::to_string(mumps_par.info[0])+errmsg,__FILE__,__LINE__);
    }

    if (this->options.matrixProperty == MatrixProperties::SYMMETRIC
        && this->options.mumps.inertia) {
      this->info_.mumps.negativePivots = mumps_par.infog[11];
      this->info_.mumps.zeroPivots = mumps_par.infog[27];
    }
  }

  /**
   * \brief Solves the system for the given right hand side @arg b.
   *
   * \param x is resized to the number of matrix columns.
   */
  void solve(std::vector<Scalar> const& b, std::vector<Scalar>& x, bool transpose=false) const
  {
    assert(b.size()>=n);
    assert(&x != &b);
    x = b;
    solve(x);
  }

  void solve(Scalar const* b, Scalar* x, bool transpose=false) const
  {
    std::copy(b,b+n,x);
    solve(x);
  }



  /**
   * \brief Solves the system for the given right hand side \arg b, which is
   * overwritten with the solution.
   */
  void solve(std::vector<Scalar>& b) const
  {
    assert(b.size()>=n);
    solve(&b[0]);
  }

  void solve(Scalar* b) const
  {
    mumps_par.job = 3;
    mumps_par.rhs = b;
    callMumps(&mumps_par);
    if (mumps_par.info[0] < 0)
      throw DirectSolverException("MUMPS solve failed with error " + std::to_string(mumps_par.info[0]),__FILE__,__LINE__);
  }

  SparseIndexInt size() const
  {
    return n;
  }



private:
  void callMumps(SMUMPS_STRUC_C* mumpsPar) const
  {
    smumps_c(mumpsPar);
  }

  void callMumps(DMUMPS_STRUC_C* mumpsPar) const
  {
    dmumps_c(mumpsPar);
  }

  size_t N;   // number of nonzero entries
  size_t n;   // matrix dimension
  using MumpsStruct = std::conditional_t<std::is_same<Scalar, float>::value, SMUMPS_STRUC_C, DMUMPS_STRUC_C>;
  mutable MumpsStruct mumps_par;
};

}  // namespace Kaskade
#endif
