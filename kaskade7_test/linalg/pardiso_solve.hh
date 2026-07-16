/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2018 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef PARDISO_SOLVE_HH
#define PARDISO_SOLVE_HH


#include <vector>
#include <memory>

#include "linalg/factorization.hh"

namespace Kaskade
{

/** \ingroup linalg
 * \brief Factorization with DirectType::PARDISO
 *
 * Link with intel mkl if you want to use it. (Version from 2013 will not work, but from 2017 and 2018 works)
 */
template <class Scalar,class SparseIndexInt=int, class DIL=int>
class PardisoFactorization: public Factorization<Scalar,SparseIndexInt>
{
public:

  /** 
   * Construction is factorization! Input matrix in triplet format.
   * @param n size of the (square) matrix, i.e. the number of rows
   * @param ridx row indices
   * @param cidx column indices
   * @param values entry values
   * @param property matrix property
   */
  PardisoFactorization(SparseIndexInt n,
                       std::vector<SparseIndexInt>const & ridx,
                       std::vector<SparseIndexInt>const & cidx,
                       std::vector<Scalar>const & values,
                       MatrixProperties property = MatrixProperties::GENERAL);

  ~PardisoFactorization();

  PardisoFactorization & operator=(PardisoFactorization && other);

  void swap(PardisoFactorization & other);

  /**
   * Solves the system for the given right hand side @arg b. @arg x is
   * resized to the number of matrix columns.
   */
  virtual void solve(std::vector<Scalar> const& b,
                     std::vector<Scalar>& x, int n, bool transposed=false) const;

  virtual void solve(std::vector<Scalar> const& b,
                     std::vector<Scalar>& x, bool transpose=false) const { solve(b,x,1,false); }

  void solve(std::vector<Scalar>& b) const { solve(b,1,false); }

  void solve(std::vector<Scalar>& b, int n, bool transposed=false) const
	{
	  std::vector<Scalar> x(b.size());
	  solve(b,x,n,transposed);
	  std::swap(b,x);
	}
private:
  int n;
  std::vector<SparseIndexInt> ap, ai;
  std::vector<Scalar> az;

  int maxfct = 1;
  int mnum = 1;
  int mtype = 0;
  mutable int phase = 11;
  int nrhs = 1;
  int msglvl = 0;
  mutable int error = 0;
  mutable int idum = 0;
  mutable Scalar ddum = 0.0;
  mutable std::array<int, 64> iparm {};
  mutable std::array<void *, 64> pt {};
};
}  // namespace Kaskade

#endif
