/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2022-2022 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef ADD_HH
#define ADD_HH

#include <cassert>
#include <vector>

#include "linalg/threadedMatrix.hh"
#include "mg/bddc.hh"
#include "utilities/timing.hh"

/**
 * \file 
 * \brief Algebraic domain decomposition based on the sparsity pattern of the matrix.
 */

namespace Kaskade::BDDC
{
  /**
   * \ingroup multigrid
   * \brief Creates a decomposition of indices into overlapping subsets that can be used for 
   *        constructing nonoverlapping domain decomposition solvers.
   * 
   * \tparam Index the integral number type to be used for indices
   * 
   * \param A the symmetric sparse matrix pattern
   * \param n the number of subdomains to create
   * 
   * \return a sequence of subdomains (index subsets)
   */
  template <class Index>
  std::vector<std::vector<Index>> algebraicDomainDecomposition(NumaCRSPattern<Index> const& A, int n);

  template <class Scalar, class Index>
  std::vector<NumaBCRSMatrix<Dune::FieldMatrix<Scalar,1,1>,Index>>
  algebraicMatrixDecomposition(NumaBCRSMatrix<Dune::FieldMatrix<Scalar,1,1>,Index> const& A,
                               std::vector<std::vector<Index>> const& subIndices,
                               Scalar tolerance, int maxIter);


  /**
   * \ingroup multigrid
   * \brief
   * \param tolerance defines the termination criterion by asking for
   *                  \f[ \max_{k,i} \sum_j \frac{|a_{ij}^k|}{a_{ii}^k} -1  \le \mathrm{tol}. \f]
   * \param maxIter defines a termination criterion by imposing an upper bound on the
   *                iteration count.
   *
   * \return
   */
  template <class Scalar, int m, class Index>
  std::tuple< std::vector<NumaBCRSMatrix<Dune::FieldMatrix<Scalar,m,m>,Index>>,
              std::vector<Dune::BlockVector<Dune::FieldVector<Scalar,m>>>,
              std::vector<std::vector<Index>>,
              std::vector<std::vector<LocalDof>> >
  matrixDomainDecommposition(NumaBCRSMatrix<Dune::FieldMatrix<Scalar,m,m>,Index> const& A,
                             Dune::BlockVector<Dune::FieldVector<Scalar,m>> const& f, int n,
                             Scalar tolerance=1000*std::numeric_limits<Scalar>::epsilon(),
                             int maxIter=100)
  {
    using Matrix = NumaBCRSMatrix<Dune::FieldMatrix<Scalar,m,m>,Index>;
    using Vector = Dune::BlockVector<Dune::FieldVector<Scalar,m>>;

    ScopedTimingSection tsec("matrixDomainDecomposition");

    auto& timer = tsec.timer();
    timer.start("partitioning");
    assert(n >= 1);
    assert(A.N() == A.M());
    auto subIndices = algebraicDomainDecomposition(*A.getPattern(),n);
    timer.stop("partitioning");

    timer.start("rhs decomposition");
    std::vector<std::vector<LocalDof>> localDofs(A.N());
    for (int k=0; k<n; ++k)
      for (int iloc=0; iloc<subIndices[k].size(); ++iloc)
      {
        Index iglob = subIndices[k][iloc];
        localDofs[iglob].push_back({k,iloc});
      }
    timer.stop("rhs decomposition");

    std::vector<Matrix> As;
    if constexpr (m==1)
      As = algebraicMatrixDecomposition(A,subIndices,tolerance,maxIter);
    else
      assert("needs to be implemented"==0);

    std::vector<Vector> Fs(n);
    for (int k=0; k<n; ++k)
    {
      Fs[k].resize(subIndices[k].size());
      for (Index i=0; i<subIndices[k].size(); ++i)
      {
        Index iglob = subIndices[k][i];
        int nsub = localDofs[iglob].size();
        Fs[k][i] = f[iglob] / nsub;
      }
    }

    // Finally remove the interior dofs belonging to only one subdomain - they're of no
    // further interest.
    localDofs.erase(std::remove_if(begin(localDofs),end(localDofs),
                                   [](auto const& lds) { return lds.size() <= 1; }),
                    end(localDofs));

    return std::tuple(As,Fs,subIndices,localDofs);
  }
  
}

#endif
