/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2011 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef FACTORIZATION_HH
#define FACTORIZATION_HH

#include <sstream>
#include <tuple>
#include <vector>

#include <dune/grid/config.h>

#include "utilities/enums.hh"
#include "dune/istl/bcrsmatrix.hh"

#include "linalg/factory.hh"
#include "linalg/triplet.hh"

namespace Kaskade
{

//---------------------------------------------------------------------

/**
 * \ingroup direct
 * \brief Abstract base class for matrix factorizations.
 *
 * \tparam Scalar the underlying field type of the matrix elements.
 * \tparam SparseIndexInt the integral type used for indices
 */
template <class Scalar,class SparseIndexInt=int>
class Factorization 
{
public:
  /**
   * \brief The type of matrix elements (a field type).
   */
  typedef Scalar field_type;

  /**
   * @brief The Options struct allows to specify further options of factorization.
   */
  struct Options
  {
    Options() = default;

    Options(MatrixProperties matrixProperty_)
    : matrixProperty(matrixProperty_) {}

    Options(MatrixProperties matrixProperty_, int verbosity_)
    : matrixProperty(matrixProperty_), verbosity(verbosity_) {}

    /// Some direct solvers (Mumps and Pardiso) have parameter specifying matrix property.
    MatrixProperties matrixProperty = MatrixProperties::GENERAL;
    int verbosity = 0;

    /**
     * @brief The Mumps struct allows to specify options which are exclusively relevant for the MUMPS solver.
     */
    struct Mumps
    {
      bool inertia = false; ///< States whether to determine inertia of given matrix, only considered if matrixProperty is SYMMETRIC.
      double zeroPivotThreshold = 0.0; ///< Threshold that defines which pivots will be considered as zero. See CNTL(3) in MUMPS documentation.
    } mumps;
  };

  /**
   * @brief The Info struct provides output information of direct solvers.
   */
  struct Info
  {
    /**
     * @brief The Mumps struct provides output information specific to MUMPS
     */
    struct Mumps
    {
      int negativePivots = 0; ///< number of negative eigenvalues (if inertia detection for symmetric matrix was enabled)
      int zeroPivots = 0; ///< number of zero eigenvalues (if inertia detection for symmetric matrix was enabled)
    } mumps;
  };

  /**
   * \name System solution
   * @{
   */

  /**
   * \brief Solves the system \f$ Ax=b \f$ for the given right hand side \f$ b \f$.
   * \param[inout] b right hand side, is overwritten with the solution \f$ x \f$.
   */
  virtual void solve(std::vector<field_type>& b) const
  {
    assert(b.size() >= size());
    solve(&b[0]);
  }

  /**
   * \brief Solves the system \f$ Ax=b \f$ for the given right hand side \f$ b \f$.
   * \param[inout] b right hand side, is overwritten with the solution \f$ x \f$.
   *                 This must point to a memory region of length of system dimension.
   */
  virtual void solve(field_type* b) const = 0;

  /**
   * \brief Solves the system \f$ Ax=b \f$ for the given right hand side \f$ b \f$.
   * \param[in] b the right hand side
   * \param[out] x the solution. x is resized if needed.
   */
  virtual void solve(std::vector<Scalar> const& b, std::vector<Scalar>& x, bool transposed=false) const
  {
    assert(b.size() >= size());
    if (x.size() < size())
      x.resize(size());
    solve(&b[0], &x[0]);
  }

  /**
   * \brief Solves the system \f$ Ax=b \f$ for the given right hand side \f$ b \f$.
   * \param[in] b the right hand side
   * \param[out] x the solution. x must point to a memory region of length of system dimension.
   */
  virtual void solve(Scalar const* b, Scalar* x, bool transposed=false) const = 0;

  /// @}
  
  Factorization() = default;
  Factorization(Options options_) : options(options_) {}
  virtual ~Factorization() {}
  void setVerbose(int verbose_) { options.verbosity = verbose_; }
  int getVerbose() const { return options.verbosity; }

  Info const& info() const { return info_; }

  /**
   * \brief reports the dimension of the system
   */
  virtual SparseIndexInt size() const = 0;

protected:
  /**
   * \brief Converts a triplet to a compressed column format using umfpack.
   */
  void tripletToCompressedColumn(SparseIndexInt nRows, SparseIndexInt nCols,
                                 size_t nNonZeros, std::vector<SparseIndexInt> const& ridx,
                                 std::vector<SparseIndexInt> const& cidx, std::vector<Scalar> const& values,
                                 std::vector<SparseIndexInt> &Ap, std::vector<SparseIndexInt>& Ai,
                                 std::vector<Scalar>& Az, std::vector<SparseIndexInt>& map);

  Options options;
  Info info_;
};

//---------------------------------------------------------------------

/**
 * \ingroup direct
 * \brief Abstract base class for factorization creators to be plugged into a factorization factory.
 *
 * Derive from this class if you need to plug a direct linear solver
 * into the factory.
 */
template <class Scalar,class SparseIndexInt>
struct Creator<Factorization<Scalar,SparseIndexInt> > 
{
  typedef std::tuple<MatrixAsTriplet<Scalar,SparseIndexInt> const&,typename Factorization<Scalar,SparseIndexInt>::Options> Argument;

  virtual std::unique_ptr<Factorization<Scalar,SparseIndexInt> > create(Argument a) const = 0;

  virtual ~Creator() {}
  
  static std::string lookupFailureHint(DirectType direct)
  {
    std::ostringstream hint;
    hint << "Solver key " << direct << " not registered.\n"
         << "Has the corresponding object file been linked in?\n"
         << "Did you specify the right index type (int/long) in your matrix?\n";
    return hint.str();
  }
};


//---------------------------------------------------------------------

/**
 * \ingroup direct
 * \brief Creates a factorization of the given triplet matrix.
 */
template <class Scalar,class SparseIndexInt>
std::unique_ptr<Factorization<Scalar,SparseIndexInt>> getFactorization(DirectType directType,
                                                                       MatrixAsTriplet<Scalar,SparseIndexInt> const& A,
                                                                       typename Factorization<Scalar,SparseIndexInt>::Options options)
{
  assert(A.nrows()==A.ncols());
  if (A.nnz()==0)
    throw SingularMatrixException("No nonzero entries.",__FILE__,__LINE__);
  auto args = std::make_tuple(A, options);
  if (directType == DirectType::ANY)
    return Factory<DirectType,Factorization<Scalar,SparseIndexInt> >::createAny(args);
  else
    return Factory<DirectType,Factorization<Scalar,SparseIndexInt> >::create(directType,args);
}

/**
 * \ingroup direct
 * \brief Creates a factorization of the given BCRSmatrix.
 */
template <class Scalar,class SparseIndexInt=int /* typename Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,1,1> >::size_type*/ >  // TODO: what is the best default?
std::unique_ptr<Factorization<Scalar,SparseIndexInt>>
getFactorization(DirectType directType, Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,1,1> > const& A, typename Factorization<Scalar,SparseIndexInt>::Options options)
{
  return getFactorization(directType,MatrixAsTriplet<Scalar,SparseIndexInt>(A),options);
}

/**
 * \ingroup direct
 * \brief Creates a factorization of the given NumaBCRS matrix.
 */
template <class Scalar, int n, int m, class Index, class SparseIndexInt=int>  // TODO: what is the best default?
std::unique_ptr<Factorization<Scalar,SparseIndexInt> >  
getFactorization(DirectType directType, NumaBCRSMatrix<Dune::FieldMatrix<Scalar,n,m>,Index> const& A, typename Factorization<Scalar,SparseIndexInt>::Options options)
{
  return getFactorization(directType,MatrixAsTriplet<Scalar,SparseIndexInt>(A),options);
}

}  // namespace Kaskade

#endif 
