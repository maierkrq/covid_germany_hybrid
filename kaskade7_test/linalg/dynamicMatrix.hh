/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2013-2022 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef DYNAMICMATRIX_HH
#define DYNAMICMATRIX_HH

#include <tuple>
#include <type_traits>
#include <vector>

#include "dune/common/dynvector.hh"
#include "dune/common/densematrix.hh"
#include "dune/common/fmatrix.hh"

#include "fem/fixdune.hh"

#include "utilities/scalar.hh"

namespace Dune
{
  // Forward declaration
  template <class B, class A>
  class BlockVector;
}


/// \internal
namespace Kaskade {
  template<class K> class DynamicMatrix;

  namespace DynamicMatrixDetail {
    
    // We support DynamicMatrix<Dune::FieldMatrix<double,n,m>> as well as DynamicMatrix<double>.
    // Occasionally we need the dimensions of the entry types. For FieldMatrix entries, these
    // are member enums, while for doubles this is implicitly 1x1. For a uniform access, we define 
    // some traits class here.

    // --------------------------------------------------------------------------------------------

    template <class EntryA, class EntryB, class Enable=void> struct ProductTraits;

    template <class Scalar>
    struct ProductTraits<Scalar,Scalar,std::enable_if_t<std::is_arithmetic_v<Scalar>>>
    {
      using type = Scalar;
    };

    template <class EntryA, class EntryB, int n, int m, int k>
    struct ProductTraits<Dune::FieldMatrix<EntryA,n,m>,Dune::FieldMatrix<EntryB,m,k>,void>
    {
      using type = Dune::FieldMatrix<typename ProductTraits<EntryA,EntryB>::type,n,k>;
    };
    
    // --------------------------------------------------------------------------------------------
    

    // The row vector class accessing DynamicMatrix rows. It simply stores a pointer to the
    // first element, the number (of columns in the matrix) and the stride (number of rows,
    // as the matrix stores elements column-major as in BLAS/LAPACK). It does not own or manage
    // the memory.
    // The interface adheres to the Dune::DenseVector requirements.
    template <class K>
    class StridedVector: public Dune::DenseVector<StridedVector<K>>
    {
      using Self = StridedVector<K>;
      typedef Dune::DenseVector<Self> Base;
      using Entry = std::remove_const_t<K>;

    public:
      using size_type = typename Base::size_type;
      using value_type = K;
      using real_type = typename EntryTraits<Entry>::real_type;

      StridedVector(StridedVector const& x): Dune::DenseVector<StridedVector<K>>(), p(x.p), n(x.n), stride(x.stride) {}

      /**
       * \brief Constructor.
       * \param p pointer to the first vector element
       * \param n number of elements in the vector
       * \param stride the offset between vector elements (a contiguously stored one has stride 1)
       */
      StridedVector(K* p_, size_type n_, size_type stride_): p(p_), n(n_), stride(stride_) {  }

      StridedVector& operator=(StridedVector const& v) = default;
      StridedVector& operator=(K const& x) { for (size_type i=0; i<n; ++i) vec_access(i) = x; return *this; }

      // for DUNE 2.4.x
      value_type&       vec_access(size_type i)       { return *(p+i*stride); }
      value_type const& vec_access(size_type i) const { return *(p+i*stride); }
      size_type         vec_size() const { return n; }
      // for DUNE 2.5
      value_type&       _access(size_type i)       { return *(p+i*stride); }
      value_type const& _access(size_type i) const { return *(p+i*stride); }
      size_type         _size() const { return n; }
      size_type         size() const { return n; } // bug in Dune 2.5 DenseVector<>::size() calling Implementation's size()
      
      value_type&       operator[] (size_type i)       { return vec_access(i); }
      value_type const& operator[] (size_type i) const { return vec_access(i); }
      
      // With Dune 2.4.1 and GCC 5.3 there is some weird thing going on that leads to an infinite loop.
      // Maybe due to wrong detection of enable_if in DenseVector::operator+=()? We provide a fix here.
      template <class X>
      Self& operator+=(StridedVector<X> const& x)
      {
        for (size_type i=0; i<n; ++i)
          (*this)[i] += x[i];
        return *this;
      }
      template <class X>
      Self& operator-=(StridedVector<X> const& x)
      {
        for (size_type i=0; i<n; ++i)
          (*this)[i] -= x[i];
        return *this;
      }
      
      // Dune::DenseVector cannot treat FieldMatrices as entries as of Dune 2.4.1
      real_type two_norm2() const 
      {
        real_type r = 0;
        for (size_type i=0; i<n; ++i)
          r += EntryTraits<Entry>::frobenius_norm2((*this)[i]);
        return r;
      }
      
      real_type two_norm() const 
      {
        return std::sqrt(two_norm2());
      }
      
    private:
      K* p;
      size_type n;
      size_type stride;
    };

    // Helper class for extracting a pointer to the first entry (or the entry itself in case of elementary types
    template <class K>
    struct GetAddress
    {
      static K      * from(K      & k) { return &k; }
      static K const* from(K const& k) { return &k; }
    };

    template <class K, int n, int m>
    struct GetAddress<Dune::FieldMatrix<K,n,m>>
    {
      using field_type = typename EntryTraits<K>::field_type;
      
      static field_type*       from(Dune::FieldMatrix<K,n,m>&       k) { return GetAddress<K>::from(k[0][0]); }
      static field_type const* from(Dune::FieldMatrix<K,n,m> const& k) { return GetAddress<K>::from(k[0][0]); }
    };

    template <class K, int n>
    struct GetAddress<Dune::FieldVector<K,n>>
    {
      using field_type = typename EntryTraits<K>::field_type;

      static field_type*       from(Dune::FieldVector<K,n>&       k) { return GetAddress<K>::from(k[0]); }
      static field_type const* from(Dune::FieldVector<K,n> const& k) { return GetAddress<K>::from(k[0]); }
    };

    template <class K, class A>
    struct GetAddress<Dune::BlockVector<K,A>>
    {
      using field_type = typename EntryTraits<K>::field_type;

      static field_type*       from(Dune::BlockVector<K,A>&       k) { return GetAddress<K>::from(k[0]); }
      static field_type const* from(Dune::BlockVector<K,A> const& k) { return GetAddress<K>::from(k[0]); }
    };
    
    template <class K>
    auto getAddress(K& a)
    {
      return GetAddress<std::remove_const_t<K>>::from(a);
    }

    // Helper class for assigning Dune::FieldMatrix of different scalar type (which leads to segmentation faults as of Dune 2.4)
    template <class A, class B>
    struct Copy
    {
      static void apply(A const& a, B& b)
      {
        b = a;
        assert(!std::isnan(b));
      }
    };

    template <class A, class B, int n, int m>
    struct Copy<Dune::FieldMatrix<A,n,m>, Dune::FieldMatrix<B,n,m>>
    {
      static void apply(Dune::FieldMatrix<A,n,m> const& a, Dune::FieldMatrix<B,n,m>& b)
      {
        for (int i=0; i<n; ++i)
          for (int j=0; j<n; ++j)
            Copy<A,B>::apply(a[i][j],b[i][j]);
      }
    };


    // matrix-vector multiplication calling BLAS
    // y <- alpha*A*x + beta*y
    void gemv(bool transpose, int n, int m, double alpha, double const* A, int lda, double const* x, double beta, double* y);
    void gemv(bool transpose, int n, int m, float  alpha, float  const* A, int lda, float  const* x, float  beta, float * y);

    /**
     * \ingroup linalgbasic
     * \brief computes \f$ C \leftarrow \alpha A B + \beta C \f$
     *
     * \f$ A\in\R^{m\times k}, B\in\R^{k\times n}, C\in\R^{m\times n} \f$
     */
    void gemm(bool transposeA, bool transposeB, int m, int n, int k, double alpha, double const* A, int lda, double const* B, int ldb,
              double beta, double* C, int ldc);
    void gemm(bool transposeA, bool transposeB, int m, int n, int k, float alpha, float const* A, int lda, float const* B, int ldb,
              float beta, float* C, int ldc);

    inline double transpose(double x) { return x; }
    inline float transpose(float x) { return x; }
    
    // For accessing flat matrices
    template <class Scalar>
    Scalar& scalarEntry(Scalar& x, int row, int col)
    { 
      return x;
    }
    
    template <class Entry, int n, int m>
    typename EntryTraits<Entry>::field_type scalarEntry(Dune::FieldMatrix<Entry,n,m> const& A, int row, int col)
    {
      using L = EntryTraits<Entry>;
      return scalarEntry(A[row/L::rows][col/L::cols],row%L::rows,col%L::cols);
    }
    
    template <class Entry, int n, int m>
    typename EntryTraits<Entry>::field_type& scalarEntry(Dune::FieldMatrix<Entry,n,m>& A, int row, int col)
    {
      using L = EntryTraits<Entry>;
      return scalarEntry(A[row/L::rows][col/L::cols],row%L::rows,col%L::cols);
    }
    
    template <class Entry>
    typename EntryTraits<Entry>::field_type scalarEntry(DynamicMatrix<Entry> const& A, int row, int col)
    {
      using L = EntryTraits<Entry>;
      return scalarEntry(A[row/L::rows][col/L::cols],row%L::rows,col%L::cols);
    }

    template <class Entry>
    typename EntryTraits<Entry>::field_type& scalarEntry(DynamicMatrix<Entry>& A, int row, int col)
    {
      using L = EntryTraits<Entry>;
      return scalarEntry(A[row/L::rows][col/L::cols],row%L::rows,col%L::cols);
    }
    
    
    
    
    
    template <class Scalar, class enable = typename std::enable_if_t<std::is_floating_point_v<Scalar>,int>>
    DynamicMatrix<Scalar> flatMatrix(DynamicMatrix<Scalar> const& A)
    {
      static_assert(EntryTraits<Scalar>::rows == 1 && EntryTraits<Scalar>::cols==1,"Default specialization only for scalars!");
      return A;
    }

    template <class Entry, int n, int m>
    DynamicMatrix<typename EntryTraits<Entry>::field_type> flatMatrix(DynamicMatrix<Dune::FieldMatrix<Entry,n,m>> const& A)
    {
      using L = EntryTraits<Dune::FieldMatrix<Entry,n,m>>;
      
      DynamicMatrix<typename EntryTraits<Entry>::field_type> B(L::rows*A.N(),L::cols*A.M());
      for (int j=0; j<B.M(); ++j)
        for (int i=0; i<B.N(); ++i)
          B[i][j] = scalarEntry(A,i,j);
      return B;
    }
    
    template <class Entry>
    void unflatten(DynamicMatrix<Entry>& A, DynamicMatrix<typename EntryTraits<Entry>::field_type> const& B)
    {
      using L = EntryTraits<Entry>;
      
      A.resize(B.N()/L::rows,B.M()/L::cols);
      
      for (int j=0; j<B.M(); ++j)
        for (int i=0; i<B.N(); ++i)
          scalarEntry(A,i,j) = B[i][j];
    }
    
  }
}

namespace Dune
{

  // type traits accessed by Dune::Dense matrix base class
  template <class K>
  struct DenseMatVecTraits<Kaskade::DynamicMatrix<K>>
  {
    typedef Kaskade::DynamicMatrix<K> derived_type;
    typedef Kaskade::DynamicMatrixDetail::StridedVector<K> row_type;
    typedef Kaskade::DynamicMatrixDetail::StridedVector<K const> const_row_type;

    typedef row_type       row_reference;
    typedef const_row_type const_row_reference;

    typedef std::vector<K>                      container_type;
    typedef K                                   value_type;
    typedef typename container_type::size_type  size_type;
  };

  template <class K>
  struct DenseMatVecTraits<Kaskade::DynamicMatrixDetail::StridedVector<K>>
  {
    typedef Kaskade::DynamicMatrixDetail::StridedVector<K> derived_type;
    typedef derived_type row_type;

    typedef K                                                                    value_type;
    typedef typename std::vector<typename std::remove_const<K>::type>::size_type size_type;
  };

  template <class K>
  struct FieldTraits<Kaskade::DynamicMatrix<K>>
  {
    typedef typename FieldTraits<K>::field_type field_type;
    typedef typename FieldTraits<K>::real_type  real_type;
  };
}
/// \endinternal

namespace Kaskade {
  /**
   * \ingroup linalgbasic
   * \brief A LAPACK-compatible dense matrix class with shape specified at runtime.
   * 
   * \tparam K the entry type, either a scalar type or Dune::FieldMatrix or Dune::FieldVector
   *
   * The memory layout is BLAS-like column major. For K = Dune::FieldMatrix<double,n,1>,
   * classical LAPACK and BLAS routines can be applied.
   *
   * In contrast to Dune::DynamicMatrix this works with contiguous memory.
   * In contrast to Dune::Matrix this does not reallocate on resize (unless the number of entries grows).
   *
   * Note that the inherited "solve" method does not work for K=Dune::FieldMatrix<??,1,1>, probably due
   * to a Dune bug around common/densematrix.hh:808 (wrong argument type - pivots are matrix entries, not
   * vector entries...).
   */
  template<class K>
  class DynamicMatrix: public Dune::DenseMatrix<DynamicMatrix<K>>
  {
    typedef Dune::DenseMatrix<DynamicMatrix<K>> Base;
    using Self = DynamicMatrix<K>;
    typedef typename Base::row_reference row_reference;
    typedef typename Base::const_row_reference const_row_reference;
    typedef DynamicMatrixDetail::StridedVector<K> row_base;
    typedef DynamicMatrixDetail::StridedVector<K const> const_row_base;
    using ET = EntryTraits<K>;

  public:
    using size_type  = typename Base::size_type;
    using value_type = typename Base::value_type;
    using field_type = typename Dune::FieldTraits<K>::field_type; // buggy in base class (defined as value_type there)
    typedef typename Base::row_type   row_type;
    
    /**
     * \name Constructors
     * @{
     */

    /**
     * \brief Creates a 0 x 0 matrix.
     */
    DynamicMatrix (): rs(0), cs(0), dat(1) {} // allocate at least one entry to prevent segmentation faults

    /**
     * \brief Copy constructor.
     */
    DynamicMatrix(Self const& a) = default;

    template <class OtherK>
    friend class DynamicMatrix;

    /**
     * \brief Copies from a matrix with different entry type.
     * 
     * If the entry types have different size, the "flattend" matrices are copied scalar by scalar.
     */
    template <class OtherK>
    DynamicMatrix(DynamicMatrix<OtherK> const& A)
    : DynamicMatrix(A.N()*EntryTraits<OtherK>::rows/ET::rows,
                    A.M()*EntryTraits<OtherK>::cols/ET::cols)
    {
      constexpr int otherRows = EntryTraits<OtherK>::rows;
      constexpr int otherCols = EntryTraits<OtherK>::cols;
      
      int n = A.N()*otherRows;      // number of scalar rows
      int m = A.M()*otherCols;      // and cols
      
      assert(n % ET::rows == 0);        // make sure the matrices can have the same scalar size
      assert(m % ET::cols == 0);
      
      // Copy scalar by scalar. Even if the entries have compatible row and column sizes,
      // there need not be an assignment operator (e.g., one is a block matrix and the other one not).
      for (int j=0; j<m; ++j)
        for (int i=0; i<n; ++i)
          DynamicMatrixDetail::scalarEntry(*this,i,j) = DynamicMatrixDetail::scalarEntry(A,i,j);
    }

    /**
     * \brief Move constructor.
     */
    DynamicMatrix (Self&& a) = default;

    /**
     * \brief Creates a r x c matrix, initialized with given value.
     */
    DynamicMatrix (size_type r, size_type c, value_type v = value_type()) :
      rs(r), cs(c), dat(r*c,v)
    {}
    
    /**
     * @}
     */

    /**
     * \todo docme
     */
    using Base::operator=;


    /**
     * \brief Copy assignment.
     * 
     * The matrix size is adjusted as needed.
     */
    DynamicMatrix& operator=(DynamicMatrix const& a) = default;

    /**
     * \brief Move assignment.
     */
    DynamicMatrix& operator=(DynamicMatrix&& a) = default;

    /**
     * \brief Assignment from matrix with other scalar type.
     */
    template <class OtherK>
    DynamicMatrix& operator=(DynamicMatrix<OtherK> const& A)
    {
      static_assert(ET::rows==EntryTraits<OtherK>::rows && ET::rows==EntryTraits<OtherK>::rows,
                    "assignment only from compatible entries");
      resize(A.N(),A.M());
      copy(A.dat,dat);
      return *this;
    }

    /**
     * \brief Resizes the matrix to r x c, leaving the entries in an undefined state.
     *
     * On resize, all row objects and iterators are invalidated. This method adheres to the Dune::DynamicMatrix interface.
     */
    void resize (size_type r, size_type c)
    {
      rs = r;
      cs = c;
      dat.resize(r*c);
    }

    /**
     * \brief Resizes the matrix to r x c, leaving the entries in an undefined state.
     *
     * On resize, all row objects and iterators are invalidated. This method adheres to the Dune::Matrix interface.
     */
    void setSize (size_type r, size_type c)
    {
      resize(r,c);
    }

    void fill (K val)
    {
      for (int i=0; i<rs*cs; ++i) 
        dat[i]=val;
    }
    
    /**
     * \brief Submatrix extraction
     * 
     * This extracts a submatrix, as in Matlab A(ridx,cidx).
     * 
     * \tparam RRange an STL random access range type
     * \tparam CRange an STL random access range type
     * \todo change the implementation to work with forward ranges
     */
    template <class RRange, class CRange>
    DynamicMatrix operator()(RRange const& ridx, CRange const& cidx) const 
    {
      DynamicMatrix B(ridx.size(),cidx.size());
      for (int c=0; c<cidx.size(); ++c)
        for (int r=0; r<ridx.size(); ++r)
          B[r][c] = (*this)[ridx[r]][cidx[c]];
      return B;
    }

    /**
     * \name Linear algebra operations
     * @{
     */
    
    /**
     * \brief matrix-vector product
     *
     * This implementation calls BLAS gemv if possible, and falls back to a simple double loop otherwise.
     * The BLAS call is possible if the entries are either scalars or Dune FieldMatrix of size n x 1, with
     * arbitrary n.
     */
    template<class X, class Y>
    void mv (const X& x, Y& y) const
    {
      using namespace DynamicMatrixDetail;

      if (ET::lapackLayout && rs>0)
        // If the scalar entries are arranged in LAPACK memory layout (and there is anything to do), 
        // we just call BLAS.
        gemv(false,rs*ET::rows,cs*ET::cols,1.0,data(),lda(),
             GetAddress<typename X::value_type>::from(x[0]),
             0.0,GetAddress<typename Y::value_type>::from(y[0]));
      else
        // As of Dune 2.4, the DenseMatrix::mv implementation is buggy and prevents the use of
        // non-scalar entry types. We provide a fixed version here. (2016-02) Apparently corrected
        // in later Dune version.
        for (size_type i=0; i<rs; ++i)
        {
          y[i] = 0; // this was initialized by value_type(0), which can be a matrix...
          for (size_type j=0; j<cs; ++j)
            y[i] += (*this)[i][j] * x[j];
        }
    }

    /**
     * \brief transpose matrix vector product
     * 
     * This computes \f$ y \leftarrow A^T x \f$. We cannot rely on the Dune implementation in DenseMatrix, as
     * it does not transpose the individual entries. This leads to errors when the entries are not square.
     */
    template <class X, class Y>
    void mtv(const X& x, Y& y) const
    {
      for (int i=0; i<this->M(); ++i)     // we could use y = 0, but the loop works for std::vector as well
        y[i] = 0;
      umtv(x,y);  
    }
    
    /**
     * \brief update transpose matrix vector product \f$ y \leftarrow y + A^T x \f$.
     * 
     * We cannot rely on the Dune implementation in DenseMatrix, as
     * it does not transpose the individual entries. This leads to errors when the entries are not square.
     */
    template <class X, class Y>
    void umtv(const X& x, Y& y) const
    {
      usmtv(1.0,x,y);
    }
    
    /**
     * \brief update and scale matrix vector product \f$ y \leftarrow y + \alpha A x \f$.
     * 
     * Dune::DenseMatrix::usmv appears to run into a compile time bug, because DenseMatrix::field_type is
     * defined as the same as value_type, which may not be scalar (as of Dune 2.5).
     */
    template <class X, class Y>
    void usmv(typename Dune::FieldTraits<Y>::field_type alpha, X const& x, Y& y) const
    {
      using namespace DynamicMatrixDetail;

      if (ET::lapackLayout && rs>0)
        gemv(false,rs*ET::rows,cs*ET::cols,alpha,data(),lda(),
             GetAddress<typename X::value_type>::from(x[0]),
             1.0,GetAddress<typename Y::value_type>::from(y[0]));
      else
        for (size_type i=0; i<rs; ++i)
        {
          typename Y::value_type tmp(0.0);
          for (size_type j=0; j<cs; ++j)
            tmp += (*this)[i][j] * x[j];
          y[i] += alpha*tmp;
        }
    }

    /**
     * \brief  update and scale transposed matrix vector product \f$ y \leftarrow y + \alpha A^T x \f$.
     *
     * Dune 2.6 has a bug (the matrix entries are not transposed for multiplication), such that we provide
     * a corrected version here.
     */
    template <class X, class Y>
    void usmtv(typename Dune::FieldTraits<Y>::field_type alpha, X const& x, Y& y) const
    {
      using namespace DynamicMatrixDetail;

      if (ET::lapackLayout && rs>0)
        gemv(true,rs*ET::rows,cs*ET::cols,alpha,data(),lda(),
             GetAddress<typename X::value_type>::from(x[0]),
             1.0,GetAddress<typename Y::value_type>::from(y[0]));
      else
        for(size_type i = 0; i<cs; ++i)
          for(size_type j=0; j<rs; ++j)
            y[i] += alpha * (transpose((*this)[j][i]) * x[j]);
    }

    /**
     * \brief Scales the matrix by a scalar: \f$ A \gets sA \f$
     *
     * This should be provided by the base class Dune::DenseMatrix, but as of Dune 2.5.2, the definition of
     * field_type as value_type in the base class prevents this operator to work for matrix-valued entries.
     * Thus we provide a working implementation.
     */
    DynamicMatrix& operator *=(field_type s)
    {
      for (K& a: dat)
        a *= s;
      return *this;
    }

    /**
     * @}
     */
    
    /**
     * \name Low level access
     * @{
     */
    
    /**
     * \brief The leading dimension (for calling BLAS/LAPACK).
     *
     * Since entries are stored gap-less, this is always the number of matrix rows times the
     * number of scalar rows in the entries.
     */
    size_type lda() const 
    { 
      return rs*EntryTraits<K>::rows;
    }

    /**
     * \brief Raw access to data (for calling BLAS/LAPACK).
     * \warning Unless the entry type is scalar or Dune::FieldMatrix<Scalar,n,1>, the order of matrix entries is NOT column-major.
     */
    field_type* data() { return DynamicMatrixDetail::GetAddress<K>::from(dat[0]); }

    /**
     * \brief Raw access to data (for calling BLAS/LAPACK).
     * \warning Unless the entry type is scalar or Dune::FieldMatrix<Scalar,n,1>, the order of matrix entries is NOT column-major.
     */
    field_type const* data() const { return DynamicMatrixDetail::GetAddress<K>::from(dat[0]); }

    // make this thing a matrix in the sense of the Dune interface
    size_type mat_rows() const { return rs; }
    size_type mat_cols() const { return cs; }
    row_reference mat_access(size_type i) { return row_reference(&dat[i],cs,rs); }
    const_row_reference mat_access(size_type i) const { return const_row_reference(&dat[i],cs,rs); }
    
    /**
     * @}
     */

  private:
    size_type rs, cs;
    std::vector<K> dat;

    // The assignment Dune::FieldMatrix<float,3,3> = Dune::FieldMatrix<double,3,3> leads to a segmentation fault
    // possibly after interminate recursion as of Dune 2.4. Thus we provide a fix. (2016-02)
    template <class OtherK>
    void copy(std::vector<OtherK> const& a, std::vector<K>& b)
    {
      for (size_t i=0; i<a.size(); ++i)
        DynamicMatrixDetail::Copy<OtherK,K>::apply(a[i],b[i]);
    }
  };
  
  // --------------------------------------------------------------------------------------------------------
  // --------------------------------------------------------------------------------------------------------

  /**
   * \ingroup io
   * \relates DynamicMatrix
   * \brief pretty output of dense matrices 
   */
  template <class Entry>
  std::ostream& operator<<(std::ostream& out, DynamicMatrix<Entry> const& A)
  {
    using namespace DynamicMatrixDetail;
    using L = EntryTraits<Entry>;
    
    int n = A.N()*L::rows;
    int m = A.M()*L::cols;
    
    for (int i=0; i<n; ++i)
    {
      for (int j=0; j<m; ++j)
        out << scalarEntry(A,i,j) << " ";
      out << '\n';
    }
    return out;
  }
  
  /**
   * \ingroup linalgbasic
   * \relates DynamicMatrix
   * \brief Computes the matrix-matrix product
   * \warning this creates a temporary
   */
  template <class EntryA, class EntryB>
  DynamicMatrix<typename DynamicMatrixDetail::ProductTraits<EntryA,EntryB>::type>
  operator*(DynamicMatrix<EntryA> const& A, DynamicMatrix<EntryB> const& B)
  {
    assert(A.M()==B.N());
    using namespace DynamicMatrixDetail;
    using EntryC = typename ProductTraits<EntryA,EntryB>::type;
    
    // TODO: use BLAS if possible
    DynamicMatrix<EntryC> C(A.N(),B.M());

    if (EntryTraits<EntryA>::lapackLayout && EntryTraits<EntryB>::lapackLayout && EntryTraits<EntryC>::lapackLayout
        && A.N()>0 && A.M()>0 && B.M()>0)
      gemm(false,false,A.N(),B.M(),A.M(),1.0,A.data(),A.lda(),B.data(),B.lda(),0.0,C.data(),C.lda());
    else
      for (int j=0; j<B.M(); ++j)
        for (int i=0; i<A.N(); ++i)
        {
          C[i][j] = 0;
          for (int l=0; l<A.M(); ++l)
            C[i][j] += A[i][l]*B[l][j];
        }

    return C;
  }

  template <class T, typename = std::enable_if_t<std::is_arithmetic<T>::value>>
  T transpose(T x)
  {
    return x;
  }
  
  /**
   * \ingroup linalgbasic
   * \relates DynamicMatrix
   * \brief Computes the transpose of a matrix 
   * \warning this creates a temporary
   */
  template <class Entry>
  DynamicMatrix<typename EntryTraits<Entry>::transpose_type>
  transpose(DynamicMatrix<Entry> const& A)
  {
    DynamicMatrix<typename EntryTraits<Entry>::transpose_type> C(A.M(),A.N());
    for (int j=0; j<C.M(); ++j)
      for (int i=0; i<C.N(); ++i)
        C[i][j] = transpose(A[j][i]);
    return C;
  }

  // ----------------------------------------------------------------------------------------------

  /**
   * \cond internals
   */
  namespace DynamicMatrixDetail
  {
    template <class Scalar>
    std::tuple<DynamicMatrix<Dune::FieldMatrix<Scalar,1,1>>,
             std::vector<typename ScalarTraits<Scalar>::Real>,
             DynamicMatrix<Dune::FieldMatrix<Scalar,1,1>>>
    svd(DynamicMatrix<Scalar>& A);
  }
  /**
   * \endcond
   */

  /**
   * \ingroup linalgbasic
   * \relates DynamicMatrix
   * \brief Computes the singular value decomposition \f$ A = U \Sigma V^T \f$.
   *
   * If \f$ A \f$ is a \f$ n \times m \f$ matrix, then the SVD gives orthonormal \f$ U, V \f$ of size
   * \f$ n\times n \f$ and \f$ m\times m\f$, respectively, and a \f$ n\times m\f$ matrix \f$ \Sigma \f$ with
   * only nonnegative diagonal entries.
   *
   * \return the tuple (U,sigma,V), where sigma is a std::vector of the \f$ \min(n,m) \f$ real-valued
   *         diagonal entries of \f$ \Sigma \f$.
   */
  template <class Entry>
  std::tuple<DynamicMatrix<Dune::FieldMatrix<typename EntryTraits<Entry>::field_type,1,1>>,
             std::vector<typename EntryTraits<Entry>::real_type>,
             DynamicMatrix<Dune::FieldMatrix<typename EntryTraits<Entry>::field_type,1,1>>>
  svd(DynamicMatrix<Entry> const& A)
  {
    DynamicMatrix<typename EntryTraits<Entry>::field_type> Acopy(A);
    return DynamicMatrixDetail::svd(Acopy);
  }

  /**
   * \ingroup linalgbasic
   * \relates DynamicMatrix
   * \brief Computes the spectral norm of the given matrix
   *
   * This computes \f$ \max_{\|v\|_2=1} \|Av\|_2 \f$.
   *
   * \todo Currently, the norm is computed inefficiently via the signular value decomposition, 
   * and takes \f$ O(n^3) \f$ time. Consider implementing a power iteration.
   */
  template <class Entry>
  typename EntryTraits<Entry>::real_type two_norm(DynamicMatrix<Entry> const& A)
  {
    // Condition number is largest  singular value.
    // TODO: there are faster ways of computing the spectral norm (at least approximately and for square matrices).
    auto sigma = std::get<1>(svd(A));
    return *std::max_element(begin(sigma),end(sigma));
  }


  /**
   * \ingroup linalgbasic
   * \relates DynamicMatrix
   * \brief Computes the condition number of the given square matrix
   *
   * For invertible matrices, the condition number is given as \f$ \kappa(A) = \|A\| \, \|A^{-1}\| \f$. For singular
   * matrices, \f$ +\infty \f$ is returned.
   *
   * Currently, the condition is computed via the signular value decomposition, and takes \f$ O(n^3) \f$ time.
   */
  template <class Entry>
  typename EntryTraits<Entry>::real_type condition(DynamicMatrix<Entry> const& A)
  {
    assert(A.N() == A.M());

    // Condition number is ratio of largest and smallest singular value. Take care if a singular value is zero.
    // TODO: there are faster ways of computing the condition number (at least in special cases).
    auto sigma = std::get<1>(svd(A));
    auto minmax = std::minmax_element(begin(sigma),end(sigma));
    if (*minmax.first <= 0)
      return std::numeric_limits<typename EntryTraits<Entry>::real_type>::infinity();
    else
      return *minmax.second / *minmax.first;
  }

  /**
   * \ingroup linalgbasic
   * \relates DynamicMatrix
   * \brief Computes the solution of \f$ Ax = b \f$ for general invertible matrices.
   * 
   * \tparam Vector a vector type, usually Dune::DynamicVector or Dune::BlockVector
   *
   * The field types of MEntry and VEntry must be the same.
   */
  template <class MEntry, class Vector>
  Vector gesv(DynamicMatrix<MEntry> const& A, Vector const& b);

  /**
   * \ingroup linalgbasic
   * \relates DynamicMatrix
   * \brief Computes the solution of \f$ Ax = b \f$ if \f$ A \f$ is symmetric positive definite.
   * 
   * \tparam Vector a vector type, usually Dune::DynamicVector or Dune::BlockVector
   *
   * \param A must be symmetric positive definite
   * 
   * The field types of MEntry and VEntry must be the same.
   */
  template <class MEntry, class Vector>
  Vector posv(DynamicMatrix<MEntry> const& A, Vector const& b);

  /**
   * \ingroup linalgbasic
   * \relates DynamicMatrix
   * \brief Computes the inverse \f$ A^{-1} \f$.
   * \param A the input matrix to be inverted, is overwritten with its inverse
   */
  template <class Scalar>
  void invert(DynamicMatrix<Scalar>& A);

  /**
   * \ingroup linalgbasic
   * \relates DynamicMatrix
   * \brief Computes the inverse \f$ A^{-1} \f$.
   * \tparam Entry the type of matrix entries, either a scalar or a FieldMatrix<scalar,n,n>
   * \param A the quadratic input matrix to be inverted (has to be spd), is overwritten with its inverse
   * \param symmetrize if true, the complete inverse is computed, otherwise only the lower triangular part
   * \return a reference to the (inverted) input matrix
   * 
   * Can throw NonpositiveMatrixException or SingularMatrixException.
   */
  template <class Entry>
  DynamicMatrix<Entry>& invertSpd(DynamicMatrix<Entry>& A, bool symmetrize=true);
  
  
  /**
   * \ingroup linalgbasic
   * \relates DynamicMatrix
   * \brief Computes the eigenvalues of a symmetric matrix \f$ A \f$.
   * \return a STL sequence of the eigenvalues in ascending order
   */
  template <class Entry>
  std::vector<typename EntryTraits<Entry>::real_type> eigenvalues(DynamicMatrix<Entry> A);
  
  
  /**
   * \ingroup linalgbasic
   * \brief Computes the eigenvalues of a symmetric matrix \f$ A \f$.
   * \return a sequence of the eigenvalues in ascending order
   */
  template <class Scalar, int d>
  Dune::FieldVector<Scalar,d> eigenvalues(Dune::FieldMatrix<Scalar,d,d> A);


  /**
   * \ingroup linalgbasic
   * \relates DynamicMatrix
   * \brief Computes a positive semi-definite approximation to \f$ A \f$ by cutting off negative eigenvalues.
   * 
   * \f$ A \f$ has to be symmetric.
   */
  template <class Scalar, int d>
  bool makePositiveSemiDefinite(Dune::FieldMatrix<Scalar,d,d>& A);
  
  
  /**
   * \ingroup linalgbasic
   * \relates DynamicMatrix
   * \brief Computes a positive semi-definite approximation to \f$ A \f$ by cutting off negative eigenvalues.
   * 
   * \f$ A \f$ has to be symmetric.
   */
  template <class Entry>
  bool makePositiveSemiDefinite(DynamicMatrix<Entry>& A,
                                typename ScalarTraits<typename EntryTraits<Entry>::field_type>::Real threshold = 0);
  
  
  /**
   * \brief reshapes a Dune::BlockVector block structure
   * 
   * Reshapes a Dune::BlockVector with entry-type Dune::FieldVector<double,dimIn> 
   * to entry type Dune::FieldVector<double,dimOut>, i.e. the size of the contained blocks 
   * are reshaped. Block sizes must be such that dimIn == 1 or dimOut == 1.
   * 
   * \tparam dimIn reshape block size from dimIn
   * \tparam dimOut to desired block size dimOut
   * 
   * \param b A Dune::BlockVector<> whose entry-type is to be reshaped
   * 
   */
  template<int dimIn, int dimOut>
  Dune::BlockVector<Dune::FieldVector<double,dimOut>> reshapeBlocks(Dune::BlockVector<Dune::FieldVector<double,dimIn>> const& b)
  {
    assert( (dimIn == 1 && b.N()%dimOut == 0) || (dimIn > 1 && dimOut  == 1) );

    Dune::BlockVector<Dune::FieldVector<double,dimOut>> c;
    
    if(dimIn == 1) // --> dimOut > 1
    {
      c.resize(b.N()/dimOut);
      for(auto i=c.begin(); i!=c.end(); ++i)
        for(int k=0; k<dimOut; ++k)
          c[i.index()][k] = b[i.index()*dimOut+k][0];
    } else      // --> dimIn > 1 and dimOut == 1
    {
      c.resize(b.N()*dimIn);
      for(auto i=b.begin(); i!=b.end(); ++i)
        for(int k=0; k<dimIn; ++k)
          c[i.index()*dimIn+k][0] = b[i.index()][k];
    }
    return c;
  }
  
}

#endif
