/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2013-2019 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef DYNAMICMATRIXOPS_HH
#define DYNAMICMATRIXOPS_HH

#include "linalg/dynamicMatrix.hh"

/**
 * \file contains linear algebra operations on dense matrices and vectors
 */

namespace Kaskade 
{
  /**
   * \ingroup linalgbasic
   * \relates DynamicMatrix
   * \brief matrix-vector product
   */
  template <class MEntry, class VEntry>
  auto operator*(DynamicMatrix<MEntry> const& A, Dune::DynamicVector<VEntry> const& x)
  {
    using YEntry = decltype(A[0][0]*x[0]);
    Dune::DynamicVector<YEntry> y(A.N());
    A.mv(x,y);
    return y;
  }
  
  /**
   * \ingroup linalgbasic
   * \relates DynamicMatrix
   * \brief matrix-vector product
   */
  template <class MEntry, class VEntry>
  auto operator*(DynamicMatrix<MEntry> const& A, Dune::BlockVector<VEntry> const& x)
  {
    using YEntry = decltype(A[0][0]*x[0]);
    Dune::BlockVector<YEntry> y(A.N());
    A.mv(x,y);
    return y;
  }
  
  /**
   * \ingroup linalgbasic
   * \relates DynamicMatrix
   * \brief multiplication with scalar
   */
  template <class MEntry>
  auto operator*(typename Dune::FieldTraits<DynamicMatrix<MEntry>>::field_type x, DynamicMatrix<MEntry> A)
  {
    // Dune 2.5.1 appears to have a bug in densematrix.hh, where DenseMatrix<M>::field_type is defined 
    // as DenseMatVecTraits<M>::value_type (which is the entry type, and maybe a matrix). Thus, A *= x
    // fails, as x is interpreted to be the entry type instead of a scalar.
    for (int j=0; j<A.M(); ++j)
      for (int i=0; i<A.N(); ++i)
        A[i][j] *= x;
    return A;
  }
  
  /**
   * \ingroup linalgbasic
   * \relates DynamicMatrix
   * \brief Multiplication with -1.
   */
  template <class MEntry>
  auto operator-(DynamicMatrix<MEntry> const& A)
  {
    return -1.0*A;
  }
  
  /**
   * \ingroup linalgbasic
   * \relates DynamicMatrix
   * \brief matrix subtraction A-B
   */
  template <class MEntry>
  auto operator-(DynamicMatrix<MEntry> A, DynamicMatrix<MEntry> const& B)
  {
    assert(A.N()==B.N() && A.M()==B.M());
//     A -= B; <--- strange, this runs into a loop (2016-09-29).
    for (int i=0; i<A.N(); ++i)
      for (int j=0; j<A.N(); ++j)
        A[i][j] -= B[i][j];
    return A;
  }

  /**
   * \ingroup linalgbasic
   * \relates DynamicMatrix
   * \brief matrix addition A+B
   */
  template <class MEntry>
  auto operator+(DynamicMatrix<MEntry> A, DynamicMatrix<MEntry> const& B)
  {
    assert(A.N()==B.N() && A.M()==B.M());
    A += B;
    return A;
  }
  
  /**
   * \ingroup linalgbasic
   * \relates DynamicMatrix
   * \brief extracts contiguous submatrices, copying the data
   * 
   * The submatrix is defined by half-open row and column ranges.
   */
  template <class MEntry>
  DynamicMatrix<MEntry> submatrix(DynamicMatrix<MEntry> const& A, int rowstart, int rowend, int colstart, int colend)
  {
    assert(rowend>=rowstart && colend>=colstart);
    assert(0<=rowstart && 0<=colstart);
    assert(rowend<=A.N() && colend<=A.M());
    DynamicMatrix<MEntry> S(rowend-rowstart,colend-colstart);
    for (int j=0; j<colend-colstart; ++j)
      for (int i=0; i<rowend-rowstart; ++i)
        S[i][j] = A[i+rowstart][j+colstart];
    return S;
  }
  
  
  /**
   * \ingroup linalgbasic
   * \relates DynamicMatrix
   * \brief Computes the Kronecker product of two matrices.
   * 
   * kron(A,B) with \f$ A\in \mathbf{R}^{n_a\times m_a} \f$ and \f$ B\in \mathbf{R}^{n_b\times m_b} \f$
   * computes a block matrix \f$ C \in\mathbf{R}^{n_a n_b \times m_a m_b}\f$ with \f$ n_b \times m_b \f$ blocks 
   * such that the \f$ (i,j)\f$-block is \f$ A_{ij} B \f$. This is compatible with Matlab's kron function.
   * 
   * The entry types of A and B must be multiplicable.
   */
  template <class EntryA, class EntryB>
  auto kron(DynamicMatrix<EntryA> const& A, DynamicMatrix<EntryB> const& B)
  {
    using EntryC = decltype(A[0][0]*B[0][0]);
    DynamicMatrix<EntryC> C(A.N()*B.N(),A.M()*B.M());
    
    for (int i=0; i<A.N(); ++i)
      for (int j=0; j<A.M(); ++j)
        for (int k=0; k<B.N(); ++k)
          for (int l=0; l<B.M(); ++l)
            C[i*B.N()+k][j*B.M()+l] = A[i][j]*B[k][l];
    return C;
  }
  
  // --------------------------------------------------------------------------------------------------
  // --------------------------------------------------------------------------------------------------

  namespace DynamicMatrixDetail
  {
    // TODO: move this into global horzcat/vertcat functions using constexpr if once we switch to C++17
    template <class EntryA, class EntryB>
    struct cat 
    {
      static auto horz(DynamicMatrix<EntryA> const& A, DynamicMatrix<EntryB> const& B)
      {
        using field_type = typename EntryA::field_type;
        static_assert(std::is_same<field_type,typename EntryB::field_type>::value,"A and B need to have the same underlying field type!");
        
        int const N = A.N()*EntryA::rows;
        assert(N==B.N()*EntryB::rows);
        
        using EntryAB = Dune::FieldMatrix<field_type,1,1>;
        
        DynamicMatrix<EntryAB> C(N,A.M()*EntryA::cols+B.M()*EntryB::cols);
        unflatten(C,horzcat(flatMatrix(A),flatMatrix(B)));  // TODO: a lot of copying around -> improve
        return C;
      }
      
      static auto vert(Dune::DynamicVector<EntryA> const& x, Dune::DynamicVector<EntryB> const& y)
      {
      }
    };
    
    template <class Entry>
    struct cat<Entry,Entry>
    {
      static auto horz(DynamicMatrix<Entry> const& A, DynamicMatrix<Entry> const& B) 
      {
        assert(A.N()==B.N());
        DynamicMatrix<Entry> C(A.N(),A.M()+B.M());
        for (int i=0; i<A.N(); ++i)
        {
          for (int j=0; j<A.M(); ++j)
            C[i][j] = A[i][j];
          for (int j=0; j<B.M(); ++j)
            C[i][j+A.M()] = B[i][j];
        }
        return C;
      }
      
      static auto vert(Dune::DynamicVector<Entry> const& x, Dune::DynamicVector<Entry> const& y)
      {
        Dune::DynamicVector<Entry> z(x.N()+y.N());
        for (int i=0; i<x.N(); ++i)
          z[i] = x[i];
        for (int i=0; i<y.N(); ++i)
          z[i+x.N()] = y[i];
        return z;
      }
    };
  }

  /**
   * \ingroup linalgbasic
   * \relates DynamicMatrix
   * \brief Concatenates two matrices horizontally
   * 
   * A and B must have the same number of scalar rows.
   */
  template <class EntryA, class EntryB>
  auto horzcat(DynamicMatrix<EntryA> const& A, DynamicMatrix<EntryB> const& B)
  {
    return DynamicMatrixDetail::cat<EntryA,EntryB>::horz(A,B);
  }

  /**
   * \ingroup linalgbasic
   * \relates DynamicMatrix 
   * \brief Concatenes two (column) vectors vertically
   */
  template <class EntryX, class EntryY>
  auto vertcat(Dune::DynamicVector<EntryX> const& x, Dune::DynamicVector<EntryY> const& y)
  {
    return DynamicMatrixDetail::cat<EntryX,EntryY>::horz(x,y);
  }

  /**
   * \ingroup linalgbasic
   * \relates DynamicMatrix
   * \brief Concatenates two matrices vertically
   * 
   * A and B must have the same number of columns.
   */
  template <class Entry>
  DynamicMatrix<Entry> vertcat(DynamicMatrix<Entry> const& A, DynamicMatrix<Entry> const& B)
  {
    assert(A.M()==B.M());
    DynamicMatrix<Entry> C(A.N()+B.N(),A.M());
    for (int j=0; j<A.M(); ++j)
    {
      for (int i=0; i<A.N(); ++i)
        C[i][j] = A[i][j];
      for (int i=0; i<B.N(); ++i)
        C[i+A.N()][j] = B[i][j];
    }
    return C;
  }
  
  /**
   * \ingroup linalgbasic
   * \relates DynamicMatrix
   * \brief Returns the diagonal matrix with the given entry replicated on the diagonal.
   * 
   * In order to obtain a unit matrix with m by m unit matrix entries on the diagonal, use 
   * \code
   * dynamicUnitMatrix(unitMatrix<Scalar,m>(),n);
   * \endcode
   */
  template <class Entry>
  auto dynamicUnitMatrix(Entry const& d, int n)
  {
    Entry zero = 0.0;
    DynamicMatrix<Entry> I(n,n,zero);  // make sure it is zero-initialized
    for (int i=0; i<n; ++i)
      I[i][i] = d;
    return I;
  }
}

#endif
