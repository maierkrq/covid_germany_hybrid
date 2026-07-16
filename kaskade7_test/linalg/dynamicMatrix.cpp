/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2014-2020 Zuse Institute Berlin                             */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <algorithm>
#include <limits>
#include <vector>

#ifdef __cplusplus
extern "C" {
#endif
#include <cblas.h>

#include "lapack/lapacke.h"
// void dgesvd_(char*,char*,int*,int*,double*,int*,double*,double*,int*,
//              double*,int*,double*,int*,int*);
#ifdef __cplusplus
}
#endif

#include "dune/grid/config.h"

#include "linalg/dynamicMatrix.hh"
#include "utilities/detailed_exception.hh"

namespace Kaskade 
{
  
  namespace DynamicMatrixDetail
  {
    int invertLAPACK(int n, double* A, int lda)
    {
      int info;
      std::vector<int> ipiv(n);
      info=LAPACKE_dgetrf(LAPACK_COL_MAJOR,n,n,A,lda,&ipiv[0]);
      assert(info>=0);
      if (info>0)
        throw SingularMatrixException("zero pivot",__FILE__,__LINE__);
      info=LAPACKE_dgetri(LAPACK_COL_MAJOR,n,A,lda,&ipiv[0]);
      assert(info>=0);
      return info;
    }
    
    int invertLAPACK(int n, float* A, int lda)
    {
      int info;
      std::vector<int> ipiv(n);
      info = LAPACKE_sgetrf(LAPACK_COL_MAJOR,n,n,A,lda,&ipiv[0]);
      assert(info>=0);
      if (info>0)
        throw SingularMatrixException("zero pivot",__FILE__,__LINE__);
      info = LAPACKE_sgetri(LAPACK_COL_MAJOR,n,A,lda,&ipiv[0]);
      return info;
    }
    
    void invertSpdLAPACK(int n, float* A, int lda)
    {
      int info = LAPACKE_spotrf(LAPACK_COL_MAJOR,'L',n,A,lda);
      assert(info>=0);
      if (info>0)
        throw NonpositiveMatrixException("nonpositive pivot in Cholesky (LAPACK spotrf)",__FILE__,__LINE__);
      info = LAPACKE_spotri(LAPACK_COL_MAJOR,'L',n,A,lda);
      if (info>0)
        throw SingularMatrixException("zero pivot (LAPACK spotri)",__FILE__,__LINE__);
    }
    
    void invertSpdLAPACK(int n, double* A, int lda)
    {
      int info = LAPACKE_dpotrf(LAPACK_COL_MAJOR,'L',n,A,lda);
      assert(info>=0);
      if (info>0)
        throw NonpositiveMatrixException("nonpositive pivot in Cholesky (LAPACK dpotrf)",__FILE__,__LINE__);
      info = LAPACKE_dpotri(LAPACK_COL_MAJOR,'L',n,A,lda);
      if (info>0)
        throw SingularMatrixException("zero pivot (LAPACK dpotri)",__FILE__,__LINE__);
    }
    
    
    
    void gemv(bool transpose, int n, int m, double alpha, double const* A, int lda, double const* x, double beta, double* y)
    {
      cblas_dgemv(CblasColMajor,transpose? CblasTrans: CblasNoTrans,n,m,alpha,const_cast<double*>(A),lda,
                  const_cast<double*>(x),1,beta,y,1);
    }
    void gemv(bool transpose, int n, int m, float  alpha, float  const* A, int lda, float  const* x, float  beta, float * y)
    {
      cblas_sgemv(CblasColMajor,transpose? CblasTrans: CblasNoTrans,n,m,alpha,const_cast<float*>(A),lda,
                  const_cast<float*>(x),1,beta,y,1); 
    }

    // C <- alpha*A*B + beta*C
    void gemm(bool transposeA, bool transposeB, int m, int n, int k, double alpha, double const* A, int lda, double const* B, int ldb,
              double beta, double* C, int ldc)
    {
      cblas_dgemm(CblasColMajor,transposeA? CblasTrans: CblasNoTrans,transposeB? CblasTrans: CblasNoTrans,
                  m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
    }
    void gemm(bool transposeA, bool transposeB, int m, int n, int k, float alpha, float const* A, int lda, float const* B, int ldb,
              float beta, float* C, int ldc)
    {
      cblas_sgemm(CblasColMajor,transposeA? CblasTrans: CblasNoTrans,transposeB? CblasTrans: CblasNoTrans,
                  m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
    }
  }
  
// ----------------------------------------------------------------------------
  
  namespace
  {
    int gesvd(int n, int m, double* a, int lda, double* sing, double* u, int ldu, double* vt, int ldvt, double* superb)
    {
      int info = LAPACKE_dgesvd(LAPACK_COL_MAJOR,'A','A', n,m,a,lda,sing,u,ldu,vt,ldvt,superb);
      return info;
    }

    int gesvd(int n, int m, float* a, int lda, float* sing, float* u, int ldu, float* vt, int ldvt, float* superb)
    {
      int info = LAPACKE_sgesvd(LAPACK_COL_MAJOR,'A','A', n,m,a,lda,sing,u,ldu,vt,ldvt,superb);
      return info;
    }

    int gesv(int n, double* a, int lda, int* ipiv, double* b)
    {
      int info =  LAPACKE_dgesv(LAPACK_COL_MAJOR,n,1,a,lda,ipiv,b,n);
      assert(info>=0);
      return info;
    }
        
    int gesv(int n, float* a, int lda, int* ipiv, float* b)
    {
      return LAPACKE_sgesv(LAPACK_COL_MAJOR,n,1,a,lda,ipiv,b,n);
    }
    
    int posv(int n, double* a, int lda, double* b)
    {
      return LAPACKE_dposv(LAPACK_COL_MAJOR,'L',n,1,a,lda,b,n);
    }
        
    int posv(int n, float* a, int lda, float* b)
    {
      return LAPACKE_sposv(LAPACK_COL_MAJOR,'L',n,1,a,lda,b,n);
    }
    
    
    template <class Scalar>
    void syev(char job, int n, Scalar* a, typename ScalarTraits<Scalar>::Real* w)
    {
      static_assert(ScalarTraits<Scalar>::isBLAStype,
                    "Eigenvalue computation implemented only for LAPACK scalar types");
      int info = 0;

      if constexpr (std::is_same_v<Scalar,double>)
        info = LAPACKE_dsyev(LAPACK_COL_MAJOR,job,'L',n,a,n,w);
      if constexpr (std::is_same_v<Scalar,float>)
        info = LAPACKE_ssyev(LAPACK_COL_MAJOR,job,'L',n,a,n,w);
      // TODO: implement complex hermitian EV computation
      assert(info >= 0);
      if (info > 0)
        throw LinearAlgebraException("LAPACK dsyev failed to converge",__FILE__,__LINE__);
    }
    
  }
  
  
  template <class MEntry, class Vector>
  Vector gesv(DynamicMatrix<MEntry> const& A_, Vector const& b_)  
  {
    using namespace DynamicMatrixDetail;
    
    auto A = flatMatrix(A_);
    auto b = b_;
    assert(A.rows()==A.cols());
    assert(A.rows()==b.size());
    std::vector<int> ipiv(A.rows());
    auto info = gesv(A.rows(),A.data(),A.lda(),&ipiv[0],getAddress(b[0]));
    if (info==-7)
      throw DirectSolverException("illegal right hand side values (probably nan or inf).",__FILE__,__LINE__);
    assert(info>=0); // info < 0 means illegal argument value
    if (info>0)
      throw SingularMatrixException("zero pivot in LAPACK gesv reported.",__FILE__,__LINE__);
    
    return b;
  }
  
  // explicit instantiation
  template Dune::DynamicVector<double> gesv(DynamicMatrix<double> const& A_,
                                            Dune::DynamicVector<double> const& b_);
  template Dune::DynamicVector<Dune::FieldVector<double,1>> gesv(DynamicMatrix<double> const& A_,
                                                                 Dune::DynamicVector<Dune::FieldVector<double,1>> const& b_);
  template Dune::DynamicVector<Dune::FieldVector<double,1>> gesv(DynamicMatrix<Dune::FieldMatrix<double,1,1>> const& A_,
                                                                 Dune::DynamicVector<Dune::FieldVector<double,1>> const& b_);
  template Dune::BlockVector<Dune::FieldVector<double,1>> gesv(DynamicMatrix<double> const& A_,
                                                               Dune::BlockVector<Dune::FieldVector<double,1>> const& b_);
  template Dune::BlockVector<Dune::FieldVector<double,1>> gesv(DynamicMatrix<Dune::FieldMatrix<double,1,1>> const& A_,
                                                               Dune::BlockVector<Dune::FieldVector<double,1>> const& b_);
  template Dune::DynamicVector<Dune::FieldVector<float,1>> gesv(DynamicMatrix<float> const& A_,
                                                                Dune::DynamicVector<Dune::FieldVector<float,1>> const& b_);
  template Dune::DynamicVector<Dune::FieldVector<float,1>> gesv(DynamicMatrix<Dune::FieldMatrix<float,1,1>> const& A_,
                                                                Dune::DynamicVector<Dune::FieldVector<float,1>> const& b_);
  template Dune::BlockVector<Dune::FieldVector<float,1>> gesv(DynamicMatrix<float> const& A_,
                                                              Dune::BlockVector<Dune::FieldVector<float,1>> const& b_);
  template Dune::BlockVector<Dune::FieldVector<float,1>> gesv(DynamicMatrix<Dune::FieldMatrix<float,1,1>> const& A_,
                                                              Dune::BlockVector<Dune::FieldVector<float,1>> const& b_);
  
  template <class MEntry, class Vector>
  Vector posv(DynamicMatrix<MEntry> const& A_, Vector const& b)
  {
    using namespace DynamicMatrixDetail;
    
    Vector x = b;
    auto A = flatMatrix(A_);

    assert(A.rows()==A.cols());
    int info = posv(A.rows(),A.data(),A.lda(),getAddress(x[0]));
    if (info>0)
      throw NonpositiveMatrixException("nonpositive pivot in LAPACK posv reported.",__FILE__,__LINE__);
    return x;
  }
  
  // explicit instantiation
  template Dune::DynamicVector<Dune::FieldVector<double,1>> posv(DynamicMatrix<double> const& A_,
                                                                 Dune::DynamicVector<Dune::FieldVector<double,1>> const& b_);
  template Dune::DynamicVector<Dune::FieldVector<double,1>> posv(DynamicMatrix<Dune::FieldMatrix<double,1,1>> const& A_,
                                                                 Dune::DynamicVector<Dune::FieldVector<double,1>> const& b_);
  template Dune::BlockVector<Dune::FieldVector<double,1>> posv(DynamicMatrix<double> const& A_,
                                                               Dune::BlockVector<Dune::FieldVector<double,1>> const& b_);
  template Dune::BlockVector<Dune::FieldVector<double,1>> posv(DynamicMatrix<Dune::FieldMatrix<double,1,1>> const& A_,
                                                               Dune::BlockVector<Dune::FieldVector<double,1>> const& b_);
  template Dune::BlockVector<Dune::FieldVector<double,2>> posv(DynamicMatrix<Dune::FieldMatrix<double,2,2>> const& A_,
                                                               Dune::BlockVector<Dune::FieldVector<double,2>> const& b_);
  template Dune::BlockVector<Dune::FieldVector<double,3>> posv(DynamicMatrix<Dune::FieldMatrix<double,3,3>> const& A_,
                                                               Dune::BlockVector<Dune::FieldVector<double,3>> const& b_);
  
  
  template Dune::DynamicVector<Dune::FieldVector<float,1>> posv(DynamicMatrix<float> const& A_,
                                                                Dune::DynamicVector<Dune::FieldVector<float,1>> const& b_);
  template Dune::DynamicVector<Dune::FieldVector<float,1>> posv(DynamicMatrix<Dune::FieldMatrix<float,1,1>> const& A_,
                                                                Dune::DynamicVector<Dune::FieldVector<float,1>> const& b_);
  template Dune::BlockVector<Dune::FieldVector<float,1>> posv(DynamicMatrix<float> const& A_,
                                                              Dune::BlockVector<Dune::FieldVector<float,1>> const& b_);
  template Dune::BlockVector<Dune::FieldVector<float,1>> posv(DynamicMatrix<Dune::FieldMatrix<float,1,1>> const& A_,
                                                              Dune::BlockVector<Dune::FieldVector<float,1>> const& b_);
  template Dune::BlockVector<Dune::FieldVector<float,2>> posv(DynamicMatrix<Dune::FieldMatrix<float,2,2>> const& A_,
                                                               Dune::BlockVector<Dune::FieldVector<float,2>> const& b_);
  template Dune::BlockVector<Dune::FieldVector<float,3>> posv(DynamicMatrix<Dune::FieldMatrix<float,3,3>> const& A_,
                                                               Dune::BlockVector<Dune::FieldVector<float,3>> const& b_);
  

  // ----------------------------------------------------------------------------
  
  namespace DynamicMatrixDetail
  {
    template <class Scalar>
    std::tuple<DynamicMatrix<Dune::FieldMatrix<Scalar,1,1>>,
               std::vector<typename ScalarTraits<Scalar>::Real>,
               DynamicMatrix<Dune::FieldMatrix<Scalar,1,1>>>
    svd(DynamicMatrix<Scalar>& A)
    {
      using namespace DynamicMatrixDetail;
      using Real = typename ScalarTraits<Scalar>::Real;
      DynamicMatrix<Dune::FieldMatrix<Scalar,1,1>> U(A.rows(),A.rows()), VT(A.cols(),A.cols());
      int minMN = std::min(A.rows(),A.cols());
      std::vector<Real> sing(minMN), superb(minMN);

      int info = gesvd(A.rows(), A.cols(), A.data(), A.lda(), &sing[0], U.data(), U.lda(),
                       VT.data(), VT.lda(), &superb[0]);
      assert(info>=0);
      if (info > 0)
        throw LinearAlgebraException("LAPACK SVD failed to converge",__FILE__,__LINE__);

      // LAPACK computes V^T, but we return V. Transpose.
      for (int i=1; i<VT.rows(); ++i)
        for (int j=0; j<i; ++j)
          std::swap(VT[i][j],VT[j][i]);

      return std::make_tuple(U,sing,VT);
    }

    // explicit instantiation
    template std::tuple<DynamicMatrix<Dune::FieldMatrix<double,1,1>>,std::vector<double>,DynamicMatrix<Dune::FieldMatrix<double,1,1>>>
    svd(DynamicMatrix<double>& A);
    template std::tuple<DynamicMatrix<Dune::FieldMatrix<float,1,1>>,std::vector<float>,DynamicMatrix<Dune::FieldMatrix<float,1,1>>>
    svd(DynamicMatrix<float>& A);
  }

// ----------------------------------------------------------------------------

  template <class Scalar>
  void invert(DynamicMatrix<Scalar>& A)
  {
    DynamicMatrixDetail::invertLAPACK(A.N(),A.data(),A.lda());
  }
  // explicit instantiation
  template void invert(DynamicMatrix<double>& A);
  template void invert(DynamicMatrix<float>& A); 
  
  
  template <class Entry>
  DynamicMatrix<Entry>& invertSpd(DynamicMatrix<Entry>& A, bool symmetrize)
  {
    using namespace DynamicMatrixDetail;
    
    // Note that the symmetrization cannot be factorized out below the if/else:
    // then diagonal blocks would not be symmetrized correctly
    if (EntryTraits<Entry>::lapackLayout)
    {
      invertSpdLAPACK(A.N(),A.data(),A.lda());
      // lower part of B is filled - symmetrize
      if (symmetrize)
        for (int i=0; i<A.N()-1; ++i)
          for (int j=i+1; j<A.M(); ++j)
            A[i][j] = A[j][i];
    }
    else
    {
      auto B = flatMatrix(A);
      invertSpdLAPACK(B.N(),B.data(),B.lda());
      // lower part of B is filled - symmetrize
      if (symmetrize)
        for (int i=0; i<B.N()-1; ++i)
          for (int j=i+1; j<B.M(); ++j)
            B[i][j] = B[j][i];
      unflatten(A,B);
    }
        
    return A;
  }
  // explicit instantiation
  template DynamicMatrix<double>& invertSpd(DynamicMatrix<double>& A, bool);
  template DynamicMatrix<float>& invertSpd(DynamicMatrix<float>& A, bool); 
  template DynamicMatrix<Dune::FieldMatrix<double,1,1>>& invertSpd(DynamicMatrix<Dune::FieldMatrix<double,1,1>>& A, bool);
  template DynamicMatrix<Dune::FieldMatrix<float,1,1>>& invertSpd(DynamicMatrix<Dune::FieldMatrix<float,1,1>>& A, bool); 
  template DynamicMatrix<Dune::FieldMatrix<double,2,2>>& invertSpd(DynamicMatrix<Dune::FieldMatrix<double,2,2>>& A, bool);
  template DynamicMatrix<Dune::FieldMatrix<float,2,2>>& invertSpd(DynamicMatrix<Dune::FieldMatrix<float,2,2>>& A, bool); 
  template DynamicMatrix<Dune::FieldMatrix<double,3,3>>& invertSpd(DynamicMatrix<Dune::FieldMatrix<double,3,3>>& A, bool);
  template DynamicMatrix<Dune::FieldMatrix<float,3,3>>& invertSpd(DynamicMatrix<Dune::FieldMatrix<float,3,3>>& A, bool); 
  

  // -----------------------------------------------------------------------------------------------
  // -----------------------------------------------------------------------------------------------

  
  template <class Scalar, int d>
  Dune::FieldVector<Scalar,d> eigenvalues(Dune::FieldMatrix<Scalar,d,d> A)
  {
    Dune::FieldVector<double,d> w;
    
    int info = LAPACKE_dsyev(LAPACK_ROW_MAJOR,'N','L', d, &A[0][0], d, &w[0]);
    // todo: check info
    assert(info>=0);
    if (info > 0)
      throw LinearAlgebraException("LAPACK dsyev failed to converge",__FILE__,__LINE__);
    return w;
  }
  
  template Dune::FieldVector<double,2> eigenvalues(Dune::FieldMatrix<double,2,2> A);
  template Dune::FieldVector<double,3> eigenvalues(Dune::FieldMatrix<double,3,3> A);
  
  // -----------------------------------------------------------------------------------------------

  template <class Entry>
  std::vector<typename EntryTraits<Entry>::real_type> eigenvalues(DynamicMatrix<Entry> A)
  {
    using namespace DynamicMatrixDetail;
    using L = EntryTraits<Entry>;
    assert(A.N()==A.M() && L::rows==L::cols);
    std::vector<typename L::real_type> ev(A.N()*L::rows);
    
    if (L::lapackLayout)
      syev('N',ev.size(),GetAddress<Entry>::from(A[0][0]),&ev[0]);
    else
    {
      DynamicMatrix<typename L::field_type> a = flatMatrix(A);
      syev('N',ev.size(),&a[0][0],&ev[0]);
    }
    
    return ev;
  }
  
  template std::vector<float> eigenvalues(DynamicMatrix<Dune::FieldMatrix<float,1,1>> A);
  template std::vector<double> eigenvalues(DynamicMatrix<Dune::FieldMatrix<double,1,1>> A);
  template std::vector<float> eigenvalues(DynamicMatrix<Dune::FieldMatrix<float,2,2>> A);
  template std::vector<double> eigenvalues(DynamicMatrix<Dune::FieldMatrix<double,2,2>> A);
  template std::vector<float> eigenvalues(DynamicMatrix<Dune::FieldMatrix<float,3,3>> A);
  template std::vector<double> eigenvalues(DynamicMatrix<Dune::FieldMatrix<double,3,3>> A);

  
  // -----------------------------------------------------------------------------------------------
  // -----------------------------------------------------------------------------------------------
  
  template <class Entry>
  bool makePositiveSemiDefinite(DynamicMatrix<Entry>& A,
                                typename ScalarTraits<typename EntryTraits<Entry>::field_type>::Real alpha)
  {
    using Scalar = typename EntryTraits<Entry>::field_type;
    using Real = typename ScalarTraits<Scalar>::Real;

    static_assert(ScalarTraits<Scalar>::isBLAStype,
                  "makePositiveSemiDefinite only defined for LAPACK-compatible real scalar types");
    // TODO: extend to complex hermitian matrices.
    
    using namespace DynamicMatrixDetail;
    DynamicMatrix<Scalar> B = flatMatrix(A);
    
    // Compute eigenvalues and eigenvectors
    int const n = B.N();
    std::vector<Real> w(n);
    syev('V', n, &B[0][0], &w[0]);

    // Check if change is necessary in the first place.
    auto [lMin,lMax] = std::minmax_element(begin(w),end(w));             // min and max eigenvalues
    Real threshold = alpha*std::max(std::abs(*lMax),std::abs(*lMin));   // accuracy threshold
    if (*lMin >= threshold)               // minimum EV positive enough
      return false;                       // -> A is already positive semidefinite

    // Now we have A = \sum_l \lambda_l B_{*,l} B_{*,l}^T as a sum of rank-1 matrices.
    // (Formally, if B would contain columns of eigenvectors - but it contains rows of EV.)
    A = 0;
    for (int l=0; l<n; ++l)
    {
      Real lambda = std::max(w[l],threshold);
      for (int j=0; j<n; ++j)
        for (int i=0; i<n; ++i)
          scalarEntry(A,i,j) += B[i][l] * lambda * B[j][l];
    }
    return true;
  }
  
  template bool makePositiveSemiDefinite(DynamicMatrix<Dune::FieldMatrix<double,1,1>>& A, double alpha);
  template bool makePositiveSemiDefinite(DynamicMatrix<Dune::FieldMatrix<double,2,2>>& A, double alpha);
  template bool makePositiveSemiDefinite(DynamicMatrix<Dune::FieldMatrix<double,3,3>>& A, double alpha);

  // -----------------------------------------------------------------------------------------------

  
  template <class Scalar, int d>
  bool makePositiveSemiDefinite(Dune::FieldMatrix<Scalar,d,d>& A)
  {
    std::array<typename ScalarTraits<Scalar>::Real,d> w;
    
    // compute eigenvalues and eigenvectors
    auto B = A;
    syev('V', d, &B[0][0], &w[0]);

    // Now LAPACK syev works in column major mode, but Dune::FieldMatrix uses
    // row major. While this was irrelevant when passing B into syev (because of
    // it being symmetric), the matrix V with columns of eigenvalues which is computed
    // is in general not symmetric. Thus, if we interpret output B as matrix, we have
    // B = V^T. Remember this below, when using B.

    // Check if change is necessary in the first place.
    if (*std::min_element(begin(w),end(w)) >= 0)
      return false;

    // Now we have A = \sum_l \lambda_l V_{*,l} V_{*,l}^T = V * L * V^T, where
    // L = diag(lambda_1,...,lambda_n).
    // We shift negative eigenvalues \lambda_l < 0 to \lambda_l = 0.
    // This implementation may be inefficient for large matrices, but 
    // for small ones as typically used for field matrices (2x2 or 3x3 mostly),
    // this should be ok.
    Dune::FieldMatrix<Scalar,d,d> L = 0;
    for (int i=0; i<d; ++i)
    {
      L[i][i] = std::max(0.0,w[i]);
      std::cerr << w[i] << ' ';
    }
    std::cerr << "\n";
    
    // Reconstruct A from its modified spectral decomposition. Remember
    // that B = V^T.
    A = transpose(B)*L*B;

    return true;
  }
  
  template bool makePositiveSemiDefinite(Dune::FieldMatrix<double,2,2>& A);
  template bool makePositiveSemiDefinite(Dune::FieldMatrix<double,3,3>& A);
  
}




// ----------------------------------------------------------------------------

#ifdef UNITTEST

#include <algorithm>
#include <iostream>

int main(void) 
{
  Kaskade::DynamicMatrix<Dune::FieldMatrix<double,1,1>> A(3,5);
  for (int i=0; i<A.rows(); ++i)
    for (int j=0; j<A.cols(); ++j) 
      A[i][j] = i+j;
  
  std::cout << "A = \n" << A << "\n";
  Kaskade::DynamicMatrix<Dune::FieldMatrix<double,1,1>> U, V;
  std::vector<double> sing;
  std::tie(U,sing,V) = svd(A);
  
  std::cout << "U = \n" << U << "\n";
  std::cout << "V = \n" << V << "\n";
  std::cout << "singular values: "; std::copy(begin(sing),end(sing),std::ostream_iterator<double>(std::cout,", ")); std::cout << "\n";
  return 0;
}

#endif
