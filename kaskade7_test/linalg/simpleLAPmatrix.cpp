/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2016-2020 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <cassert>

#ifdef __cplusplus
extern "C" {
#endif
#include <cblas.h>

#include "lapack/lapacke.h"
// void dsyevd_(char*,char*,int*,double*,int*,double*,double*,int*,int*,int*,int*);
// void dgelsy_(int*,int*,int*,double*,int*,double*,int*,int*,double*,int*,double*,int*,int*);
// void dgesvd_(char*,char*,int*,int*,double*,int*,double*,double*,int*,double*,int*,double*,int*,int*);
#ifdef __cplusplus
}
#endif

#include "simpleLAPmatrix.hh"


// prevent unused variable from causing warnings (or errors if checking with -Werror)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"


namespace Kaskade
{

void MatrixKT_A_K(SLAPMatrix<double> A, SLAPMatrix<double>& K, SLAPMatrix<double>& out)
{
  assert(K.rows()==A.cols());
  assert(K.cols()==out.cols());
  assert(K.cols()==out.rows());
  SLAPMatrix<double> AK(A.rows(),K.cols());
  cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,AK.rows(),AK.cols(), K.rows(), 1.0,
              A.ptrToData(), A.LDA(), K.ptrToData(), K.LDA(), 0.0, AK.ptrToData(), AK.LDA());
  cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,out.rows(),out.cols(), AK.rows(), 1.0,
              K.ptrToData(), K.LDA(), AK.ptrToData(), AK.LDA(), 0.0, out.ptrToData(), out.LDA());
}

void MatrixMultiplication(SLAPMatrix<double> A, SLAPMatrix<double>& B, SLAPMatrix<double>& AB)
{
  assert(B.rows()==A.cols());
  assert(B.cols()==AB.cols());
  assert(A.rows()==AB.rows());
  cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,AB.rows(),AB.cols(), B.rows(), 1.0,
              A.ptrToData(), A.LDA(), B.ptrToData(), B.LDA(), 0.0, AB.ptrToData(), AB.LDA());
}

void MatrixMultiply(SLAPMatrix<double> A, std::vector<double>& in, std::vector<double>& out)
{
  if(in.size() !=A.cols())
  {
    std::cout << "Size Mismatch:" << in.size() << " " << A.cols() << std::endl;
  }
  assert(in.size()==A.cols());
  out.resize(A.rows());
  cblas_dgemv(CblasColMajor,CblasNoTrans,A.rows(),A.cols(), 1.0, A.ptrToData(), A.LDA(), &in[0],
              1, 0.0, &out[0], 1);
}

void TransposedMatrixMultiply(SLAPMatrix<double> A, std::vector<double>& in, std::vector<double>& out)
{
  assert(in.size()==A.rows());
  out.resize(A.cols());
  cblas_dgemv(CblasColMajor,CblasTrans,A.rows(),A.cols(), 1.0, A.ptrToData(), A.LDA(), &in[0], 1, 0.0, &out[0], 1);
}

void LeastSquares(SLAPMatrix<double> a, std::vector<double> const& b, std::vector<double>& x)
{
  x.reserve(std::max(a.cols(),a.rows()));
  x=b;
  int m=a.rows(), n=a.cols();
  int info = LAPACKE_dgels(LAPACK_COL_MAJOR,'N',m,n,1,a.ptrToData(),a.LDA(),&x[0],m);
  assert(info==0);
}

void SymmetricEigenvalues(SLAPMatrix<double> a, std::vector<double>& eig)
{
  eig.resize(a.rows());
  LAPACKE_dsyevd(LAPACK_COL_MAJOR,'N','U',a.rows(),a.ptrToData(),a.LDA(),&eig[0]);
}

//---------------------------------------------------------------------

template <>
void pseudoinverse(SLAPMatrix<double> A, SLAPMatrix<double>& X)
{
  assert(A.rows()==X.cols() && A.cols()==X.rows());

  int m=A.rows(), n = A.cols();
  SLAPMatrix<double> B(m,m,std::max(A.rows(),A.cols()));
  for (int i=0; i<m; ++i) B(i,i) = 1;

  //int jpvt[A.cols()]; // probably this causes a valgrind error/warning message: use of not inititalized value
  std::vector<int> jpvt(A.cols(),0);
  int rank, info;
  double cond=1.0e-10;
  info=LAPACKE_dgelsy(LAPACK_COL_MAJOR,m,n,B.cols(),A.ptrToData(),A.LDA(),B.ptrToData(),B.LDA(),&jpvt[0],cond,&rank);

  if(rank < std::min(A.rows(),A.cols()))
  {
    std::cout << "Warning: rank deficiency detected:" << std::endl;
    std::cout << "A is " << A.rows() << "x" << A.cols() << " but Rank(A) = " << rank << std::endl;
    std::cout << "Info:" << info << std::endl;
    for(int i=0; i<A.rows(); ++i)
    {
      for(int j=0; j<A.cols(); ++j)
        std::cout << A(i,j) << " ";
      std::cout << std::endl;
    }
  }

  for (int j=0; j<X.cols(); ++j)
    for (int i=0; i<X.rows(); ++i)
      X(i,j) = B(i,j);
}

void basisOfKernel(SLAPMatrix<double> A, SLAPMatrix<double>& Basis) 
{
  assert(A.cols()==Basis.rows() && (A.cols()-A.rows())==Basis.cols());

  int m=A.rows(), n=A.cols();
  double zero=0.0;
  std::vector<double> sing(std::max(m,n)), superb(std::min(m,n));
  SLAPMatrix<double> B(n,n);
  int info = LAPACKE_dgesvd(LAPACK_COL_MAJOR,'N','A',A.rows(),A.cols(),A.ptrToData(),A.LDA(),&sing[0],&zero,1,B.ptrToData(),A.cols(),
          &superb[0]);
  assert(info==0);
  
  for (int j=0; j<Basis.cols(); ++j)
    for (int i=0; i<Basis.rows(); ++i)
      Basis(i,j) = B(j+A.rows(),i);
}


}  // namespace Kaskade

#pragma GCC diagnostic pop


// explicit instantiation
