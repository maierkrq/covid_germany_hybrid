/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2016 Zuse Institute Berlin                                 */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef __CLAPACK_TO_STDLAPACK__
#define __CLAPACK_TO_STDLAPACK__

#include <algorithm>
#include <cstdio>
#include <vector>

#if defined(__CYGWIN__) || defined(StdLAPACK)
#define __CLPK_integer int
extern "C" {
int ilaenv_(int*,char*,int,char*,int,int*,int*,int*,int*);
}
#endif

#ifdef DEBUG
#define VDEBUG 1
#else
#define VDEBUG 0
#endif

int clapack_dgetrf(const enum CBLAS_ORDER Order, int mA, int nA, double *a, int ldaA, int *ipiv)
{
  __CLPK_integer m=mA, n=nA, lda=ldaA, info;
  dgetrf_( &m, &n, a, &lda, ipiv, &info );
  return info;
};

int clapack_dgetri(const enum CBLAS_ORDER Order, int nA, double *a, int ldaA, int *ipiv)
{
  __CLPK_integer n=nA, lda=ldaA, info;
  int minusOne=-1, one=1;
#if defined(__CYGWIN__) || defined(StdLAPACK)
  __CLPK_integer lwork = n* (__CLPK_integer) ilaenv_( &one, "dgetri",6, "",0, (int*) &n, &minusOne, &minusOne, &minusOne );
#else
  __CLPK_integer lwork = n* (__CLPK_integer) ilaenv_( &one, "dgetri", "", (int*) &n, &minusOne, &minusOne, &minusOne );
#endif
  std::vector<double> work(lwork);
  if (VDEBUG) printf("*** entering dgetri ***\n");
  dgetri_( &n, a, &lda, ipiv, &work[0], &lwork, &info );
  if (VDEBUG) printf("***dgetri finished***\n");
  return info;
};

int clapack_sgetrf(const enum CBLAS_ORDER Order, int mA, int nA, float *a, int ldaA, int *ipiv)
{
  __CLPK_integer m=mA, n=nA, lda=ldaA,info;
  sgetrf_( &m, &n, a, &lda, ipiv, &info );
  return info;
};

int clapack_sgetri(const enum CBLAS_ORDER Order, int nA, float *a, int ldaA, int *ipiv)
{
  __CLPK_integer n=nA, lda=ldaA, info;
  int minusOne=-1, one=1;
#if defined(__CYGWIN__) || defined(StdLAPACK)
  __CLPK_integer lwork = n* (__CLPK_integer) ilaenv_( &one, "sgetri",6, "",0, (int*) &n, &minusOne, &minusOne, &minusOne );
#else
  __CLPK_integer lwork = n* (__CLPK_integer) ilaenv_( &one, "sgetri", "", (int*) &n, &minusOne, &minusOne, &minusOne );
#endif
  std::vector<float> work(lwork);
  if (VDEBUG) printf("*** entering sgetri ***\n");
  sgetri_( &n, a, &lda, ipiv, &work[0], &lwork, &info );
  if (VDEBUG) printf("***sgetri finished***\n");
  return info;
};

int clapack_dgesv(const enum CBLAS_ORDER Order, int nA, int nrhsA, double *a, int ldaA, int *pivA, double *b, int ldbA)
{
  __CLPK_integer n=nA, nrhs=nrhsA, lda=ldaA, *ipiv=pivA, ldb=ldbA;
  __CLPK_integer info;
  if (VDEBUG) printf("*** entering dgesv ***\n");
  dgesv_(&n,&nrhs,a,&lda,ipiv,b,&ldb,&info);
  if (VDEBUG) printf("***dgesv finished***\n");
  return info;
};

int clapack_dgels(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE trans, int mA, int nA, int nrhsA, double *a, int ldaA, double *b, int ldbA) 
{
  __CLPK_integer m=mA, n=nA, nrhs=nrhsA, lda=ldaA, ldb=ldbA, info;
  __CLPK_integer mn = std::min(m,n);
  __CLPK_integer lwork = std::max( 1, mn + std::max( mn, nrhs )*n );
  char transc;
  if ( trans == CblasNoTrans ) transc = 'N';
  if ( trans == CblasTrans )   transc = 'T';
  std::vector<double> work(lwork);
  if (VDEBUG) printf("*** entering dgels ***\n");
  dgels_(&transc,&m,&n,&nrhs,a,&lda,b,&ldb,&work[0],&lwork,&info);
  if (VDEBUG) printf("***dgels finished***\n");
  return info;
};

// example for a new placeholder routine - use copy/paste, then adapt the name and parameter list,
// and uncomment the code
//
// int clapack_sgetrf(const enum CBLAS_ORDER Order, int mA, int nA, float *a, int ldaA, int *ipiv)
// {
//   printf("*** clapack_sgetrf: Please supply mapping routine to standard LAPACK routine sgetrf_ ***\n");
//   printf("*** by editing this routine in header file linalg/clapack_to_mac.hh\n");
//   printf("*** Execution of program aborted! ***\n");
//   exit(10);
// };
// 

#endif
