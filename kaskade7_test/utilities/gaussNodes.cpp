/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2018-2018 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <cassert>
#include <iostream>
#include <vector>

extern "C" {
#include "linalg/lapack/lapacke.h"
}

#include "utilities/detailed_exception.hh"
#include "utilities/gaussNodes.hh"

namespace Kaskade 
{
  namespace
  {
    void computeTridiagonalEigenvalues(std::vector<double>& diag, std::vector<double>& subdiag)
    {
      auto const n = diag.size();
      assert(subdiag.size()+1>=n);
      int info = LAPACKE_dstev(LAPACK_COL_MAJOR,'N',n,&diag[0],&subdiag[0],nullptr,1);
      assert(info>=0);
      if (info>0)
        throw LinearAlgebraException("Lapack routine dstev failed to converge",__FILE__,__LINE__);
    }
  }
  
  // Computes the n Lobatto nodes on [-1,1]
  std::vector<double> lobattoNodes(int const n)
  {
    assert(n >= 1);
    std::vector<double> t(n);
    
    // The interior Lobatto points are eigenvalues of the symmetric tridiagonal matrix
    // with A_i,i+1 = (i+1)*(i+3)/(2*i+3)/(2*i+5) and every other uppper diagonal
    // entries zero. See, e.g., J. Shen, T. Tang and L. Wang, Spectral Methods:
    // Algorithms, Analysis and Applications, Springer Series in Compuational
    // Mathematics, 41, Springer, 2011.
    //
    // Of course, this works only for n at least 4.
    t[0] = -1;
    t[n-1] = 1;
    
    if (n==1)
      t[0] = 0;
    if (n==3)
      t[1] = 0;
    
    if (n<=3)
      return t;
    
    // Build up matrix. Symmetric tridiagonal matrices are stored in two arrays, one 
    // for diagonal and one for subdiagonal.
    std::vector<double> diag(n-2,0.0), subdiag(n-3);
    for (int i=0; i<n-3; ++i)
      subdiag[i] = sqrt((i+1)*(i+3)/(2.*i+3)/(2.*i+5));
    
    // compute eigenvalues
    computeTridiagonalEigenvalues(diag,subdiag);
    
    // write result
    std::copy(diag.begin(),diag.end(),begin(t)+1);
    
    return t;
  }
  
  // Computes the n Radau nodes on ]-1,1]
  std::vector<double> radauNodes(int const n)
  {
    assert(n >= 1);
    std::vector<double> t(n);
    
    // The interior Radau points are the zeros of the Jacobi polynomials P^(1,0) and are obtained as 
    // eigenvalues of the symmetric tridiagonal matrix with 
    // A_00 = -1/3, A_jj = -1/(2j+1)/(2j+3), A_j(j-1) = 2j(j+1)/(2j+1)/sqrt(2j(2j+2))
    // See Deuflhard/Bornemann 6.3.2 and Wikipedia on Jacobi polynomials.
    
    std::vector<double> diag(n-1,-1.0/3), subdiag(n-1); // we provide larger subdiagonal since address of subdiag[0] is taken
    for (int i=1; i<n-1; ++i)
    {
      diag[i]      = -1.0/(2.0*i+1)/(2.0*i+3);
      subdiag[i-1] = 2.*i*(i+1)/(2.0*i+1)/std::sqrt(2*i*(2.0*i+2));
    }
    
    // compute eigenvalues
    computeTridiagonalEigenvalues(diag,subdiag);
    
    std::copy(diag.begin(),diag.end(),begin(t));
    t[n-1] = 1;
    
    return t;
  }
  
  // Computes the n Gauss nodes on ]-1,1[
  std::vector<double> gaussNodes(int const n)
  {
    assert(n >= 1);
    
    
    // trivial case
    if (n == 1)
      return std::vector<double>{0.0};
    
    // The Gauss points are obtained by the Golub-Welsch algorithm as eigenvalues of the symmetric 
    // tridiagonal matrix with A_jj = 0, A_j(j-1) = .5 / sqrt(1-(2*j)^(-2)).
    // See Deuflhard/Bornemann 6.3.2 and Wikipedia on Jacobi polynomials.
    
    std::vector<double> diag(n,0.0), subdiag(n); // we provide larger subdiagonal since address of subdiag[0] is taken
    for (int i=1; i<n; ++i)
      subdiag[i-1] = .5 / sqrt(1-1.0/(4*i*i));
    
    // compute eigenvalues
    computeTridiagonalEigenvalues(diag,subdiag);
    
    return diag;
  }
  

}

