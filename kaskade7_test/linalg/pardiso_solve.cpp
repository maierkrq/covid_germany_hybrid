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

#include <cassert>
#include <memory>
#include <cstdio>

#include "linalg/pardiso_solve.hh"
#include "linalg/triplet.hh"

// alternatively one could include mkl_pardiso.h and mkl_spblas.h but this would require to always provide these headers for installation of Kaskade
extern "C" {
  void pardiso(void *, const int *, const int *, const int *,  const int *, const int *, const void *, const int *,
               const int *, int *, const int *, int *, const int *, void *, void *, int * );
  void mkl_dcsrcoo(const int *, const int *, double *, int *, int *, int *, double *,  int *,  int *,  int *);
  void mkl_scsrcoo(const int *, const int *, float *, int *, int *, int *, float *,  int *,  int *,  int *);
}

namespace Kaskade
{

void mklTripletToCompressedColumn(int n, std::vector<double> & az, std::vector<int> & ai,
                                  std::vector<int> & ap, int nnz, std::vector<double> & values,
                                  std::vector<int> & ridx, std::vector<int> & cidx)
{
  std::vector<int> job {2, 0, 0, 0, 0, 0};
  int idum = 0;
  mkl_dcsrcoo(job.data(), &n, az.data(), ai.data(), ap.data(), &nnz, values.data(), ridx.data(), cidx.data(), &idum);
}
void mklTripletToCompressedColumn(int n, std::vector<float> & az, std::vector<int> & ai,
                                  std::vector<int> & ap, int nnz, std::vector<float> & values,
                                  std::vector<int> & ridx, std::vector<int> & cidx)
{
  std::vector<int> job {2, 0, 0, 0, 0, 0};
  int idum = 0;
  mkl_scsrcoo(job.data(), &n, az.data(), ai.data(), ap.data(), &nnz, values.data(), ridx.data(), cidx.data(), &idum);
}

//---------------------------------------------------------------------

template <class Scalar,class SparseIndexInt, class DIL>
PardisoFactorization<Scalar,SparseIndexInt,DIL>::~PardisoFactorization()
{
  phase = -1; // Release internal memory.
  pardiso(pt.data(), &maxfct, &mnum, &mtype, &phase,
          &n, &ddum, ap.data(), ai.data(), &idum, &nrhs,
          iparm.data(), &msglvl, &ddum, &ddum, &error);

  assert(!error);
}

template <class Scalar,class SparseIndexInt, class DIL>
void PardisoFactorization<Scalar,SparseIndexInt,DIL>::swap(PardisoFactorization & other)
{
  std::swap(n, other.n);
  std::swap(ap, other.ap);
  std::swap(ai, other.ai);
  std::swap(az, other.az);
  std::swap(phase, other.phase);
  std::swap(maxfct, other.maxfct);
  std::swap(mnum, other.mnum);
  std::swap(mtype, other.mtype);
  std::swap(nrhs, other.nrhs);
  std::swap(msglvl, other.msglvl);
  std::swap(error, other.error);
  std::swap(iparm, other.iparm);
  std::swap(pt, other.pt);
}

template <class Scalar,class SparseIndexInt, class DIL>
PardisoFactorization<Scalar,SparseIndexInt,DIL> &
PardisoFactorization<Scalar,SparseIndexInt,DIL>::operator=(PardisoFactorization && other)
{
  swap(other);
  return *this;
}

template <class Scalar,class SparseIndexInt, class DIL>
PardisoFactorization<Scalar,SparseIndexInt,DIL>::PardisoFactorization(SparseIndexInt n_,
                                                                 std::vector<SparseIndexInt>const & rid,
                                                                 std::vector<SparseIndexInt>const & cid,
                                                                 std::vector<Scalar>const & val,
                                                                 MatrixProperties property)
  : n(n_), ap(n+1), ai(rid.size()), az(rid.size()), maxfct(1),
    mnum(1), mtype(0), nrhs(1), msglvl(0), error(0), idum(0), ddum(0.0)
{
  int nnz(rid.size());
  std::vector<int> ridx(rid), cidx(cid);
  std::vector<Scalar> values(val);
  if (cidx.size()!=nnz || values.size()!=nnz)
  {
    throw DirectSolverException("wrong size of column indices or values!", __FILE__, __LINE__);
  }

  std::string s = "UNKNOWN";
  switch (property)
  {
    case MatrixProperties::GENERAL:
      s = "GENERAL";
      mtype = 11;
      break;
    case MatrixProperties::SYMMETRIC:
      s = "SYMMETRIC";
      mtype = -2;
      // swap rows and columns, since maybe only lower triangle is given but pardiso needs upper
      std::swap(ridx, cidx);
      deleteLowerTriangle(ridx,cidx,values);
      nnz = values.size();
      break;
    case MatrixProperties::SYMMETRICSTRUCTURE:
      s = "SYMMETRICSTRUCTURE";
      mtype = 1;
      break;
    case MatrixProperties::POSITIVEDEFINITE:
      s = "POSITIVEDEFINITE";
      mtype = 2;
      // swap rows and columns, since maybe only lower triangle is given but pardiso needs upper
      std::swap(ridx, cidx);
      deleteLowerTriangle(ridx,cidx,values);
      nnz = values.size();
      break;
  }

  if (this->getVerbose()>=2)
  {
    std::cout << "Pardiso" << " solver, n=" << n << ", nnz=" << nnz << ", matrix is " << s << std::endl;
  }

  iparm[0] = 1;         // No solver default
  iparm[1] = 0;         // Fill-in reducing reordering (0=minimum degree, 2=METIS nested dissection, 3=parallel nested dissection)
  iparm[7] = 2;         // Max numbers of iterative refinement steps
  iparm[9] = 13;        // Perturb the pivot elements with 1E-13
  iparm[10] = 1;        // Use nonsymmetric permutation and scaling MPS
  iparm[12] = 1;        // Maximum weighted matching algorithm is switched-on
  iparm[13] = 0;        // Output: Number of perturbed pivots
  iparm[17] = -1;       // Output: Number of nonzeros in the factor LU
  iparm[18] = -1;       // Output: Mflops for LU factorization
  iparm[24] = 0;        // use 1 for parallel forward/backward solve
  iparm[27] = std::is_same<Scalar, double>::value ? 0 : 1; // single or double precision
  iparm[36] = -80;      // convert matrix internally to variable bsr format
  iparm[34] = 1;        // zero based indexing

  mklTripletToCompressedColumn(n, az, ai, ap, nnz, values, ridx, cidx);

  phase = 11; 
  
  pardiso(pt.data(), &maxfct, &mnum, &mtype, &phase,
          &n, az.data(), ap.data(), ai.data(), &idum, &nrhs,
          iparm.data(), &msglvl, &ddum, &ddum, &error);

  if (error)
  {
    throw DirectSolverException(std::string("phase 11, error=") + std::to_string(error), __FILE__, __LINE__);
  }

  phase = 22; 
  
  pardiso(pt.data(), &maxfct, &mnum, &mtype, &phase,
          &n, az.data(), ap.data(), ai.data(), &idum, &nrhs,
          iparm.data(), &msglvl, &ddum, &ddum, &error);

  if (error)
  {
    throw DirectSolverException(std::string("phase 22, error=") + std::to_string(error), __FILE__, __LINE__);
  }
}

template <class Scalar,class SparseIndexInt, class DIL>
void PardisoFactorization<Scalar,SparseIndexInt,DIL>::solve(std::vector<Scalar> const& b,
                                                               std::vector<Scalar>& x, int nr, bool transposed) const
{
  if (b.size()<n*nr)
  {
    throw DirectSolverException("b too small!", __FILE__, __LINE__);
  }
  if (&x == &b)
  {
    throw DirectSolverException("x must be different object than b!", __FILE__, __LINE__);
  }
  x.resize(n*nr);
  phase = 33;

  if(transposed) iparm[11]=1;
  else iparm[11]=0;

  pardiso(pt.data(), &maxfct, &mnum, &mtype, &phase,
          &n, az.data(), ap.data(), ai.data(), &idum, &nr,
          iparm.data(), &msglvl, const_cast<Scalar*>(b.data()), x.data(), &error);

  if (error)
  {
    throw DirectSolverException(std::string("phase 33, error=") + std::to_string(error), __FILE__, __LINE__);
  }
}


// Creator for DirectType::PARDISO factorizations to be plugged into the factory.

template <class Scalar>
struct PARDISOCreator: public Creator<Factorization<Scalar,int>>
{
  PARDISOCreator(bool plugin)
  {
    if (plugin)
      Factory<DirectType,Factorization<Scalar,int> >::plugin(DirectType::PARDISO,this);
  }
  
  virtual std::unique_ptr<Factorization<Scalar,int>>
  create(typename Creator<Factorization<Scalar,int>>::Argument a) const
  {
    MatrixAsTriplet<Scalar,int> const& A = std::get<0>(a);
    assert(A.nrows()==A.ncols());
    MatrixProperties property = std::get<1>(a).matrixProperty;

    return std::make_unique<PardisoFactorization<Scalar,int>> (A.nrows(),A.ridx,A.cidx,A.data,property);
  }
};

PARDISOCreator<double> pardisoCreator(true);
PARDISOCreator<float> pardisoSingleCreator(true);

} // End of kaskade namespace ---------------------------------------------------------------------


#ifdef UNITTEST

#include <iostream>
#include <vector>

int main(void)
{
  using namespace Kaskade;

  using Scalar = double;

  {
    std::vector<int> r {0, 1, 3, 4, 1, 0, 4, 2, 1, 2, 0, 2};
    std::vector<int> c {1, 2, 2, 4, 0, 0, 1, 3, 4, 1, 2, 2};
    std::vector<Scalar> v {3.0, -3.0, 2.0, 1.0, 3.0, 2.0, 4.0, 2.0, 6.0, -1.0, 4.0, 1.0};

    std::vector<Scalar> x {20.0, 24.0, 9.0, 6.0, 13.0};

    PardisoFactorization<Scalar> f(5, r, c, v);
    f.solve(x);

    std::cout << "Supposed solution: {1, 2, 3, 4, 5}" << std::endl;
    std::cout << "Actual solution:   {" << x[0] << ", " << x[1] << ", " << x[2] << ", " << x[3] << ", " << x[4] << "}" << std::endl;

    r = {0, 0, 1, 1, 1, 2, 2};
    c = {0, 1, 0, 1, 2, 1, 2};
    v = {2.0, -1.0, -1.0, 2.0, -1.0, -1.0, 2.0};

    x = {1.0, 0.0, 1.0};

    f = PardisoFactorization<Scalar>(3, r, c, v, MatrixProperties::POSITIVEDEFINITE);
    f.solve(x);

    std::cout << "\nSupposed solution: {1, 1, 1}" << std::endl;
    std::cout << "Actual solution:   {" << x[0] << ", " << x[1] << ", " << x[2] << "}" << std::endl;

    x = {1.0, 0.0, 1.0};

    f = PardisoFactorization<Scalar>(3, r, c, v, MatrixProperties::SYMMETRIC);
    f.solve(x);

    std::cout << "\nSupposed solution: {1, 1, 1}" << std::endl;
    std::cout << "Actual solution:   {" << x[0] << ", " << x[1] << ", " << x[2] << "}" << std::endl;

    x = {1.0, 0.0, 1.0};

    f = PardisoFactorization<Scalar>(3, r, c, v, MatrixProperties::SYMMETRICSTRUCTURE);
    f.solve(x);

    std::cout << "\nSupposed solution: {1, 1, 1}" << std::endl;
    std::cout << "Actual solution:   {" << x[0] << ", " << x[1] << ", " << x[2] << "}" << std::endl;
  }

  {
//    std::vector<int> r {0, 0, 0, 0, 1, 1, 1, 2, 2, 3, 3, 4, 4, 4, 5, 5, 6, 7};
//    std::vector<int> c {0, 2, 5, 6, 1, 2, 4, 2, 7, 3, 6, 4, 5, 6, 5, 7, 6, 7};
    // provide lower triangle, will be transposed internally
    std::vector<int> r {0, 2, 5, 6, 1, 2, 4, 2, 7, 3, 6, 4, 5, 6, 5, 7, 6, 7};
    std::vector<int> c {0, 0, 0, 0, 1, 1, 1, 2, 2, 3, 3, 4, 4, 4, 5, 5, 6, 7};
    std::vector<Scalar> v {7.0, 1.0, 2.0, 7.0, -4.0, 8.0, 2.0, 1.0, 5.0, 7.0, 9.0, 5.0, 1.0, 5.0, -1.0, 5.0, 11.0, 5.0};

    std::vector<Scalar> x(8, 1.0);

    PardisoFactorization<Scalar> f(8, r, c, v, MatrixProperties::SYMMETRIC);
    f.solve(x);

    std::cout << "\nSupposed solution: {-0.041860, -0.003413, 0.117250, -0.112640, 0.024172, -0.107633, 0.198720, 0.190383}" << std::endl;
    std::cout << "Actual solution:   {";
    std::cout << x[0] << ", " << x[1] << ", " << x[2] << ", " << x[3] << ", ";
    std::cout << x[4] << ", " << x[5] << ", " << x[6] << ", " << x[7] << "}" << std::endl;

    r = {0, 0, 0, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4};
    c = {0, 1, 3, 0, 1, 2, 3, 4, 0, 2, 3, 1, 4};
    v = { 1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0, 2.0, 7.0, 8.0, -5.0};

    x = std::vector<Scalar>(5, 1.0);

    PardisoFactorization<Scalar> g(5, r, c, v);
    g.solve(x);

    std::cout << "\nSupposed solution: {-0.522321, -0.008929, 1.220982, -0.504464, -0.214286}" << std::endl;
    std::cout << "Actual solution:   {";
    std::cout << x[0] << ", " << x[1] << ", " << x[2] << ", " << x[3] << ", " << x[4] << "}" << std::endl;
  }

  return 0;
}
#endif
