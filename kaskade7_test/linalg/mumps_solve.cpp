/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2012 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <cassert>
#include <memory>

#include "linalg/mumps_solve.hh"

namespace Kaskade
{

// Creator for DirectType::MUMPS factorizations to be plugged into the factory.

struct MUMPSCreator: public Creator<Factorization<double,int> > 
{
  MUMPSCreator(bool plugin) 
  {
    if (plugin)
      Factory<DirectType,Factorization<double,int>>::plugin(DirectType::MUMPS,this);
  }
  
  virtual std::unique_ptr<Factorization<double,int>>
  create(Creator<Factorization<double,int>>::Argument a) const
  {
    // Extract the matrix. We convert it directly to a mutable matrix as we
    // need to modify it for the mumps factorization.
    // TODO: this appears to thwart const correctness! Change interfaces!
    MatrixAsTriplet<double,int>& A = const_cast<MatrixAsTriplet<double,int>&>(std::get<0>(a));
    assert(A.nrows()==A.ncols());
    Factorization<double,int>::Options options = std::get<1>(a);

    return std::make_unique<MUMPSFactorization<double,int>>(A.nrows(),A.ridx,A.cidx,A.data,options);
  }
};

MUMPSCreator mumpsCreator(true);

#undef ICNTL
#undef JOB_INIT
#undef JOB_END
#undef USE_COMM_WORLD
} // namespace Kaskade 

#ifdef UNITTEST

#include <iostream>
#include <vector>

int main(void)
{
  using namespace Kaskade;

  using Scalar = double;

  std::vector<int> r {0, 1, 3, 4, 1, 0, 4, 2, 1, 2, 0, 2};
  std::vector<int> c {1, 2, 2, 4, 0, 0, 1, 3, 4, 1, 2, 2};
  std::vector<Scalar> v {3.0, -3.0, 2.0, 1.0, 3.0, 2.0, 4.0, 2.0, 6.0, -1.0, 4.0, 1.0};

  std::vector<Scalar> x {20.0, 24.0, 9.0, 6.0, 13.0};

  MUMPSFactorization<Scalar> f(5, r, c, v);
  f.solve(x);

  std::cout << "Supposed solution: {1, 2, 3, 4, 5}" << std::endl;
  std::cout << "Actual solution:   {" << x[0] << ", " << x[1] << ", " << x[2] << ", " << x[3] << ", " << x[4] << "}" << std::endl; 

  return 0;
}
#endif
