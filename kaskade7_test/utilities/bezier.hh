/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2019-2020 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef BEZIER_HH
#define BEZIER_HH

#include "dune/common/fvector.hh"

#include "fem/barycentric.hh"
#include "utilities/combinatorics.hh"
#include "utilities/power.hh"


namespace Kaskade
{
  /**
   * \ingroup utilities
   * \brief Computes the i-th Bezier function of order p at the given point in the unit simplex.
   *
   * The Bezier functions are the summands in \f[ 1 = 1^p = \left(\sum_{j=0}^d \lambda_j\right)^p, \f]
   * where \f$ \lambda_j \f$ are the barycentric coordinates.
   *
   * \param p the order of Bezier functions
   * \param i the number of the Bezier function of order p. 0 <= i < numberOfMultiindices(d,p) has to hold.
   * \param x the simplex coordinates, i.e. \f$ \lambda = (x,1-\|x\|_1) \f$.
   *
   * \see numberOfMultiindices
   */
  template <int d, class Real>
  Real bezier(int p, int i, Dune::FieldVector<Real,d> const& x)
  {
    auto b = barycentric2(x);
    if (d==1)
      return binomial(p,i)*power(b[0],p-i)*power(b[1],i);
    else if (d==2)
    {
      // First get the i-th multiindex of integers such that k[0]+...+k[d-1] <= p.
      auto ks = multiIndex<d>(p,i);
      // Then compute k[d] such that k[0]+...+k[d] == p
      std::array<size_t,d+1> Ks;
      std::copy(begin(ks),end(ks),begin(Ks));
      Ks[d] = p - std::accumulate(begin(ks),end(ks),0);

      // Now the Bezier function is multinomial(p,Ks) * lambda_0^Ks[0] * ... * lambda_d^Ks[d].

      Real r = multinomial(Ks);
      for (int j=0; j<=d; ++j)
        r *= power(b[j],Ks[j]);
      return r;
    }

    static_assert(d<=2,"not yet implemented");
    return -1; // never get here
  }



}

#endif
