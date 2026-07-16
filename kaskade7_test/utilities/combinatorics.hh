/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2019 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef COMBINATORICS_HH
#define COMBINATORICS_HH

#include <algorithm>
#include <numeric>
#include <cassert>
#include <functional>

namespace Kaskade
{
  /**
   * \ingroup utilities
   * \brief Computes the binomial coefficient \f$ \binom{n}{k} \f$.
   *
   * Precondition: n >= k >= 0
   */
  constexpr int binomial(int n, int k)
  {
    assert(0<=k && k<=n);

    if (2*k>n)
      k = n-k;

    if (k==0)
      return 1;

    // Use the product form for efficiency and overflow avoidance.
    int r = n;
    for (int i=2; i<=k; ++i)
    {
      r *= n+1-i;
      r /= i;
    }

    return r;
  }

  // ----------------------------------------------------------------------------------------------

  /**
   * \ingroup utilities
   * \brief Computes the multinomial coefficient \f$ \binom{n}{k_1,\dots,k_m}\f$, where \f$ n = \sum_{i=1}^m k_i \f$.
   * \tparam m the number of multinomials
   * \param ks an array of nonegative integers
   *
   * The multinomial coefficient is \f[ \binom{n}{k_1,\dots,k_m} = \frac{n!}{k_1! \cdots k_m!}. \f]
   */
  template <size_t m>
  size_t multinomial(std::array<size_t,m> const& ks)
  {
    size_t r = 1;
    size_t n = std::accumulate(begin(ks),end(ks),0);

    // Compute the multinomial factor as the binomial one above.
    for (size_t k: ks)
      for (size_t i=1; i<=k; ++i)
      {
        r *= n;
        r /= i;
        --n;
      }

    return r;
  }

  // ----------------------------------------------------------------------------------------------

  /**
   * \ingroup utilities
   * \brief Computes the number of multiindices of order m and dimension d.
   *
   * A multiindex of dimension d is a d-tuple \f$ (i_1,\dots,i_d)\f$ of nonnegative integers that sums
   * up to at most m: \f[ \sum_{j=1}^d i_j \le m \f]
   * Multiindices correspond directly to points on a cartesian grid covering the
   * d-dimensional unit simplex and to Bezier functions.
   */
  constexpr int numberOfMultiindices(int d, int m)
  {
    int n = 0;

    if (d==1)
      return m+1;

    if (d==2)
      return (m+2)*(m+1)/2;

    for (int i=0; i<=m; ++i)
      n += numberOfMultiindices(d-1,i);
    return n;
  }

  // ----------------------------------------------------------------------------------------------

  /**
   * \ingroup utilities
   * \brief Random access to multiindices.
   *
   * This yields multiindices of length d that add up to at most m:
   * \f[ \sum_{i=0}^{d-1} k_i \le m \f]
   *
   * \tparam d the number of indices
   * \param m the maximum sum of the indices
   * \param n compute the n-th index (0 <= n < numberOfMultiindices(d,m) has to hold).
   */
  template <int d>
  std::array<size_t,d> multiIndex(size_t m, size_t n)
  {
    static_assert(d>=1,"multiindices only defined for d>0");
    assert(n < numberOfMultiindices(d,m));

    std::array<size_t,d> ks;

    if constexpr (d==1)
    {
      ks[0] = n;
    }
    else
    {
      size_t i = 0;
      size_t l = 0;
      while (n >= l+numberOfMultiindices(d-1,m-i))
      {
        l += numberOfMultiindices(d-1,m-i);
        ++i;
      }
      assert(n >= l);
      ks[d-1] = i;
      auto hs = multiIndex<d-1>(m-i,n-l);
      std::copy(begin(hs),end(hs),begin(ks));
    }

    return ks;
  }

  // ----------------------------------------------------------------------------------------------

  /**
   * \ingroup utilities
   * \brief Computes the next nonnegative multiindex of order m. 
   * 
   * Using this function one can cycle through all integer tuples that sum up to at
   * most m.
   */
  template <class It>
  void next_multiindex(It first, It last, int const m)
  {
    // If the sum limit is reached, we need to "make room" for an increment.
    if (std::accumulate(first,last,0)==m)
    {
      // Find the first (i.e. least significant) nonzero entry and set it to zero.
      first = std::find_if(first,last,std::bind2nd(std::not_equal_to<int>(),0));
      assert(first!=last);
      *first = 0;
      // Move to the next higher (i.e. more significant) entry.
      ++first;
    }

    // Increment the relevant index.
    if (first!=last)
      ++(*first);
  }
} /* end of namespace Kaskade */


#endif
