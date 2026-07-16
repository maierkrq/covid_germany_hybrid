/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2020 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef UTILITIES_POWER_HH
#define UTILITIES_POWER_HH


namespace Kaskade
{
  /**
   * \ingroup utilities
   * \brief Computes integral powers of values in a multiplicative half-group.
   */
  // unfortunately not part of the standard C++99, but only GCC extension... here it is: <cmath>
  // template <class R> R power(R,int) __attribute__((deprecated));
  // 
  // No: std::pow works for built-in floating point values and std::complex only 
  // and calls transzendental functions, which may be slower than simple multiplication.

  template <class R>
  R power(R base, int exp)
  {
    if (exp>1 && exp%2==0) {
      R p = power(base,exp/2);
      return p*p;
    }
    else if (exp > 0)
      return base * power(base,exp-1);
    else if (exp < 0)
      return static_cast<R>(1)/power(base,-exp);
    else
      return static_cast<R>(1);
  }
  
  /**
   * \ingroup utilities
   * \brief returns the square of its argument.
   */
  template <class R>
  R square(R x)
  {
    return x*x;
  }
  
  /**
   * \ingroup utilities
   * \brief Computes faculties. Is this somewhere in the standard? no.
   */
  template <class R>
  R faculty(int n)
  {
    R result = 1;
    for (int i=2; i<=n; ++i)
      result *= i;
    return result;
  }
  
  /**
   * \cond internals
   */
  namespace PowerDetails
  {
    // look for x^2 = n by bisection
    constexpr int sqrti(int lo, int mid, int hi, int n)
    {
      return   lo == hi      ? lo
             : n / mid < mid ? sqrti(lo,(lo+mid)/2,mid-1,n)
             :                 sqrti(mid,(mid+hi+1)/2,hi,n);
    }
  }
  /**
   * \endcond 
   */
  
  /**
   * \ingroup utilities
   * \brief Computes floor(sqrt(n)) for integers n at compile time.
   */
  constexpr int sqrti(int n)
  {
    return PowerDetails::sqrti(0,(n+1)/2,n,n);
  }

}

#endif
