/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2012-2019 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef SCALAR_HH_
#define SCALAR_HH_

#include <complex>
#include <type_traits>

/// \internal 
// forward declarations
namespace Dune
{
 template <class Scalar, int n>        class FieldVector; 
 template <class Scalar, int n, int m> class FieldMatrix; 
}
/// \endinternal

namespace Kaskade {

  /**
   * \brief Helper class for working with scalar field types
   */
  template <class Scalar>
  class ScalarTraits
  {
  public:
    /**
     * \brief Whether the scalar is BLAS/LAPACK compatible or not.
     */
    static bool const isBLAStype = std::is_same_v<Scalar,float> || std::is_same_v<Scalar,double>;

    /**
     * \brief The real type on which the scalar field is based.
     *
     * For real types such as \a float or \a double, this is just \a Scalar. For complex types
     * such as \a std::complex<T>, this is \a T.
     */
    typedef Scalar Real;

    /**
     * \brief Conversion to the real type, ignoring the imaginary part if nonzero.
     */
    static Real real(Scalar const& s) { return s; }
  };

  // Specialization for standard complex types
  template <class T>
  class ScalarTraits<std::complex<T> > {
  public:
    static bool const isBLAStype = std::is_floating_point_v<T> && ScalarTraits<T>::isBLAStype;

    typedef T Real;

    static Real real(std::complex<T> const& z) { return z.real(); }
  };

  //-----------------------------------------------------------------------------------------------

  template <class Entry>      // Default implementation for scalars (float, double, complex<float>,...)
  struct EntryTraits
  {
    static int const rows = 1;    // number of scalar rows
    static int const cols = 1;    // number of scalar columns
    using transpose_type = Entry;
    using field_type = Entry;
    using real_type = typename ScalarTraits<field_type>::Real;
    static int const lapackLayout = std::is_floating_point<Entry>::value;

    static real_type frobenius_norm2(Entry const& x) { return x*x; } // TODO: needs specialization for complex
  };

  template <class Entry, int n, int m>
  struct EntryTraits<Dune::FieldMatrix<Entry,n,m>>
  {
    static int const rows = n*EntryTraits<Entry>::rows;
    static int const cols = m*EntryTraits<Entry>::cols;
    using transpose_type = Dune::FieldMatrix<typename EntryTraits<Entry>::transpose_type,m,n>;
    using field_type = typename EntryTraits<Entry>::field_type;
    using real_type = typename ScalarTraits<field_type>::Real;
    static int const lapackLayout = m==1 && EntryTraits<Entry>::lapackLayout;   // FieldMatrices are stored row by row

    static real_type frobenius_norm2(Dune::FieldMatrix<Entry,n,m> const& x) { return x.frobenius_norm2(); }
  };

  template <class Entry, int n>
  struct EntryTraits<Dune::FieldVector<Entry,n>>
  {
    static int const rows = n*EntryTraits<Entry>::rows;
    static int const cols = EntryTraits<Entry>::cols;
    using transpose_type = Dune::FieldMatrix<typename EntryTraits<Entry>::transpose_type,1,n>;
    using field_type = typename EntryTraits<Entry>::field_type;
    using real_type = typename ScalarTraits<field_type>::Real;
    static int const lapackLayout = EntryTraits<Entry>::lapackLayout;

    static real_type frobenius_norm2(Dune::FieldVector<Entry,n> const& x) { return x.two_norm2(); }
  };

  template <class Entry>
  typename EntryTraits<Entry>::real_type frobenius_norm2(Entry const& x)
  {
    return EntryTraits<Entry>::frobenius_norm2(x);
  }

  // ----------------------------------------------------------------------------------------------

  /**
   * \cond internals
   */
  namespace ScalarDetail
  {
    template <class K>
    struct Rank
    {
      static constexpr int value = 0;
    };

    template <class K, int n>
    struct Rank<Dune::FieldVector<K,n>>
    {
      static constexpr int value = 1;
    };

    template <class K, int n, int m>
    struct Rank<Dune::FieldMatrix<K,n,m>>
    {
      static constexpr int value = 2;
    };
  }
  /**
   * \endcond
   */

  /**
   * \ingroup linalgbasic
   * \brief Reports the rank of vector, matrix, and tensor types of static size.
   */
  template <class T>
  constexpr int rank = ScalarDetail::Rank<T>::value;


  // ----------------------------------------------------------------------------------------------

  /**
   * \ingroup utilities
   * \brief Reports the converted type.
   */
  template <class T, class Real>
  struct ConvertTo 
  {
    using type = std::conditional_t<std::is_convertible<T,Real>::value,Real,T>;
  };
  
  template <class T, int n, int m, class Real>
  struct ConvertTo<Dune::FieldMatrix<T,n,m>,Real>
  {
    using type = Dune::FieldMatrix<typename ConvertTo<T,Real>::type,n,m>;
  };
}

#endif
