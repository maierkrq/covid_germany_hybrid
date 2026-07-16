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
#ifndef TENSOR_HH
#define TENSOR_HH


#include "dune/common/fvector.hh"
#include "dune/common/fmatrix.hh"

#include "utilities/scalar.hh"

namespace Kaskade
{

  /**
   * \ingroup linalgbasic
   * \brief A compile-time index range.
   */
  template <int first, int last>
  struct StaticIndexRange {};

  // forward declaration
  template <class Entry, int Size0, int ... Sizes>
  class Tensor;

  /**
   * \cond internals
   */
  namespace TensorDetail
  {
    template <class Entry, int Size>
    struct Prepend { using type = Tensor<Entry,Size>; };

    template <class Entry, int SizePre, int Size0, int ... Sizes>
    struct Prepend<Tensor<Entry,Size0,Sizes...>,SizePre> { using type = Tensor<Entry,SizePre,Size0,Sizes...>; };
  }

  // Provide the rank of tensors with the same interface as for vectors and matrices.
  namespace ScalarDetail
  {
      template <class K, int Size0, int ... Sizes>
      struct Rank<Tensor<K,Size0,Sizes...>>
      {
        static constexpr int value = Tensor<K,Size0,Sizes...>::rank;
      };
  }
  /**
   * \endcond
   */

  namespace TensorAccess
  {
    /**
     * \ingroup linalgbasic
     * \brief Convenience static range denoting the full available range, analogous to Matlab's :.
     *
     * Use this to take complete dimensions when slicing tensors:
     * \code
     *     Tensor<double,2,2,2> t(0.0);
     *     using TensorAccess::_;
     *     std::cout << "vector is " <<  t(_,_,1);
     * \endcode
     * \relates Tensor
     */
    constexpr StaticIndexRange<0,std::numeric_limits<int>::max()> _;
  }

  /**
   * \ingroup linalgbasic
   * \brief A class for representing tensors of arbitrary static rank and extents.
   * \tparam Entry the type of tensor entries
   * \tparam Size0 the extent of the first dimension (must be >0)
   * \tparam Sizes the extents of the remaining dimensions
   */
  template <class Entry, int Size0, int ... Sizes>
  class Tensor: public Dune::FieldVector<Tensor<Entry,Sizes...>,Size0>
  {
    using Base = Dune::FieldVector<Tensor<Entry,Sizes...>,Size0>;
    static_assert(Size0>0);
    static constexpr int Size1 = Tensor<Entry,Sizes...>::dimension;
  public:

    /**
     * \brief The rank of the tensor.
     */
    static constexpr int rank = Tensor<Entry,Sizes...>::rank + 1;

    /**
     * \name Constructors
     * @{
     */
    Tensor() = default;

    Tensor(Entry const& e)
    {
      (*this) = e;
    }

    /**
     * @}
     */

    /**
     * \name Assignment
     * @{
     */
    /**
     * \brief Assignment from a scalar.
     */
    Tensor& operator=(Entry const& e)
    {
      for (auto& t: *this)
        t = e;
      return *this;
    }

    Tensor& operator=(Dune::FieldMatrix<Entry,Size0,Size1> const& a)
    {
      static_assert(rank==2);
      if constexpr (rank==2)
        for (int i=0; i<Size0; ++i)
            (*this)[i] = a[i];
      return *this;
    }
    /**
     * @}
     */

    /**
     * \name Element and subtensor access
     *
     * The tensor class supports element access via the subscript operator [], indexed by integers, and
     * the access via the call operator (). Both have different semantics. While subscript provides read
     * and write access by returning (mutable or const) references to entries (or sub-tensors ordered by
     * dimension), the call operator returns subtensors only by value, but is not restricted to ordered
     * dimension access or complete ranges along dimensions.
     *
     * #### Examples
     *
     * Let \f$ a_{ijk} \f$ be a rank-3 tensor.
     *
     * - `a[0]` yields a reference to the rank-2 subtensor \f$ a_{0jk} \f$
     * - `a[0][1][2]` yields a reference to the tensor entry \f$ a_{012} \f$
     * - `a(0,1,2)` yields a copy of the tensor entry \f$ a_{012} \f$
     * - `a(StaticIndexedRange<0,2>(),1,2)` yields the rank-1 tensor (vector) \f$ (a_{i12})_{i\in\{0,1\}} \f$ by value
     * - `a(_,_,1)` yields the rank-2 tensor (convertible to matrix) \f$ (a_{ij1})_{i,j} \f$ by value, taking the
     *    complete dimensions 0 and 1.
     *
     * \see _
     *
     * @{
     */
    template <class ... Access>
    auto operator()(int i0, Access... is) const
    {
      return (*this)[i0](is...);
    }

    template <int first, int last, class ... Access>
    auto operator()(StaticIndexRange<first,last> , Access... is) const
    {
      constexpr int s = std::min(Size0,last-first);
      using SubT = decltype( (*this)[0](is...) );
      using R = typename TensorDetail::Prepend<SubT,s>::type;
      R r;
      for (int i=0; i<s; ++i)
        r[i] = (*this)[i+first](is...);
      return r;
    }

    /**
     * @}
     */

    /**
     * \brief Implicit conversion of rank 2 tensors to Dune::FieldMatrix
     */
    operator std::conditional_t<rank==2,
                                Dune::FieldMatrix<Entry,Size0,Base::block_type::dimension>,
                                Base&> () const
    {
      if constexpr(rank==2)
      {
        Dune::FieldMatrix<Entry,Size0,Base::value_type::dimension> A;
        // We want that A[i][j] = (*this)[i][j]. A[i] are the rows of A (of type FieldVector),
        // and (*this)[i] are Tensors of rank 1, and thus derived from FieldVector. So we can
        // assign the matrix row by row.
        for (int i=0; i<Size0; ++i)
          A[i] = (*this)[i];
        return A;
      }
      else
        return *this;
    }

  };

  // ----------------------------------------------------------------------------------------------

  // Terminal of parameter pack recursion
  template <class Entry, int Size0>
  class Tensor<Entry,Size0>: public Dune::FieldVector<Entry,Size0>
  {
    using Base = Dune::FieldVector<Entry,Size0>;
  public:

    using field_type = typename EntryTraits<Entry>::field_type;

    static constexpr int rank = 1;

    Tensor() = default;

    Tensor(Entry const& e)
    {
      (*this) = e;
    }

    Tensor& operator=(Entry const& e)
    {
      for (auto& t: *this)
        t = e;
      return *this;
    }

    Tensor& operator=(Base const& v)
    {
      static_cast<Base&>(*this) = v;
      return *this;
    }

    Entry operator()(int i) const
    {
      return (*this)[i];
    }

    template <int first, int last>
    auto operator()(StaticIndexRange<first,last>) const
    {
      constexpr int s = std::min(Size0,last-first);

      Tensor<Entry,s> r;
      for (int i=0; i<s; ++i)
        r[i] = (*this)[i+first];
      return r;
    }

    Tensor& operator*=(field_type const& x)
    {
      static_cast<Base&>(*this) *= x;
      return *this;
    }
  };


  // Template specialization for accessing type information.
  template <class Entry, int ... Sizes>
  struct EntryTraits<Tensor<Entry,Sizes...>>
  {
    using field_type = typename EntryTraits<Entry>::field_type;
    using real_type = typename ScalarTraits<field_type>::Real;
    static int const lapackLayout = false;
  };

  // ----------------------------------------------------------------------------------------------

  /**
   * \brief Scalar times tensor multiplication.
   * \relates Tensor
   */
  template <class Entry, int Size0, int ...Sizes>
  Tensor<Entry,Size0,Sizes...> operator * (typename EntryTraits<Entry>::field_type a, Tensor<Entry,Size0,Sizes...> A)
  {
    A *= a;
    return A;
  }

  // ----------------------------------------------------------------------------------------------

  // Tensor output
  template <class Entry, int ... Sizes>
  std::ostream& operator<<(std::ostream& out, Tensor<Entry,Sizes...> const& t)
  {
    out << "[" << t[0];
    for (int i=1; i<std::size(t); ++i)
      out << "," << t[i];
    out << "]";
    return out;
  }

}

// ------------------------------------------------------------------------------------------------

/**
 * \cond internals
 */
namespace Dune
{
  // For interoperability with Dune (and as we put Tensors as entries into FieldVectors,
  // we need to specify the underlying field type.
  template <class K, int Size0, int ... Sizes>
  struct FieldTraits<Kaskade::Tensor<K,Size0,Sizes...>>
  {
    using field_type = typename FieldTraits<K>::field_type;
    using real_type = typename FieldTraits<K>::real_type;
  };
}
/**
 * \endcond
 */

#endif
