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

#if !defined(FIXDUNE_HH)
#define FIXDUNE_HH

#include <type_traits>

#include "dune/common/fmatrix.hh"
#include "dune/common/fvector.hh"
#include "dune/common/version.hh"
#include "dune/geometry/type.hh"
#include "dune/geometry/referenceelements.hh"

// Here we'd rather include dune/common/config.h, but this messes around with
// including "FC.h" without specifying the dune/... path. Fortunately,
// the grid config header has the common module version defined as well.
#include "dune/grid/config.h"

#include "dune/istl/matrix.hh"
#include "dune/istl/bvector.hh"
#include "dune/istl/bcrsmatrix.hh"

#include "utilities/deprecation.hh"
#include "utilities/scalar.hh"
/**
 * \file fixdune.hh
 * \brief This file contains various utility functions that augment the basic functionality of Dune.
 */


namespace Dune
{

  /**
   * \ingroup linalgbasic
   * \brief  Scalar-vector multiplication \f$ (s,x) \mapsto sx \f$
   */
  template <class T, int n, class S,
            class enable = typename std::enable_if_t<std::is_arithmetic_v<S>>>
  Dune::FieldVector<T,n> operator*(Dune::FieldVector<T,n> x, S s)
  {
    x *= static_cast<T>(s);
    return x;
  }
  
  /**
   * \ingroup linalgbasic
   * \brief  Scalar-vector multiplication \f$ (x,s) \mapsto sx \f$
   */
  template <class T, int n, class S,
            class enable = typename std::enable_if_t<std::is_arithmetic_v<S>>>
  Dune::FieldVector<T,n> operator*(S s, Dune::FieldVector<T,n> x)
  {
    x *= static_cast<typename Kaskade::EntryTraits<T>::field_type>(s);
    return x;
  }
  
  /**
   * \ingroup linalgbasic
   * \brief Vector addition \f$ (x,y) \mapsto x+y \f$
   */
  template <class T, int n>
  Dune::FieldVector<T,n> operator+(Dune::FieldVector<T,n> x, Dune::FieldVector<T,n> const& y)
  {
    x += y;
    return x;
  }
  
  /**
   * \ingroup linalgbasic
   * \brief Vector addition \f$ (x,y) \mapsto x+y \f$
   */
  template <class T, int n>
  Dune::FieldVector<T,n> operator+(Dune::FieldVector<T,n> x, T const& y)
  {
    for (int i=0; i<n; ++i)
      x[i] += y;
    return x;
  }
  
  
  /**
   * \ingroup linalgbasic
   * \brief Vector negation \f$ x \mapsto - x \f$
   */
  template <class T, int n>
  Dune::FieldVector<T,n> operator-(Dune::FieldVector<T,n> x)
  {
    for (int i=0; i<n; ++i)
      x[i] = -x[i];
    return x;
  }
  
  /**
   * \ingroup linalgbasic
   * \brief Componentwise maximum.
   */
  template <class T, int n>
  Dune::FieldVector<T,n> max(Dune::FieldVector<T,n> x, Dune::FieldVector<T,n> const& y)
  {
    for (int i=0; i<n; ++i)
      x[i] = std::max(x[i],y[i]);
    return x;
  }
  
  /**
   * \ingroup linalgbasic
   * \brief Componentwise minimum.
   */
  template <class T, int n>
  Dune::FieldVector<T,n> min(Dune::FieldVector<T,n> x, Dune::FieldVector<T,n> const& y)
  {
    for (int i=0; i<n; ++i)
      x[i] = std::min(x[i],y[i]);
    return x;
  }
  
  
  /**
   * \ingroup linalgbasic
   * \brief  Converts a matrix of size n x m to a vector of size n*m by concatenating all the columns.
   */
  template <class T, int n, int m>
  Dune::FieldVector<T,n*m> asVector(Dune::FieldMatrix<T,n,m> const& x)
  {
    Dune::FieldVector<T,n*m> y;
    for (int j=0; j<m; ++j)
      for (int i=0; i<n; ++i)
        y[j*n+i] = x[i][j];
    return y;
  }
  
  /**
   * \ingroup linalgbasic
   * \brief Converts a vector of size nm to a matrix of size n x m by filling columns successively.
   */
  template <int n, int m, class T>
  Dune::FieldMatrix<T,n,m> asMatrix(Dune::FieldVector<T,n*m> const& x)
  {
    Dune::FieldMatrix<T,n,m> y;
    for (int j=0; j<m; ++j)
      for (int i=0; i<n; ++i)
        y[i][j] = x[j*n+i];
    return y;
  }
  
  
  /**
   * \ingroup linalgbasic
   * \brief outer vector product \f$ (x,y) \mapsto x y^T \f$.
   */
  template <class T, int n, int m>
  Dune::FieldMatrix<T,n,m> outerProduct(Dune::FieldVector<T,n> const& x, Dune::FieldVector<T,m> const& y)
  {
    Dune::FieldMatrix<T,n,m> A;
    for (int i=0; i<n; ++i)
      for (int j=0; j<m; ++j)
        A[i][j] = x[i]*y[j];
    return A;
  }
  
  /**
   * \ingroup linalgbasic
   * \brief returns the normalized vector
   */
  template <class T, int n> 
  Dune::FieldVector<T,n> normalize(Dune::FieldVector<T,n> x)
  {
    T nm = x.two_norm();
    if (nm>0)
      x /= nm;
    return x;
  }

  /**
   * \ingroup linalgbasic
   * \brief Concatenation of vectors.
   */
  template <class T, int n, int m>
  Dune::FieldVector<T,n+m> vertcat(Dune::FieldVector<T,n> const& x,
                                   Dune::FieldVector<T,m> const& y)
  {
    Dune::FieldVector<T,n+m> z;
    for (int i=0; i<n; ++i)
      z[i] = x[i];
    for (int i=0; i<m; ++i)
      z[i+n] = y[i];
    return z;
  }
  
  
#ifndef DUNE_ALBERTA_ALGEBRA_HH
  /**
   * \ingroup linalgbasic
   * \brief vector product \f$ (x,y) \mapsto x \times y \in \mathbb{R} \f$.
   */
  template <class T>
  T vectorProduct(Dune::FieldVector<T,2> const& x, Dune::FieldVector<T,2> const& y)
  {
    return x[0]*y[1] - x[1]*y[0];
  }
  
  /**
   * \ingroup linalgbasic
   * \brief the matrix \f$ A(x) \in \mathbb{R}^{1\times 2}\f$ satisfying \f$ A(x) b = x\times b \f$ for all \f$ b \f$
   */
  template <class T>
  Dune::FieldMatrix<T,1,2> vectorProductMatrix(Dune::FieldVector<T,2> const& x)
  {
    Dune::FieldMatrix<T,1,2> A;
    A[0][0] = -x[1];
    A[0][1] = x[0];
  }

  /**
   * \ingroup linalgbasic
   * \brief vector product \f$ (x,y) \mapsto x \times y \f$.
   */
  template <class T>
  Dune::FieldVector<T,3> vectorProduct(Dune::FieldVector<T,3> const& x, Dune::FieldVector<T,3> const& y)
  {
    Dune::FieldVector<T,3> t;
    for (int i=0; i<3; ++i)
      t[i] = x[(i+1)%3]*y[(i+2)%3]-x[(i+2)%3]*y[(i+1)%3];
    return t;
  }
  
  /**
   * \ingroup linalgbasic
   * \brief the skew-symmetric matrix \f$ A(x)\in \mathbb{R}^{3\times 3} \f$ satisfying \f$ A(x) b = x\times b \f$ for all \f$ b \f$
   */
  template <class T>
  Dune::FieldMatrix<T,3,3> vectorProductMatrix(Dune::FieldVector<T,3> const& x)
  {
    Dune::FieldMatrix<T,3,3> A(0.0);
    for (int i=0; i<3; ++i)
    {
      A[i][(i+2)%3] = x[(i+1)%3];
      A[i][(i+1)%3] = -x[(i+2)%3];
    }
    return A;
  }
  
#endif
  


  /**
   * \ingroup linalgbasic
   * \brief Matrix-vector multiplication \f$ (A,x) \mapsto Ax \f$
   */
  template <class T, class S, int n, int m>
  Dune::FieldVector<S,n> operator*(Dune::FieldMatrix<T,n,m> const& A,
                                   Dune::FieldVector<S,m> const& x)
  {
    Dune::FieldVector<S,n> b(0);
    //A.umv(x,b);
    for (int i=0; i<n; ++i)
      for (int j=0; j<m; ++j)
        b[i] += A[i][j]*x[j];
    return b;
  }
  
  /**
   * \ingroup linalgbasic
   * \brief Matrix negation \f$ A \mapsto - A \f$
   */
  template <class T, int n, int m>
  Dune::FieldMatrix<T,n,m> operator-(Dune::FieldMatrix<T,n,m> const& A)
  {
    Dune::FieldMatrix<T,n,m> C;
    for (int i=0; i<n; ++i)
      C[i] = -A[i]; // vector negation of rows
    return C;
  }

  template <class T, int n, int m>
  auto normalForm(Dune::FieldMatrix<T,n,m> const& A)
  {
    if constexpr (n*m==1)
      return A[0][0];
    else
      return A;
  }


  template <class T, typename = std::enable_if_t<std::is_arithmetic<T>::value>>
  T transpose(T x)
  {
    return x;
  }

  /**
   * \ingroup linalgbasic
   * \brief Matrix transposition \f$ A \mapsto A^T \f$
   */
  template <class T, int n, int m>
  Dune::FieldMatrix<T,m,n> transpose(Dune::FieldMatrix<T,n,m> const& A)
  {
    Dune::FieldMatrix<T,m,n> At;
    for (int i=0; i<n; ++i)
      for (int j=0; j<m; ++j)
        At[j][i] = transpose(A[i][j]);
    return At;
  }
} // end of namespace Dune

template <class Scalar, int n>
auto horzcat(Dune::FieldVector<Scalar,n> const& c1)
{
  return c1;
}

/**
 * \ingroup linalgbasic
 * \brief Concatenates an arbitrary number of column vectors and matrices with coinciding number of rows to a larger matrix.
 *
 * This works as Matlab's [x,y,z] construction.
 */
template <class Scalar, int n, class ... Rest>
auto horzcat(Dune::FieldVector<Scalar,n> const& a, Rest ... rest)
{
  // Works, but I think it's way too slow unless the compiler
  // does some magic...

  using Kaskade::rank;

  auto b = horzcat(rest...);
  using B = decltype(b);
  static_assert(Kaskade::EntryTraits<B>::rows == n*Kaskade::EntryTraits<Scalar>::rows);

// leads to compiler error as of GCC 8.2 - why?
//   static_assert(rank<B>==1 || rank<B>==2);
  if constexpr (rank<B> == 1)
  {
    Dune::FieldMatrix<Scalar,n,2> result;
    for (int i=0; i<n; ++i)
    {
      result[i][0] = a[i];
      result[i][1] = b[i];
    }
    return result;
  }
  else if constexpr (rank<B> == 2)
  {
    Dune::FieldMatrix<Scalar,n,1+B::cols> result;
    for (int i=0; i<n; ++i)
    {
      result[i][0] = a[i];
      for (int j=0; j<B::cols; ++j)
        result[i][1+j] = b[i][j];
    }
    return result;
  }
  else
  {
    // never get here
    return 0;
  }
}

template <class Scalar, int n, int m>
auto horzcat(Dune::FieldMatrix<Scalar,n,m> const& a)
{
  return a;
}

template <class Scalar, int n, int m, class ... Rest>
auto horzcat(Dune::FieldMatrix<Scalar,n,m> const& a, Rest ... rest)
{
  // Works, but I think it's way too slow unless the compiler
  // does some magic...

  using Kaskade::rank;

  auto b = horzcat(rest...);
  using B = decltype(b);
  static_assert(Kaskade::EntryTraits<B>::rows == n*Kaskade::EntryTraits<Scalar>::rows);

  if constexpr (rank<B> == 1)
  {
    Dune::FieldMatrix<Scalar,n,m+1> result;
    for (int i=0; i<n; ++i)
    {
      for (int j=0; j<m; ++j)
        result[i][j] = a[i][j];
      result[i][m] = b[i];
    }
    return result;
  }
  else if constexpr (rank<B> == 2)
  {
    Dune::FieldMatrix<Scalar,n,m+B::cols> result;
    for (int i=0; i<n; ++i)
    {
      for (int j=0; j<m; ++j)
        result[i][j] = a[i][j];
      for (int j=0; j<B::cols; ++j)
        result[i][m+j] = b[i][j];
    }
    return result;
  }
  else
  {
    // never get here
    return 0;
  }
}


// TODO: generic version via parameter packs?
template <class Scalar, int n>
Dune::FieldMatrix<Scalar,n,3> horzcat(Dune::FieldVector<Scalar,n> const& col1,
                                      Dune::FieldVector<Scalar,n> const& col2,
                                      Dune::FieldVector<Scalar,n> const& col3)
{
  Dune::FieldMatrix<Scalar,n,3> A;
  for (int i=0; i<n; ++i)
  {
    A[i][0] = col1[i]; A[i][1] = col2[i]; A[i][2] = col3[i];
  }
  return A;
}

/**
 * \brief Extracts the j-th column of the given matrix A.
 */
template <class Scalar, int n, int m>
Dune::FieldVector<Scalar,n> column(Dune::FieldMatrix<Scalar,n,m> const& A, int j)
{
  assert(0<=j && j<m);
  Dune::FieldVector<Scalar,n> x;
  for (int i=0; i<n; ++i)
    x[i] = A[i][j];
  return x;
}

/**
 * \ingroup linalgbasic
 * \brief matrix entries as vector, concatenating columns
 */
template <class Scalar, int n>
Dune::FieldVector<Scalar,n*n> vectorize(Dune::FieldMatrix<Scalar,n,n> const& A)
{
  Dune::FieldVector<Scalar,n*n> result(0.);
  
  for(int i=0; i<n; ++i)
    for(int j=0; j<n; ++j)
      result[i+j*n] = A[i][j];
    
  return result;
}

/**
 * \ingroup linalgbasic
 * \brief forming quadratic matrix columns of vector segments
 */
template <int n, class Scalar>
Dune::FieldMatrix<Scalar,n,n> unvectorize(Dune::FieldVector<Scalar,n*n> const& v)
{
  Dune::FieldMatrix<Scalar,n,n> result(0.);
  
  for(int i=0; i<n; ++i)
    for(int j=0; j<n; ++j)
      result[i][j] = v[i+j*n];
    
  return result;
}




/**
 * \ingroup linalgbasic
 * \brief Matrix contraction (Frobenius product) \f$ (A,B) \mapsto A:B = \sum_{i,j} A_{ij} B_{ij} \f$
 */
template <class T, int n, int m>
T contraction(Dune::FieldMatrix<T,n,m> const& A, Dune::FieldMatrix<T,n,m> const& B) 
{
  T x = 0;
  for (int i=0; i<n; ++i)
    for (int j=0; j<m; ++j)
      x += A[i][j] * B[i][j];
  return x;
}

/**
 * \ingroup linalgbasic
 * \brief Matrix trace \f$ A \mapsto \mathrm{tr}\,A = \sum_{i} A_{ii} \f$
 *
 * Note: Instead of trace(A*B) use the equivalent, but faster contraction(A,B).
 */
template <class T, int n>
T trace(Dune::FieldMatrix<T,n,n> const& A) 
{
  T x = 0;
  for (int i=0; i<n; ++i)
    x += A[i][i];
  return x;
}

/**
 * \ingroup linalgbasic
 * \brief Matrix determinant \f$ A \mapsto \mathrm{det}\,A  \f$
 * \warning DEPRECATED, use A.determinant() directly, as provided by Dune.
 */
template <class T, int n>
T determinant(Dune::FieldMatrix<T,n,n> const& A) 
{
  static Kaskade::Deprecated dummy("detrminant(A) is deprecated, use A.determinant() instead",2020);
  return A.determinant();
}

/**
 * \ingroup linalgbasic
 * \brief Matrix diagonal as a vector.
 */
template <class T, int n>
Dune::FieldVector<T,n> diag(Dune::FieldMatrix<T,n,n> const& A)
{
  Dune::FieldVector<T,n> d;
  for (int i=0; i<n; ++i)
    d[i] = A[i][i];
  return d;
}

/**
 * \ingroup linalgbasic
 * \brief Returns the identity matrix of size n \f$ \mapsto I \f$
 */
template <class T, int n>
Dune::FieldMatrix<T,n,n> unitMatrix() 
{
  Dune::FieldMatrix<T,n,n> I(0);
  for (int i=0; i<n; ++i)
    I[i][i] = 1;
  return I;
}

/**
 * \ingroup linalgbasic
 * \brief Returns the unit vector \f$ e_i \f$ with \f$ (e_i)_k = \delta_{ik} \f$.
 */
template <class T, int n>
Dune::FieldVector<T,n> unitVector(int i)
{
  assert(0<=i && i<n);
  Dune::FieldVector<T,n> ei(0.0);
  ei[i] = 1;
  return ei;
}

/**
 * \brief Returns the matrix type or its transposed type, depending on the given transpose flag.
 */
template <class Matrix, bool transposed>
struct Transpose {};

template <class Scalar, int n, int m, bool transposed>
struct Transpose<Dune::FieldMatrix<Scalar,n,m>,transposed>
{
  typedef Dune::FieldMatrix<Scalar,transposed? m:n,transposed? n:m> type;
};


/**
 * \ingroup linalgbasic
 * \brief Computes Z = X*Y.
 * X and Y need to be compatible, i.e. X.M()==Y.N()
 * has to hold. No aliasing may occur, i.e. &Z != &X and &Z != &Y must
 * hold. Z is reshaped as needed.
 */
template <class MatrixZ, class Matrix>
void MatMult(MatrixZ& z, Matrix const& x, Matrix const& y) 
{
  assert(((void*)&z != (void*)&x) && ((void*)&z != (void*)&y));
  assert(x.M() == y.N());
  
  z.setSize(x.N(),y.M());
  
  for (int i=0; i<z.N(); i++ ) {
      for (int j=0; j<z.M(); j++ ) {
          typename MatrixZ::block_type temp = 0.0;
          for (int k=0; k<x.M(); k++)
            temp += x[i][k]*y[k][j];
          z[i][j] = temp;
        }
    }
}


//---------------------------------------------------------------------

/**
 * \cond internals
 *
 * This namespace contains helper functions for implementing Hierarchic mappers.
 */
namespace FixDuneDetail
{

  //TODO: in Dune 2.0 this is not needed anymore!
  template <int Codim>
  struct GetIndexOfSubEntity
  {
    template <class Cell, class IS>
    static int value(Cell const& cell, int codim, int entity, IS const& is)
    {
      if (codim==Codim)
        //       return is.template subIndex<Codim>(cell,entity);
        return is.subIndex(cell,entity,Codim);
      else
        return GetIndexOfSubEntity<Codim-1>::value(cell,codim,entity,is);
    }
  };

  template <>
  struct GetIndexOfSubEntity<0>
  {
    template <class Cell, class IS>
    static int value(Cell const& cell, int codim, int entity, IS const& is)
    {
      assert (codim==0);
      //     return is.template subIndex<0>(cell,entity);
      return is.subIndex(cell,entity,0);
    }
  };

  template <int Codim>
  struct GetNumberOfSubEntities
  {
    template <class Cell>
    static int value( Cell const& cell, int codim )
    {
      if (codim==Codim)
        return cell.template count<Codim>() ;
      else
        return GetNumberOfSubEntities<Codim-1>::value(cell, codim) ;
    }
  } ;

  template <>
  struct GetNumberOfSubEntities<0>
  {
    template <class Cell>
    static int value( Cell const& cell, int codim )
    {
      assert (codim==0);
      return cell.template count<0>();
    }
  };

} // End of namespace FixDuneDetail
/// \endcond ---------------------------------------------------------------------

/**
 * Computes the index of the k-th subentity of codimension c in the
 * index set is.  This is available in Dune only with c given
 * statically.
 *
 * TODO: use IndexSet::IndexType instead of size_t for Dune 1.1
 */
template <class IndexSet, class Cell>
size_t subIndex(IndexSet const& is, Cell const& cell, int codim, int subentity) 
{
  return FixDuneDetail::GetIndexOfSubEntity<Cell::dimension>::value(cell,codim,subentity,is);
}

/**
 *   \brief Computes the number of subentities of codimension c of a given entity.
 *
 *  This is available in Dune only with c given statically.
 *
 */
template <class Cell>
int count(Cell const& cell, int codim)
{
  return FixDuneDetail::GetNumberOfSubEntities<Cell::dimension>::value(cell, codim) ;
}


//  ---------------------------------------------------------------------

/**
 * \ingroup utilities
 * \brief Checks whether a point is inside or outside the reference element, and how much
 *
 * Computational geometry predicates are notoriously difficult to get right. In Dune, it may happen that
 * a point inside a cell is reported not to be contained in any of its children. Often this is due to
 * rounding errors. This function returns not just a binary decision subject to rounding errors, but
 * a floating point value that gives the magnitude of being inside or outside, such that the result can
 * be tested by arbitrary tolerances or sorted according to magnitude.
 *
 * \tparam LocalCoordinate the type of local coordinate vectors. Usually a Dune::FieldVector<ctype,dim>.
 *
 * \todo The current implementation is just for cubes and simplices. Pyramids and prisms are not covered.
 *
 * \return a value less or equal to zero if the point is inside, a value greater zero if outside
 */
template <class LocalCoordinate>
typename LocalCoordinate::field_type checkInside(Dune::GeometryType const& gt, LocalCoordinate const& x)
{
  if (gt.isSimplex())
    // The reference simplex is the set defined by componentwise x_i>=0 and 1-sum x_i>=0. Thus we return
    // the negative of the minimum of these conditions.
    return -std::min(*std::min_element(x.begin(),x.end()),1-x.one_norm());
  else if (gt.isCube())
    // The reference cube is the set defined by componentwise xi>=0 and xi_<=1.
    return -std::min(*std::min_element(x.begin(),x.end()),1-*std::max_element(x.begin(),x.end()));
  else
    // Other elements are not (yet) covered. Fallback to standard Dune implementation
    return Dune::ReferenceElements<typename LocalCoordinate::field_type,
                                   LocalCoordinate::dimension>::general(gt).checkInside(x)? 0: 1;
}

/**
 * \ingroup utilities
 * \brief Computes the bounding box of a cell
 *
 * \tparam Cell the type of cell considered. Usually a typename GridView::template Codim<0>::Entity.
 *
 * 
 * \return a pair of Dune::FieldVector<double,dim> which are the minimal and the maximal coordiantes of the bounding box.
 */
template <class Cell>
std::pair<Dune::FieldVector<double,Cell::dimension>,Dune::FieldVector<double,Cell::dimension>> boundingBox(Cell const& cell)
{
  Dune::FieldVector<double,Cell::dimension> boxmin(std::numeric_limits<double>::max()),boxmax(std::numeric_limits<double>::min());
  
  auto geo = cell.geometry();
  for(int i=0; i<geo.corners(); ++i)
  {
    auto corner = geo.corner(i);
    for(int d=0; d<Cell::dimension; ++d)
    {
      boxmin[d] = std::min(boxmin[d],corner[d]);
      boxmax[d] = std::max(boxmax[d],corner[d]);
    }
  }
  return std::make_pair(boxmin,boxmax);
}


/// copy between two block vectors of different element lengths
template <class Scalar, int n, int m, class Allocator1, class Allocator2>
void transferBlockVector(Dune::BlockVector<Dune::FieldVector<Scalar,n>,Allocator1> const& from, Dune::BlockVector<Dune::FieldVector<Scalar,m>,Allocator2>& to)
{
  to.resize(from.dim()*m);
  for(size_t i=0; i<from.N(); ++i)
    for(int j=0; j<n; ++j)
      to[(i*n+j)/m][(i*n+j)%m] = from[i][j];
}

/**
 * \ingroup IO
 * \brief Pretty prints the sparse matrix A to the given output stream.
 */
template <class Scalar, int n, int m, class Allocator>
void bcrsPrint(Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,n,m>,Allocator> const& A, std::ostream& os=std::cout)
{
  for(size_t i0=0; i0<A.N(); ++i0)
    for(size_t i1=0; i1<n; ++i1)
      {
        for(size_t j0=0; j0<A.M(); ++j0)
          for(size_t j1=0; j1<m; ++j1)
            if(A.exists(i0,j0)) os << A[i0][j0][i1][j1] << " ";
        os << std::endl;
      }
  os << std::endl;
}

/**
 * \ingroup IO
 * \brief Pretty prints the sparse matrix A to the given output stream.
 */
template <class Scalar, int n, int m, class Allocator>
std::ostream& operator<<(std::ostream& os, Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,n,m>,Allocator> const& A)
{
  bcrsPrint(A,os);
  return os;
}

#endif
