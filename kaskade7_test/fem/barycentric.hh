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

#ifndef BARYCENTRIC_HH
#define BARYCENTRIC_HH

#include <cassert>
#include <numeric>   // std::accumulate, with clang++

#include "dune/common/fvector.hh"
#include "dune/common/fmatrix.hh"

namespace Kaskade
{
  /**
   * \ingroup grid
   * \brief Computes the barycentric coordinates of a point in the unit simplex.
   * 
   * The barycentric coordinates of \f$ x \f$ are just the cartesian ones with \f$ 1-\|x\|_{\infty}\f$ appended.
   * This corresponds to the cartesian origin being the vertex with index dim.
   */
  template <class CoordType, int dim>
  Dune::FieldVector<CoordType,dim+1> barycentric(Dune::FieldVector<CoordType,dim> const& x)
  {
    Dune::FieldVector<CoordType,dim+1> zeta;
    zeta[dim] = 1;
    for (int i=0; i<dim; ++i)
    {
      zeta[i] = x[i];
      zeta[dim] -= x[i];
    }

    // Shape functions vanishing on the boundary will probably rely on
    // one barycentric coordinate being exactly zero, in order not to be
    // affected by large penalty values from Dirichlet boundary
    // conditions. Hence we here enforce exactly zero values for very
    // small barcyentric coordinates.
    constexpr CoordType eps = std::numeric_limits<CoordType>::epsilon();
    for (int i=0; i<=dim; ++i)
      if (std::abs(zeta[i]) < 100*eps) // well, that's probably *exactly* on the boundary
        zeta[i] = 0;

    return zeta;
  }

  // ----------------------------------------------------------------------------------------------

  /**
   * \ingroup linalgbasic
   * \brief Projects a vector onto the unit simplex.
   */
  template <class CoordType, int dim,
            std::enable_if_t<std::is_floating_point_v<CoordType>,int> = 0>
  Dune::FieldVector<CoordType,dim> projectToUnitSimplex(Dune::FieldVector<CoordType,dim> x)
  {
    constexpr CoordType eps = std::numeric_limits<CoordType>::epsilon();

    auto b = barycentric(x);    // This is already modified such that b>=0.
    b /= (b.one_norm()+eps);    // Now sum(b) = 1

    std::copy(b.begin(),b.begin()+dim,x.begin());
    return x;
  }

  // ----------------------------------------------------------------------------------------------

  /**
   * \ingroup grid
   * \brief Computes barycentric coordinates of a point in the unit simplex.
   * 
   * Here, the barycentric coordinates of \f$ x \f$ are \f$ (1-\sum_{i=0}^d xi,x_0,\dots,x_{d-1}) \f$.
   * This corresponds to vertex zero being the cartesian origin.
   *
   * Points extremely close to the unit simplex boundary (less than \f$ 100 \epsilon \f$ distance)
   * are snapped to the boundary. This is intended to prevent Dirichlet penalty factors
   * with huge penalty factors from affecting the computation due to roundoff.
   *
   * \see barycentric2Derivative
   */
  template <class CoordType, int dim>
  Dune::FieldVector<CoordType,dim+1> barycentric2(Dune::FieldVector<CoordType,dim> const& x)
  {
    Dune::FieldVector<CoordType,dim+1> zeta;
    zeta[0] = 1;
    for (int i=0; i<dim; ++i) {
      zeta[i+1] = x[i];
      zeta[0] -= x[i];
    }

    // Shape functions vanishing on the boundary will probably rely on
    // one barycentric coordinate being exactly zero, in order not to be
    // affected by large penalty values from Dirichlet boundary
    // conditions. Hence we here enforce exactly zero values for very
    // small barcyentric coordinates.
    constexpr CoordType eps = std::numeric_limits<CoordType>::epsilon();
    for (int i=0; i<=dim; ++i)
      if (std::abs(zeta[i]) < 100*eps) // well, that's probably *exactly* on the boundary
        zeta[i] = 0;

    return zeta;
  }

  /**
   * \ingroup grid
   * \brief Computes the derivative of the barycentric coordinates of a point in the unit simplex.
   *
   * Here, the barycentric coordinates of \f$ x \f$ are \f$ (1-\sum_{i=0}^d xi,x_0,\dots,x_{d-1}) \f$.
   * This corresponds to vertex zero being the cartesian origin.
   *
   * \see barycentric2
   */
  template <class CoordType, int dim>
  Dune::FieldMatrix<CoordType,dim+1,dim> barycentric2Derivative(Dune::FieldVector<CoordType,dim> const& x)
  {
    Dune::FieldMatrix<CoordType,dim+1,dim> B(0.0);
    for (int i=0; i<dim; ++i)
    {
      B[0][i] = -1.0;
      B[i+1][i] = 1.0;
    }
    return B;
  }
  
  // ----------------------------------------------------------------------------------------------

  
  /**
   * \ingroup fem
   * \brief Converts integer coordinate index to barycentric indices.
   * \tparam dim the spatial dimension
   * \param x the coordinate index
   * \param bsum the barycentric sum
   * 
   * The entries of x need to be between 0 and bsum, and their sum 
   * shall not exceed bsum. Then the returned index contains in the 
   * first dim entries just x, and in the last one bsum-sum.
   */
  template <size_t dim>
  std::array<int,dim+1> barycentric(std::array<int,dim> const& x, int bsum)
  {
    std::array<int,dim+1> zeta;
    std::copy(begin(x),end(x),begin(zeta));
    zeta[dim] = bsum - std::accumulate(begin(x),end(x),0);
    return zeta;
  }
  
} /* end of namespace Kaskade */

#endif
