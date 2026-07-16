/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2013-2022 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*
 *  Created on: 25.06.2013
 *      Author: bzflubko
 */

#ifndef KASKADE_GRID_GENERATION_HH_
#define KASKADE_GRID_GENERATION_HH_

#include <boost/math/constants/constants.hpp>

#include <algorithm>
#include <cmath>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/common/fvector.hh>
#include "dune/grid/uggrid.hh"

#include "fem/fixdune.hh"
#include "GridGeneration/cube.hh"
#include "utilities/detailed_exception.hh"
#include "utilities/timing.hh"

namespace Kaskade
{
  /**
   * \cond internals
   */
  namespace GridGeneration_Detail
  {
    /**
     * \brief Extract vector of vertices and tetrahedra's vertex ids.
     *
     * Taking a list of cubes, which is each represented as union of a set of tetrahedra, 
     * this function collects all the local vertices, removes duplicates,
     * and stores them in the first entry of the result-pair. The second contains the 
     * connectivity of the tetrahedra wrt. the returned global vertex list.
     *
     * \param cubes list of cubes
     * \return pair of vectors, vertices in the first, connectivity in the second entry
     */
    template <class Type,
               int dim = std::is_same_v<Type,Cube<typename Type::Real>>? 3: 2>
    std::pair<std::vector<Dune::FieldVector<typename Type::Real,dim>>, std::vector<std::vector<unsigned int>>> extractSimplices(std::vector<Type> const& cubes)
    {
      using Vec = Dune::FieldVector<typename Type::Real,dim>;

      // A lexical ordering on the vertices (with tolerance).
      typename Type::Real tol = 1e-6;
      auto vLess = [&](Vec const& x, Vec const& y)
      {
        for (int i=0; i<dim; ++i)
          if (x[i] < y[i]-tol)
            return true;
          else if (x[i] > y[i]+tol)
            return false;

        return false;                   // equal, i.e. all components within tolerance, and thus not less
      };

      auto vEq = [&](Vec const& x, Vec const& y)
      {
        return !(vLess(x,y) || vLess(y,x));
      };

      // First we create a list of unique vertices, eliminating any duplicates,
      // such that we get unique and stable vertex ids.
      std::vector<Vec> vertices;
      for (auto const& cube: cubes)
      {
        auto const& localVertices = cube.getVertices();
        vertices.insert(vertices.end(),localVertices.begin(),localVertices.end());
      }
      std::sort(vertices.begin(),vertices.end(),vLess);
      vertices.erase(std::unique(vertices.begin(),vertices.end(),vEq),vertices.end());



      std::vector<std::vector<unsigned int>> tetrahedra; // TODO: use std::array here, the number of corners is fixed

      for (Type const& cube: cubes)
      {
        auto const& localVertices = cube.getVertices();
        auto const& localTetrahedra = cube.getSimplices();
        size_t const nVertices = localVertices.size();

        // maps local vertex ids to global ones (i.e. the index in 'vertices')
        std::vector<size_t> vertexIdMap(nVertices);
        for(size_t i=0; i<nVertices; ++i) 
        {
          auto it = std::lower_bound(vertices.begin(),vertices.end(),localVertices[i],vLess);
          assert(it != vertices.end());
          vertexIdMap[i] = it - vertices.begin();
        }

        // map local connectivity onto global one (wrt. numbering in vector 'vertices')
        std::vector<unsigned int> newTet;
        for(std::vector<unsigned int> const& tet: localTetrahedra)
        {
          newTet.resize(tet.size());
          for(size_t i=0; i<newTet.size(); ++i)
            newTet[i] = vertexIdMap[tet[i]];
          tetrahedra.push_back(newTet);
        }
      }

      return std::make_pair(vertices,tetrahedra);
    }


    template <class Grid, bool isUGGrid>
    struct SetDefaultHeapSizeImpl
    {
      static void apply(size_t){}
    };

    template <class Grid>
    struct SetDefaultHeapSizeImpl<Grid,true>
    {
      static void apply(size_t heapSize)
      {
        Grid::setDefaultHeapSize(heapSize);
      }
    };


    /**
     * \brief Fill rectangle with squares.
     *
     * \param x0 left, lower, front corner of the cuboid
     * \param dx side lengths
     * \param dh side length of cubes that are used to fill the cuboid
     * \param symmetric true: each cube provides a symmetric tetrahedal decomposition
     */
    template <class Scalar>
    std::vector<Square<Scalar> > fillRectangle(Dune::FieldVector<Scalar,2> const& x0, Dune::FieldVector<Scalar,2> const& dx, Dune::FieldVector<Scalar,2> const& dh, bool symmetric=true)
    {
      std::vector<Square<Scalar> > squares;
      Dune::FieldVector<Scalar,2> dh2(dh);
      for(int i=0; i<2; ++i)
        if(dh[i]>dx[i])     // if dh larger than sidelength, set dh = dx
          dh2[i] = dx[i];
        
        
      size_t xi_end = (size_t)round(fabs(dx[0]/dh2[0])),
          yi_end = (size_t)round(fabs(dx[1]/dh2[1]));
      Dune::FieldVector<Scalar,2> dsquare(dh2), x(x0);

      for(size_t xi=0; xi<xi_end; ++xi)
      {
        x[0] = x0[0] + xi*dh2[0];
        for(size_t yi=0; yi<yi_end; ++yi)
        {
          x[1] = x0[1] + yi*dh2[1];
          squares.push_back(Square<Scalar>(x,dsquare,symmetric));
        }
      }

      return squares;
    }

    /**
     * \brief Fill rectangle with squares.
     *
     * \param x0 left, lower, front corner of the cuboid
     * \param dx side lengths
     * \param dh side length of cubes that are used to fill the cuboid
     * \param symmetric true: each cube provides a symmetric tetrahedal decomposition
     */
    template <class Scalar>
    std::vector<Square<Scalar> > fillRectangle(Dune::FieldVector<Scalar,2> const& x0, Dune::FieldVector<Scalar,2> const& dx, Scalar dh, bool symmetric=true)
    {
      return fillRectangle(x0,dx,Dune::FieldVector<Scalar,2>(dh),symmetric);
    }

    /**
     * \brief Fill cuboid with smaller cuboids.
     *
     * \param x0 left, lower, front corner of the cuboid
     * \param dx side lengths
     * \param dh side lengths of cuboids that are used to fill the cuboid
     * \param symmetric true: each cuboid provides a symmetric tetrahedal decomposition
     */
    template <class Scalar>
    std::vector<Cube<Scalar>> fillCuboid(Dune::FieldVector<Scalar,3> const& x0, Dune::FieldVector<Scalar,3> const& dx,
                                         Dune::FieldVector<Scalar,3> const& dh, bool symmetric=true)
    {
      std::vector<Cube<Scalar>> cubes;
      
      Dune::FieldVector<Scalar,3> dh2(dh);
      for(int i=0; i<3; ++i)
        if(dh[i]>dx[i])           // if dh larger than sidelength, set dh = dx
          dh2[i] = dx[i];
        
      size_t xi_end = (size_t)round(fabs(dx[0]/dh2[0])),
             yi_end = (size_t)round(fabs(dx[1]/dh2[1])),
             zi_end = (size_t)round(fabs(dx[2]/dh2[2]));
             
      if (xi_end==0 || yi_end==0 || zi_end==0)
        throw GridException("empty cubes list",__FILE__,__LINE__);
             
      Dune::FieldVector<Scalar,3> x(x0);

      for(size_t xi=0; xi<xi_end; ++xi)
      {
        x[0] = x0[0] + xi*dh2[0];
        for(size_t yi=0; yi<yi_end; ++yi)
        {
          x[1] = x0[1] + yi*dh2[1];
          for(size_t zi=0; zi<zi_end; ++zi)
          {
            x[2] = x0[2] + zi*dh2[2];
            cubes.push_back(Cube<Scalar>(x,dh2,symmetric));
//             std::cout << "insert cube at " << x << " width " << dh2 << "\n";
          }
        }
      }

      return cubes;
    }


    /**
     * \brief Fill cuboid with cubes
     *
     * \param x0 left, lower, front corner of the cuboid
     * \param dx side lengths
     * \param dh side length of cubes that are used to fill the cuboid
     * \param symmetric true: each cube provides a symmetric tetrahedal decomposition
     */
    template <class Scalar>
    std::vector<Cube<Scalar>> fillCuboid(Dune::FieldVector<Scalar,3> const& x0, Dune::FieldVector<Scalar,3> const& dx, 
                                         Scalar dh, bool symmetric=true)
    {
      return fillCuboid(x0,dx,Dune::FieldVector<Scalar,3>(dh),symmetric);
    }



    template<class Scalar, int dim>
    std::vector<Cube<Scalar>> fillTwoLayerCuboid(Dune::FieldVector<Scalar,dim> x0, Dune::FieldVector<Scalar,dim> const& dx1, Dune::FieldVector<Scalar,dim> const& dx2,
                                                    Dune::FieldVector<Scalar,dim> const& dh1, Scalar dh2, bool symmetric=true)
    {
      auto cubes = fillCuboid(x0,dx1,dh1,symmetric);

      x0[dim-1] += dx1[dim-1];
      auto cubes2 = fillCuboid(x0,dx2,dh2,symmetric);
      for( Cube<Scalar> const& cube : cubes2 ) cubes.push_back(cube);

      return cubes;
    }

    template<class Scalar, int dim>
    std::vector<Cube<Scalar> > fillTwoLayerCuboid(Dune::FieldVector<Scalar,dim> x0, Dune::FieldVector<Scalar,dim> const& dx1, Dune::FieldVector<Scalar,dim> const& dx2,
                                                    Dune::FieldVector<Scalar,dim> const& dh1, Dune::FieldVector<Scalar,dim> const& dh2, bool symmetric=true)
    {
      auto cubes = fillCuboid(x0,dx1,dh1,symmetric);

      x0[dim-1] += dx1[dim-1];
      auto cubes2 = fillCuboid(x0,dx2,dh2,symmetric);
      for( Cube<Scalar> const& cube : cubes2 ) cubes.push_back(cube);

      return cubes;
    }

    template <int dim>
    struct FillCuboid;

    template <>
    struct FillCuboid<3>
    {
      template <class Scalar, class EdgeLength>
      static std::vector<Cube<Scalar> > apply(Dune::FieldVector<Scalar,3> const& x0, Dune::FieldVector<Scalar,3> const& dx, 
                                              EdgeLength const& dh, bool symmetric=true)
      {
        return fillCuboid(x0, dx, dh, symmetric);
      }
    };

    template <>
    struct FillCuboid<2>
    {
      template <class Scalar, class EdgeLength>
      static std::vector<Square<Scalar> > apply(Dune::FieldVector<Scalar,2> const& x0, Dune::FieldVector<Scalar,2> const& dx, EdgeLength const& dh, bool symmetric=true)
      {
        return fillRectangle(x0, dx, dh, symmetric);
      }
    };


    /**
     * \brief Fill L-shaped domain with cubes
     *
     * \param x0 left, lower, front corner of the L-shaped domain
     * \param dx length of horizontal line
     * \param dy length of vertical line
     * \param thickness
     * \param dh side length of cubes that are used to fill the cuboid
     * \param symmetric true: each cube provides a symmetric tetrahedal decomposition
     */
    template <class Scalar, int dim>
    std::vector< typename std::conditional<dim==3,Cube<Scalar>,Square<Scalar> >::type > 
    fillLShape(Dune::FieldVector<Scalar,dim> x0, Scalar dx, Scalar dy, Dune::FieldVector<Scalar,dim> const& thickness, Scalar dh, bool symmetric=true)
    {
      // fill horizontal line
      Dune::FieldVector<Scalar,dim> d0(thickness);    d0[0] = dx;
      auto cubes = FillCuboid<dim>::apply(x0,d0,dh,symmetric);

      // fill vertical line
      d0 = thickness;   d0[1] = dy;
      x0[1] += thickness[1];
      auto cubes2 = FillCuboid<dim>::apply(x0,d0,dh,symmetric);

      // merge
      for(auto const& cube : cubes2)
        cubes.push_back(cube);

      return cubes;
    }

    /**
     * \brief Fill L-shaped domain with cubes
     *
     * \param x0 center on lower boundary of the T-shaped domain
     * \param dx length of horizontal line
     * \param dy length of vertical line
     * \param thickness
     * \param dh side length of cubes that are used to fill the cuboid
     * \param symmetric true: each cube provides a symmetric tetrahedal decomposition
     */
    template <class Scalar, int dim>
    std::vector< typename std::conditional<dim==3,Cube<Scalar>,Square<Scalar>>::type> 
    fillTShape(Dune::FieldVector<Scalar,dim> x0, Scalar dx, Scalar dy, Dune::FieldVector<Scalar,dim> const& thickness, Scalar dh, bool symmetric=true)
    {
      // fill vertical line
      x0[0] -= 0.5*thickness[0];
      Dune::FieldVector<Scalar,dim> d0(thickness);    d0[1] = dy - thickness[1];
      auto cubes = FillCuboid<dim>::apply(x0,d0,dh,symmetric);

      // fill vertical line
      x0[0] += 0.5*thickness[0] - 0.5*dx;
      x0[1] += d0[1];
      d0[1] = thickness[1];   d0[0] = dx;
      auto cubes2 = FillCuboid<dim>::apply(x0,d0,dh,symmetric);

      // merge
      for(auto const& cube : cubes2) cubes.push_back(cube);

      return cubes;
    }


//#ifdef HAVE_UG
    template <class Grid>
    using SetDefaultHeapSize = SetDefaultHeapSizeImpl<Grid,std::is_same<Grid,Dune::UGGrid<Grid::dimension> >::value>;
/*#else
    template <class Grid>
    using SetDefaultHeapSize = SetDefaultHeapSizeImpl<Grid,false>;
#endif*/
  } // end of namespace GridGeneration_Detail
  /**
   * \endcond
   */


  /**
   * \ingroup gridInput
   * \brief Extract simplicial grid from list of cubes.
   * 
   * From a given list of (conforming) cubes or rectangles, this creates a simplicial grid.
   * 
   * \tparam Grid the type of grid to create 
   * \tparam Type
   * 
   * Usage: 
   * \code
   * grid = createGrid<Grid>(GridGeneration_Detail::fillCuboid(x0,dx,dh,true));
   * \endcode
   */
  template <class Grid, class Type>
  std::unique_ptr<Grid> createGrid(std::vector<Type> const& cubes, size_t heapSize=0)
  {
    Timings& timer = Timings::instance();
    timer.start("extract simplices");
    auto simplices = GridGeneration_Detail::extractSimplices(cubes);
    timer.stop("extract simplices");
    if (heapSize > 0)
      GridGeneration_Detail::SetDefaultHeapSize<Grid>::apply(heapSize);
    Dune::GridFactory<Grid> factory;
    Dune::GeometryType gt(Dune::GeometryType::simplex,Type::dim);

    timer.start("insert vertex");
    for (Dune::FieldVector<double,Type::dim> const& vertex: simplices.first)
      factory.insertVertex(vertex);
    timer.stop("insert vertex");
    timer.start("insert element");
    for (std::vector<unsigned int> const& elem: simplices.second)
      factory.insertElement(gt,elem);
    timer.stop("insert element");

    return std::unique_ptr<Grid>(factory.createGrid());
  }

  /**
   * \ingroup gridInput
   * \brief Creates a uniform simplicial grid on a rectangular or cuboid domain.
   * 
   * \tparam Grid type of grid to create
   * \tparam Scalar
   * \tparam dim spatial dimension (usually 2 or 3)
   * \tparam EdgeLength either a scalar or a Dune::FieldVector type
   *
   * \param x0 left, lower, front corner of the cuboid
   * \param dx desired side lengths of the hexahedron to fill
   * \param dh side length(s) of hexahedra that are used to fill the cuboid
   * \param symmetric true: each cube provides a symmetric tetrahedal decomposition
   * 
   * If the edge length dh is a scalar, exact cubes are created. If it is a vector (of size dim), 
   * hexahedra with the given edge length are created.
   * 
   * Note that the desired size dx of the domain to be filled is just an indication - the resulting
   * grid will approximately have this size. In contrast, the hexahedra used for filling the domain 
   * will always have the shape given by dh.
   */
  template <class Grid, class Scalar, int dim, class EdgeLength>
  std::unique_ptr<Grid> createCuboid(Dune::FieldVector<Scalar,dim> const& x0, Dune::FieldVector<Scalar,dim> const& dx, 
                                     EdgeLength const& dh, bool symmetric=true, size_t heapSize=0)
  {
    Timings& timer = Timings::instance();
    timer.start("create cubes");
    auto const& cubes = GridGeneration_Detail::FillCuboid<dim>::apply(x0,dx,dh,symmetric);
    timer.stop("create cubes");
    return createGrid<Grid>(cubes,heapSize);
  }

  /**
   * \ingroup gridInput
   * \brief Creates a unit square, filled regularly.
   * \param dh side length of the subsquares.
   * \param symmetric if true, each subsquare is filled symmetrically with four triangles
   */
  template <class Grid>
  std::unique_ptr<Grid> createUnitSquare(double dh = 1., bool symmetric=true, size_t heapSize=0)
  {
    using Vec = Dune::FieldVector<double,2>;
    return createCuboid<Grid>(Vec(0.), Vec(1.), dh, symmetric, heapSize);
  }

  /**
   * \ingroup gridInput
   * \brief Creates a unit cube.
   * \param dh desired side length of the subcubes.
   */
  template <class Grid>
  std::unique_ptr<Grid> createUnitCube(double dh = 1., bool symmetric=true, size_t heapSize=0)
  {
    typedef Dune::FieldVector<double,3> Vec;
    return createCuboid<Grid>(Vec(0.), Vec(1.), dh, symmetric, heapSize);
  }


  /**
   * \ingroup gridInput
   * \brief Fill rectangle with squares.
   *
   * \param x0 left, lower, front corner of the cuboid
   * \param dx side lengths
   * \param dh side length of cubes that are used to fill the cuboid
   * \param symmetric true: each cube provides a symmetric tetrahedal decomposition
   */
  template <class Grid, class Scalar>
  std::unique_ptr<Grid> createRectangle(Dune::FieldVector<Scalar,2> const& x0, Dune::FieldVector<Scalar,2> const& dx, Scalar dh, bool symmetric=true, size_t heapSize=0)
  {
    return createGrid<Grid>(GridGeneration_Detail::fillRectangle(x0,dx,dh,symmetric),heapSize);
  }

  /**
   * \ingroup gridInput
   * \brief Fill L-shaped domain with cubes
   *
   * \param x0 left, lower, front corner of the L-shaped domain
   * \param dx length of horizontal line
   * \param dy length of vertical line
   * \param thickness
   * \param dh side length of cubes that are used to fill the cuboid
   * \param symmetric true: each cube provides a symmetric tetrahedal decomposition
   */
  template <class Grid, class Scalar, int dim>
  std::unique_ptr<Grid> createLShape(Dune::FieldVector<Scalar,dim> const& x0, Scalar dx, Scalar dy,
                                     Dune::FieldVector<Scalar,dim> const& thickness, Scalar dh,
                                     bool symmetric=true, size_t heapSize = 0)
  {
    return createGrid<Grid>(GridGeneration_Detail::fillLShape(x0,dx,dy,thickness,dh,symmetric),heapSize);
  }

  /**
   * \ingroup gridInput
   * \brief Fill L-shaped domain with cubes
   *
   * \param x0 center on lower boundary of the T-shaped domain
   * \param dx length of horizontal line
   * \param dy length of vertical line
   * \param thickness
   * \param dh side length of cubes that are used to fill the cuboid
   * \param symmetric true: each cube provides a symmetric tetrahedal decomposition
   */
  template <class Grid, class Scalar, int dim>
  std::unique_ptr<Grid> createTShape(Dune::FieldVector<Scalar,dim> const& x0, Scalar dx, Scalar dy, Dune::FieldVector<Scalar,dim> const& thickness, Scalar dh, bool symmetric=true, size_t heapSize = 0)
  {
    return createGrid<Grid>(GridGeneration_Detail::fillTShape(x0,dx,dy,thickness,dh,symmetric),heapSize);
  }

  /**
   * \ingroup gridInput
   */
  template <class Grid, class Scalar, int dim>
  std::unique_ptr<Grid> createTwoLayerCuboid(Dune::FieldVector<Scalar,dim> const& x0, Dune::FieldVector<Scalar,dim> const& dx1, Dune::FieldVector<Scalar,dim> const& dx2,
                                             Dune::FieldVector<Scalar,dim> const& dh1, Scalar dh2, bool symmetric=true, size_t heapSize = 0)
  {
    return createGrid<Grid>(GridGeneration_Detail::fillTwoLayerCuboid(x0,dx1,dx2,dh1,dh2,symmetric),heapSize);
  }

  /**
   * \ingroup gridInput
   */
  template <class Grid, class Scalar, int dim>
  std::unique_ptr<Grid> createTwoLayerCuboid(Dune::FieldVector<Scalar,dim> const& x0, Dune::FieldVector<Scalar,dim> const& dx1, Dune::FieldVector<Scalar,dim> const& dx2,
                                             Dune::FieldVector<Scalar,dim> const& dh1, Dune::FieldVector<Scalar,dim> const& dh2, bool symmetric=true, size_t heapSize = 0)
  {
    return createGrid<Grid>(GridGeneration_Detail::fillTwoLayerCuboid(x0,dx1,dx2,dh1,dh2,symmetric),heapSize);
  }
  
  // --------------------------------------------------------------------------------------------------------------------
  // --------------------------------------------------------------------------------------------------------------------
  
  /**
   * \ingroup gridInput
   * \defgroup platonian Platonian solids
   * Functions for creating basic grids for Platonian solids.
   */
  
  /**
   * \ingroup platonian
   * \brief Creates a hexahedron consisting of five tetrahedra.
   * 
   * The hexahedron is then the unit cube \f$ [0,1]^3 \f$.
   * 
   * \tparam Grid the type of grid to create 
   */
  template <class Grid>
  std::unique_ptr<Grid> createHexahedron(size_t heapSize=0)
  {
    if(heapSize > 0) 
      GridGeneration_Detail::SetDefaultHeapSize<Grid>::apply(heapSize);
    Dune::GridFactory<Grid> factory;
    Dune::GeometryType gt(Dune::GeometryType::simplex,3);

    for (int i=0; i<2; ++i)
      for (int j=0; j<2; ++j)
        for (int k=0; k<2; ++k)
        {
          Dune::FieldVector<double,3> x;
          x[0] = i;
          x[1] = j; 
          x[2] = k;
          factory.insertVertex(x);
        }
    
    factory.insertElement(gt,std::vector<unsigned int>{7,1,4,2});
    factory.insertElement(gt,std::vector<unsigned int>{0,4,1,2});
    factory.insertElement(gt,std::vector<unsigned int>{5,1,4,7});
    factory.insertElement(gt,std::vector<unsigned int>{6,2,4,7});
    factory.insertElement(gt,std::vector<unsigned int>{3,1,7,2});

    return std::unique_ptr<Grid>(factory.createGrid());
  }

  /**
   * \ingroup platonian
   * \brief Creates an octahedron consisting of eight tetrahedra.
   * 
   * The octahedron is centered at the origin,
   * with the vertices located at given positions on the cartesian axes.
   * 
   * \tparam Grid the type of grid to create 
   */
  template <class Grid>
  std::unique_ptr<Grid> createOctahedron(Dune::FieldVector<typename Grid::ctype,3> extent = Dune::FieldVector<typename Grid::ctype,3>(1.0),
                                         size_t heapSize=0)
  {
    if(heapSize > 0) 
      GridGeneration_Detail::SetDefaultHeapSize<Grid>::apply(heapSize);
    Dune::GridFactory<Grid> factory;
    Dune::GeometryType gt(Dune::GeometryType::simplex,3);

    for (int i=0; i<3; ++i)
    {
      Dune::FieldVector<double,3> x;
      x[i] = extent[i];
      factory.insertVertex(x);
      factory.insertVertex(-x);
    }
    factory.insertVertex(Dune::FieldVector<double,3>());
    
    factory.insertElement(gt,std::vector<unsigned int>{0,2,4,6});
    factory.insertElement(gt,std::vector<unsigned int>{0,3,4,6});
    factory.insertElement(gt,std::vector<unsigned int>{0,2,5,6});
    factory.insertElement(gt,std::vector<unsigned int>{0,3,5,6});
    factory.insertElement(gt,std::vector<unsigned int>{1,2,4,6});
    factory.insertElement(gt,std::vector<unsigned int>{1,3,4,6});
    factory.insertElement(gt,std::vector<unsigned int>{1,2,5,6});
    factory.insertElement(gt,std::vector<unsigned int>{1,3,5,6});

    return std::unique_ptr<Grid>(factory.createGrid());
  }

  /**
   * \ingroup platonian
   * \brief Creates an icosahedron consisting of 20 tetrahedra.
   * 
   * The icosahedron is centered around the origin and has a diameter of two.
   * 
   * \tparam Grid the type of grid to create 
   */
  template <class Grid>
  std::unique_ptr<Grid> createIcosahedron(size_t heapSize=0)
  {
    if(heapSize > 0) 
      GridGeneration_Detail::SetDefaultHeapSize<Grid>::apply(heapSize);
    Dune::GridFactory<Grid> factory;
    Dune::GeometryType gt(Dune::GeometryType::simplex,3);
    
    auto vector = [=](double x, double y, double z) 
    { 
      Dune::FieldVector<double,3> p;
      p[0] = x; p[1] = y; p[2] = z;
      return p;
    };

    double s = (1+std::sqrt(5))/2;
    
    factory.insertVertex(vector( 0, 1, s));    // 0
    factory.insertVertex(vector( 0, 1,-s));    // 1
    factory.insertVertex(vector( 0,-1, s));    // 2
    factory.insertVertex(vector( 0,-1,-s));    // 3
    factory.insertVertex(vector( 1, s, 0));    // 4
    factory.insertVertex(vector( 1,-s, 0));    // 5
    factory.insertVertex(vector(-1, s, 0));    // 6
    factory.insertVertex(vector(-1,-s, 0));    // 7
    factory.insertVertex(vector( s, 0, 1));    // 8
    factory.insertVertex(vector( s, 0,-1));    // 9
    factory.insertVertex(vector(-s, 0, 1));    // 10
    factory.insertVertex(vector(-s, 0,-1));    // 11
    factory.insertVertex(vector( 0, 0, 0));    // 12
    
    factory.insertElement(gt,std::vector<unsigned int>{0,2, 8,12});
    factory.insertElement(gt,std::vector<unsigned int>{0,2,10,12});
    factory.insertElement(gt,std::vector<unsigned int>{1,3, 9,12});
    factory.insertElement(gt,std::vector<unsigned int>{1,3,11,12});
    
    factory.insertElement(gt,std::vector<unsigned int>{4,6, 0,12});
    factory.insertElement(gt,std::vector<unsigned int>{4,6, 1,12});
    factory.insertElement(gt,std::vector<unsigned int>{5,7, 2,12});
    factory.insertElement(gt,std::vector<unsigned int>{5,7, 3,12});
    
    factory.insertElement(gt,std::vector<unsigned int>{8,9, 4,12});
    factory.insertElement(gt,std::vector<unsigned int>{8,9, 5,12});
    factory.insertElement(gt,std::vector<unsigned int>{10,11,6,12});
    factory.insertElement(gt,std::vector<unsigned int>{10,11,7,12});
    
    factory.insertElement(gt,std::vector<unsigned int>{0,4, 8,12});
    factory.insertElement(gt,std::vector<unsigned int>{0,6,10,12});
    factory.insertElement(gt,std::vector<unsigned int>{3,7,11,12});
    factory.insertElement(gt,std::vector<unsigned int>{3,5,9,12});
    factory.insertElement(gt,std::vector<unsigned int>{1,4,9,12});
    factory.insertElement(gt,std::vector<unsigned int>{1,6,11,12});
    factory.insertElement(gt,std::vector<unsigned int>{2,7,10,12});
    factory.insertElement(gt,std::vector<unsigned int>{2,5,8,12});

    return std::unique_ptr<Grid>(factory.createGrid());
  }
  
  //-----------------------------------------------------------------------------------------------
  
  /**
   * \ingroup gridInput
   * \brief Creates a cylinder.
   * 
   * This creates a grid for the cyclinder domain \f$ \{ (x,y,z) \mid x^2 + y^2 \le r, 0\le z \le H \} \f$,
   * i.e. the cylinder's axis is the z-axis. The grid consists of a series of prisms along the cylinder's axis.
   * 
   * \param radius the cylinder's radius \f$ r \f$
   * \param height the cylinder's height \f$ H \f$
   * \param layers the number of anti-prsimatic disks making up the grid
   * \param circumEdges the number of edges (or vertices) around the circumference
   * \param 
   * 
   * Cell quality (in terms of aspect ratio) depends on the parameter values. Reasonable cell quality
   * is achieved for values approximately satisfying radius*layers = height and circumEdges = 6.
   */
  template <class Grid>
  std::unique_ptr<Grid> createCylinder(double radius, double height, int layers, int circumEdges=6)
  {
    return createCylinder<Grid>(radius,height,layers,circumEdges,[](auto const& x) { return x; });
  }
  
  /**
   * \ingroup gridInput
   * \brief Creates a deformed cylinder.
   * 
   * This creates a grid for the cyclinder domain \f$ \{ (x,y,z) \mid x^2 + y^2 \le r, 0\le z \le H \} \f$,
   * and deforms it as specified. The grid consists of a series of prisms along the cylinder's axis.
   * 
   * \param radius the cylinder's radius \f$ r \f$
   * \param height the cylinder's height \f$ H \f$
   * \param layers the number of anti-prsimatic disks making up the grid
   * \param circumEdges the number of edges (or vertices) around the circumference
   * \param deformation a callable object `Dune::FieldVector<double,3> deformation(Dune::FieldVector<double,3>)`
   * 
   * Cell quality (in terms of aspect ratio) depends on the parameter values (and the deformation, of course). 
   * Reasonable cell quality is achieved for values approximately satisfying radius*layers = height and circumEdges = 6.
   */
  template <class Grid, class Deformation>
  std::unique_ptr<Grid> createCylinder(double radius, double height, int layers, int circumEdges, 
                                       Deformation const& deformation)
  {
    using boost::math::double_constants::pi;
    
    Dune::GridFactory<Grid> factory;
    Dune::GeometryType gt(Dune::GeometryType::simplex,3);
    
    auto vector = [](double x, double y, double z) 
    { 
      Dune::FieldVector<double,3> p;
      p[0] = x; p[1] = y; p[2] = z;
      return p;
    };
    
    double h = height/layers;
    // Create vertices. These are circumEdges+1 on each of the layers+1 circular disks
    for (int i=0; i<=layers; ++i)
    {
      double z = i*h;
      factory.insertVertex(deformation(vector(0,0,z)));
      for (int j=0; j<circumEdges; ++j)
      {
//         double angle = 2*pi*(j-i/6.0)/circumEdges;
        double angle = 2*pi*j/circumEdges;
        auto v = vector(radius*std::cos(angle),radius*std::sin(angle),z);
        factory.insertVertex(deformation(v));
      }
    }
    
    // Create the cells. Each prismatic disc consists of circumEdges prisms
    // connecting the central axis with a surface quadrilateral. These prisms are split
    // into three tetrahedra.
    for (int i=0; i<layers; ++i)
    {
      // The offset of the vertex indices: circular disks below current anti-prismatic disk 
      // times the number of vertices in each circular disk.
      int n = circumEdges+1;
      int idxOffset = i*n;
      
      unsigned int c0 = idxOffset;    // lower center vertex
      unsigned int c1 = idxOffset+n;  // top center vertex
      
      for (int j=0; j<circumEdges; ++j)
      {
        unsigned int a0 = c0+1+j;                 // bottom left corner of boundary quadrilateral
        unsigned int b0 = c0+1+(j+1)%circumEdges; // bottom right corner
        unsigned int a1 = c1+1+j;                 // top left corner
        unsigned int b1 = c1+1+(j+1)%circumEdges; // top right corner
        factory.insertElement(gt,std::vector<unsigned int>{a0,b0,c0,c1});
        factory.insertElement(gt,std::vector<unsigned int>{a0,b1,a1,c1});
        factory.insertElement(gt,std::vector<unsigned int>{a0,b0,b1,c1});
      }
    }
    
    return std::unique_ptr<Grid>(factory.createGrid());
  }
  
  //-----------------------------------------------------------------------------------------------

  /**
   * \ingroup gridInput
   * \brief Creates an regular pentagon consisting of five triangles, centered at the origin and with radius r.
   * 
   * \param radius of the circumcircle of the pentagon
   * 
   * \tparam Grid the type of grid to create 
   * 
   */
  template <class Grid>
  std::unique_ptr<Grid> createPentagon(double radius=1, size_t heapSize=0)
  {
    
    if(heapSize > 0) 
      GridGeneration_Detail::SetDefaultHeapSize<Grid>::apply(heapSize);
    Dune::GridFactory<Grid> factory;
    Dune::GeometryType gt(Dune::GeometryType::simplex,2);
    
    auto vector = [=](double x, double y) 
    { 
      Dune::FieldVector<double,2> p;
      p[0] = x; p[1] = y;
      return p;
    };
    
    double c1, c2, s1, s2;
    c1 = 0.25*(std::sqrt(5)-1);
    c2 = 0.25*(std::sqrt(5)+1);
    s1 = 0.25*(std::sqrt(10+2*std::sqrt(5)));
    s2 = 0.25*(std::sqrt(10-2*std::sqrt(5)));
    
    factory.insertVertex(radius*vector( 1., 0.));          // 0
    factory.insertVertex(radius*vector( c1, s1));          // 1
    factory.insertVertex(radius*vector(-1.*c2,s2));        // 2
    factory.insertVertex(radius*vector(-1.*c2,-1.*s2));    // 3
    factory.insertVertex(radius*vector( c1, -1.0*s1));     // 4
    factory.insertVertex(radius*vector( 0., 0.));          // 5

    
    factory.insertElement(gt,std::vector<unsigned int>{0,1,5});
    factory.insertElement(gt,std::vector<unsigned int>{1,2,5});
    factory.insertElement(gt,std::vector<unsigned int>{2,3,5});
    factory.insertElement(gt,std::vector<unsigned int>{3,4,5});
    factory.insertElement(gt,std::vector<unsigned int>{4,0,5});
   
    return std::unique_ptr<Grid>(factory.createGrid());
  }
  
  
  /**
   * \ingroup gridInput
   * \brief Creates an U-shaped domain.
   * 
   *  =======================|
   * ||
   * ||
   * ||
   *  =======================|
   * 
   * \tparam Grid the type of grid to create 
   * 
   * \param l   side length
   * \param d   thickness of gemetry
   * \param eps slit width of geometry
   * 
   */
  template <class Grid>
  std::unique_ptr<Grid> createUshape(double l, double d, double eps, size_t heapSize=0)
  {
    assert(eps > 0);
    assert(d > 0);
    assert(l > d);
    
    if(heapSize > 0) 
      GridGeneration_Detail::SetDefaultHeapSize<Grid>::apply(heapSize);
    Dune::GridFactory<Grid> factory;
    Dune::GeometryType gt(Dune::GeometryType::simplex,2);
    
    auto vector = [=](double x, double y) 
    { 
      Dune::FieldVector<double,2> p;
      p[0] = x; p[1] = y;
      return p;
    };
    
    factory.insertVertex(vector( 0, 0));              // 0
    factory.insertVertex(vector( (l+d)/2, 0));        // 1
    factory.insertVertex(vector( l, 0));              // 2
    factory.insertVertex(vector( l, d));              // 3
    factory.insertVertex(vector( (l+d)/2, d));        // 4
    factory.insertVertex(vector( d, d));              // 5
    factory.insertVertex(vector( d, d+eps));          // 6
    factory.insertVertex(vector( (l+d)/2, d+eps));    // 7
    factory.insertVertex(vector( l, d+eps));          // 8
    factory.insertVertex(vector( l, 2*d+eps));        // 9
    factory.insertVertex(vector( (l+d)/2, 2*d+eps));  // 10
    factory.insertVertex(vector( 0, 2*d+eps));        // 11
    factory.insertVertex(vector( 0, d+eps/2));        // 12
    
    factory.insertElement(gt,std::vector<unsigned int>{0,1,5});
    factory.insertElement(gt,std::vector<unsigned int>{1,4,5});
    factory.insertElement(gt,std::vector<unsigned int>{1,2,4});
    factory.insertElement(gt,std::vector<unsigned int>{2,3,4});
    factory.insertElement(gt,std::vector<unsigned int>{5,6,12});
    factory.insertElement(gt,std::vector<unsigned int>{6,7,10});
    factory.insertElement(gt,std::vector<unsigned int>{7,8,9});
    factory.insertElement(gt,std::vector<unsigned int>{7,9,10});
    factory.insertElement(gt,std::vector<unsigned int>{6,10,11});
    factory.insertElement(gt,std::vector<unsigned int>{6,11,12});
    factory.insertElement(gt,std::vector<unsigned int>{0,5,12});

    return std::unique_ptr<Grid>(factory.createGrid());
  }

} // end of namespace Kaskade

#endif /* GRIDGENERATION_HH_ */
