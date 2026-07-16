/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2017-2020 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CELLLOCATOR_HH
#define CELLLOCATOR_HH

#include <memory>
#include <optional>
#include <tuple>

#include "dune/common/fvector.hh"
#include "dune/geometry/type.hh"

#include "utilities/timing.hh"

namespace Kaskade
{
  /// \internal 
  // forward declaration for PIMPL idiome
  namespace CellLocatorDetail
  {
    template <class GridView, int dimw>
    class SpatialIndex;    
  }
  /// \endinternal
  
  // ----------------------------------------------------------------------------------------------



  /**
   * \ingroup fem
   * \brief Supports finding cell intersections by geometric queries.
   * \tparam GridView
   * 
   * An instance of CellLocator builds an (boost) RTree with bounding boxes of every cell of a grid and therefore can
   * efficiently compute geometric queries. In particular, this is to be used to find the cell in a grid, which contains 
   * a particular point (global coordinate).
   * 
   * Explicit instantiations for Dune::UGGrid<2> and Dune::UGGrid<3> leaf grid views are provided in
   * the file cellLocator.cpp, and are compiled into the library. For other grids views, include
   * cellLocator.hpp, which contains the implementations.
   */
  template <class GridView, int dimw=GridView::Grid::dimension>
  class CellLocator
  {
    using Index = CellLocatorDetail::SpatialIndex<GridView,dimw>;
    using Grid = typename GridView::Grid;

  public:
    static int const dimension = GridView::Grid::dimension;
    static int const dimension_world = dimw;
    using ctype    = typename GridView::ctype;
    using Position = Dune::FieldVector<ctype,dimension_world>;
    using Cell = typename GridView::template Codim<0>::Entity;

    
    /**
     * \brief Constructs a locator based on the grids view.
     */
    CellLocator(GridView const& gridView);

    // we need this because the default destructor doesn't know about the size of SpatialIndex when deallocating via unique_ptr
    ~CellLocator();

    /**
     * \brief The diagonal length of the grid's cartesian bounding box.
     * 
     */
    ctype diameter() const 
    {
      return diam;
    }
    /**
     * \brief Returns a sequence of all cells.
     *
     * The sequence of cells is extracted from the RTree and returned by value, which is a considerable computational effort.
     * The order of cells is unspecified.
     */
    std::vector<Cell> cells() const;
    
    /**
    * \brief Returns a pair consisting of the closest cell to pos and the distance to this cell.
    *       
    * If the point is contained in any of the grid cells, the distance is 0 or negative, a positive distance 
    * means that the point is not contained in any of the grid cells. The closest cell and the distance to this
    * cell are returned.
    * 
    * \return (cell,distance) tuple
    */
    std::pair<Cell,double> closestCell(Position const& pos) const;

  private:
    std::unique_ptr<CellLocatorDetail::SpatialIndex<GridView,dimension_world>> spatialIndex;
    ctype diam;

    void init(GridView const& gridView);
    
  };
}

#endif
