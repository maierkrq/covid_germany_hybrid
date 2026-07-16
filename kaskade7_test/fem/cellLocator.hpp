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

#ifndef CELLLOCATOR_HPP
#define CELLLOCATOR_HPP

#include <iterator>
#include <utility>
#include <vector>

#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/geometry/geometries/point.hpp> 
#include <boost/geometry/geometries/linestring.hpp>

#include "dune/common/fvector.hh"
#include "dune/grid/config.h"
#include "dune/grid/uggrid.hh"
 
#include "fem/cellLocator.hh"
#include "fem/fixdune.hh"
#include "fem/forEach.hh"
#include "fem/lagrangespace.hh"
#include "fem/spaces.hh"
#include "utilities/detailed_exception.hh"

namespace Kaskade
{
  namespace bg = boost::geometry;

  namespace CellLocatorDetail
  {
    // convert Dune points to bg points. Adapting Dune::FieldVector as point type in bg
    // didn't work (for unknown reason)
    bg::model::point<double,2,bg::cs::cartesian> bgPoint(Dune::FieldVector<double,2> const& p)
    {
      return bg::model::point<double,2,bg::cs::cartesian>(p[0],p[1]);
    }
    bg::model::point<double,3,bg::cs::cartesian> bgPoint(Dune::FieldVector<double,3> const& p)
    {
      return bg::model::point<double,3,bg::cs::cartesian>(p[0],p[1],p[2]);
    }
    
    // ---------------------------------------------------------------------------------------------

    template <class GridView, int dimw=GridView::Grid::dimension>
    struct SpatialIndex
    {
      static int const dimension = GridView::Grid::dimension;
      static int const dimension_world = dimw;
      using Grid = typename GridView::Grid;
      using Cell = typename GridView::template Codim<0>::Entity;
      using ctype    = typename GridView::ctype;
      using Position = Dune::FieldVector<ctype,dimension_world>;
      using Point = bg::model::point<double,dimension_world,bg::cs::cartesian>;
      using Box = bg::model::box<Point>;
      using Value = std::pair<Box,Cell>;
      using RTree = bg::index::rtree<Value,bg::index::rstar<8>>; // is 8 a good value?
      
      RTree tree;      
      
      ctype diameter() const 
      {
        auto box = tree.bounds();
        return bg::distance(box.max_corner(),box.min_corner());
      }
    };
  }
  
  // ---------------------------------------------------------------------------------------------------
  // ---------------------------------------------------------------------------------------------------

  template <class GridView, int dimw>
  CellLocator<GridView, dimw>::CellLocator(GridView const& gridView)
  : spatialIndex(new Index)
  {
    init(gridView);
  }

  // ---------------------------------------------------------------------------------------------------

  template <class GridView, int dimw>
  void CellLocator<GridView,dimw>::init(GridView const& gridView)
  {
    using namespace CellLocatorDetail;
    using Box   = typename Index::Box;
    using Value = typename Index::Value;
    
    ScopedTimingSection constructorTiming("cell locator construction");
    
    auto& timer = Timings::instance();
    
    spatialIndex->tree.clear();       // remove old displaced faces from spatial index
    
    std::vector<Value> boxes;         // std::vector of bounding box and cell

    timer.start("cell bounding boxes");
    // TODO: parallelize this using cell ranges (around 10% of constructor time)
    // enter all cells
    for (auto const& cell: Dune::elements(gridView))
    {
      auto [box_min,box_max] = boundingBox(cell);
      boxes.push_back(Value(Box(bgPoint(box_min),bgPoint(box_max)),cell));
    }
    timer.stop("cell bounding boxes");

    // insert all bounding boxes together with their cells into the spatial index
    timer.start("constructing spatial index");
    spatialIndex->tree.insert(boxes);
    diam = spatialIndex->diameter();
    timer.stop("constructing spatial index");
  }

  // ---------------------------------------------------------------------------------------------------

  template <class GridView, int dimw>
  CellLocator<GridView,dimw>::~CellLocator()
  {}
  

  // ----------------------------------------------------------------------------------------------

  template <class GridView, int dimw>
  std::pair<typename GridView::template Codim<0>::Entity,double>
  CellLocator<GridView,dimw>::closestCell(Position const& pos) const
  {
    using namespace CellLocatorDetail;
    using Point = bg::model::point<ctype,dimension,bg::cs::cartesian>;
    Point point = bgPoint(pos);

    // vector of Values(Box,Cell) for the rtree query
    std::vector<typename Index::Value> queryCells;
    
    // closest cells with distance to point
    std::vector<std::pair<Cell,double>> nearestCells;

    // query the rtree for bounding boxes containing the point
    spatialIndex->tree.query(bg::index::contains(point),std::back_inserter(queryCells));
    
    // point is NOT contained in any of the bounding boxes
    if(queryCells.empty()) 
    {
      // query the rtree for boxes closest to the point
      spatialIndex->tree.query(bg::index::nearest(point,3),std::back_inserter(queryCells));
    }
     
    for (auto const& value: queryCells)
    {
      auto cell = value.second;
      // checkInside computes the distance to the cell (negative, if point contained)
      double d = checkInside(cell.type(),cell.geometry().local(pos));
      nearestCells.push_back(std::make_pair(cell,d));
    }
    
    // return minimal distance cell (d=0, if point contained in the cell)
    return *std::min_element(nearestCells.cbegin(), nearestCells.cend(), 
                            [](const auto& v1, const auto& v2) { return v1.second < v2.second; });
  }
  
  // ----------------------------------------------------------------------------------------------

  template <class GridView, int dimw>
  std::vector<typename CellLocator<GridView,dimw>::Cell>
  CellLocator<GridView,dimw>::cells() const
  {
    std::vector<Cell> cs;
    cs.reserve(spatialIndex->tree.size());

    std::transform(spatialIndex->tree.begin(),spatialIndex->tree.end(),back_inserter(cs),
                   [&](auto const& value) { return value.second; });
    return cs;
  }

}

#endif
