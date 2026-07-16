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

#ifndef BOUNDARYLOCATOR_HPP
#define BOUNDARYLOCATOR_HPP

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
 
#include "fem/boundaryLocator.hh"
#include "fem/fixdune.hh"
#include "fem/forEach.hh"
#include "fem/lagrangespace.hh"
#include "fem/spaces.hh"
#include "utilities/detailed_exception.hh"

namespace Kaskade
{
  namespace bg = boost::geometry;

  namespace BoundaryLocatorDetail
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

    template <class GridView, class Function, int dimw=GridView::Grid::dimension>
    struct SpatialIndex
    {
      static int const dimension = GridView::Grid::dimension;
      static int const dimension_world = dimw;
      using Grid = typename GridView::Grid;
      using Face = std::conditional_t<dimension==dimension_world, 
                                    typename GridView::Intersection,
                                    typename GridView::template Codim<0>::Entity>;
      using DisplacedFace = BoundaryFace<Grid,Face,Function,dimw>;
      using ctype    = typename GridView::ctype;
      using Position = Dune::FieldVector<ctype,dimension_world>;
      using Point = bg::model::point<double,dimension_world,bg::cs::cartesian>;
      using Box = bg::model::box<Point>;
      using Value = std::pair<Box,DisplacedFace>;
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

  template <class GridView, class Function, int dimw>
  BoundaryLocator<GridView,Function, dimw>::BoundaryLocator(GridView const& gridView)
  : spatialIndex(new Index)
  {
    init(gridView,std::vector<bool>());
  }

  // ---------------------------------------------------------------------------------------------------

  template <class GridView, class Function, int dimw>
  BoundaryLocator<GridView,Function, dimw>::BoundaryLocator(GridView const& gridView,
                                                            std::vector<bool> const& boundarySegments)
  : spatialIndex(new Index)
  {
    init(gridView,boundarySegments);
  }

  // ---------------------------------------------------------------------------------------------------

  template <class GridView, class Function, int dimw>
  void BoundaryLocator<GridView,Function, dimw>::init(GridView const& gridView, std::vector<bool> const& boundarySegments)
  {
    using namespace BoundaryLocatorDetail;
    
    
    ScopedTimingSection constructorTiming("boundary locator construction");

    if(!boundarySegments.empty() && (boundarySegments.size() != gridView.grid().numBoundarySegments()))
      throw Kaskade::DetailedException("Boundary segments vector has wrong length.", __FILE__, __LINE__);
    
    auto& timer = Timings::instance();
    timer.start("boundary faces");

    // TODO: parallelize this using cell ranges (around 10% of constructor time)
    std::vector<Face> fs;
    if constexpr(dimension==dimension_world)
    {
      forEachBoundaryFace(gridView,[&](auto const& face)
      {
        if (face.neighbor())   // on periodic boundary
          return;             // - ignore
        if(!boundarySegments.empty() && !boundarySegments[face.boundarySegmentIndex()])
          return;
          
        fs.push_back(face);
      });
    }
    else
    {
      // nothing to skip here, enter all cells
      forEachCell(gridView,[&](auto const& face)
      {          
        fs.push_back(face);
      });
    }
    timer.stop("boundary faces");

    enterFaces(fs);
  }

  // ---------------------------------------------------------------------------------------------------

  template <class GridView, class Function, int dimw>
  BoundaryLocator<GridView,Function, dimw>::~BoundaryLocator()
  {}

  // ----------------------------------------------------------------------------------------------

  template <class GridView, class Function, int dimw>
  void BoundaryLocator<GridView,Function, dimw>::enterFaces(std::vector<Face> const& fs)
  {
    using namespace BoundaryLocatorDetail;
    using Box   = typename Index::Box;
    using Value = typename Index::Value;

    auto& timer = Timings::instance();

    spatialIndex->tree.clear();       // remove old displaced faces from spatial index

    std::vector<Value> boxes;
    
    timer.start("boundary face bounding boxes");
    // TODO: parallelize this. Amounts to 50-60% of construction time.
    for (Face const& face: fs)
    {
      BoundaryFace<Grid,Face,Function,dimw> bf(face,f.get());
      auto bbox = bf.boundingBox();
      boxes.push_back(Value(Box(bgPoint(bbox.first),bgPoint(bbox.second)),bf));
    }
    timer.stop("boundary face bounding boxes");
    
    // insert all bounding boxes together with their faces into the spatial index
    timer.start("constructing spatial index");
    spatialIndex->tree.insert(boxes);
    diam = spatialIndex->diameter();
    timer.stop("constructing spatial index");
  }

  // ----------------------------------------------------------------------------------------------

  template <class GridView, class Function, int dimw>
  typename BoundaryLocator<GridView,Function, dimw>::IntersectionSet
  BoundaryLocator<GridView,Function,dimw>::getIntersectingFaces(Position const& from, Position const& to, ctype minAngle) const
  {
    using namespace BoundaryLocatorDetail;
    using Point = bg::model::point<ctype,dimension,bg::cs::cartesian>;
    
    // extract all faces from the spatial index whose bounding boxes intersect the given line segment
    bg::model::linestring<Point> segment{bgPoint(from),bgPoint(to)};
    std::vector<typename Index::Value> intersectingFaces;

    spatialIndex->tree.query(bg::index::intersects(segment),std::back_inserter(intersectingFaces));
    
    // Note that the found faces are just those whose bounding box 
    // intersects the line segment. We filter those which actually intersect the line parametrized over [0,1].
    IntersectionSet result;
    Position direction = to-from;
    auto length = direction.two_norm();
    direction /= length;


    // TODO: avoid computing useless intersection points and return only the nearest one.
    // For accomplishing that: (i) sort bounding boxes according to their intersection point
    //                         (ii) maintain the currently first encountered intersection
    //                         (iii) drop bounding box if possible intersection is *behind* current best intersection
    
    for (auto const& value: intersectingFaces)
    {
      // Find the intersection of the line segment with the tangent plane at the face center.
      // The intersection point is then translated into local coordinates.
      auto const& face = value.second;
     
      constexpr int facedimension = Index::DisplacedFace::facedimension;
      auto const& refElem = Dune::referenceElement<ctype,facedimension>(face.type());
      LocalPosition xiCenter = refElem.position(0,0);
      Position x = face.global(xiCenter);
      Position n = face.unitOuterNormal(xiCenter);

      double q = n*direction;                     // >0 if the line leaves the domain
      double tau = (n*(x-from)) / (q*length);     // intersection parameter where the face tangent plane intersects

      if (0<=tau && tau<=1                                 // consider only faces that have a chance to intersect
          && (dimension<dimension_world || q<=minAngle) )  // or (in volume meshes) have the correct orientation
      {
        auto xitau = face.intersection(from,to);

        if (xitau)                             // only if an intersection point has been found...
          result.push_back(std::make_tuple(face,xitau->first,xitau->second,q));
      }
    }
    
    return result;
  }
  
  // ----------------------------------------------------------------------------------------------

  template <class GridView, class Function, int dimw>
  std::optional<std::tuple<typename BoundaryLocator<GridView,Function,dimw>::DisplacedFace,
                           typename BoundaryLocator<GridView,Function,dimw>::LocalPosition,
                           typename GridView::ctype>>
  BoundaryLocator<GridView,Function, dimw>::byLineSegment(Position const& from,
                                                          Position const& to,
                                                          ctype minAngle) const
  {
    // extract all faces from the spatial index which intersect the given line segment and the line enters the domain
    auto intersectingFaces = getIntersectingFaces(from,to,minAngle);
    
    // Find the closest intersection. 
    std::optional<std::tuple<DisplacedFace,LocalPosition,ctype>> retval;
    
    auto it = std::min_element(begin(intersectingFaces),end(intersectingFaces),
                               [](auto const& fa, auto const& fb) { return std::get<2>(fa) < std::get<2>(fb); });
    if (it != end(intersectingFaces))
      retval = std::make_tuple(std::get<0>(*it),std::get<1>(*it),std::get<3>(*it));

    return retval;
  }
  
  // ----------------------------------------------------------------------------------------------

  template <class GridView, class Function, int dimw>
  typename BoundaryLocator<GridView,Function, dimw>::Position
  BoundaryLocator<GridView,Function,dimw>::global(Face const& face, LocalPosition const& xi) const
  {    
    if constexpr(dimension==dimension_world)
    {
      auto x = face.geometry().global(xi);
      if (auto displacement = f.get())
      {
        auto cellXi = face.geometryInInside().global(xi);
        x += displacement->value(face.inside(),cellXi);
      }
      return x;
    }
    else
    {
      auto x = face.geometry().global(xi);
      if (f.get())
        x += displacement->value(face,xi);
      return x;
    }
  }

  // ----------------------------------------------------------------------------------------------

  template <class GridView, class Function, int dimw>
  typename BoundaryLocator<GridView,Function,dimw>::Position
  BoundaryLocator<GridView,Function,dimw>::unitOuterNormal(Face const& face, LocalPosition const& xi) const
  {
    DisplacedFace df(face,f.get());
    return df.unitOuterNormal(xi);
  }
  
  // ----------------------------------------------------------------------------------------------

  template <class GridView, class Function, int dimw>
  std::vector<typename BoundaryLocator<GridView,Function,dimw>::Face>
  BoundaryLocator<GridView,Function,dimw>::faces() const
  {
    std::vector<Face> fs;
    fs.reserve(spatialIndex->tree.size());

    std::transform(spatialIndex->tree.begin(),spatialIndex->tree.end(),back_inserter(fs),
                   [&](auto const& value) { return value.second.gridFace(); });
    return fs;
  }

  template <class GridView, class Function, int dimw>
  std::vector<typename BoundaryLocator<GridView,Function,dimw>::DisplacedFace>
  BoundaryLocator<GridView,Function,dimw>::displacedFaces() const
  {
    std::vector<DisplacedFace> fs;
    fs.reserve(spatialIndex->tree.size());

    std::transform(spatialIndex->tree.begin(),spatialIndex->tree.end(),back_inserter(fs),
                   [&](auto const& value) { return value.second; });
    return fs;
  }

}

#endif
