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

#ifndef BOUNDARYLOCATOR_HH
#define BOUNDARYLOCATOR_HH

#include <memory>
#include <optional>
#include <tuple>

#include "dune/common/fvector.hh"
#include "dune/geometry/type.hh"

#include "fem/boundaryFace.hh"
#include "utilities/timing.hh"

namespace Kaskade
{
  /// \internal 
  // forward declaration for PIMPL idiome
  namespace BoundaryLocatorDetail
  {
    template <class GridView, class Function, int dimw>
    class SpatialIndex;    
    
    // Dummy implementation for FE function, just to be syntactically correct
    // TODO: replace this by void type and use constexpr if.
    struct ZeroDisplacement
    {
      template <class Cell>
      int order(Cell const&) const { return 0; }

      template <class Cell, class B>
      B value(Cell, B b) const { return B(0.0); }

      template <class Cell, class B>
      Dune::FieldMatrix<typename B::field_type,B::dimension,B::dimension> derivative(Cell, B) const { return 0.0; }
    };

  }
  /// \endinternal
  
  // ----------------------------------------------------------------------------------------------



  /**
   * \ingroup fem
   * \brief Supports finding boundary intersections by geometric queries.
   * \tparam GridView
   * \tparam Function a vector valued function defined on the grid view. It has to support a subset of the
   *                  FunctionSpaceElement interface, namely the methods `order(cell)`, `value(cell,xi)`,
   *                  and `derivative(cell,pos)`.
   * 
   * A displacement function defined on the given grid view can be provided. In that case, instead
   * of the original grid the deformed grid, defined by the displacement, is considered.
   * 
   * Design rationale: Dune::GeometryGrid provides similar functionality, but incurs the (small) overhead of
   *                   defining a complete grid structe and is restricted to piecewise linear deformations.
   * 
   * Explicit instantiations for Dune::UGGrid<2> and Dune::UGGrid<3> leaf grid views are provided in
   * the file boundaryLocator.cpp, and are compiled into the library. For other grids views, include
   * boundaryLocator.hpp, which contains the implementations.
   */
  template <class GridView, class Function=BoundaryLocatorDetail::ZeroDisplacement, int dimw=GridView::Grid::dimension>
  class BoundaryLocator
  {
    using Index = BoundaryLocatorDetail::SpatialIndex<GridView,Function,dimw>;
    using Grid = typename GridView::Grid;

  public:
    static int const dimension = GridView::Grid::dimension;
    static int const dimension_world = dimw;
    using ctype    = typename GridView::ctype;
    using Position = Dune::FieldVector<ctype,dimension_world>;
    using LocalPosition = Dune::FieldVector<ctype,dimension_world-1>;    
    using Face = std::conditional_t<dimension==dimension_world, 
                                    typename GridView::Intersection,
                                    typename GridView::template Codim<0>::Entity>;
    using DisplacedFace = BoundaryFace<Grid,Face,Function,dimw>;

    
    /**
     * \brief Constructs a locator, optionally providing a displacement.
     */
    BoundaryLocator(GridView const& gridView);

    /**
     * \brief Constructor.
     *
     * \param boundarySegments states which boundary faces (by their boundary segment index) are to be considered;
     * if empty all are considered
     */
    BoundaryLocator(GridView const& gridView, std::vector<bool> const& boundarySegments);

    /**
     * \brief Constructs a locator, constructing the displacement function with the given argument.
     *
     * \tparam Args the argument type for the displacement function. A constructor `Function(Args)` must exist.
     */
    template <class... Args>
    BoundaryLocator(GridView const& gridView, Args&&... args)
    : BoundaryLocator(gridView)
    {
      f = std::make_unique<Function>(std::forward<Args>(args)...);
    }

    // we need this because the default destructor doesn't know about the size of SpatialIndex when deallocating via unique_ptr
    ~BoundaryLocator();
    
    /**
     * \brief Creates or updates the displacement function.
     *
     * \tparam Update the argument type for constructing or updating the displacement. `Function` must have both a
     *                constructor `Function(Update)` and an assignment operator `operator=(Update)`.
     * 
     * Calling this method recreates the spatial search tree, which can lead to a new tree structure. This affects
     * the order of faces in the sequence returned by \ref faces().
     */
    template <class Update>
    void updateDisplacement(Update const& update)
    {
      // Extract the faces from the spatial index. This is much faster than running through the whole
      // grid (O(n^(1-1/d)) vs O(n)). Moreover, faces on periodic boundaries and on not considered
      // boundary segments are already filtered out.
      auto& timer = Timings::instance();
      timer.start("boundary faces from rtree");
      auto fs = faces();
      timer.stop("boundary faces from rtree");

      if (f)
        *f = update;              // register new displacement
      else
        f = std::make_unique<Function>(update);

      enterFaces(fs);                   // recreate the spatial tree with newly displaced faces
    }


    /**
     * \brief Finds the first boundary face that is intersected by the given line segment from outside towards inside.
     * 
     * "First" means that the intersection is closest to the start point of the line.
     * 
     * \param from the origin of the line segment
     * \param to the end point of the line segment (shall not coincide with from)
     * \param minAngle the minimum angle for which an intersection is accepted. Only faces with unit outer
     *                 normal \f$ n \f$ and \f$ n^T (\mathrm{to}-\mathrm{from}) / \| \mathrm{to}-\mathrm{from} \|
     *                 \le \mathrm{minAngle} \f$ are considered.
     * 
     * An (optional) triple (i,c,q) is returned. If valid, i contains the boundary face and
     * c is the face-local coordinate of the intersection point. q is the scalar product between
     * normalized ray direction and intersecting face unit outer normal.
     */
    std::optional<std::tuple<typename BoundaryLocator<GridView,Function,dimw>::DisplacedFace,
                             typename BoundaryLocator<GridView,Function,dimw>::LocalPosition,
                             typename GridView::ctype>>
    byLineSegment(Position const& from, Position const& to, ctype minAngle=-0.5) const;
    
    /**
     * \brief Finds the first boundary face that is intersected by the given ray from outside towards inside.
     * 
     * "First" means that the intersection is closest to the start point of the ray.
     * 
     * \param from the origin of the line segment
     * \param direction the direction of the ray (shall not be zero)
     * \param minAngle the minimum angle for which an intersection is accepted. Only faces with unit outer
     *                 normal \f$ n \f$ and \f$ n^T \mathrm{dir} / \| \mathrm{dir} \|
     *                 \le \mathrm{minAngle} \f$ are considered.
     *
     * An (optional) triple (i,c,q) is returned. If valid, i contains the boundary face and
     * c is the face-local coordinate of the intersection point. q is the scalar product between
     * normalized ray direction and intersecting face unit outer normal.
     */
    std::optional<std::tuple<typename BoundaryLocator<GridView,Function,dimw>::DisplacedFace,
                             typename BoundaryLocator<GridView,Function,dimw>::LocalPosition,
                             typename GridView::ctype>>
    byRay(Position const& from, Position const& direction, ctype minAngle=-0.5) const
    {
      return byLineSegment(from,from+(diam/direction.two_norm())*direction, minAngle);
    }
    
    /**
     * \brief The diagonal length of the grid's cartesian bounding box.
     * 
     * The diameter of the geometry comes handy if intersections with a ray instead of a line segment
     * is desired. Taking an end point on the ray in diameter distance to the origin defines a line
     * segment with all possible intersections.
     */
    ctype diameter() const 
    {
      return diam;
    }
    
    /**
     * \brief Returns the global position of a boundary point.
     * 
     * This respects a potentially given nontrivial displacement of the boundary.
     */
    Position global(Face const& face, LocalPosition const& xi) const;
    
    /**
     * \brief DEPRECATED use BoundaryFace::unitOuterNormal() instead
     */
    Position unitOuterNormal(Face const& face, LocalPosition const& xi) const;
    
    /**
     * \brief Returns the current displacement of the mesh.
     * 
     * The returned pointer can be null, which signals a trivial zero displacement.
     */
    Function const* displacement() const
    {
      return f.get();
    }
    
    /**
     * \brief Returns a sequence of all faces.
     *
     * The sequence of faces is extracted and returned by value, which is a considerable computational effort.
     * The order of faces is unspecified, and subject to change on calling \ref updateDisplacement().
     */
    std::vector<Face> faces() const;

    /**
     * \brief Returns a sequence of all deformed faces.
     *
     * The sequence of faces is extracted and returned by value, which is a considerable computational effort.
     * The order of faces is unspecified, and subject to change on calling \ref updateDisplacement().
     */
    std::vector<DisplacedFace> displacedFaces() const;

  private:
    std::unique_ptr<BoundaryLocatorDetail::SpatialIndex<GridView,Function,dimension_world>> spatialIndex;
    std::unique_ptr<Function> f;
    ctype diam;
    
    using IntersectionSet = std::vector<std::tuple<BoundaryFace<Grid,Face,Function>,
                                                   LocalPosition,ctype,ctype>>;

    void init(GridView const& gridView, std::vector<bool> const& boundarySegments);
    
    /**
     * \brief Returns a sorted sequence of faces intersecting the given line segment.
     * 
     * \return (face,localpos,t,in/out) tuple
     */
    IntersectionSet getIntersectingFaces(Position const& from, Position const& to, ctype minAngle) const;
    
    /**
     * \brief Recreates the spatial index with given faces.
     *
     * The faces are displaced by the displacement function (if given).
     */
    void enterFaces(std::vector<Face> const& fs);
  };
}

#endif
