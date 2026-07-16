/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2020-2020 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef BOUNDARYFACE_HH
#define BOUNDARYFACE_HH

#include <memory>
#include <optional>
#include <utility>

#include "dune/common/fvector.hh"
#include "dune/geometry/type.hh"


namespace Kaskade
{
  /**
   * \ingroup fem
   * \brief A class for representing displaced/deformed boundary faces.
   *
   * We distinguish three cases:
   * - boundary faces of volume grids: these are characterized by facedimension==celldimension+1
   *   and celldimension==griddimension==worlddimension.
   * - cells of a d-dimensional grid in a d+1-dimensional space: these are characterized by
   *   dimension==dimensionworld-1 and mydimension==dimension. Examples are surface meshes
   *   embedded in3D space.
   * - cells of a d-dimensional grid in a d-dimensional space that are deformed into a d+1-dimensional
   *   world: these are characterized by dimensionworld==dimension+1. Examples are planar grids
   *   deformed into 3D space, e.g., for shell or membrane mechanics.
   *
   * In the first case, Face denotes a face in the grid, i.e. a codimension 1 entity. In the
   * other two cases, Face denotes a grid cell (a codimension 0 entity).
   */
  template <class Grid, class Face, class Displacement, int dimw=Face::dimensionworld>
  class BoundaryFace
  {
  public:
    static int const facedimension = Face::mydimension;
    static int const celldimension = Face::dimension;
    static int const griddimension = Face::dimensionworld;
    static int const worlddimension = dimw;

    using ctype = typename Face::ctype;
    using Cell = typename Grid::Codim<0>::Entity;

  private:
    enum Kind { VOLUMEMESH, SURFACEMESH, PLANARMESH };
    static Kind const kind = facedimension==celldimension-1? VOLUMEMESH
                           : griddimension==worlddimension?  SURFACEMESH
                           :                                 PLANARMESH;

    using PositionExtension = Dune::FieldVector<ctype,worlddimension-celldimension>;

  public:
    /**
     * \brief Local position vector in the face's reference simplex.
     */
    using LocalPosition = Dune::FieldVector<ctype,facedimension>;

    /**
     * \brief Global position in the world.
     */
    using GlobalPosition = Dune::FieldVector<ctype,worlddimension>;

    /**
     * \brief Local position in the grid cell to which this face belongs.
     */
    using CellLocalPosition = Dune::FieldVector<ctype,celldimension>;

    /**
     * \brief Global position of the grid.
     */
    using CellGlobalPosition = Dune::FieldVector<ctype,griddimension>;

    /**
     * \brief Constructor.
     */
    BoundaryFace(Face const& f, Displacement const* d = nullptr);

    ~BoundaryFace();

    Face const& gridFace() const
    {
      return *face;
    }

    /**
     * \brief Computes the axis-aligned bounding box.
     */
    std::pair<GlobalPosition,GlobalPosition> boundingBox() const;

    /**
     * \brief Returns the type of reference face.
     */
    Dune::GeometryType type() const
    {
      return gt;
    }

    /**
     * \brief Computes the area of the deformed face.
     */
    ctype volume() const
    {
      return vol;
    }

    /**
     * \brief Computes the global world position of the given local coordinate.
     */
    GlobalPosition global(LocalPosition const& xi) const;


    /**
     * \brief Computes the unit outer normal of the deformed face.
     */
    GlobalPosition unitOuterNormal(LocalPosition const& xi) const;

    /**
     * \brief Computes the intersection of the boundary face with the given line segment
     *        \f$ [a,b] \f$, if there is any.
     *
     * \return a pair \f$ (\xi,t) \f$ of local position and line parameter, such that
     *         \f$ \mathrm{global}(\xi) \approx (1-t)a+tb \f$ holds.
     */
    std::optional<std::pair<LocalPosition,ctype>> intersection(GlobalPosition const& a,
                                                               GlobalPosition const& b) const;


  private:
    // In planar meshes, the grid-global coordinates as returned by mesh functions are "too short".
    // We extend the missing entries by the given values in z.
    static GlobalPosition extend(CellGlobalPosition const& x,
                                 PositionExtension const& z);

    // Get the cell corresponding to the face. In surface and embedded meshes, this is
    // the face itself.
    Cell cell() const;

    GlobalPosition globalUncached(LocalPosition const& xi) const;
    GlobalPosition unitOuterNormalUncached(LocalPosition const& xi) const;


    using DeformationJacobian = Dune::FieldMatrix<ctype,worlddimension,facedimension>;
    using GridJacobian = Dune::FieldMatrix<ctype,griddimension,facedimension>;

    // Compute the Jacobian of the deformation at given local position.
    std::pair<DeformationJacobian,GridJacobian> jacobians(LocalPosition const& xi) const;
    std::pair<DeformationJacobian,GridJacobian> jacobiansUncached(LocalPosition const& xi) const;

    std::shared_ptr<Face> face;
    Dune::GeometryType gt;    // geometry type
    bool faceAffine;

    Displacement const* disp;
    bool displacementAffine;
    ctype vol;                // deformed face's area

    // Precomputed (because of being constant) values in case of affine displacement and affine face
    // (implying that the deformed face is affine).
    DeformationJacobian deformationJacobian;    // grid -> world
    GridJacobian gridJacobian;                  // reference -> grid
    GlobalPosition origin;                      // world position of reference origin
    GlobalPosition normal;
  };

  // ----------------------------------------------------------------------------------------------
}


#endif
