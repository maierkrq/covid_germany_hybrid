/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2017-2017 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef DEFORMINGGRIDMANAGER_HH_
#define DEFORMINGGRIDMANAGER_HH_

#include "dune/common/fvector.hh"
#include "dune/grid/geometrygrid/coordfunction.hh"
#include "dune/grid/geometrygrid/grid.hh"
#include "dune/grid/uggrid.hh"

#include "fem/gridmanager.hh"
#include "fem/boundaryInterpolation.hh"
#include "fem/spaces.hh"

namespace Kaskade
{
  /// \cond internals
  namespace DeformingGridManagerDetail
  {
    template <class Grid>
    class Deformation: public Dune::DiscreteCoordFunction<typename Grid::ctype,
                                                          Grid::dimension,
                                                          Deformation<Grid>>
    {
      using ctype = typename Grid::ctype;
      constexpr static int dim = Grid::dimension;
      using Vector = Dune::FieldVector<ctype,dim>;
      using LevelView = typename Grid::LevelGridView;
      using Index = typename LevelView::IndexSet::IndexType;
      
    public:
      /**
       * \brief Constructor 
       * \param featureCurveAngle the angle (in radians) between two edges (feature edges in 3D, boundary edges in 2D)
       *                          that determines whether they form a differentiable
       *                          feature curve. Defaults to roughly 57°.
       * \param featureEdgeAngle the angle (in radians) between adjacent facets' normals that determines whether they share a
       *                         feature edge or a normal (\f$ G^1 \f$ continuous) edge. Defaults to roughly 29°. This is
       *                         only of effect in 3D grids.
       */
      Deformation(Grid const& grid_, ctype featureEdgeAngle=0.5, ctype featureCurveAngle=1.0)
      : grid(grid_), tangents_(grid,GridAngleFeatureDetector<LevelView>(featureCurveAngle,featureEdgeAngle))
      {
        adapt();
      }

      Deformation(Grid const& grid_, GeometricFeatureDetector<Index,ctype,dim> const& featureDetector)
      : grid(grid_)
      , tangents_(grid,featureDetector)
      {
        adapt();
      }
      
      void adapt()
      {
        vertexPositions_ = interpolateVertexPositions(grid.leafGridView(),tangents_);
      }
      
      void setFeatureAngles(ctype featureCurveAngle=1.0, ctype featureEdgeAngle=0.5)
      {
        GridAngleFeatureDetector<LevelView> features(featureCurveAngle,featureEdgeAngle);
        tangents_ = GridBoundaryTangents<LevelView>(grid,features);
        adapt();
      }
      
      void evaluate(typename Grid::template Codim<dim>::Entity const& vertex, unsigned int /* corner */,
                    Vector& y) const
      {
        y = vertexPositions_[grid.leafGridView().indexSet().index(vertex)];
      }
      
      void evaluate(typename Grid::template Codim<0>::Entity const& cell, unsigned int corner,
                    Vector& y) const
      {
        y = vertexPositions_[grid.leafGridView().indexSet().subIndex(cell,corner,dim)];
      }

      GridBoundaryTangents<typename Grid::LevelGridView> const& tangents() const
      {
        return tangents_;
      }

      std::vector<Vector> const& vertexPositions() const
      {
        return vertexPositions_;
      }
      
    private:
      Grid const& grid;
      GridBoundaryTangents<typename Grid::LevelGridView> tangents_;
      std::vector<Vector> vertexPositions_;
    };
  }
  /// \endcond
  
  
  
  /**
   * \ingroup grid 
   * \brief A grid manager for asymptotically \f$ G^1 \f$-continuous domains.
   * 
   * On construction, a deformation of the grid domain is computed, such that the deformed domain has a
   * \f$ G^1 \f$-continuous boundary. On mesh refinement, the new vertices are moved such that the 
   * smooth boundary is interpolated.
   */
  template <class Grid>
  class DeformingGridManager 
  : public GridManagerBase<Dune::GeometryGrid<Grid,DeformingGridManagerDetail::Deformation<Grid>>>
  {
    using Deformation = DeformingGridManagerDetail::Deformation<Grid>;
    using DeformedGrid = Dune::GeometryGrid<Grid,Deformation>;
    using Base = GridManagerBase<DeformedGrid>;
    using ctype = typename Grid::ctype;
    constexpr static int dim = Grid::dimension;
    using Vector = Dune::FieldVector<ctype,dim>;
    using LevelView = typename Grid::LevelGridView;
    using Index = typename LevelView::IndexSet::IndexType;

  public:
    /**
     * \brief Constructor 
     *
     * See BoundaryTangents for a detailed description of feature edge definition.
     *
     * \param featureCurveAngle the angle (in radians) between two feature edges that determines
     *                          whether they form a differentiable
     *                          feature curve. Defaults to roughly 57°.
     * \param featureEdgeAngle the angle (in radians) between adjacent facets' normals that determines whether they share a
     *                         feature edge or a normal (\f$ G^1 \f$ continuous) edge. Defaults to roughly 29°.
     */
    DeformingGridManager(std::unique_ptr<Grid>&& grid,
                         ctype featureCurveAngle=1.0, ctype featureEdgeAngle=0.5,
                         bool verbose=false, bool enforceConcurrentReads=false)
    : Base(new DeformedGrid(grid.release()),                         // geometry grid takes ownership of grid, and
           verbose,enforceConcurrentReads)                           // creates its own deformation function
    {
      this->gridptr->coordFunction().setFeatureAngles(featureCurveAngle,featureEdgeAngle);
    }

    DeformingGridManager(std::unique_ptr<Grid>&& grid,
                         GeometricFeatureDetector<Index,ctype,dim> const& featureDetector,
                         bool verbose=false, bool enforceConcurrentReads=false)
    : Base(makeDeformedGrid(std::move(grid),featureDetector), // geometry grid takes ownership of grid and
           verbose,enforceConcurrentReads)                                          // deformation function
    {
    }

    virtual void update()
    {}
    
    virtual void refineGrid(int refCount)
    {
      this->gridptr->globalRefine(refCount);
    }
    
    virtual bool adaptGrid()
    {
      return this->gridptr->adapt();
    }
    
    /**
     * \brief correctingDisplacement provides a displacement function which displaces the points on the boundary of the grid
     * on to the actual domain boundary.
     *
     * This function is zero on all grid vertices, so it is only useful if requested for FE order greater than one.
     * We leave it undefined what values the provided correction has in the interior of the grid.
     *
     * \param correction FE function of any order >1, will be the computed boundary correction on exit.
     */
    template <class Scalar>
    void correctingDisplacement(FunctionSpaceElement<H1Space<DeformedGrid,Scalar>,dim>& correction)
    {
      if(&correction.space().gridManager() != this)
        throw DetailedException("Correction has wrong underlying gridManager.",__FILE__,__LINE__);


      if (correction.space().mapper().maxOrder()==1)      // The boundary displacement function of order 1
      {                                                   // just interpolates the vertex displacement (which is zero).
        correction = static_cast<Scalar>(0.0);
        return;
      }

      BoundaryInterpolationDisplacement<typename DeformedGrid::LevelGridView> u(this->gridptr->levelGridView(0),
                                                                                this->gridptr->coordFunction().tangents());
      // Todo: This is quite inefficient, since expensive surface interpolation is performed for many points
      //       (those which lie on the cell boundary) several times
      //       Instead: Consider iterating over the shape function nodes in each cell and writing into the coefficient vector directly
      //       with computing surface interpolation only if coefficient was not already set previously
      interpolateGloballyWeak(correction,u);

      // here we have correction as displacement with respect to the undeformed grid
      // but we want it with respect to the current deformed grid, therefore we do the following

      H1Space<DeformedGrid, Scalar> firstOrderSpace(*this, this->gridptr->leafGridView(), 1);
      FunctionSpaceElement<H1Space<DeformedGrid, Scalar>, dim> vertexAdjustment(firstOrderSpace);
      vertexAdjustment = correction;

      // vertexAdjustment is now exactly the displacement from the undefomed to the current deformed grid
      // therefore we subtract it from the computed correction
      correction -= vertexAdjustment;
    }

  private:
    static DeformedGrid* makeDeformedGrid(std::unique_ptr<Grid>&& grid,
                                          GeometricFeatureDetector<Index,ctype,dim> const& featureDetector)
    {
      Deformation* deformation = new Deformation(*grid,featureDetector);
      Grid* g = grid.release();
      return new DeformedGrid(g,deformation);
    }
  };
}

#endif /* DEFORMINGGRIDMANAGER_HH_ */
