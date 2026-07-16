/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2013 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef IO_TOOLS_HH
#define IO_TOOLS_HH

#if HAVE_UG
#include <dune/grid/uggrid/uggridfactory.hh>
#endif

#if ENABLE_ALUGRID
#include <dune/alugrid/3d/gridfactory.hh>
#endif

namespace Kaskade
{
  namespace IOTools
  {
    /**
     * \cond internals
     */

    /// Wrapper hiding the differences in the grids constructors;
    template <typename Grid>
    struct FactoryGenerator;

#ifdef HAVE_UG
    template <int dim>
    struct FactoryGenerator<Dune::UGGrid<dim>>
    {
      using Grid = Dune::UGGrid<dim>;
      using Factory = Dune::GridFactory<Grid>;

      static Factory createFactory(unsigned int initialHeapSize=0)
      {
        if(initialHeapSize > 0) Grid::setDefaultHeapSize(initialHeapSize);

        return Factory();
      }
    };
#endif /* HAVE_UG */

#ifdef ENABLE_ALUGRID
    template <int dim, Dune::ALUGridElementType et, Dune::ALUGridRefinementType rt>
    struct FactoryGenerator<Dune::ALUGrid<dim, dim, et, rt>>
    {
      using Factory = Dune::GridFactory<Dune::ALUGrid<dim, dim, et, rt>>;

      static Factory createFactory(unsigned int) { return Factory(); }
    };

#endif /* ENABLE_ALUGRID */

    template <class Grid>
    struct BoundarySegmentIndexWrapper
    {
      template <class Factory, class VertexVector, class BoundaryVector, class TmpIdVector, class IdVector>
      static void readBoundarySegmentIndices(Grid const& grid, Factory const& factory, VertexVector const& vertices, BoundaryVector const& boundary, TmpIdVector const& boundaryIdsFromFile, IdVector& boundaryIds)
      {
        boundaryIds.resize(boundaryIdsFromFile.size());
        for (auto const& cell : elements(grid.leafGridView()))
        {
          for (auto const& intersection : intersections(grid.leafGridView(), cell))
          {
            if (intersection.boundary())
            {
              boundaryIds[intersection.boundarySegmentIndex()] = boundaryIdsFromFile[factory.insertionIndex(intersection)][0];
            }
          }
        }
      }
    };

    /**
     * \endcond
     */
  } // namespace IOTools
} // namespace Kaskade

#endif // IO_TOOLS_HH
