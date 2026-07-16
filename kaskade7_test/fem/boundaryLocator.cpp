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

#include "dune/grid/config.h"
#include "dune/grid/uggrid.hh"

#include "fem/boundaryLocator.hpp"
#include "fem/deforminggridmanager.hh"
#include "fem/boundaryInterpolation.hh"

namespace Kaskade
{
  
  using namespace BoundaryLocatorDetail;
  
  
  // explicit instantiation
  using UG2 = Dune::UGGrid<2>;
  using UG3 = Dune::UGGrid<3>;
  template class BoundaryLocator<UG2::LeafGridView,ZeroDisplacement>;
  template class BoundaryLocator<UG2::LeafGridView,H1Space<UG2>::Element_t<2>>;
  template class BoundaryLocator<UG3::LeafGridView,ZeroDisplacement>;
  template class BoundaryLocator<UG3::LeafGridView,H1Space<UG3>::Element_t<3>>;
  template class BoundaryLocator<UG2::LeafGridView,BoundaryDisplacementByInterpolation<UG2::LevelGridView>>;

// NOT YET IMPLEMENTED FOR 3D: cellInterpolationDerivative
//   template class BoundaryLocator<UG3::LeafGridView,BoundaryDisplacementByInterpolation<UG3::LevelGridView>>;

  using DUG2 = DeformingGridManager<UG2>::Grid;
  using DUG3 = DeformingGridManager<UG3>::Grid;
  template class BoundaryLocator<DUG2::LeafGridView,ZeroDisplacement>;
  template class BoundaryLocator<DUG2::LeafGridView,H1Space<DUG2>::Element_t<2>>;
  template class BoundaryLocator<DUG3::LeafGridView,ZeroDisplacement>;
  template class BoundaryLocator<DUG3::LeafGridView,H1Space<DUG3>::Element_t<3>>;
  template class BoundaryLocator<DUG2::LeafGridView,BoundaryDisplacementByInterpolation<DUG2::LevelGridView>>;

// NOT YET IMPLEMENTED FOR 3D: cellInterpolationDerivative
//   template class BoundaryLocator<DUG3::LeafGridView,BoundaryDisplacementByInterpolation<DUG3::LevelGridView>>;
}


#ifdef UNITTEST

using namespace Kaskade;

#include <iostream>

#include "fem/gridmanager.hh"
#include "io/readPoly.hh"

int main(void)
{
  using Grid = Dune::UGGrid<2>;
  GridManager<Grid> gridManager(readPolyData("testgrid.1"));
  auto const& view = gridManager.grid().leafGridView();
  
  BoundaryLocator<Grid::LeafGridView> boundaryLocator(view);
  
  for (auto const& cell: elements(view))
    for (auto fi=view.ibegin(cell); fi!=view.iend(cell); ++fi)
      if (!fi->neighbor() && fi->boundary())
      {
        auto center = fi->geometry().center();
        auto n = fi->centerUnitOuterNormal();
        auto intersection = boundaryLocator.byLineSegment(center,center+1e6*n);
        std::cout << "From " << center << " in direction " << n << ": ";
        if (intersection)
          std::cout << "intersection at " << intersection->first.geometry().global(intersection->second) << "\n";
        else
          std::cout << "no intersection\n";
      }
  return 0;
}

#endif
