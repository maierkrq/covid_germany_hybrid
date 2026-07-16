/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2020-2020 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "dune/grid/config.h"
#include "dune/grid/uggrid.hh"

#include "fem/boundaryFace.hpp"
#include "fem/deforminggridmanager.hh"
#include "fem/boundaryInterpolation.hh"
#include "fem/boundaryLocator.hh"

namespace Kaskade
{
  using namespace BoundaryLocatorDetail;

  // explicit instantiation
  using UG2 = Dune::UGGrid<2>;
  using UG3 = Dune::UGGrid<3>;
  template class BoundaryFace<UG2,UG2::LeafGridView::Intersection,ZeroDisplacement>;
  template class BoundaryFace<UG2,UG2::LeafGridView::Intersection,H1Space<UG2>::Element_t<2>>;
  template class BoundaryFace<UG3,UG3::LeafGridView::Intersection,ZeroDisplacement>;
  template class BoundaryFace<UG3,UG3::LeafGridView::Intersection,H1Space<UG3>::Element_t<3>>;
  template class BoundaryFace<UG2,UG2::LeafGridView::Intersection,BoundaryDisplacementByInterpolation<UG2::LevelGridView>>;

// NOT YET IMPLEMENTED FOR 3D: cellInterpolationDerivative
//   template class BoundaryLocator<UG3::LeafGridView,BoundaryDisplacementByInterpolation<UG3::LevelGridView>>;

  using DUG2 = DeformingGridManager<UG2>::Grid;
  using DUG3 = DeformingGridManager<UG3>::Grid;
  template class BoundaryFace<DUG2,DUG2::LeafGridView::Intersection,ZeroDisplacement>;
  template class BoundaryFace<DUG2,DUG2::LeafGridView::Intersection,H1Space<DUG2>::Element_t<2>>;
  template class BoundaryFace<DUG3,DUG3::LeafGridView::Intersection,ZeroDisplacement>;
  template class BoundaryFace<DUG3,DUG3::LeafGridView::Intersection,H1Space<DUG3>::Element_t<3>>;
  template class BoundaryFace<DUG2,DUG2::LeafGridView::Intersection,BoundaryDisplacementByInterpolation<DUG2::LevelGridView>>;
}
