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

#include "fem/cellLocator.hpp"
#include "fem/deforminggridmanager.hh"
#include "fem/boundaryInterpolation.hh"

namespace Kaskade
{
  
  using namespace CellLocatorDetail;
  
  // explicit instantiation
  using UG2 = Dune::UGGrid<2>;
  using UG3 = Dune::UGGrid<3>;
  template class CellLocator<UG2::LeafGridView>;
  template class CellLocator<UG3::LeafGridView>;

  using DUG2 = DeformingGridManager<UG2>::Grid;
  using DUG3 = DeformingGridManager<UG3>::Grid;
  template class CellLocator<DUG2::LeafGridView>;
  template class CellLocator<DUG3::LeafGridView>;
}


#ifdef UNITTEST

using namespace Kaskade;

#include <iostream>

#include "fem/gridmanager.hh"
#include "io/readPoly.hh"
#include "dune/grid/config.h"
#include "dune/grid/uggrid.hh"
#include "utilities/gridGeneration.hh"
#include "utilities/kaskopt.hh"
#include "utilities/timing.hh"

int main(void)
{

  
  return 0;
}

#endif
