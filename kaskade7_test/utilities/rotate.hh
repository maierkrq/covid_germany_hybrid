/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2015 Zuse Institute Berlin                                 */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef UTILITIES_ROTATE_HH
#define UTILITIES_ROTATE_HH

#include <array>
#include <vector>
#include <cmath>
#include <dune/common/fvector.hh>
//#include "fem/forEach.hh"
//#include "fem/gridmanager.hh"

namespace Kaskade 
{

/**
   * \brief Rotates a 2D object around it's center (and restores the upper boundary)
   *
   */

void rotate(double alpha, std::array<double,2> *barycenter, std::vector<Dune::FieldVector<double,2> > & nodes);

void find2DBarycenter(std::vector<Dune::FieldVector<double,2> > const & nodes, std::array<double,2>* barycenter);

void translate(double deviation, std::vector<Dune::FieldVector<double,2> > & nodes);

void translateHorizontally(double shift, std::vector<Dune::FieldVector<double,2> > & nodes);

void translateVertically(double shift, std::vector<Dune::FieldVector<double,2> > & nodes);

}

#endif

