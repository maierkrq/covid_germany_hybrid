/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2018-2018 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


#ifndef GAUSSNODES
#define GAUSSNODES

namespace Kaskade
{
  /**
   * \ingroup utilities
   * \brief Computes the n Lobatto nodes on [-1,1] including the interval end points.
   */
  std::vector<double> lobattoNodes(int n);
  
  /**
   * \ingroup utilities
   * \brief Computes the n Radau nodes on [-1,1] including the right end point.
   */
  std::vector<double> radauNodes(int n);
  
  /**
   * \ingroup utilities
   * \brief Computes the n Gauß nodes on [-1,1] excluding the interval end points.
   */
  std::vector<double> gaussNodes(int n);
}

#endif