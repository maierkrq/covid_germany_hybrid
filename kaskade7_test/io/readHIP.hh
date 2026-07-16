/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/en/numerik/software/kaskade-7.html               */
/*                                                                           */
/*  Copyright (C) 2015 Zuse Institute Berlin                                 */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef READHIP_HH
#define READHIP_HH

#include <fstream>
//#include <memory>
#include <string>
#include "boost/multi_array.hpp"

/**
 * \brief Read out motion gate data from .HIP files
 *
 * \param name the base name of the file.
 */
//pointer to what container? vector, array, ...?
boost::multi_array<double,2> readHIPData(std::string const& name);

boost::multi_array<double,2> readHIPData2(std::string const& name);

#endif // READHIP_HH
