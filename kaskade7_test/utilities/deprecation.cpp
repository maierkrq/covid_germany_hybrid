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

#include <chrono>
#include <cstdlib>
#include <iostream>

#include "utilities/deprecation.hh"

namespace Kaskade 
{
  Deprecated::Deprecated(std::string const& comment, int year)
  {
      
    
    std::cerr << "WARNING: deprecated function called.\n"
              << comment << std::endl;

    // If the grace period passed, abort the program in order to increase the motivation for 
    // replacing the deprecated function call.
    auto now = std::time(nullptr);
    auto nowptr = std::localtime(&now);
    if (nowptr && nowptr->tm_year+1900 <= year)
      std::cerr << "The deprecated function will cease working after the end of year " << year << "\n";
    else
    {
      std::cerr << "The grace period ended after year " << year << ".";
      std::abort();
    }
  }
}

