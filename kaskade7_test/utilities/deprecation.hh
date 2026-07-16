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

#ifndef DEPRECATION_HH
#define DEPRECATION_HH

#include <string>

namespace Kaskade 
{
  /**
   * \ingroup utilities
   * \brief Convenience class for marking deprecated functions.
   * A static object of this class can be created inside deprecated functions. It will 
   * output a warning message to stderr if the deprecated function is called, but only 
   * once.
   * \code
   * void f()
   * {
   *   static Deprecated dummy("f() is deprecated, use g() instead",2018);
   * }
   * \endcode
   * After the given grace period is passed, the program is aborted.
   */
  class Deprecated
  {
  public:
    /**
     * \brief Constructor.
     * \param comment The comment to be written to stderr along with a generic warning
     * \param year the grace period
     */
    Deprecated(std::string const& comment, int year);
  };
  
}

#endif
