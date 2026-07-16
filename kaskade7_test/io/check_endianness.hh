/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2018 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CHECK_ENDIANNESS_HH_
#define CHECK_ENDIANNESS_HH_

namespace Kaskade{

  // endian types
  enum class Endian {
    Unknown,
    Big,
    Little,
    BigWord,   /* Middle-endian, Honeywell 316 style */
    LittleWord /* Middle-endian, PDP-11 style */
  };

  /// Return endianness
  Endian endianness();

  /// Return true if byte order is big endian
  bool bigEndian();

  /// Return true if byte order is little endian
  bool littleEndian();
}

#endif
