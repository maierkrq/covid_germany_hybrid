/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/Numerik/numsoft/kaskade7/                        */
/*                                                                           */
/*  Copyright (C) 2002-2020 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


#include <iostream>
#include <memory>
#include <string>

#define HAVE_UG 1

#include "dune/grid/config.h"
#include "dune/grid/uggrid.hh"

#include "fem/functionspace.hh"
#include "fem/lagrangespace.hh"
#include "fem/spaces.hh"

#include "io/amirameshreader.hh"
#include "io/vtk.hh"

#include "utilities/kaskopt.hh"
#include "utilities/enums.hh"


int main(int argc, char *argv[])
{
  using namespace boost::fusion;
  using namespace Kaskade;

  int verbosity = 1;
  bool dump = false;
  int const dim = 3;

  std::string gridfile, location, component;
  if (getKaskadeOptions(argc,argv,Options
    ("grid",gridfile, "", "file name of Amira grid to be converted")
    ("location",location,"Tetrahedra","location of material parameters in the Amira grid file")
    ("component",component,"Materials","name of material parameter data")
      ))
    return 0;

  // Read the Amira grid
  typedef Dune::UGGrid<dim> Grid;
  GridManager<Grid> gridManager(AmiraMeshReader::readGrid<Grid,float>(gridfile));

  // define piecewise constant function space for material values
  L2Space<Grid> materialSpace(gridManager,gridManager.grid().leafGridView(),0);
  auto material = materialSpace.element<1>();

  // Read the material types from the grid file
  AmiraMeshReader::readData<unsigned char>(gridfile,location,component,material);


  // Write the grid and material data to VTK.
  IoOptions options;
  options.outputType = IoOptions::binary;

  // Remove the Amira file name suffix
  std::string const ending_1(".grid"), ending_2(".am");
  std::string name(gridfile);

  if (name.find(ending_1) != std::string::npos)
    name.replace(name.find(ending_1), ending_1.length(), "");
  if (name.find(ending_2) != std::string::npos)
    name.replace(name.find(ending_2), ending_2.length(), "");

  writeVTK(material,name,options,component);

  std::cout << "Conversion finished." << std::endl;

  return 0;
}
