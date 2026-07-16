/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2011 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>

#include "dune/grid/config.h"
#include "dune/grid/uggrid.hh"
#include "dune/grid/io/file/amirameshreader.hh"
#include "dune/grid/io/file/amirameshwriter.hh"
#include "dune/grid/common/gridinfo.hh"

#include "io/vtk.hh"
#include "io/matlab.hh"
#include "io/amirameshreader.hh"
#include "io/amira.hh"
#include "amiramesh.hh"

#include "utilities/kaskopt.hh"

bool check(std::string const& file1, std::string const& file2)
{
  bool pass = true;
  
  std::ifstream fileRead(file1);
  std::ifstream fileWritten(file2);
  std::fstream  tmpFile;
  std::string find1("nNodes ");
  std::string find2("nTetrahedra ");
  std::string find3("TetrahedronData { byte Materials }");
  std::string matId("@@@@@@@@@@");
  std::string line, line2;
  int c1, c2;
  int nNodes1, nNodes2, nTet1, nTet2;
  bool matSection = false;
  
  tmpFile.open ("tmp.txt", std::fstream::in | std::fstream::out | std::fstream::app);
  
  while( getline(fileRead, line) ){
    
    if( line.length() >= find1.length() && line.substr(0,find1.length()).compare(find1) == 0) nNodes1 = std::stoi(line.substr(find1.length()));
    if( line.length() >= find2.length() && line.substr(0,find2.length()).compare(find2) == 0) nTet1 = std::stoi(line.substr(find2.length()));
    if( line.length() >= find3.length() && line.substr(0,find3.length()).compare(find3) == 0) matId = line.substr(line.find("@"));
    
    if( line.substr(0,1).compare("@") == 0|| line.empty()) matSection = false;
    if( line.substr(0,matId.length()).compare(matId) == 0 ) {
      matSection = true;
      getline(fileRead, line);
    }
    
    if(matSection) {
      tmpFile << line << std::endl;
    }
  }
  fileRead.close();
  tmpFile.close();
  
  matSection = false;
  matId = "@@@@@@@@@@";
  
  tmpFile.open ("tmp.txt", std::fstream::in | std::fstream::out | std::fstream::app);
  while( getline(fileWritten, line) ){
    
    if( line.length() >= find1.length() && line.substr(0,find1.length()).compare(find1) == 0) nNodes2 = std::stoi(line.substr(find1.length()));
    if( line.length() >= find2.length() && line.substr(0,find2.length()).compare(find2) == 0) nTet2 = std::stoi(line.substr(find2.length()));
    if( line.length() >= find3.length() && line.substr(0,find3.length()).compare(find3) == 0) matId = line.substr(line.find("@"));

    
    if( line.substr(0,1).compare("@") == 0 || line.empty() ) matSection = false;
    if( line.substr(0,matId.length()).compare(matId) == 0) {
      matSection = true;
      getline(fileWritten, line);
    }
    
    if(matSection) {
      getline(tmpFile, line2);
      c1 = std::stoi(line);
      c2 = std::stoi(line2);
      if( c1 != c2 ) {
        pass = false;
      }
    }
  }
  fileWritten.close();
  tmpFile.close();
  std::remove("tmp.txt");
  
  if( nNodes1 != nNodes2 || nTet1 != nTet2 ) {
    pass = false;
  }
  
  return pass;
}

int main(int argc, char *argv[]) 
{
  
  
  int verbosityOpt = 0;
  bool dump = false; 
  std::unique_ptr<boost::property_tree::ptree> pt = getKaskadeOptions(argc, argv, verbosityOpt, dump);
  
  int verbosity   = getParameter(pt, "verbosity", 0);
  // if true, then the test result will be written in a file
  bool result = getParameter(pt, "result",0);
  
  std::stringstream message("Test succeeded", std::stringstream::out);
  
  bool pass1 = true;
  std::string outStr1;
  bool pass2 = true;
  std::string outStr2;
  
  using namespace Kaskade;
  using namespace boost::fusion;
  
  constexpr int DIM = 3;
  using Grid = Dune::UGGrid<DIM>;
  using LeafView = Grid::LeafGridView;
  using Spaces = boost::fusion::vector<H1Space<Grid> const*>;
  using VariableDescriptions = boost::fusion::vector<VariableDescription<0,DIM,0>>;
  using VarSetDesc = VariableSetDescription<Spaces,VariableDescriptions>;
  using MaterialSpace = FEFunctionSpace<DiscontinuousLagrangeMapper<double,Grid::LeafGridView>>;
  using Material = MaterialSpace::Element<1>::type;
  
  std::string fileIn("randomCube.am");
  Dune::GridFactory<Grid> factory;
  DuneAmiraMesh<DIM,DIM,Grid> mesh(fileIn.c_str());
  mesh.InsertUGGrid(factory);
  GridManager<Grid> gridManager( factory.createGrid() );
    
  // material of the geometry
  MaterialSpace materialSpace(gridManager,gridManager.grid().leafGridView(), 0);
  Material material(materialSpace);
  Kaskade::AmiraMeshReader::readData<unsigned char,Material>(fileIn, "Tetrahedra", "Materials", material);

  // construction of finite element space for the scalar solution T.
  H1Space<Grid> h1Space(gridManager,gridManager.grid().leafGridView(),1);
  Spaces spaces(&h1Space);
  std::string varNames[1] = { "u" };
  VarSetDesc varSetDesc(spaces,varNames);
  VarSetDesc::VariableSet x(varSetDesc);

  //write Amira-file
  IoOptions options;
  options.outputType = IoOptions::ascii;
  LeafView leafGridView = gridManager.grid().leafGridView();
  std::string fileOut("randomCube_out.am");
  writeAMIRAFile(leafGridView, x, fileOut, options, material.coefficients(), fileIn);
  
  std::string str("checkDefault.am");
  writeAMIRAFile(leafGridView, x, str, options, material.coefficients());
  
  std::string str2("checkDefault2.am");
  writeAMIRAFile(leafGridView, x, str2, options);

  pass1 = check(fileIn, fileOut);
  pass2 = check(str, fileOut);

  if(!(pass1 && pass2)) message << "Test failed!";
  std::cout << message.str() << std::endl;
  
  if(result) {
    std::string description = "Test reading/writing Amira-files with material data.";
    std::ofstream outfile("../testResult.txt", std::ofstream::out | std::ofstream::app);
    outfile << description << std::endl << message.str() << std::endl << std::endl;
    outfile.close();
  }

  return !(pass1 && pass2);
}
