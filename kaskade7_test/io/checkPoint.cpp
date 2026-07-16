/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2014-2017 Zuse Institute Berlin                             */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// ----------------------------------------------------------------------------

#include <dune/grid/config.h>
#include "checkPoint.hh"

namespace Kaskade {

  CheckPoint::CheckPoint() {
    setTimeStamp();
    directory = timeStamp_;
  }

  CheckPoint::CheckPoint(std::string const& pathName) : directory(pathName) {
    if(!boost::filesystem::exists(directory)) {
      setTimeStamp();
      directory += std::string("_") + timeStamp_;
      boost::filesystem::create_directories(directory);
    }
  }

  void CheckPoint::setTimeStamp() {
    using namespace std::chrono;
    std::time_t time1 = system_clock::to_time_t(system_clock::now());
    std::tm time2 = *std::localtime(&time1);
    std::ostringstream timeStream;
    timeStream << std::put_time(&time2, "%F_%H-%M-%S");
    timeStamp_ = timeStream.str();
  }
}

#ifdef UNITTEST

#include <iostream>

#include <dune/alugrid/grid.hh>
#include <dune/grid/uggrid.hh>

#include <fem/spaces.hh>
#include <fem/norms.hh>

template <class GridView>
void writeGridView(GridView const& gridView, std::ostream & out) {
  constexpr int dim = GridView::dimension;
  auto const& indexSet = gridView.indexSet();
  for(auto const& element : elements(gridView)) {
    out << indexSet.index(element);
    for(int cd = 1; cd <= dim; ++cd) {
      for(int i=0; i<element.subEntities(cd); ++i) {
        out << indexSet.subIndex(element, i, cd);
      }
    }
  }
  // The following neither works for UGGrid nor for ALUGrid.
  // Even just obtaining the face index seems to be not working for ALUGrid sometimes (at least on some grid levels).
//   for(auto const& face : facets(gridView)) {
//     out << indexSet.index(face);
//     for(int cd = 2; cd <= dim; ++cd) {
//       for(int i=0; i<face.subEntities(cd); ++i) {
//         out << indexSet.subIndex(face, i, cd);
//       }
//     }
//   }
//   for(auto const& edge : edges(gridView)) {
//     out << indexSet.index(edge);
//     for(int i=0; i<edge.subEntities(dim); ++i) {
//       out << indexSet.subIndex(edge, i, 0);
//     }
//   }
  for(auto const& vertex : vertices(gridView)) {
    out << indexSet.index(vertex);
  }
}

template <class Grid>
void writeAllViews(Grid const& grid, std::string const& fileName) {
  std::ofstream fs(fileName);
  for(int i=0; i<=grid.maxLevel(); ++i) {
    writeGridView(grid.levelGridView(i), fs);
  }
  writeGridView(grid.leafGridView(), fs);
}

int main(void) {
  using namespace Kaskade;

  constexpr int dim = 3;
  using AluGrid = Dune::ALUGrid<dim,dim,Dune::ALUGridElementType::simplex,Dune::ALUGridRefinementType::conforming>;
  using UgGrid = Dune::UGGrid<dim>;
  using Grid = UgGrid;
  using Spaces = boost::fusion::vector<H1Space<Grid> const*, L2HierarchicSpace<Grid> const*>;
  using VariableDescriptions = boost::fusion::vector<Variable<SpaceIndex<0>, Components<3>>, Variable<SpaceIndex<1>, Components<7>>>;
  using VariableSetDesc = VariableSetDescription<Spaces,VariableDescriptions>;
  using CoefficientVectors = VariableSetDesc::CoefficientVectorRepresentation<0,2>::type;

  std::string pathName = "checkPoint";
  std::string timeStamp;

  double referenceH1Norm = 0.0;
  double referenceL2Norm = 0.0;

  {// write check points
    // refinement parameter
    int globalRefinements1 = 1;
    int adaptRefinements1 = 1;
    int globalRefinements2 = 2;
    int adaptRefinements2 = 1;
    double adaptFraction = 0.1;

    // random numbers
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    auto seed = rd();
    std::mt19937 gen(seed); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0, 1);

    // create check point in write mode, store time stamp for reading later
    CheckPoint checkPoint(pathName);
    timeStamp = checkPoint.timeStamp();
    std::cout << timeStamp << std::endl;

    // create initial grid (which is just one element, the reference element)
    Dune::GridFactory<Grid> factory;
    // insert corners of standard tetrahedron/triangle by vertex coordinates v[0], v[1], (v[2])
    Dune::FieldVector<double,dim> v;
    factory.insertVertex(v);
    v[0] = 1;  factory.insertVertex(v); v[0] = 0;
    v[1] = 1;  factory.insertVertex(v); v[1] = 0;
    if(dim==3) {
      v[2] = 1;  factory.insertVertex(v);
    }
    // insert tetrahedron/triangle by vertex indices
    std::vector<unsigned int> vid(dim+1);
    vid[0]= 0; vid[1]= 1; vid[2]= 2;
    if(dim==3) vid[3] = 3;
    factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex, dim), vid);
    GridManager<Grid> gridManager(factory.createGrid());

    // refine grid and store refinment steps in refHist
    std::vector<std::vector<std::pair<int, typename Grid::LeafIndexSet::IndexType>>> refHist(globalRefinements1);
    for(int i=0; i<globalRefinements1; ++i) {
      checkPoint.appendRefinementHistory(refHist[i], "grid2");
    }
    gridManager.globalRefine(globalRefinements1);
    for(int i=0; i<adaptRefinements1; ++i) {
      refHist.emplace_back();
      auto const& indexSet = gridManager.grid().leafIndexSet();
      for (const auto& element : elements(gridManager.grid().leafGridView())) {
        if(dis(gen) < adaptFraction) {
          gridManager.mark(1, element);
          refHist.back().push_back({1, indexSet.index(element)}); // Todo: make this part of the gridManager
        }
      }
      gridManager.adaptAtOnce();
      if(refHist.back().empty()) { // Todo: make this part of the gridManager
        refHist.pop_back();
      } else {
        checkPoint.appendRefinementHistory(refHist.back(), "grid2");
      }
    }
    refHist.resize(refHist.size()+globalRefinements2);
    for(int i=0; i<globalRefinements2; ++i) {
      checkPoint.appendRefinementHistory(refHist[refHist.size()-1-i], "grid2");
    }
    gridManager.globalRefine(globalRefinements2);
    for(int i=0; i<adaptRefinements2; ++i) {
      refHist.emplace_back();
      auto const& indexSet = gridManager.grid().leafIndexSet();
      for (const auto& element : elements(gridManager.grid().leafGridView())) {
        if(dis(gen) < adaptFraction) {
          gridManager.mark(1, element);
          refHist.back().push_back({1, indexSet.index(element)});
        }
      }
      gridManager.adaptAtOnce();
      if(refHist.back().empty()) {
        refHist.pop_back();
      } else {
        checkPoint.appendRefinementHistory(refHist.back(), "grid2");
      }
    }

    std::cout << "Grid has " << gridManager.grid().maxLevel() << " levels, "
              << gridManager.grid().size(0) << " elements and "
              << gridManager.grid().size(dim) << " vertices." << std::endl;

    // write refined grid
    writeAllViews(gridManager.grid(), "gridO");
    if(Dune::Capabilities::hasBackupRestoreFacilities<Grid>::v) {
      checkPoint.writeGrid(gridManager, "grid");
    } else {
      checkPoint.writeRefinementHistory(refHist, "grid");
    }

    // create data on grid
    H1Space<Grid> h1Space(gridManager, gridManager.grid().leafGridView(), 2);
    L2HierarchicSpace<Grid> l2Space(gridManager, gridManager.grid().leafGridView(), 1);
    Spaces spaces(&h1Space, &l2Space);
    VariableSetDesc variableSetDesc(spaces, {"u", "v"});
    VariableSetDesc::VariableSet u(variableSetDesc);
    H1Space<Grid>::Element_t<3> & feFunction = component<0>(u);
    for(auto & coefficient : feFunction.coefficients()) {
      for(auto & entry : coefficient) {
        entry = dis(gen);
      }
    }
    for(auto & coefficient : component<1>(u).coefficients()) {
      for(auto & entry : coefficient) {
        entry = dis(gen);
      }
    }
    referenceH1Norm = H1Norm{}(feFunction);
    std::cout << "reference H1-norm: " << referenceH1Norm << std::endl;
    referenceL2Norm = L2Norm{}(component<1>(u));
    std::cout << "reference L2-norm: " << referenceL2Norm << std::endl;
    CoefficientVectors uc(u);

    // write data
    checkPoint.writeData(u, "data1");
    checkPoint.writeData(feFunction, "data2");
    checkPoint.writeData(uc, "data3");
    checkPoint.writeData(component<0>(uc), "data4");
  }

  {// read check points
    // create check point in read mode
    CheckPoint checkPoint(pathName + "_" + timeStamp);

    // read grid
//    GridManager<Grid> gridManager(checkPoint.readGrid<Grid>("grid"));
    std::unique_ptr<Grid> grid;
    if(Dune::Capabilities::hasBackupRestoreFacilities<Grid>::v) {
      grid = checkPoint.readGrid<Grid>("grid");
    } else {
      // for UGGrid we create the same intial grid, read the refinement history in the check point, and then refine accordingly
      // create grid
      Dune::GridFactory<Grid> factory;
      // insert corners of standard tetrahedron/triangle by vertex coordinates v[0], v[1], (v[2])
      Dune::FieldVector<double,dim> v;
      factory.insertVertex(v);
      v[0] = 1;  factory.insertVertex(v); v[0] = 0;
      v[1] = 1;  factory.insertVertex(v); v[1] = 0;
      if(dim==3) {
        v[2] = 1;  factory.insertVertex(v);
      }
      // insert tetrahedron/triangle by vertex indices
      std::vector<unsigned int> vid(dim+1);
      vid[0]= 0; vid[1]= 1; vid[2]= 2;
      if(dim==3) vid[3] = 3;
      factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex, dim), vid);
      grid = std::unique_ptr<Grid>(factory.createGrid());

      // refine grid
      auto refHist = checkPoint.readRefinementHistory<typename Grid::LeafIndexSet::IndexType>("grid");
      for(auto const& refStep : refHist) { // Todo: make this part of the gridManager
        if(refStep.empty()) grid->globalRefine(1);
        else {
          auto const& indexSet = grid->leafIndexSet();
          int refCellIdx = 0;
          // the following only works if grid elements are visited in the same order as above. This is probably not guaranteed in theory but should work in praxis.
          for (const auto& element : elements(grid->leafGridView())) {
            if(indexSet.index(element) == refStep[refCellIdx].second) {
              grid->mark(refStep[refCellIdx].first, element);
              ++refCellIdx;
            }
          }
          if(refCellIdx != refStep.size()) {
            throw std::logic_error("Grid traversion order not as expectet.");
          }
          grid->preAdapt();
          grid->adapt();
          grid->postAdapt();
        }
      }
    }
    GridManager<Grid> gridManager(std::move(grid));

    // write some grid information for comparison
    writeAllViews(gridManager.grid(), "gridC");

    std::cout << "Grid has " << gridManager.grid().maxLevel() << " levels, "
              << gridManager.grid().size(0) << " elements and "
              << gridManager.grid().size(dim) << " vertices." << std::endl;

    // read in data and check
    H1Space<Grid> h1Space(gridManager, gridManager.grid().leafGridView(), 2);
    L2HierarchicSpace<Grid> l2Space(gridManager, gridManager.grid().leafGridView(), 1);
    Spaces spaces(&h1Space, &l2Space);
    VariableSetDesc variableSetDesc(spaces, {"u", "v"});
    VariableSetDesc::VariableSet u(variableSetDesc);
    H1Space<Grid>::Element_t<3> & feFunction = component<0>(u);

    checkPoint.readData(feFunction.coefficients(), "data2");
    double norm = H1Norm{}(feFunction);
    std::cout << "Relative error in H1-norm after checkpoint write and read is: " << (referenceH1Norm-norm)/referenceH1Norm << std::endl;
    feFunction = 0.0;

    CoefficientVectors uc(u);
    checkPoint.readData(uc, "data1");
    u = uc;
    norm = H1Norm{}(feFunction);
    std::cout << "Relative error in H1-norm after checkpoint write and read is: " << (referenceH1Norm-norm)/referenceH1Norm << std::endl;
    norm = L2Norm{}(component<1>(u));
    std::cout << "Relative error in L2-norm after checkpoint write and read is: " << (referenceL2Norm-norm)/referenceL2Norm << std::endl;
    u = 0.0;

    checkPoint.readData(feFunction, "data4");
    norm = H1Norm{}(feFunction);
    std::cout << "Relative error in H1-norm after checkpoint write and read is: " << (referenceH1Norm-norm)/referenceH1Norm << std::endl;
    feFunction = 0.0;

    checkPoint.readData(u, "data3");
    norm = H1Norm{}(feFunction);
    std::cout << "Relative error in H1-norm after checkpoint write and read is: " << (referenceH1Norm-norm)/referenceH1Norm << std::endl;
    norm = L2Norm{}(component<1>(u));
    std::cout << "Relative error in L2-norm after checkpoint write and read is: " << (referenceL2Norm-norm)/referenceL2Norm << std::endl;
  }
  return 0;
}

#endif
