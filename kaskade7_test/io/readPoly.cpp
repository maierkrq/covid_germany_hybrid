/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2019 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

///#include <cctype>
#include <fstream>
#include <sstream>

#include "dune/grid/config.h"

#include "io/readPoly.hh"
#include "utilities/detailed_exception.hh"

namespace  {

  // reads the next nonempty line from the given stream into the given string. 
  // Comments starting with '#' extend up to the next newline and are considered as "empty".
  void readLine(std::istream& in, std::string& line) 
  {
    do {
      std::getline(in,line);                                          // read complete line
      line.erase(std::min(line.find('#'),line.size()),line.size());   // erase everything from # on
      
      size_t i=0;                                                     // Find first non-space character
      for (; i<line.size(); ++i)
        if (!std::isspace(line[i]))
          break;
      line.erase(0,i);                                                // erase leading whitespace
    } while(line.size()==0);                                          // until we found a nonempty line
  }
  
  std::vector<Dune::FieldVector<double,2>> readNodeData(std::string const& name) 
  {
    using namespace std;
    
    //open file with nodes
    ifstream in(name.c_str());
    if (!in) // check for successful opening
      throw Kaskade::FileIOException("Cannot open for reading.",name,__FILE__,__LINE__);
    
    string line;
    readLine(in,line);
    istringstream is(line);
    
    int nVertices, dimension, nAttributes, nBoundaryMarkers;
    is >> nVertices >> dimension >> nAttributes >> nBoundaryMarkers;
    assert(dimension==2); // TODO: throw exception
    assert(0<=nBoundaryMarkers && nBoundaryMarkers<=1);
    
    std::vector<Dune::FieldVector<double,2>> nodes(nVertices);
    for (int i=0; i<nVertices; ++i) {
      readLine(in,line);
      istringstream is(line);
      int id; // the node number
      is >> id;
      
      if (id<0 || id>nVertices)
        throw Kaskade::FileIOException("Node index out of range.",name,__FILE__,__LINE__); 
      
      // Note that .poly files may be 0-based or 1-based. We just wrap around to cover both cases.
      id = id % nVertices;
      // Read in the coordinates
      is >> nodes[id][0] >> nodes[id][1];
    }
    
    return nodes;
  }
  
  // Reads a sequence of edges (as pairs of vertex list and boundary marker ([v1,v2],bm)) from a .edge file
  std::vector<std::pair<std::vector<unsigned int>,int>> readEdgeData(std::string const& name, int nNodes)
  {
    using namespace std;
    
    ifstream in(name+".edge");
    if (!in)
      throw Kaskade::FileIOException("Cannot open for reading.",name+".edge",__FILE__,__LINE__);
    
    // read header
    string line;
    readLine(in,line);
    istringstream is(line);
    int nEdges, nBoundaryMarkers;
    is >> nEdges >> nBoundaryMarkers;
    if (nEdges<3 || nBoundaryMarkers<0 || nBoundaryMarkers>1)
      throw Kaskade::FileIOException("Edge file header corrupt.",name+".edge",__FILE__,__LINE__);
    if (nBoundaryMarkers==0)    // no boundary markers specified -> clear the map
      throw Kaskade::FileIOException("No boundary markers present.",name+".edge",__FILE__,__LINE__);

    // read boundary segments
    vector<pair<vector<unsigned int>,int>> edges;
    for (int i=0; i<nEdges; ++i)
    {
      readLine(in,line);
      istringstream is(line);
      int id, v1, v2, bm;
      is >> id >> v1 >> v2 >> bm;
      if (v1<0 || v1>nNodes || v2<0 || v2>nNodes)
        throw Kaskade::FileIOException("Vertex numbers in edge file corrupt.",name+".edge",__FILE__,__LINE__);
      
      if (bm != 0) // Boundary marker 0 is assigned to *interior* edges by triangle. Skip them.
        edges.push_back(make_pair(vector<unsigned int>{static_cast<unsigned int>(v1%nNodes),
                                                       static_cast<unsigned int>(v2%nNodes)},
                                  bm));
    }
    
    return edges;
  }
  
  
} // End of namespace 

std::unique_ptr<Dune::UGGrid<2>> readPolyData(std::string const& name, std::vector<int>* boundaryMarker, 
                                              int heapSize, int envSize) 
{
  using namespace std;
  
  Dune::UGGrid<2>::setDefaultHeapSize(heapSize);
  Dune::UGGrid<2>* grid(new Dune::UGGrid<2>);
  Dune::GridFactory<Dune::UGGrid<2>> factory(grid);
  
  
   // read node file
  auto nodes = readNodeData(name+".node"); 
  for (auto node: nodes) 
    factory.insertVertex(node);

  // read elements file
  ifstream in(name+".ele");
  if (!in)
    throw Kaskade::FileIOException("Cannot open for reading.",name+".ele",__FILE__,__LINE__);
  
  string line;
  readLine(in,line);
  istringstream is(line);

  // read header
  int nEle, nNode, nAttributes;
  is >> nEle >> nNode >> nAttributes;

  assert(nNode>=3); // TODO throw

  std::vector<unsigned int> elemNodes(3);
  for (int i=0; i<nEle; ++i) {
    readLine(in,line);
    istringstream is(line);
    int id;
    is >> id >> elemNodes[0] >> elemNodes[1] >> elemNodes[2];
    
    if (!(0<=elemNodes[0] && elemNodes[0]<=nodes.size() && 
          0<=elemNodes[1] && elemNodes[1]<=nodes.size() && 
          0<=elemNodes[2] && elemNodes[2]<=nodes.size()) )
      throw Kaskade::FileIOException("Node index of element out of range.",name,__FILE__,__LINE__);
    
    // TODO: read attributes
    
    // Remember that node numbers can be 0-based or 1-based. As in readNodeData, we just wrap around.
    for (auto& n: elemNodes) 
      n = n % nodes.size();
    
    factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2),elemNodes);
  }
  
  // read edge file if requested
  if (boundaryMarker)
  {
    boundaryMarker->clear();
    
    auto edges = readEdgeData(name,nodes.size());                 // read boundary edges and their markers
    for (auto const& e: edges)
    {
      factory.insertBoundarySegment(e.first);                     // creates a new boundary segment for this edge
      boundaryMarker->push_back(e.second);                        // and stores its boundary marker
    }
  }
   
  std::unique_ptr<Dune::UGGrid<2>> gridptr(factory.createGrid());
  return gridptr;
}

  
