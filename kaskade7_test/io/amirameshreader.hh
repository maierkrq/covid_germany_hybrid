/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2020 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef AMIRAMESHREADER_HH
#define AMIRAMESHREADER_HH

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <memory>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <cassert>
#include <stdexcept>

#include <boost/timer/timer.hpp>
#include <boost/lexical_cast.hpp>

#include <dune/common/fvector.hh>

#include <amiramesh/AmiraMesh.h>
#include "fem/lagrangespace.hh"
#include "io/ioTools.hh"
#include "io/vtk.hh"
#include "utilities/detailed_exception.hh"

// forward declaration
template <> class std::complex<float>;

namespace Kaskade
{
  /**
   * \ingroup gridInput
   * \brief Reader for Amira meshes.
   * 
   * Note that the Amira mesh file must be in ASCII format. (I think this is not true, 
   * can also read binary Amira mesh, but only Big Endian.)
   * Currently UG and ALUSimplex grids are supported.
   *
   */
  namespace AmiraMeshReader
  {
    namespace ImplementationDetails
    {
      std::string exceptionMessage(std::string const& function, std::string const& gridfile, int line)
      {
        return std::string("In AmiraMeshReader::") + function + ", line " + boost::lexical_cast<std::string>(line) + ": Error reading " + gridfile + "\n";
      }

      std::string exceptionMessage(std::string const& function, std::string const& location, std::string const& component, int line)
      {
        return std::string("In AmiraMeshReader::") + function + ", line " + boost::lexical_cast<std::string>(line)  + std::string(": Error reading ") + component + " of " + location + "\n";
      }


      template <class Source, class value_type, int components>
      void extractData(size_t const n, Source const* source,
                       std::vector<Dune::FieldVector<value_type,components> >& data,
                       bool verbose=false)
      {
        if(verbose)
        {
          std::cout << n << " entries...\n"; 
          std::cout.flush();
        }
        data.clear();
        data.reserve(n);
        data.insert(data.begin(),n,Dune::FieldVector<value_type,components>(0));

        for (size_t i=0; i<n; ++i)
          for (int j=0; j<components; ++j)
            data[i][j] = static_cast<value_type>(source[i*components+j]);
      }

      template <class value_type, int components>
      bool readData(AmiraMesh& am, char const* location, char const* name,
                    std::vector<Dune::FieldVector<value_type,components>>& out,
                    bool verbose=false)
      {
        if(verbose)
          std::cout << "Looking for " << name << " at " << location << " with " << components << " components\n";

        // Scan the contained data for the desired field
        for (int i=0; i<am.dataList.size(); ++i)
          {
            auto const& data = *am.dataList[i];
            if (   components==data.dim()                                   // correct number of components
                && std::strcmp(name,    data.name())==0                     // correct name
                && std::strcmp(location,data.location()->name())==0 )       // correct location
            {
              // Number of nodes contained in the data set
              size_t const n = data.location()->nNodes();

              // Check which basic type actually is contained within the
              // data set and convert the pointer type accordingly.
              if (data.primType()==MC_FLOAT)
                extractData(n,static_cast<float const*>(data.dataPtr()),out);
              else if (data.primType()==MC_DOUBLE)
                extractData(n,static_cast<double const*>(data.dataPtr()),out);
              else if (data.primType()==MC_UINT8)
                extractData(n,static_cast<std::uint8_t const*>(data.dataPtr()),out);
              else if (data.primType()==MC_INT16)
                extractData(n,static_cast<std::int16_t const*>(data.dataPtr()),out);
              else if (data.primType()==MC_UINT16)
                extractData(n,static_cast<std::uint16_t const*>(data.dataPtr()),out);
              else if (data.primType()==MC_INT32)
                extractData(n,static_cast<std::int32_t const*>(data.dataPtr()),out);
              else
                continue; // didn't find a known data type - look on

              return true;      // found
            }
          }

          // nothing appropriate found
          if(verbose)
            std::cout << "Did not find specified data.\n";
          for (int i=0; i<am.dataList.size(); ++i)
          {
            auto const& data = *am.dataList[i];
            if(verbose)
              std::cout << "data " << data.name() << " at " << data.location()->name()
                      << " with " << data.dim() << " components of type " << (int)data.primType() << "\n";
          }
        return false;
      }


      // Checks whether BoundaryTriangleData in the Amira mesh file is of type "byte" (true) or "int" (false).
      bool boundaryIdByByte(std::string const& gridFile) 
      {
        using namespace std;
        
        ifstream infile(gridFile);
        string line;
        while(getline(infile, line)) {
          istringstream iss(line);
          string firstWord;
          if((iss >> firstWord) && (firstWord == "BoundaryTriangleData")) {
            string bracket, type;
            if(iss >> bracket >> type) {
              if(type == "byte") return true;
              if(type == "int") return false;
            }
            throw Kaskade::LookupException("Wrong data type for BoundaryTriangleData. Has to be either byte or int.",__FILE__,__LINE__);
          }
        }
        return false;
      }

      template <class Scalar, class Grid, class ScalarEd>
      void readEdgePositions(Grid const& grid, Dune::GridFactory<Grid> const& factory, AmiraMesh & am, std::vector<Dune::FieldVector<ScalarEd, Grid::dimension>> & edgeDisplacement)
      {
        int const dim = Grid::dimension;
        std::vector<Dune::FieldVector<int, 2>> edges;
        std::vector<Dune::FieldVector<Scalar, dim>> edgePos;

        bool exists = AmiraMeshReader::ImplementationDetails::readData(am,"Edges","fromTo",edges);
        if (!exists)
        {
          throw Kaskade::LookupException("Could not read edges.",__FILE__,__LINE__);
        }
        exists = AmiraMeshReader::ImplementationDetails::readData(am,"Edges","Coordinates",edgePos);
        if (!exists)
        {
          throw Kaskade::LookupException("Could not read edge positions.",__FILE__,__LINE__);
        }

        edgeDisplacement.resize(edges.size());

        std::vector<std::pair<std::array<int, 2>, Dune::FieldVector<Scalar, dim>>> edgeTable (edges.size());
        for (int i=0; i<edges.size(); ++i)
        {
          edgeTable[i].first[0] = edges[i][0]-1;
          edgeTable[i].first[1] = edges[i][1]-1;
          if(edgeTable[i].first[0] > edgeTable[i].first[1]) std::swap(edgeTable[i].first[0], edgeTable[i].first[1]);

          edgeTable[i].second = edgePos[i];
        }
        std::sort(edgeTable.begin(), edgeTable.end(), [](auto const& a, auto const& b) {return a.first < b.first;});

        std::vector<bool> edgeRead(grid.size(dim-1), false);

        for (auto const& element : elements(grid.leafGridView()))
        {
          auto const& refElement = Dune::ReferenceElements<double, dim>::general(element.type());
          for (int j=0; j<refElement.size(dim-1); j++)
          {
            size_t edgeIndex = grid.leafGridView().indexSet().subIndex(element,j,dim-1);
            if (!edgeRead[edgeIndex])
            {
              auto c0 = element.template subEntity<dim>(refElement.subEntity(j, dim-1, 0, dim));
              auto c1 = element.template subEntity<dim>(refElement.subEntity(j, dim-1, 1, dim));

              std::array<int, 2> wantedEdge {factory.insertionIndex(c0), factory.insertionIndex(c1)};
              if(wantedEdge[0] > wantedEdge[1]) std::swap(wantedEdge[0], wantedEdge[1]);

              auto edge = std::lower_bound(edgeTable.begin(), edgeTable.end(), wantedEdge,
                                      [] (auto const& a, auto const& b) { return a.first < b;});

              if (edge == edgeTable.end() || edge->first != wantedEdge)
              {
                throw Kaskade::LookupException("Could not find edge.",__FILE__,__LINE__);
              }

              auto edgePosition = c0.geometry().corner(0);
              edgePosition += c1.geometry().corner(0);
              edgePosition /= 2.0;

              edgeDisplacement[edgeIndex] = edge->second;
              edgeDisplacement[edgeIndex] -= edgePosition;

              edgeRead[edgeIndex] = true;
            }
          }
        }
      }

    } // End of ImplementationDetails
    
    // --------------------------------------------------------------------------------------------
    // --------------------------------------------------------------------------------------------

    /// function reading an Amira mesh file in ASCII format.
    /**
     * This function reads vertices, tetrahedra and, if defined, boundary triangles
     * of an Amira mesh file into a grid and returns this grid.
     *
     * \tparam Grid the type of grid to be constructed
     * \tparam Scalar scalar type of the vertex coordinates, i.e. how it is specified in the amira file, 
     *                it is important to choose it correctly
     * \param gridfile the name of the Amira mesh file
     * \param initialGridSize optional parameter for reserving memory for the grid
     * \param boundaryIds if ids (e.g. membership to patches) for boundary segments are specified in the grid file, 
     *                    and this pointer is not null, then on return the pointed to vector will contain 
     *                    for each boundary segment (by its index, see Dune::Intersection) the corresponding id
     * \param edgeDisplacement if edges (by node indices) and their positions (of midpoints) are specified 
     *                         in the grid file and this pointer is not null, then on return the pointed to vector 
     *                         will contain for each edge (by its global index in the returned grid) 
     *                         the corresponding displacement
     * \param scale a factor with which the geometry is scaled. Use this to convert the geometry to standard SI units
     *              if it is given in the grid file in a different scale (e.g., millimeter). This is overwritten if 
     *              the grid file contains a UnitLength parameter.
     *  \param factory a Dune::GridFactory for the mesh to be created in. If non is passed, a new factory is created.
     *                 You can pass multiple grids to one single factor by calling this function more than one on different 
     *                 gridfiles. 
     * \return pointer on the created grid
     */
    template<class Grid, class Scalar, class ScalarEd = Scalar>
    std::unique_ptr<Grid> readGrid(std::string const& gridfile, int initialGridSize = 0,
                                   std::vector<int>* boundaryIds = nullptr, 
                                   std::vector<Dune::FieldVector<ScalarEd, Grid::dimension>>* edgeDisplacement = nullptr,
                                   double scale = 1.0)
    {
      // we apply no transformation, i.e. the identity to
      // each vertex of the grid
      auto identity = [] (auto x) {return x;};
      return readGrid<Grid,Scalar>(gridfile,identity,initialGridSize,boundaryIds,edgeDisplacement,scale);
    }

    /// function reading an Amira mesh file in ASCII format.
    /**
     * This function reads vertices, tetrahedra and, if defined, boundary triangles
     * of an Amira mesh file into a grid and returns this grid.
     *
     * \tparam Grid the type of grid to be constructed
     * \tparam Scalar scalar type of the vertex coordinates, i.e. how it is specified in the amira file, 
     *                it is important to choose it correctly
     * \param gridfile the name of the Amira mesh file
     * \param deformation performs a deformation/transformation to each vertex of the mesh to be read (typically a lambda function)
     * \param initialGridSize optional parameter for reserving memory for the grid
     * \param boundaryIds if ids (e.g. membership to patches) for boundary segments are specified in the grid file, 
     *                    and this pointer is not null, then on return the pointed to vector will contain 
     *                    for each boundary segment (by its index, see Dune::Intersection) the corresponding id
     * \param edgeDisplacement if edges (by node indices) and their positions (of midpoints) are specified 
     *                         in the grid file and this pointer is not null, then on return the pointed to vector 
     *                         will contain for each edge (by its global index in the returned grid) 
     *                         the corresponding displacement
     * \param scale a factor with which the geometry is scaled. Use this to convert the geometry to standard SI units
     *              if it is given in the grid file in a different scale (e.g., millimeter). This is overwritten if 
     *              the grid file contains a UnitLength parameter.
     *  \param factory a Dune::GridFactory for the mesh to be created in. If non is passed, a new factory is created.
     *                 You can pass multiple grids to one single factor by calling this function more than one on different 
     *                 gridfiles. 
     * \return pointer on the created grid
     */
    template<class Grid, class Scalar, class ScalarEd = Scalar, class Deformation>
    std::unique_ptr<Grid> readGrid(std::string const& gridfile, 
                                   Deformation const& deformation,
                                   int initialGridSize = 0,
                                   std::vector<int>* boundaryIds = nullptr, 
                                   std::vector<Dune::FieldVector<ScalarEd, Grid::dimension>>* edgeDisplacement = nullptr,
                                   double scale = 1.0)
    {
      int const dim = Grid::dimension;
      std::vector<Dune::FieldVector<Scalar,dim>>    vertices;
      std::vector<Dune::FieldVector<int,dim+1>>     cells;
      std::vector<Dune::FieldVector<int,dim>>       boundary;

      // open grid file.
      std::cout << "Opening " << gridfile << " ...\n";
      std::unique_ptr<AmiraMesh> am(AmiraMesh::read(gridfile.c_str()));
      if(!am)
        throw DetailedException(ImplementationDetails::exceptionMessage("readGrid(std::string const&, std::vector<int>&, int)", gridfile, __LINE__),__FILE__,__LINE__);

      // read vertices
      if(!ImplementationDetails::readData(*am,"Nodes","Coordinates",vertices))
        throw std::runtime_error(ImplementationDetails::exceptionMessage("readGrid(std::string const&, std::vector<int>&, int)", "Nodes", "Coordinates", __LINE__));

      // read cells
      std::string cellName;
      if(dim == 2) cellName = "Triangles";
      if(dim == 3) cellName = "Tetrahedra";
      if(dim != 2 && dim != 3)
        throw std::runtime_error("Sorry!\nOnly 2D and 3D grids are supported.");

      if(!ImplementationDetails::readData(*am,cellName.c_str(),"Nodes",cells))
        throw std::runtime_error(ImplementationDetails::exceptionMessage("readGrid(std::string const&, std::vector<int>&, int)", cellName, "Nodes", __LINE__));

      // read boundary segments if defined
      bool boundaryExists = ImplementationDetails::readData(*am,"BoundaryTriangles","Nodes",boundary);
      if(boundaryIds && !boundaryExists)
        throw Kaskade::LookupException("No boundary triangle nodes found in AmiraMesh file.",__FILE__,__LINE__);

      std::vector<Dune::FieldVector<int,1>> boundaryID;
      if (boundaryIds) 
      {
        bool exists;
        if(ImplementationDetails::boundaryIdByByte(gridfile)) {
          std::vector<Dune::FieldVector<unsigned char,1>> boundaryIDchar;
          exists = ImplementationDetails::readData(*am,"BoundaryTriangles","Id",boundaryIDchar);
          if (exists)
          {
            std::transform(begin(boundaryIDchar),end(boundaryIDchar),back_inserter(boundaryID),
                           [](auto id) { return Dune::FieldVector<int,1>(id[0]); });
          }
        } 
        else 
          exists = ImplementationDetails::readData(*am,"BoundaryTriangles","Id",boundaryID);
        
        if(!exists)
          throw Kaskade::LookupException("No boundary triangle ids found in AmiraMesh file [neither of type char nor int].",__FILE__,__LINE__);
      }
      
      // Look for length scale parameter.
      if (auto* unitLength = am->parameters.find("UnitLength"))
      {
        scale = unitLength->getReal();
        std::cout << "Found unit length of " << scale << "m.\n";
        if (scale <= 0)
          throw Kaskade::FileIOException("Nonpositive unit length specified.",gridfile,__FILE__,__LINE__);
      }
      

      // We use start indices 0. Amira uses 1. Correct indices.
      for (size_t i=0; i<cells.size(); ++i) 
      {
        for (int j=0; j<dim+1; ++j) 
        {
          cells[i][j] -= 1;
          if (cells[i][j]>=vertices.size())
          {
            std::string message = std::string("In ") + __FILE__ + " line " + boost::lexical_cast<std::string>(__LINE__) + ": cell node " + boost::lexical_cast<std::string>(i) + "[" + boost::lexical_cast<std::string>(j) + "]=" + boost::lexical_cast<std::string>(cells[i][j]) + " exceeds vertex numbers!\n";
            throw std::runtime_error(message);
          }
        }
      }

      if (boundaryExists)
      {
        for (size_t i=0; i<boundary.size(); ++i)
          for (int j=0; j<dim; ++j)
            boundary[i][j] -= 1;
      }


      // Create grid.
      Dune::GridFactory<Grid> factory = IOTools::FactoryGenerator<Grid>::createFactory(initialGridSize);
      std::cout << "inserting vertices...";
      for (size_t i=0; i<vertices.size(); ++i) 
      {
        Dune::FieldVector<double,dim> pos;
        for(int j=0; j<dim; ++j)
          pos[j] = scale*(double)vertices[i][j];
        factory.insertVertex(deformation(pos));   // apply deformation to each vertex
      }
      std::cout << "done.      ";

      if (boundaryExists && boundaryIds) 
      {
        std::cout << "inserting boundary segments...";
        for(size_t i=0; i<boundary.size(); ++i){
          std::vector<unsigned int> tmp(dim);
          for(size_t j=0; j<dim; ++j){
            assert(boundary[i][j]>=0 && boundary[i][j]<vertices.size());
            tmp[j] = boundary[i][j];
          }
          factory.insertBoundarySegment(tmp);
        }
        std::cout << "done.      ";
      }

      std::cout << "inserting elements...";
      for (size_t i=0; i<cells.size(); ++i) {
        std::vector<unsigned int> tmp(dim+1);
        for (int j=0; j<dim+1; ++j) {
          assert(cells[i][j]>=0 && cells[i][j]<vertices.size());
          tmp[j] = cells[i][j];
        }
        factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,dim),tmp);
      }
      std::cout << "done. " << std::endl;

      boost::timer::cpu_timer timer;
      std::unique_ptr<Grid> grid(factory.createGrid());
      std::cout << "grid creation finished in " << boost::timer::format(timer.elapsed()) << std::endl;

      if (boundaryIds) {
        // get boundary segment indices
        std::cout << "sorting boundary ids...";
        IOTools::BoundarySegmentIndexWrapper<Grid>::readBoundarySegmentIndices(*grid, factory, vertices, boundary, boundaryID, *boundaryIds);
        std::cout << "done." << std::endl;
      }

      if (edgeDisplacement)
        ImplementationDetails::readEdgePositions<Scalar>(*grid, factory, *am, *edgeDisplacement);

      return grid;
    }
    
    
    /// function joining two Amira mesh files in ASCII format into one grid.
    /**
     * This function reads vertices, tetrahedra and, if defined, boundary triangles
     * of each of two Amira mesh files. After renumbering the two are joined into
     * one grid.
     *
     * \tparam Grid the type of grid to be constructed
     * \tparam Scalar scalar type of the vertex coordinates, i.e. how it is specified in the amira file, 
     *                it is important to choose it correctly
     * \param gridfile1 the name of the first Amira mesh file
     * \param gridfile2 the name of the second Amira mesh file
     * \param initialGridSize optional parameter for reserving memory for the grid
     * \param boundaryIds if ids (e.g. membership to patches) for boundary segments are specified in the grid file, 
     *                    and this pointer is not null, then on return the pointed to vector will contain 
     *                    for each boundary segment (by its index, see Dune::Intersection) the corresponding id
     * \param edgeDisplacement if edges (by node indices) and their positions (of midpoints) are specified 
     *                         in the grid file and this pointer is not null, then on return the pointed to vector 
     *                         will contain for each edge (by its global index in the returned grid) 
     *                         the corresponding displacement
     * \param scale a factor with which the geometry is scaled. Use this to convert the geometry to standard SI units
     *              if it is given in the grid file in a different scale (e.g., millimeter). This is overwritten if 
     *              the grid file contains a UnitLength parameter.
     * 
     * 
     * \return pointer on the created grid
     */
    template<class Grid, class Scalar, class ScalarEd = Scalar>
    std::unique_ptr<Grid> mergeGrids(std::string const& gridfile1, std::string const& gridfile2,
                                     int initialGridSize = 0,
                                     std::vector<int>* boundaryIds = nullptr, 
                                     //  std::vector<Dune::FieldVector<ScalarEd, Grid::dimension>>* edgeDisplacement = nullptr,
                                     double scale1 = 1.0, double scale2 = 1.0)
    {
      int const dim = Grid::dimension;
      std::vector<Dune::FieldVector<Scalar,dim>>    vertices1, vertices2;
      std::vector<Dune::FieldVector<int,dim+1>>     cells1, cells2;
      std::vector<Dune::FieldVector<int,dim>>       boundary1, boundary2;

      // open and read the two files grid file.
      std::cout << "Reading first amira mesh\n";
      std::cout << "opening " << gridfile1 << '\n';
      std::unique_ptr<AmiraMesh> am1(AmiraMesh::read(gridfile1.c_str()));
      if(!am1)
        throw DetailedException(ImplementationDetails::exceptionMessage("readGrid(std::string const&, std::vector<int>&, int)", gridfile1, __LINE__),__FILE__,__LINE__);
      std::cout << "Reading second amira mesh\n";
      std::cout << "opening " << gridfile2 << '\n';
      std::unique_ptr<AmiraMesh> am2(AmiraMesh::read(gridfile2.c_str()));
      if(!am2)
        throw DetailedException(ImplementationDetails::exceptionMessage("readGrid(std::string const&, std::vector<int>&, int)", gridfile2, __LINE__),__FILE__,__LINE__);
      
      // read vertices in both Amira files
      // and store in vertices1 and vertices2
      if(!ImplementationDetails::readData(*am1,"Nodes","Coordinates",vertices1))
        throw std::runtime_error(ImplementationDetails::exceptionMessage("readGrid(std::string const&, std::vector<int>&, int)", "Nodes", "Coordinates", __LINE__));
      if(!ImplementationDetails::readData(*am2,"Nodes","Coordinates",vertices2))
        throw std::runtime_error(ImplementationDetails::exceptionMessage("readGrid(std::string const&, std::vector<int>&, int)", "Nodes", "Coordinates", __LINE__));
      
      // read cells in both Amira files
      // and store in cells1 and cells2
      std::string cellName;
      if(dim == 2) cellName = "Triangles";
      if(dim == 3) cellName = "Tetrahedra";
      if(dim != 2 && dim != 3)
        throw std::runtime_error("Sorry!\nOnly 2D and 3D grids are supported.");

      if(!ImplementationDetails::readData(*am1,cellName.c_str(),"Nodes",cells1))
        throw std::runtime_error(ImplementationDetails::exceptionMessage("readGrid(std::string const&, std::vector<int>&, int)", cellName, "Nodes", __LINE__));
      if(!ImplementationDetails::readData(*am2,cellName.c_str(),"Nodes",cells2))
        throw std::runtime_error(ImplementationDetails::exceptionMessage("readGrid(std::string const&, std::vector<int>&, int)", cellName, "Nodes", __LINE__));


      // read boundary segments if defined
      bool boundaryExists1 = ImplementationDetails::readData(*am1,"BoundaryTriangles","Nodes",boundary1);
      if(boundaryIds && !boundaryExists1)
        throw Kaskade::LookupException("No boundary triangle nodes found in first AmiraMesh file.",__FILE__,__LINE__);
      bool boundaryExists2 = ImplementationDetails::readData(*am2,"BoundaryTriangles","Nodes",boundary2);
      if(boundaryIds && !boundaryExists2)
        throw Kaskade::LookupException("No boundary triangle nodes found in second AmiraMesh file.",__FILE__,__LINE__);

      //  std::vector<Dune::FieldVector<int,1>> boundaryID;
      //  if (boundaryIds) 
      //  {
        //  bool exists;
        //  if(ImplementationDetails::boundaryIdByByte(gridfile)) {
          //  std::vector<Dune::FieldVector<unsigned char,1>> boundaryIDchar;
          //  exists = ImplementationDetails::readData(*am,"BoundaryTriangles","Id",boundaryIDchar);
          //  if (exists)
          //  {
            //  std::transform(begin(boundaryIDchar),end(boundaryIDchar),back_inserter(boundaryID),
                           //  [](auto id) { return Dune::FieldVector<int,1>(id[0]); });
          //  }
        //  } 
        //  else 
          //  exists = ImplementationDetails::readData(*am,"BoundaryTriangles","Id",boundaryID);
        
        //  if(!exists)
          //  throw Kaskade::LookupException("No boundary triangle ids found in AmiraMesh file [neither of type char nor int].",__FILE__,__LINE__);
      //  }
      
      // Look for length scale parameter in both files.
      if (auto* unitLength = am1->parameters.find("UnitLength"))
      {
        scale1 = unitLength->getReal();
        std::cout << "Found unit length of " << scale1 << "m in first Amira file.\n";
        if (scale1 <= 0)
          throw Kaskade::FileIOException("Nonpositive unit length specified.",gridfile1,__FILE__,__LINE__);
      }
      if (auto* unitLength = am2->parameters.find("UnitLength"))
      {
        scale2 = unitLength->getReal();
        std::cout << "Found unit length of " << scale2 << "m in second Amira file.\n";
        if (scale2 <= 0)
          throw Kaskade::FileIOException("Nonpositive unit length specified.",gridfile2,__FILE__,__LINE__);
      }
      

      // We use start indices 0. Amira uses 1. Correct indices.
      for (size_t i=0; i<cells1.size(); ++i) 
      {
        for (int j=0; j<dim+1; ++j) 
        {
          cells1[i][j] -= 1;
          if (cells1[i][j]>=vertices1.size())
          {
            std::string message = std::string("In ") + __FILE__ + " line " + boost::lexical_cast<std::string>(__LINE__) + ": cell node " + boost::lexical_cast<std::string>(i) + "[" + boost::lexical_cast<std::string>(j) + "]=" + boost::lexical_cast<std::string>(cells1[i][j]) + " exceeds vertex numbers!\n";
            throw std::runtime_error(message);
          }
        }
      }
      for (size_t i=0; i<cells2.size(); ++i) 
      {
        for (int j=0; j<dim+1; ++j) 
        {
          cells2[i][j] -= 1;
          if (cells2[i][j]>=vertices2.size())
          {
            std::string message = std::string("In ") + __FILE__ + " line " + boost::lexical_cast<std::string>(__LINE__) + ": cell node " + boost::lexical_cast<std::string>(i) + "[" + boost::lexical_cast<std::string>(j) + "]=" + boost::lexical_cast<std::string>(cells2[i][j]) + " exceeds vertex numbers!\n";
            throw std::runtime_error(message);
          }
        }
      }

      //  if (boundaryExists)
      //  {
        //  for (size_t i=0; i<boundary.size(); ++i)
          //  for (int j=0; j<dim; ++j)
            //  boundary[i][j] -= 1;
      //  }


      // Create grid.
      std::cout << std::endl;
      std::cout << "preparing grid creation\n";

      Dune::GridFactory<Grid> factory = IOTools::FactoryGenerator<Grid>::createFactory(initialGridSize);
      std::cout << "inserting vertices of first mesh...";
      for (size_t i=0; i<vertices1.size(); ++i) 
      {
        Dune::FieldVector<double,dim> pos;
        for(int j=0; j<dim; ++j)
          pos[j] = scale1*(double)vertices1[i][j];
        factory.insertVertex(pos);
      }
      std::cout << "done.\n";
      std::cout << "inserting vertices of second mesh...";
      for (size_t i=0; i<vertices2.size(); ++i) 
      {
        Dune::FieldVector<double,dim> pos;
        for(int j=0; j<dim; ++j)
          pos[j] = scale2*(double)vertices2[i][j];
        factory.insertVertex(pos);
      }
      std::cout << "done." << std::endl;

      //  if (boundaryExists && boundaryIds) 
      //  {
        //  std::cout << "inserting boundary segments...";
        //  for(size_t i=0; i<boundary.size(); ++i){
          //  std::vector<unsigned int> tmp(dim);
          //  for(size_t j=0; j<dim; ++j){
            //  assert(boundary[i][j]>=0 && boundary[i][j]<vertices.size());
            //  tmp[j] = boundary[i][j];
          //  }
          //  factory.insertBoundarySegment(tmp);
        //  }
        //  std::cout << "done." << std::endl;
      //  }

      std::cout << "inserting elements of first mesh...";
      for (size_t i=0; i<cells1.size(); ++i) {
        std::vector<unsigned int> tmp(dim+1);
        for (int j=0; j<dim+1; ++j) {
          assert(cells1[i][j]>=0 && cells1[i][j]<vertices1.size());
          tmp[j] = cells1[i][j];
        }
        factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,dim),tmp);
      }
      std::cout << "done.\n";
      std::cout << "inserting elements of second mesh...";
      for (size_t i=0; i<cells2.size(); ++i) {
        std::vector<unsigned int> tmp(dim+1);
        for (int j=0; j<dim+1; ++j) {
          assert(cells2[i][j]>=0 && cells2[i][j]<vertices2.size());
          tmp[j] = cells2[i][j] + vertices1.size();     // shift indices by number of vertices of the first mesh
        }
        factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,dim),tmp);
      }
      std::cout << "done." << std::endl;

      boost::timer::cpu_timer timer;
      std::unique_ptr<Grid> grid(factory.createGrid());
      std::cout << "grid creation finished in " << boost::timer::format(timer.elapsed()) << std::endl;

      //  if (boundaryIds) {
        //  // get boundary segment indices
        //  std::cout << "sorting boundary ids...";
        //  IOTools::BoundarySegmentIndexWrapper<Grid>::readBoundarySegmentIndices(*grid, factory, vertices, boundary, boundaryID, *boundaryIds);
        //  std::cout << "done." << std::endl;
      //  }

      //  if (edgeDisplacement)
        //  ImplementationDetails::readEdgePositions<Scalar>(*grid, factory, *am, *edgeDisplacement);

      return grid;
    }
    
    

    /// function reading additional data from an Amira mesh file in ASCII format.
    /**
     * \param DataType scalar type of the data in the Amira mesh file (e.g., unsigned char)
     * \param FunctionSpaceElement type of FE-Function to store the additional data in
     * \param datafile the name of the Amira mesh file containing additional data
     * \param dataname the location of the data in the Amira mesh file (e.g., "Tetrahedra")
     * \param componentname the name of the data's components in the Amira mesh file  (e.g., "Materials")
     * \param data FE function object for storing the additional data
     */
    template<class DataType, class FunctionSpaceElement>
    void readData(std::string const& datafile, std::string const& dataname, std::string const& componentname,
                  FunctionSpaceElement &data, bool verbose=false)
    {
      constexpr int numberOfComponents = FunctionSpaceElement::StorageValueType::dimension;
      // FunctionSpaceElement::StorageValueType in general uses floating point precision.
      // Thus, in order to read integer values or characters, the data vector is defined manually.
      // When reading the information into the function space element the data is
      // casted to the desired type.
      std::vector<Dune::FieldVector<DataType, numberOfComponents>> dataVector;

      // open data file
      std::cout << "Reading additional data." << std::endl;
      std::cout << "Opening " << datafile << " ...\n";
      std::unique_ptr<AmiraMesh> am(AmiraMesh::read(datafile.c_str()));

      if(!am) 
        throw FileIOException("could not open Amira mesh file for reading",datafile,__FILE__,__LINE__);


      bool exists = ImplementationDetails::readData(*am, dataname.c_str(), componentname.c_str(), dataVector);
      if(!exists)
        throw LookupException("data " + componentname + " not found at location " + dataname + " in Amira mesh file " + datafile,
                              __FILE__,__LINE__);

      if(verbose)
      {
        std::cout << "number of entries: " << dataVector.size() << std::endl;
        std::cout << "number of components: " << numberOfComponents << std::endl;
        std::cout << "data object: " << data.coefficients().size() << "x" << data.coefficients()[0].size() << std::endl;
      }

      // reading prescribed data into function space element
      for(size_t i=0; i<dataVector.size(); ++i)
        for(size_t j=0; j<numberOfComponents; ++j)
          data.coefficients()[i][j] = static_cast<typename FunctionSpaceElement::Scalar>(dataVector[i][j]);
          
      std::cout << "\n";
    }

    /// function reading additional data for each vertex from an Amira mesh file in ASCII format.
    /**
     * \param Scalar scalar type of the data in the Amira mesh file
     * \param numberOfComponents number of components of each entry of the data
     * \param datafile the name of the Amira mesh file containing additional node data
     * \param dataname the location of the data in the Amira mesh file (e.g., "Tetrahedra")
     * \param componentname the name of the data's components in the Amira mesh file (e.g., "Materials")
     * \param data vector containing the data
     */
    template<class Scalar, int numberOfComponents>
    void readData(std::string const& datafile, std::string const& dataname, std::string const& componentname,
                  std::vector<Dune::FieldVector<Scalar,numberOfComponents> > &data)
    {
      // open data file
      std::cout << std::endl;
      std::cout << "Reading additional vertex data" << std::endl;
      std::cout << "opening " << datafile << std::endl;
      AmiraMesh* am = AmiraMesh::read(datafile.c_str());
      if(!am)
        throw std::runtime_error(ImplementationDetails::exceptionMessage("readData(std::string const&, std::string const&, std::string const&, std::vector<Dune::FieldVector>&)", datafile, __LINE__));

      bool exists = ImplementationDetails::readData(*am, dataname.c_str(), componentname.c_str(), data);
      if(!exists)
      {
        delete(am);
        throw std::runtime_error(ImplementationDetails::exceptionMessage("readData(std::string const&, std::string const&, std::string const&, std::vector<Dune::FieldVector>&)", dataname, componentname, __LINE__));
      }

      // close data file
      delete(am);
    }

    /// function reading return a FE-function describing the deformation.
    /**
     * Function reading a deformed and an undeformed geometry from Amira mesh files
     * in ASCII format. These geometries must possess the same number of nodes. For each
     * vertex the difference between the geometries is calculated. This difference is stored
     * into the FE-function at index 0 of the variable set data.
     *
     * \param datafile the name of the Amira mesh file containing the deformed geometry
     * \param gridfile the name of the Amira mesh file containing the undeformed geometry
     * \param data fe-function object for storing the deformation of type VariableSet::VariableSet
     */
    template <class VarSetDesc>
    void readDeformationIntoVariableSetRepresentation(std::string const& gridfile, std::string const& datafile, 
						      typename VarSetDesc::VariableSet& data){
      using namespace boost::fusion;
      int const numberOfComponents = result_of::at_c<typename VarSetDesc::VariableSet::Functions, 0>::type::Components;
      typedef std::vector<Dune::FieldVector<float,numberOfComponents> > DataVector;
      std::vector<Dune::FieldVector<float,numberOfComponents> > vertices;
      DataVector nodeData;

      // open grid file.
      std::cout << "Reading undeformed geometry\n";
      std::cout << "opening " << gridfile << '\n';
      AmiraMesh* am = AmiraMesh::read(gridfile.c_str());
      if(!am) throw std::runtime_error(ImplementationDetails::exceptionMessage("readDeformationIntoVariableSetRepresentation(std::string const&, std::string const&,  VariableSet::VariableSet&)", gridfile, __LINE__));

      // read vertices
      bool exists = ImplementationDetails::readData(*am,"Nodes","Coordinates",vertices);
      if(!exists){
        delete(am);
        throw std::runtime_error(ImplementationDetails::exceptionMessage("readDeformationIntoVariableSetRepresentation(std::string const&, std::string const&,  VariableSet::VariableSet&)","Nodes","Coordinates",__LINE__));
      }

      // close grid file
      delete(am);

      // open data file
      std::cout << std::endl;
      std::cout << "Reading deformed geometry" << std::endl;
      std::cout << "opening " << datafile << std::endl;
      am = AmiraMesh::read(datafile.c_str());
      if(!am) throw std::runtime_error(ImplementationDetails::exceptionMessage("readDeformationIntoVariableSetRepresentation(std::string const&, std::string const&,  VariableSet::VariableSet&)", datafile, __LINE__));

      exists = ImplementationDetails::readData(*am,"Nodes","Coordinates",nodeData);
      if(!exists){
        delete(am);
        throw std::runtime_error("readDeformationIntoVariableSetRepresentation(std::string const&, std::string const&,  VariableSet::VariableSet&)", "Nodes", "Coordinates", __LINE__);
      }

      // close data file
      delete(am);

      // check for correct sizes
      assert(nodeData.size()==vertices.size());

      // get deformation
      for(size_t i=0; i<nodeData.size(); ++i)
        nodeData[i]=nodeData[i]-vertices[i];

      // reading prescribed data into function space
      std::vector<double> tmpVec(nodeData.size()*numberOfComponents);
      for(size_t i=0; i<nodeData.size(); ++i)
        for(size_t j=0; j<numberOfComponents; ++j)
          tmpVec[numberOfComponents*i+j] = (double)nodeData[i][j];

      data.read(tmpVec.begin());
    }

    /// function reading return a fe-function describing the deformation.
    /**
     * Function reading a deformed and a undeformed geometry from Amira mesh files
     * in ASCII format. These geometry must possess the same number of nodes. For each
     * vertex the difference between the geometries is calculated. This difference is read
     * into the fe-function object data.
     *
     * \param datafile the name of the Amira mesh file containing the deformed geometry
     * \param gridfile the name of the Amira mesh file containing the undeformed geometry
     * \param data fe-function object for storing the deformation
     */
    template<class FunctionSpaceElement, class FileScalar=float>
    void readDeformation(std::string const& gridfile, std::string const& datafile, FunctionSpaceElement &data) {
      int const numberOfComponents = FunctionSpaceElement::Components;
      typedef std::vector<Dune::FieldVector<FileScalar,numberOfComponents> > DataVector;
      typedef typename FunctionSpaceElement::Scalar Scalar;
      std::vector<Dune::FieldVector<FileScalar,numberOfComponents> > vertices;
      DataVector nodeData;

      // open grid file.
      std::cout << "Reading undeformed geometry\n";
      std::cout << "opening " << gridfile << '\n';
      AmiraMesh* am = AmiraMesh::read(gridfile.c_str());
      if(!am) throw std::runtime_error(ImplementationDetails::exceptionMessage("readDeformation(std::string const&, std::string const&, FunctionSpaceElement&)", gridfile, __LINE__));

      // read vertices
      bool exists = ImplementationDetails::readData(*am,"Nodes","Coordinates",vertices);
      if(!exists){
        delete(am);
        throw std::runtime_error(ImplementationDetails::exceptionMessage("readDeformation(std::string const&, std::string const&, FunctionSpaceElement&)", "Nodes", "Coordinates", __LINE__));
      }

      // close grid file
      delete(am);

      // open data file
      std::cout << std::endl;
      std::cout << "Reading deformed geometry" << std::endl;
      std::cout << "opening " << datafile << std::endl;
      am = AmiraMesh::read(datafile.c_str());
      if(!am) throw std::runtime_error(ImplementationDetails::exceptionMessage("readDeformation(std::string const&, std::string const&, FunctionSpaceElement&)", datafile, __LINE__));

      exists = ImplementationDetails::readData(*am,"Nodes","Coordinates",nodeData);
      if(!exists){
        delete(am);
        throw std::runtime_error(ImplementationDetails::exceptionMessage("readDeformation(std::string const&, std::string const&, FunctionSpaceElement&)", "Nodes", "Coordinates", __LINE__));
      }

      // close data file
      delete(am);

      // check for correct sizes
      assert(nodeData.size()==vertices.size());

      // get deformation
      for(size_t i=0; i<nodeData.size(); ++i)
        nodeData[i]=nodeData[i]-vertices[i];

      // reading prescribed data into function space
      std::vector<Scalar> tmpVec(nodeData.size()*numberOfComponents);
      for(size_t i=0; i<nodeData.size(); ++i)
        for(size_t j=0; j<FunctionSpaceElement::Components; ++j)
          data.coefficients()[i][j] = (Scalar) nodeData[i][j];
    }


    template<class GridView, class FunctionSpaceElement, class FileScalar=float>
    void readDeformation2(GridView const& gridView, std::string const& datafile, FunctionSpaceElement &data){
      int const numberOfComponents = FunctionSpaceElement::Components;
      typedef std::vector<Dune::FieldVector<FileScalar,numberOfComponents> > DataVector;
      typedef typename FunctionSpaceElement::Scalar Scalar;
      std::vector<Dune::FieldVector<FileScalar,numberOfComponents> > vertices;
      DataVector nodeData;

      // open data file
      std::cout << "Reading deformed geometry" << std::endl;
      std::cout << "opening " << datafile << std::endl;
      AmiraMesh *am = AmiraMesh::read(datafile.c_str());
      if(!am) throw std::runtime_error(ImplementationDetails::exceptionMessage("readDeformation(std::string const&, std::string const&, FunctionSpaceElement&)", datafile, __LINE__));

      bool exists = ImplementationDetails::readData(*am,"Nodes","Coordinates",nodeData);
      if(!exists){
        delete(am);
        throw std::runtime_error(ImplementationDetails::exceptionMessage("readDeformation(std::string const&, std::string const&, FunctionSpaceElement&)", "Nodes", "Coordinates", __LINE__));
      }

      // close data file
      delete(am);

      // check for correct sizes
      assert(nodeData.size()==vertices.size());

      // get deformation
      auto vIter = gridView.template begin<GridView::dimension>();
      auto vend = gridView.template end<GridView::dimension>();
      size_t i=0;
      for(;vIter!=vend;++vIter)
      {
        nodeData[i] -= vIter->geometry().corner(0);
        ++i;
      }
      //      for(size_t i=0; i<nodeData.size(); ++i)
//        nodeData[i]=nodeData[i]-vertices[i];

      // reading prescribed data into function space
      std::vector<Scalar> tmpVec(nodeData.size()*numberOfComponents);
      for(size_t i=0; i<nodeData.size(); ++i)
        for(size_t j=0; j<FunctionSpaceElement::Components; ++j)
          data.coefficients()[i][j] = (Scalar) nodeData[i][j];
    }


    /// Read boundary indices and save as scalar field in .vtu-file
    /**
     * \param gridfile name of the Amira mesh file
     * \param savefilename name of the output file (optional), default: "boundaryConditions.vtu"
     * \param initialGridSize initial grid size in mb(only for UGGrid)
     */
    template <class Grid>
    void boundaryConditionsToScalarField(std::string const& gridfile, std::string const& savefilename=std::string("boundaryConditions"), int initialGridSize = 0) {

      int const dim = Grid::dimension;
      // read grid
      typedef typename Grid::LeafGridView LeafView;
      std::vector<int> boundaryIndices;
      std::unique_ptr<Grid> grid = readGrid<Grid,float>(gridfile, boundaryIndices);
      GridManager<Grid> gridManager(grid);

      // create variable set
      typedef FEFunctionSpace<ContinuousLagrangeMapper<double,LeafView> > Space;
      Space space(gridManager, gridManager.grid().leafGridView(), 1);
      typedef boost::fusion::vector<Space const*> Spaces;
      Spaces spaces(&space);
      typedef boost::fusion::vector< VariableDescription<0,1,0> > VariableDescriptions;
      typedef VariableSetDescription<Spaces,VariableDescriptions> VarSetDesc;
      std::string name[1] = { "boundary ids" };
      VarSetDesc varSetDesc(spaces, name);
      typename VarSetDesc::VariableSet indices(varSetDesc);

      // read boundary indices
      std::vector<Dune::FieldVector<int,dim> > boundaryVertices;

      std::string const comp_name("BoundaryTriangles");
      std::string const node_name("Nodes");

      ImplementationDetails::readData<int,dim>(gridfile, comp_name, node_name, boundaryVertices);
      std::vector<double> tmp(gridManager.grid().size(dim),0);
      for(int i=0; i<boundaryIndices.size(); ++i){
        tmp[boundaryVertices[i][0]-1] = boundaryIndices[i];
        tmp[boundaryVertices[i][1]-1] = boundaryIndices[i];
        tmp[boundaryVertices[i][2]-1] = boundaryIndices[i];
      }
      indices.read(tmp.begin());

      // write file
      IoOptions options;
      options.outputType = IoOptions::ascii;
      writeVTKFile(gridManager.grid().leafGridView(), varSetDesc, indices, savefilename,options);
    }
  }
} // namespace Kaskade
#endif
