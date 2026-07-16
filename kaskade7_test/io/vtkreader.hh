/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2014-2021 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef VTKREADER_HH
#define VTKREADER_HH


#include <vector>
#include <map>
#include <string>
#include <memory>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <cassert>

#include <boost/fusion/container/vector.hpp>
#include <boost/fusion/include/zip_view.hpp>

#include <dune/grid/config.h>
#include <dune/grid/common/gridfactory.hh>

#include <fem/functionspace.hh>
#include <fem/lagrangespace.hh>
#include <utilities/detailed_exception.hh>

namespace Kaskade
{
  namespace VTKReaderDetail
  {
    enum CoefficientDataType
    {
      pointData,
      cellData
    };
    
    struct CoefficientData
    {
      CoefficientDataType type;
      std::size_t numComponents = 0;
      std::vector<double> data;
    };
    
    struct VTKData
    {
      bool reverseByteOrder = false;      // if loading binary
      std::string headerType = "UInt32";  // data type for binary data size header
      std::string compressor;             // compression method for binary data (not yet supported)
      std::string filename;
      std::size_t dimVertex = 0, numPoints = 0, numCells = 0;
      std::vector<double> coordinates;
      std::vector<std::size_t> connectivity, offsets, types;
      std::map<std::string, CoefficientData> coefficientData; // this includes PointData and CellData
    };
    
    struct VTKCell
    {
      int vtkGeomType;
      std::vector<std::size_t> corners, edges;
    };
    
  } // namespace VTKReaderDetail
  
  /**
   * \ingroup gridInput
   * @brief A class to create a grid and load coefficients from a VTK file.
   *
   * A low-level interface to load VTK files.
   * VTKReader will create a grid for you. From this grid you need to construct 
   * an appropriate space and variables yourself. Then you can copy the data
   * from the VTK file into your variables. Example:
   * \code
   *    VTKReader vtk("myfile.vtu");
   *    GridManager<Grid> gm(vtk.createGrid<Grid>());
   *    // create varSet
   *    vtk.getCoefficients(varSet);
   * \endcode
   * 
   * Currently supported: 2D and 3D simplical grids, Orders 0, 1, 2, continuous and discontinuous data.
   * DataArrays must be in ascii or binary format.
   */
  class VTKReader
  {
  public:
    /**
     * @brief Options for `getCoefficients()`.
     */
    enum Flags
    {
      noInterpolation = 0,
      allowOrderMismatch = 0x1,
      allowContinuityMismatch = 0x2,
      allowCoefficientDimensionMismatch = 0x4,
      
      allowEverything = allowOrderMismatch | allowContinuityMismatch | allowCoefficientDimensionMismatch
    };
    
    /**
     * @brief Information about a data set loaded from a file.
     */
    struct FunctionInfo
    {
      bool dataProvided = false, continuous = false;
      int order = 0, numComponents = 0;
    };
    
    /**
     * @brief Default constructor.
     */
    VTKReader() = default;
    
    /**
     * @brief Constructor that loads a VTK file.
     * @param filename The filename including the ending.
     */
    VTKReader(std::string const& filename);
    
    /**
     * @brief Loads a VTK file.
     * @param filename The filename including the ending.
     */
    void load(std::string const& filename);
    
    /**
     * @brief Loads only coefficient data, omitting the grid.
     * This can be used for time series.
     * @param filename The filename including the ending.
     */
    void loadOnlyData(std::string const& filename);
    
    /**
     * @brief Create a grid from a previously loaded VTK file.
     * @tparam Grid The type of grid to be created. The grid's dimension 
     * must match the dimension of shapes in the file. Note that the right 
     * dimension is not necessarily the dimension of points in the VTK file. 
     * Only cell types matter: triangles mean 2D, tetrahedrons mean 3D.
     * @return A `unique_ptr` that can be used to construct a `GridManager`.
     */
    template<class Grid>
    std::unique_ptr<Grid> createGrid();
    
    template<class Grid, class Deformation>
    std::unique_ptr<Grid> createGrid(Deformation const& deformation);
    
    /**
     * @brief Copy the data from the current VTK file into the coefficients of your FEFunctions.
     * 
     * The following heuristics are used to determine if coefficients can be copied:
     * + If the data is given as CellData, the function must be discontinuous, order 0.
     * + If the data is given as PointData, the type of function depends on the grid:
     *    + If there are multiple points with the same coordinates, the function must be discontinuous (otherwise continuous).
     *    + If Cell types are quadratic triangles or quadratic tetrahedrons, order must be 2 (otherwise order 1).
     * 
     * 
     * You can circumvent these restrictions by passing @param interpolationFlags.
     * If you are trying to load a file written by `writeVTK()`, specify order and set data mode
     * to inferred when saving the file.
     *
     * @param varSet A `VariableSet` to receive the data. Variables are assumed to live on the grid created by `createGrid()`. This is not checked.
     */
    template<class VariableSet>
    void getCoefficients(VariableSet& varSet, int interpolationFlags = noInterpolation) const;
    
    /**
     * @brief Copy to the coefficients of a single FSE.
     * @param functionName The name of the data set to be copied (as in the VTK file).
     */
    template<class FSE>
    void getCoefficients(std::string const& functionName, FSE& fse, int interpolationFlags = noInterpolation) const;
    
    /**
     * @brief Information about the expected type of function.
     * Meet these criteria to copy values into your variable without interpolation.
     * @param functionName The name of the data set as in the VTK file.
     * @return A `FunctionInfo` object that can be printed by `cout`. 
     */
    FunctionInfo getFunctionInfo(std::string const& functionName) const;
    
  private:
    
    enum LoadOption
    {
      loadOptionEverything,
      loadOptionOnlyData
    };
    
    enum CoefficientMapping
    {
      continuous,
      discontinuous
    };
    
    template<class VariableSet, class FSEVarDescPair>
    void copyCoefficientsByName(VariableSet const& varSet, FSEVarDescPair& pair, int interpolationFlags) const;
    
    template<class FSE>
    void copyCoefficients(std::string const& functionName, VTKReaderDetail::CoefficientData const& coefficientData, FSE& fse, int interpolationFlags) const;
    
    template<class FSE>
    void copyCoefficientsInterpolated(std::string const& functionName, VTKReaderDetail::CoefficientData const& coefficientData, FSE& fse) const;
    
    template<class FSE>
    void copyCoefficientsDirectly(std::string const& functionName, VTKReaderDetail::CoefficientData const& coefficientData, FSE& fse) const;
    
    void loadVTKInternal(std::string const& filename, LoadOption option);
    void prepareGridCreation(std::vector<VTKReaderDetail::VTKCell>& vtkCells, std::vector<char>& isCorner);
    void reset();
    void clearCoefficientData();
    
    
    VTKReaderDetail::VTKData vtk;   // stores global information about the VTK file
    bool gridWasLoaded = false, gridWasCreated = false;
    
    // coefficientMapping and polynomialOrder refer only to PointData. CellData is handled separately.
    CoefficientMapping coefficientMapping = continuous;
    int polynomialOrder = 0;
    
    // index mapping: for a coefficient-index i, dataIndices[i] is the corresponding vtk-index
    std::vector<std::size_t> pointDataIndices, cellDataIndices;
  };
  
  std::ostream& operator<<(std::ostream& stream, VTKReader::FunctionInfo const& info);
  
  /// \cond internals
  namespace VTKReaderDetail
  {
    // used to find duplicate points
    struct IndexedPoint
    {
      std::size_t vtkIdx = 0, numComponents = 0;
      std::vector<double>::const_iterator coordinates;
      
      bool operator<(IndexedPoint const& other) const;
      bool operator==(IndexedPoint const& other) const; // for std::adjacent_find
    };
    
    // used to keep track of indices during grid creation
    template<class Grid>
    class VTKFactoryHelper
    {
    public:
      template<class Deformation> void insertVertices(std::vector<char> const& isCorner, VTKData const& vtk, Deformation const& deformation);
      void insertVertices(std::vector<char> const& isCorner, VTKData const& vtk);
      void insertCells(std::vector<VTKCell> const& vtkCells);
      std::unique_ptr<Grid> createGrid();
      void setIndicesForDuplicatePoints(std::vector<IndexedPoint> const& indexedPoints, std::vector<std::size_t> const& vertexVTKToMerged);
      void getDataIndicesContinuous(typename Grid::LeafGridView const& leafView, std::vector<VTKCell> const& vtkCells, std::vector<std::size_t>& dataIndicesOut) const;
      void getDataIndicesDiscontinuous(typename Grid::LeafGridView const& leafView, std::vector<VTKCell> const& vtkCells, std::vector<std::size_t>& dataIndicesOut) const;
      void getCellDataIndices(typename Grid::LeafGridView const& leafView, std::vector<std::size_t>& dataIndicesOut) const;
      
    private:
    
      template<class ElementEntity>
      void getIndicesLocalToVTK(VTKCell const& vtkCell, ElementEntity const& element, std::vector<std::size_t>& vertexLocalToVTK, std::vector<std::size_t>& edgeLocalToVTK) const;
      
      Dune::GridFactory<Grid> factory;
      std::vector<std::size_t> vertexVTKToInsertion, vertexInsertionToVTK;
      std::size_t totalCorners = 0, totalEdges = 0;
    };
    
    std::string formatWarning(std::string const& filename, std::string const& message);
    
    Dune::GeometryType getDuneGeomType(int vtkGeomType);
    std::size_t getDuneCorner(int vtkGeomType, std::size_t corner);
    void countCornersAndEdges(int vtkGeomType, std::size_t& numCorners, std::size_t& numEdges);
    int countCellDim(int vtkGeomType);
    // compare index-pairs as sets
    struct CompareTwoSet
    {
      bool operator()(std::pair<std::size_t, std::size_t> const& p1, std::pair<std::size_t, std::size_t> const& p2) const;
    };
    using EdgeTable = std::map<std::pair<std::size_t, std::size_t>, std::size_t, CompareTwoSet>;
    
    std::size_t findEdgeIndex(EdgeTable const& edgeTable, std::size_t a, std::size_t b);
    std::pair<std::size_t, std::size_t> getVTKLocalCorners(int vtkGeomType, std::size_t edge);
    std::vector<std::size_t> getIndicesPerCell(std::vector<std::size_t> const& corners, std::vector<std::size_t> const& edges, int vtkGeomType);
    
    template<class ElementEntity, class IndexSet>
    EdgeTable createEdgeTable(ElementEntity const& element, IndexSet const& indexSet, std::vector<std::size_t> const& vertexGridToVTK)
    {
      // for an element, create a map {c0, c1} -> e.
      // {c0, c1} and e define the same edge. {c0, c1} are the vtk-indices of its corners and e is the grid-index of the edge-entity.
      
      constexpr int dim = IndexSet::dimension;
      constexpr int codimEdge = dim - 1, codimVertex = dim;
      EdgeTable edgeTable;
      auto const& refElement = Dune::ReferenceElements<double, dim>::general(element.type());
      
      unsigned numEdges = refElement.size(codimEdge);
      for(unsigned i = 0; i < numEdges; ++i)
      {
        // get local corner indices for local edge i
        std::size_t c0 = refElement.subEntity(i, codimEdge, 0, codimVertex);
        std::size_t c1 = refElement.subEntity(i, codimEdge, 1, codimVertex);
        
        // transform {c0, c1} and i to grid
        c0 = indexSet.subIndex(element, c0, codimVertex);
        c1 = indexSet.subIndex(element, c1, codimVertex);
        std::size_t e = indexSet.subIndex(element, i, codimEdge);
        
        // transform {c0, c1} to vtk
        c0 = vertexGridToVTK[c0];
        c1 = vertexGridToVTK[c1];
        
        // insert into edgeTable
        edgeTable[std::make_pair(c0, c1)] = e;
      }
      return edgeTable;
    }
    
    // match FSEs with their names
    template<class DataSequence, class VariableDescriptionSequence>
    boost::fusion::zip_view<boost::fusion::vector<DataSequence&, VariableDescriptionSequence const&>>
    makeZipView(DataSequence& data, VariableDescriptionSequence const& descriptions)
    {
      using Sequences = boost::fusion::vector<DataSequence&, VariableDescriptionSequence const&>;
      return boost::fusion::zip_view<Sequences>(Sequences(data, descriptions));
    }
    
    //                VTKFactoryHelper
    
    template<class Grid>
    template<class Deformation>
    void VTKFactoryHelper<Grid>::insertVertices(std::vector<char> const& isCorner, VTKData const& vtk, Deformation const& deformation)
    {
      // Insert vertices from VTK's coordinate vector into grid and keep track of indices.
      // If too many or too few coordinates are provided, excess values are ignored or set to 0 respectively.
      constexpr int dim = Grid::dimension;
      using Vertex = Dune::FieldVector<double, dim>;
      std::size_t numComponentsToRead = std::min((std::size_t)dim, vtk.dimVertex);
      vertexVTKToInsertion.resize(vtk.numPoints, 0);
      std::size_t countInsertion = 0;
      
      for(std::size_t i = 0; i < vtk.numPoints; ++i)
      {
        // only insert marked corners
        if(isCorner[i])
        {
          std::size_t offs = i * vtk.dimVertex;
          Vertex vertex(0);
          auto it = vtk.coordinates.begin() + offs;
          std::copy(it, it + numComponentsToRead, vertex.begin());
          factory.insertVertex(deformation(vertex));
          
          vertexInsertionToVTK.push_back(i);
          vertexVTKToInsertion[i] = countInsertion;
          ++countInsertion;
        }
      }
      totalCorners = countInsertion;
      // totalEdges is used only for continuous coefficient mapping
      // otherwise the number is wrong due to duplicate points
      totalEdges = vtk.numPoints - totalCorners;
    }
    
    template<class Grid>
    void VTKFactoryHelper<Grid>::insertVertices(std::vector<char> const& isCorner, VTKData const& vtk)
    {
      auto identity = [] (auto x) {return x;};
      insertVertices(isCorner, vtk, identity);
    }
    
    template<class Grid>
    void VTKFactoryHelper<Grid>::insertCells(std::vector<VTKCell> const& vtkCells)
    {
      int dim = Grid::dimension;
      // insert cells using the insertion-indices of their corners
      for(VTKCell const& vtkCell : vtkCells)
      {
        int geoType = vtkCell.vtkGeomType;
        int cellDim = countCellDim(geoType);
        
        if(dim !=cellDim){
          continue;
        }

        std::vector<unsigned int> corners(vtkCell.corners.size());
        for(std::size_t k = 0; k < corners.size(); ++k)
        {
          std::size_t corner = getDuneCorner(vtkCell.vtkGeomType, k);
          std::size_t vtkCornerIdx = vtkCell.corners[corner];
          corners[k] = (unsigned int) vertexVTKToInsertion[vtkCornerIdx];
        }
        factory.insertElement(getDuneGeomType(vtkCell.vtkGeomType), corners);
      }
    }
    
    template<class Grid>
    std::unique_ptr<Grid> VTKFactoryHelper<Grid>::createGrid()
    {
      return std::unique_ptr<Grid>(factory.createGrid());
    }
    
    template<class Grid>
    void VTKFactoryHelper<Grid>::setIndicesForDuplicatePoints(std::vector<IndexedPoint> const& indexedPoints, std::vector<std::size_t> const& vertexVTKToMerged)
    {
      // duplicate points get the insertion index of their representative
      for(IndexedPoint const& p : indexedPoints)
      {
        std::size_t mergedIdx = vertexVTKToMerged[p.vtkIdx];
        std::size_t insertionIdx = vertexVTKToInsertion[mergedIdx];
        vertexVTKToInsertion[p.vtkIdx] = insertionIdx;
      }
    }
    
    template<class Grid>
    void VTKFactoryHelper<Grid>::getDataIndicesContinuous(typename Grid::LeafGridView const& leafView, std::vector<VTKCell> const& vtkCells, std::vector<std::size_t>& dataIndicesOut) const
    {
      std::vector<std::size_t> vertexGridToVTK(totalCorners), edgeGridToVTK(totalEdges);
      
      // vertices
      auto const& indexSet = leafView.indexSet();
      for(auto const& vertexEntity : vertices(leafView))
      {
        std::size_t insertionIdx = factory.insertionIndex(vertexEntity);
        std::size_t gridIdx = indexSet.index(vertexEntity);

        if(insertionIdx >= vertexInsertionToVTK.size() || gridIdx >= vertexGridToVTK.size())
        {
          // Grid returned invalid index?
          throw Kaskade::GridException("Vertex index is out of bounds.", __FILE__, __LINE__);
        }
        
        std::size_t vtkIdx = vertexInsertionToVTK[insertionIdx];
        vertexGridToVTK[gridIdx] = vtkIdx;
      }
      
      // edges
      for(auto const& elementEntity : elements(leafView))
      {
        // get the vtk-definition used to create this cell
        std::size_t insertionIdx = factory.insertionIndex(elementEntity);
        // out of bounds check
        if(insertionIdx >= vtkCells.size())
        {
          // Grid returned invalid index?
          throw Kaskade::GridException("Cell insertion index is out of bounds.", __FILE__, __LINE__);
        }
        
        VTKCell const& vtkCell = vtkCells[insertionIdx];
        
        // the EdgeTable for this cell maps {c0, c1} -> e,
        // where c0, c1 are the vtk-indices of two of its corners and e is the grid-index of the edge between c0 and c1
        EdgeTable edgeTable = createEdgeTable(elementEntity, indexSet, vertexGridToVTK);
        for(std::size_t i = 0; i < vtkCell.edges.size(); ++i)
        {
          // get the vtk-index of an edge
          std::size_t vtkIdx = vtkCell.edges[i];
          
          // get vtk's local indices of the corresponding corners
          std::pair<int, int> vtkLocalCorners = getVTKLocalCorners(vtkCell.vtkGeomType, i);
          
          // use the local indices to look up the vtk-indices of these corners
          // finally find the grid-index of the edge via the EdgeTable
          std::size_t gridIdx = findEdgeIndex(edgeTable, vtkCell.corners[vtkLocalCorners.first], vtkCell.corners[vtkLocalCorners.second]);
          
          if(gridIdx >= edgeGridToVTK.size())
          {
            // Grid returned invalid index?
            throw Kaskade::GridException("Edge index is out of bounds.", __FILE__, __LINE__);
          }
          edgeGridToVTK[gridIdx]= vtkIdx;
        }
      }
      
      dataIndicesOut.resize(totalCorners + totalEdges);
      std::copy(edgeGridToVTK.begin(), edgeGridToVTK.end(), dataIndicesOut.begin());
      std::copy(vertexGridToVTK.begin(), vertexGridToVTK.end(), dataIndicesOut.begin() + totalEdges);
    }
    
    template<class Grid>
    void VTKFactoryHelper<Grid>::getDataIndicesDiscontinuous(typename Grid::LeafGridView const& leafView, std::vector<VTKCell> const& vtkCells, std::vector<std::size_t>& dataIndicesOut) const
    {
      // per cell indices are first stored as blocks
      std::vector<std::vector<std::size_t>> dataIndicesPerCell(vtkCells.size());
      
      for(auto const& elementEntity : elements(leafView))
      {
        // get the corresponding vtk-cell
        std::size_t insertionIdx = factory.insertionIndex(elementEntity);
        if(insertionIdx >= vtkCells.size())
        {
          // Grid returned invalid index?
          throw Kaskade::GridException("Cell insertion index is out of bounds.", __FILE__, __LINE__);
        }
        VTKCell const& vtkCell = vtkCells[insertionIdx];
        
        // get the data-indices for this cell in Kaskade discontinuous order
        std::vector<std::size_t> vertexLocalToVTK, edgeLocalToVTK;
        getIndicesLocalToVTK(vtkCell, elementEntity, vertexLocalToVTK, edgeLocalToVTK);
        std::vector<std::size_t> cellIndices = getIndicesPerCell(vertexLocalToVTK, edgeLocalToVTK, vtkCell.vtkGeomType);
        
        // insert in the right place
        std::size_t gridIdx = leafView.indexSet().index(elementEntity);
        // out of bounds check
        if(gridIdx >= dataIndicesPerCell.size())
        {
          // Grid returned invalid index?
          throw Kaskade::GridException("Cell index is out of bounds.", __FILE__, __LINE__);
        }
        dataIndicesPerCell[gridIdx] = cellIndices;
      }
      
      // make a sequence out of the blocks
      for(std::vector<std::size_t> const& cellIndices : dataIndicesPerCell)
      {
        for(std::size_t idx : cellIndices)
        {
          dataIndicesOut.push_back(idx);
        }
      }
    }
    
    template<class Grid>
    void VTKFactoryHelper<Grid>::getCellDataIndices(typename Grid::LeafGridView const& leafView, std::vector<std::size_t>& dataIndicesOut) const
    {
      std::size_t numCells = leafView.size(0);
      dataIndicesOut.resize(numCells);
      
      for(auto const& elementEntity : elements(leafView))
      {
        // cells were inserted in vtk order, so vtk-index == insertion-index
        std::size_t vtkIdx = factory.insertionIndex(elementEntity);
        std::size_t gridIdx = leafView.indexSet().index(elementEntity);
        if(gridIdx >= dataIndicesOut.size())
        {
          // Grid returned invalid index?
          throw Kaskade::GridException("Cell index is out of bounds.", __FILE__, __LINE__);
        }
        dataIndicesOut[gridIdx] = vtkIdx;
      }
    }
    
    template<class Grid>
    template<class ElementEntity>
    void VTKFactoryHelper<Grid>::getIndicesLocalToVTK(VTKCell const& vtkCell, ElementEntity const& element, std::vector<std::size_t>& vertexLocalToVTK, std::vector<std::size_t>& edgeLocalToVTK) const
    {
      // this function takes one grid-element and maps its local indices (corners and edges) to the right vtk data indices
      // edges are identified via the insertion-indices of their corners
      // vtkCell and element must define the same cell for this to work
      
      using IndexPair = std::pair<std::size_t, std::size_t>;
      
      // abusing std::map, the desired mapping will be in the values
      std::map<std::size_t, IndexPair> cornerLookup;              // c(ins) ---> (c(grid_local), c(vtk-data))
      std::map<IndexPair, IndexPair, CompareTwoSet> edgeLookup;   // {a, b}(ins) ---> (e(grid_local), e(vtk-data))
      
      std::size_t numCorners = vtkCell.corners.size(), numEdges = vtkCell.edges.size();
      constexpr int dim = Grid::dimension;
      constexpr int codimVertex = dim, codimEdge = dim - 1;
      
      assert(element.subEntities(codimVertex) == numCorners);
      assert(numEdges == 0 || element.subEntities(codimEdge) == numEdges);
      
      // corners
      std::vector<std::size_t> vertexLocalToInsertion(numCorners);
      
      for(std::size_t i = 0; i < numCorners; ++i)
      {
        // grid        
        auto vertexEntity = element.template subEntity<codimVertex>(i);
        std::size_t insertionIndex = factory.insertionIndex(vertexEntity);
        cornerLookup[insertionIndex].first = i;
        vertexLocalToInsertion[i] = insertionIndex;
        
        // vtk
        std::size_t vtkIndex = vtkCell.corners[i];
        insertionIndex = vertexVTKToInsertion[vtkIndex];
        cornerLookup[insertionIndex].second = vtkIndex;
      }
      
      // edges
      auto const& refElement = Dune::ReferenceElements<double, dim>::general(element.type());
      
      for(std::size_t i = 0; i < numEdges; ++i)
      {
        // grid
        std::size_t c0 = refElement.subEntity(i, codimEdge, 0, codimVertex); // local
        std::size_t c1 = refElement.subEntity(i, codimEdge, 1, codimVertex); // local
        c0 = vertexLocalToInsertion[c0];  // insertion
        c1 = vertexLocalToInsertion[c1];  // insertion
        edgeLookup[std::make_pair(c0, c1)].first = i; // local
        
        // vtk
        IndexPair p = getVTKLocalCorners(vtkCell.vtkGeomType, i); // local
        c0 = vtkCell.corners[p.first];  // vtk
        c1 = vtkCell.corners[p.second]; // vtk
        c0 = vertexVTKToInsertion[c0]; // insertion
        c1 = vertexVTKToInsertion[c1]; // insertion
        std::size_t edgeDataIdx = vtkCell.edges[i]; // vtk
        edgeLookup[std::make_pair(c0, c1)].second = edgeDataIdx;
      }
      
      assert(cornerLookup.size() == numCorners);
      assert(edgeLookup.size() == numEdges);
      
      // fill vectors to return
      vertexLocalToVTK.resize(numCorners);
      
      for(auto const& p : cornerLookup)
      {
        IndexPair const& ip = p.second;
        vertexLocalToVTK[ip.first] = ip.second;
      }
      
      edgeLocalToVTK.resize(numEdges);
      
      for(auto const& p : edgeLookup)
      {
        IndexPair const& ip = p.second;
        edgeLocalToVTK[ip.first] = ip.second;
      }
    }
    
  } // namespace VTKReaderDetail
  /// \endcond
  
  //            VTKReader
  
  
  template<class Grid>
  std::unique_ptr<Grid> VTKReader::createGrid()
  {
    // identity means grid vertices are not transformed
    auto identity = [] (auto x) {return x;};
    return createGrid<Grid>(identity);
  }
  
  template<class Grid, class Deformation>
  std::unique_ptr<Grid> VTKReader::createGrid(Deformation const& deformation)
  {
    using namespace VTKReaderDetail;
    
    if(!gridWasLoaded)
    {
      throw Kaskade::FileIOException("No grid data has been loaded.", vtk.filename, __FILE__, __LINE__);
    }
    
    std::vector<char> isCorner(vtk.numPoints, 0);
    std::vector<VTKCell> vtkCells(vtk.numCells);
    
    // reset indices, group cells, mark corners, set polynomialOrder
    prepareGridCreation(vtkCells, isCorner);
    
    std::vector<IndexedPoint> indexedPoints;
    std::vector<std::size_t> vertexVTKToMerged;
    
    // collect non-edge points
    for(std::size_t i = 0; i < vtk.numPoints; ++i)
    {
      if(isCorner[i])
      {
        std::size_t offs = i * vtk.dimVertex;
        IndexedPoint p;
        p.vtkIdx = i;
        p.numComponents = vtk.dimVertex;
        p.coordinates = vtk.coordinates.begin() + offs;
        indexedPoints.push_back(p);
      }
    }
    
    std::sort(indexedPoints.begin(), indexedPoints.end());
    
    // detect duplicates and set coefficientMapping.
    // coefficientMapping is set to discontinuous if there is at least one duplicate point.
    bool hasDuplicatePoints = std::adjacent_find(indexedPoints.begin(), indexedPoints.end()) != indexedPoints.end();
    coefficientMapping = hasDuplicatePoints? discontinuous : continuous;
    
    if(hasDuplicatePoints)
    {
      // The isCorner vector will be used to prevent insertion of duplicate points.
      vertexVTKToMerged.resize(vtk.numPoints, 0);
      // unmark all corners
      std::fill(isCorner.begin(), isCorner.end(), 0);
      
      auto P = indexedPoints.begin();
      while(P != indexedPoints.end())
      {
        // mark P as corner
        isCorner[P->vtkIdx] = 1;
        // map duplicates' vtk-indices to P's vtk-index
        std::size_t currentRepresentativeIdx = P->vtkIdx;
        auto range = std::equal_range(P, indexedPoints.end(), *P);
        for(auto it = range.first; it != range.second; ++it)
        {
          vertexVTKToMerged[it->vtkIdx] = currentRepresentativeIdx;
        }
        P = range.second;
      }
    }
    
    VTKFactoryHelper<Grid> factoryHelper;
    std::unique_ptr<Grid> grid;
    
    // create the grid
    try
    {
      factoryHelper.insertVertices(isCorner, vtk, deformation);
      
      // If points have been merged, map vtk-indices of duplicate points to the same insertion index.
      if(coefficientMapping == discontinuous)
      {
        factoryHelper.setIndicesForDuplicatePoints(indexedPoints, vertexVTKToMerged);
      }
      
      factoryHelper.insertCells(vtkCells);
      grid = factoryHelper.createGrid();
    }
    catch(Dune::GridError const& e)
    {
      throw Kaskade::FileIOException(e.what(), vtk.filename, __FILE__, __LINE__);
    }
    
    // Get the indices. We always obtain index-mappings for PointData and CellData.
    auto leafView = grid->leafGridView();
    
    if(coefficientMapping == continuous)
    {
      factoryHelper.getDataIndicesContinuous(leafView, vtkCells, pointDataIndices);
    }
    else if(coefficientMapping == discontinuous)
    {
      factoryHelper.getDataIndicesDiscontinuous(leafView, vtkCells, pointDataIndices);
    }
    
    factoryHelper.getCellDataIndices(leafView, cellDataIndices);
    
    gridWasCreated = true;
    return grid;
  }
    
  template<class VariableSet>
  void VTKReader::getCoefficients(VariableSet& varSet, int interpolationFlags) const
  {
    using namespace VTKReaderDetail;
    
    if(!gridWasCreated)
    {
      throw Kaskade::FileIOException("Grid must be created before accessing data.", vtk.filename, __FILE__, __LINE__);
    }
    
    // Loop through FSEs and their descriptions simultaneously to match input data by name.
    auto functions = makeZipView(varSet.data, typename VariableSet::Descriptions::Variables());
    boost::fusion::for_each(functions, [&](auto const& f){this->copyCoefficientsByName(varSet, f, interpolationFlags);});
  }
  
  template<class FSE>
  void VTKReader::getCoefficients(std::string const& functionName, FSE& fse, int interpolationFlags) const
  {
    using namespace VTKReaderDetail;
    
    if(!gridWasCreated)
    {
      throw Kaskade::FileIOException("Grid must be created before accessing data.", vtk.filename, __FILE__, __LINE__);
    }
    
    auto it = vtk.coefficientData.find(functionName);
    if(it != vtk.coefficientData.end())
    {
      copyCoefficients(functionName, it->second, fse, interpolationFlags);
    }
    else
    {
      // warning
      std::ostringstream out;
      out << "No coefficient data provided for " << functionName << ".";
      std::cout << formatWarning(vtk.filename, out.str()) << std::endl;
    }
  }
    
  template<class VariableSet, class FSEVarDescPair>
  void VTKReader::copyCoefficientsByName(VariableSet const& varSet, FSEVarDescPair& pair, int interpolationFlags) const
  {
    using namespace boost::fusion;
    using namespace VTKReaderDetail;
    
    auto& fse = at_c<0>(pair);
    auto const& desc = at_c<1>(pair);
    std::string functionName = varSet.descriptions.names[desc.id];
    
    // find data-values associated with fse's name
    auto it = vtk.coefficientData.find(functionName);
    if(it != vtk.coefficientData.end())
    {
      copyCoefficients(functionName, it->second, fse, interpolationFlags);
    }
    else
    {
      // warning
      std::ostringstream out;
      out << "No data provided for " << functionName << ".";
      std::cout << formatWarning(vtk.filename, out.str()) << std::endl;
    }
  }
    
  template<class FSE>
  void VTKReader::copyCoefficients(std::string const& functionName, VTKReaderDetail::CoefficientData const& coefficientData, FSE& fse, int interpolationFlags) const
  {
    using namespace VTKReaderDetail;
    
    constexpr int fseCoefficientDim = FSE::StorageValueType::dimension;
    int fseMaxOrder = fse.space().mapper().maxOrder(), fseContinuity = FSE::Space::continuity;
    
    bool interpolationNeeded = false;
    
    if((std::size_t)fseCoefficientDim != coefficientData.numComponents)
    {
      if((interpolationFlags & allowCoefficientDimensionMismatch) == 0)
      {
        std::ostringstream out;
        out << "\"" << functionName << "\": Coefficient dimension mismatch.\nfunction info: " << getFunctionInfo(functionName);
        throw Kaskade::FileIOException(out.str(), vtk.filename, __FILE__, __LINE__);
      }
    }
    
    // order is always 0 for CellData
    int order = coefficientData.type == cellData? 0 : polynomialOrder;
    
    if(fseMaxOrder != order)
    {
      if(interpolationFlags & allowOrderMismatch)
      {
        interpolationNeeded = true;
      }
      else
      {
        std::ostringstream out;
        out << "\"" << functionName << "\": Order mismatch.\nfunction info: " << getFunctionInfo(functionName);
        throw Kaskade::FileIOException(out.str(), vtk.filename, __FILE__, __LINE__);
      }
    }
    
    // continuity is always discontinuous for CellData
    CoefficientMapping continuity = coefficientData.type == cellData? discontinuous : coefficientMapping;
    
    if(
    (fseContinuity < 0 && continuity == continuous)
    || (fseContinuity >= 0 && continuity == discontinuous))
    {
      if(interpolationFlags & allowContinuityMismatch)
      {
        interpolationNeeded = true;
      }
      else
      {
        std::ostringstream out;
        out << "\"" << functionName << "\": Continuity mismatch.\nfunction info: " << getFunctionInfo(functionName);
        throw Kaskade::FileIOException(out.str(), vtk.filename, __FILE__, __LINE__);
      }
    }
    
    if(interpolationNeeded)
    {
      copyCoefficientsInterpolated(functionName, coefficientData, fse);
    }
    else
    {
      copyCoefficientsDirectly(functionName, coefficientData, fse);
    }
  } 
    
  template<class FSE>
  void VTKReader::copyCoefficientsInterpolated(std::string const& functionName, VTKReaderDetail::CoefficientData const& coefficientData, FSE& fse) const
  {
    using namespace Kaskade;
    using namespace VTKReaderDetail;
    
    using Mapper = typename FSE::Space::Mapper;
    using Scalar = typename FSE::Space::Scalar;
    using LeafView = typename Mapper::Grid::LeafGridView;
    constexpr int fseCoefficientDim = FSE::StorageValueType::dimension;
    
    if(coefficientData.type == pointData)
    {
      if(coefficientMapping == continuous)
      {
        using Space = FEFunctionSpace<ContinuousLagrangeMapper<Scalar, LeafView>>;
          
        Space space(fse.space().gridManager(), fse.space().grid().leafGridView(), polynomialOrder);
        auto otherFSE = space.template element<fseCoefficientDim>();
        copyCoefficientsDirectly(functionName, coefficientData, otherFSE);
        fse = otherFSE;
      }
      else if(coefficientMapping == discontinuous)
      {
        using Space = FEFunctionSpace<DiscontinuousLagrangeMapper<Scalar, LeafView>>;
          
        Space space(fse.space().gridManager(), fse.space().grid().leafGridView(), polynomialOrder);
        auto otherFSE = space.template element<fseCoefficientDim>();
        copyCoefficientsDirectly(functionName, coefficientData, otherFSE);
        fse = otherFSE;
      }
    }
    else if(coefficientData.type == cellData)
    {
      // CellData means discontinuous, order 0
      using Space = FEFunctionSpace<DiscontinuousLagrangeMapper<Scalar, LeafView>>;
      
      Space space(fse.space().gridManager(), fse.space().grid().leafGridView(), 0);
      auto otherFSE = space.template element<fseCoefficientDim>();
      copyCoefficientsDirectly(functionName, coefficientData, otherFSE);
      fse = otherFSE;
    }
  }
  
  template<class FSE>
  void VTKReader::copyCoefficientsDirectly(std::string const& functionName, VTKReaderDetail::CoefficientData const& coefficientData, FSE& fse) const
  {
    using namespace VTKReaderDetail;
    
    constexpr int fseCoefficientDim = FSE::StorageValueType::dimension;
    std::size_t numComponentsToCopy = std::min(std::size_t(fseCoefficientDim), coefficientData.numComponents);
    std::vector<std::size_t> const& indices = coefficientData.type == pointData? pointDataIndices : cellDataIndices;
    
    if(indices.size() != fse.coefficients().size())
    {
      std::ostringstream out;
      out << "\"" << functionName << "\": Number of coefficients doesn't match number of data points.\nfunction info: " << getFunctionInfo(functionName);
      throw Kaskade::FileIOException(out.str(), vtk.filename, __FILE__, __LINE__);
    }
    for(std::size_t i = 0; i < indices.size(); ++i)
    {
      std::size_t dataIdx = indices[i];
      auto& block = fse.coefficients()[i];
      auto itData = coefficientData.data.begin() + dataIdx * coefficientData.numComponents;
      std::copy(itData, itData + numComponentsToCopy, block.begin());
    }
  }
} // namespace Kaskade

#endif
