/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2014-2022 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "vtkreader.hh"

#include <stdexcept>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <io/check_endianness.hh>

namespace Kaskade
{
  namespace VTKReaderDetail
  {
    // internal exception
    class DataArrayError : public std::runtime_error
    {
    public:
      DataArrayError(std::string const& message) : std::runtime_error(message){}
    };
    
    namespace DecodeBase64
    {
      using UByte = unsigned char;
      
      static UByte const
      mask0 = 0x3f,   // 00111111
      mask1 = 0x30,   // 00110000
      mask2 = 0xf,    // 00001111
      mask3 = 0x3c,   // 00111100
      mask4 = 0x3;    // 00000011
      
      UByte base64CharToByte(UByte c)
      {
        if ('A' <= c && c <= 'Z')
          return c - 'A';
        else if ('a' <= c && c <= 'z')
          return 26 + (c - 'a');
        else if ('0' <= c && c <= '9')
          return 52 + (c - '0');
        else if (c == '+')
          return 62;
        else if (c == '/')
          return 63;
        return 0;
      }
      
      struct isNotBase64
      {
        bool operator()(UByte c)
        {
          return !(
          ('A' <= c && c <= 'Z')
          || ('a' <= c && c <= 'z')
          || ('0' <= c && c <= '9')
          || (c == '+')
          || (c == '/')
          || (c == '='));
        }
      };

      std::vector<UByte> bytesToBlock(UByte b0, UByte b1, UByte b2, UByte b3, std::size_t size)
      {
        std::vector<UByte> block(size);
        
        block[0] = ((b0 & mask0) << 2) | ((b1 & mask1) >> 4);
        if (size > 1)
        {
          block[1] = ((b1 & mask2) << 4) | ((b2 & mask3) >> 2);
        }
        if (size > 2)
        {
          block[2] = ((b2 & mask4) << 6) | (b3 & mask0);
        }
        return block;
      }

      std::vector<UByte> decodeBase64(UByte a, UByte b, UByte c, UByte d)
      {
        std::size_t blockSize = 3;
        
        if (d == '=')
        {
          blockSize = 2;
          if (c == '=')
          {
            blockSize = 1;
          }
        }
        
        a = base64CharToByte(a);
        b = base64CharToByte(b);
        c = base64CharToByte(c);
        d = base64CharToByte(d);
        
        return bytesToBlock(a, b, c, d, blockSize);
      }
      
      void base64StringToBytes(std::string const& src, std::vector<UByte>& bytes)
      {
        if (src.size() % 4 != 0)
          throw DataArrayError("Number of Base64 characters (" + std::to_string(src.size()) 
                               + ") is not divisible by 4.");
        
        for(std::size_t i=0; i<src.size(); i+=4)
        {
          std::vector<UByte> block = decodeBase64(src[i], src[i + 1], src[i + 2], src[i + 3]);
          for(std::size_t j=0; j<block.size(); ++j)
            bytes.push_back(block[j]);
        }
      }
      
      template<class T>
      void readDataFromByteVector(std::vector<UByte> const& bytes, std::vector<T>& vecOut, 
                                  bool reverseByteOrder, const T* upperBound)
      {
        constexpr std::size_t sizeOfT = sizeof(T);
        
        if (bytes.size() % sizeOfT != 0)
        {
          throw DataArrayError("Wrong number of bytes (" + std::to_string(bytes.size())
                               + ") for the given datatype of size " + std::to_string(sizeOfT) + ".");
        }
        
        if (reverseByteOrder)
        {
          for(std::size_t i = 0; i < bytes.size(); i += sizeOfT)
          {
            T t;
            UByte* bytePtr = reinterpret_cast<UByte*>(&t);
            for(std::size_t j = 0; j < sizeOfT; ++j)
            {
              bytePtr[sizeOfT - (j + 1)] = bytes[i + j];
            }
            if (upperBound && t >= *upperBound)
            {
              std::ostringstream out;
              out << "Value " << t << " is out of bounds (has to be less than " << *upperBound << ").";
              throw DataArrayError(out.str());
            }
            vecOut.push_back(t);
          }
        }
        else
        {
          for(auto it = bytes.begin(); it != bytes.end(); it += sizeOfT)
          {
            T t;
            UByte* bytePtr = reinterpret_cast<UByte*>(&t);
            std::copy(it, it + sizeOfT, bytePtr);
            vecOut.push_back(t);
          }
        }
      }
    } // namespace DecodeBase64
    
    
    void loadVTKFromPTree(VTKData& vtk, boost::property_tree::ptree const& root, bool onlyData = false);
    void loadPointsAndCells(VTKData& vtk, boost::property_tree::ptree const& points, boost::property_tree::ptree const& cells);
    void loadCoefficientData(VTKData& vtk, boost::property_tree::ptree const& data, CoefficientDataType type);
    bool vtkIsConsistent(VTKData const& vtk, std::ostream& problems);
    
    std::string formatWarning(std::string const& filename, std::string const& message)
    {
      std::ostringstream out;
      out << "VTKReader (\"" << filename << "\") Warning: " << message;
      return out.str();
    }
    
    bool IndexedPoint::operator<(IndexedPoint const& other) const
    {
      for(std::size_t i = 0; i < numComponents; ++i)
      {
        if (coordinates[i] < other.coordinates[i])
          return true;
        else if (other.coordinates[i] < coordinates[i])
          return false;
      }
      return false;
    }
    
    bool IndexedPoint::operator==(IndexedPoint const& other) const
    {
      return !(*this < other || other < *this);
    }
    
    // --------------------------------------------------------------------------------------------

    template<class T>
    void readDataFromAscii(std::string const& src, std::vector<T>& dest, const T* upperBound)
    {
      std::istringstream stream(src);
      while(stream)
      {
        if (stream.eof())
        {
          break;
        }
        T t;
        stream >> t;
        if (stream.fail())
        {
          throw DataArrayError("Value extraction failed.");
        }
        if (upperBound && t >= *upperBound)
        {
          std::ostringstream out;
          out << "Value " << t << " is out of bounds (has to be less than " << *upperBound << ").";
          throw DataArrayError(out.str());
        }
          
        dest.push_back(t);
      }
    }
    
    template<class RawType, class DestType>
    void readDataFromBinary(std::vector<DecodeBase64::UByte> const& bytes, std::vector<DestType>& dest, 
                            bool reverseByteOrder, const DestType* upperBound)
    {
      using namespace DecodeBase64;
      
      std::vector<RawType> raw;
      
      if (upperBound)
      {
        RawType ub = *upperBound;
        readDataFromByteVector(bytes, raw, reverseByteOrder, &ub);
      }
      else
      {
        readDataFromByteVector(bytes, raw, reverseByteOrder, (RawType*) 0);
      }
      
      // Copy to destination, thereby performing type conversion if necessary.
      dest.resize(raw.size());
      std::copy(raw.begin(), raw.end(), dest.begin());
    }
    
    // --------------------------------------------------------------------------------------------
    
    template<class T>
    void readDataArray(boost::property_tree::ptree const& dataArray, std::vector<T>& dest, 
                       bool reverseByteOrder, std::string const& headerType, const T* upperBound = 0)
    {
      std::string format = dataArray.get<std::string>("<xmlattr>.format");
      if (format == "ascii")
      {
        readDataFromAscii(dataArray.get_value<std::string>(), dest, upperBound);
      }
      else if (format == "binary")
      {
        using namespace DecodeBase64;
        
        std::string type = dataArray.get<std::string>("<xmlattr>.type");
        std::string rawStr = dataArray.get_value<std::string>();
        // remove non-base64 characters
        rawStr.erase(std::remove_if (rawStr.begin(), rawStr.end(), isNotBase64()), rawStr.end());
        
        // VTK binary base64 is structured as |header|data| where header is a UInt32 or UInt64 (since file format 1.0),
        // and gives the number of data bytes before base64 encoding. data follows 
        // immediately after the header. No white space is expected within header or data, as VTK reads the data 
        // apparently without parsing as a plain block of data. The header shall start immediately after a newline 
        // following the closing '>' of the DataArray tag.
        // See http://public.kitware.com/pipermail/paraview/2005-April/001391.html and 
        // https://www.paraview.org/Wiki/VTK_XML_Formats
        
        // We don't make use of the header information yet, and simply skip it. The header may or may not be
        // base64-encoded separately (probably separately encoded in file format 0.1 and jointly in 1.0).
        // Separate base64-encoding introduces at least one "=" character at position 7 or 11, for UInt32
        // and UInt64 header size, respectively.
        int headerChars = headerType=="UInt32"? 8: 12;    // header separately encoded in this many characters
        int headerBytes = headerType=="UInt32"? 4: 8;     // header jointly encoded in this many bytes
        
        std::vector<DecodeBase64::UByte> data;
        if (rawStr[headerChars-1]=='=')       // base64 end marker found => separately encoded 
        {
          rawStr.erase(0,headerChars);        // skip header
          base64StringToBytes(rawStr,data);   // then decode the rest
        }
        else                                  // jointly encoded in first 
        {
          base64StringToBytes(rawStr,data);   // decode everything
          data.erase(data.begin(),data.begin()+headerBytes);  // then skip the decoded header
        }
        // TODO: use header for consistency check
        
        if (type == "Int8")
        {
          readDataFromBinary<char, T>(data, dest, reverseByteOrder, upperBound);
        }
        else if (type == "UInt8")
        {
          readDataFromBinary<unsigned char, T>(data, dest, reverseByteOrder, upperBound);
        }
        else if (type == "Int16")
        {
          readDataFromBinary<short, T>(data, dest, reverseByteOrder, upperBound);
        }
        else if (type == "UInt16")
        {
          readDataFromBinary<unsigned short, T>(data, dest, reverseByteOrder, upperBound);
        }
        else if (type == "Int32")
        {
          readDataFromBinary<int, T>(data, dest, reverseByteOrder, upperBound);
        }
        else if (type == "UInt32")
        {
          readDataFromBinary<std::uint32_t, T>(data, dest, reverseByteOrder, upperBound);
        }
        else if (type == "Float32")
        {
          readDataFromBinary<float, T>(data, dest, reverseByteOrder, upperBound);
        }
        else if (type == "Int64")
        {
          readDataFromBinary<std::int64_t, T>(data, dest, reverseByteOrder, upperBound);
        }
        else if (type == "UInt64")
        {
          readDataFromBinary<std::uint64_t, T>(data, dest, reverseByteOrder, upperBound);
        }
        else if (type == "Float64")
        {
          readDataFromBinary<double, T>(data, dest, reverseByteOrder, upperBound);
        }
        else
        {
          std::ostringstream out;
          out << "Datatype not supported (" << type << ").";
          throw DataArrayError(out.str());
        }
      }
      else
      {
        std::ostringstream out;
        out << "Format not supported: " << format << ".";
        throw DataArrayError(out.str());
      }
    }
    
    // copy-paste from vtkwriter.hh
    enum VTKGeometryType
    {
      vtkLine = 3,
      vtkTriangle = 5,
      vtkQuadrilateral = 9,
      vtkTetrahedron = 10,
      vtkHexahedron = 12,
      vtkPrism = 13,
      vtkPyramid = 14,
      vtkquadLine = 21,
      vtkquadTriangle = 22,
      vtkquadQuadrilateral = 23,
      vtkquadTetrahedron = 24,
      vtkquadHexahedron = 25,
    };
    
    bool geometryTypeIsSupported(int vtkGeomType)
    {
      switch(vtkGeomType)
      { 
        case vtkLine:
        case vtkTriangle:
        case vtkquadTriangle:
        case vtkTetrahedron:
        case vtkquadTetrahedron:
          return true;
        
        default:
          break;
      }
      return false;
    }
    
    Dune::GeometryType getDuneGeomType(int vtkGeomType)
    {
      assert(geometryTypeIsSupported(vtkGeomType));
      
      Dune::GeometryType gt;
      
      switch(vtkGeomType)
      {
        case vtkLine:
          break;
        case vtkquadTriangle:
        case vtkTriangle:
          gt.makeTriangle();
          break;
        case vtkquadTetrahedron:
        case vtkTetrahedron:
          gt.makeTetrahedron();
          break;
          
        default:
          break;
      }
      return gt;
    }
    
    // Reorder corners of VTK cell definition when inserting into factory.
    // (unnecessary for simplices, but quadrilaterals etc would need such a function)
    std::size_t getDuneCorner(int vtkGeomType, std::size_t corner)
    {
      assert(geometryTypeIsSupported(vtkGeomType));
      
      static std::vector<std::size_t>
      triangleCorners = {0, 1, 2},
      tetrahedronCorners = {0, 1, 3, 2};
      
      switch(vtkGeomType)
      {
        case vtkLine:
          break;
        case vtkquadTriangle:
        case vtkTriangle:
          return triangleCorners[corner];
        case vtkquadTetrahedron:
        case vtkTetrahedron:
          return tetrahedronCorners[corner];
        
        default:
          break;
      }
      
      return 0;
    }
    
    void countCornersAndEdges(int vtkGeomType, std::size_t& numCorners, std::size_t& numEdges)
    {
      assert(geometryTypeIsSupported(vtkGeomType));
      
      switch(vtkGeomType)
      {
        case vtkLine:
          numCorners = 2;
          numEdges = 0;
          break;
        case vtkTriangle:
          numCorners = 3;
          numEdges = 0;
          break;
        case vtkquadTriangle:
          numCorners = 3;
          numEdges = 3;
          break;
        case vtkTetrahedron:
          numCorners = 4;
          numEdges = 0;
          break;
        case vtkquadTetrahedron:
          numCorners = 4;
          numEdges = 6;
          break;
        default:
          numCorners = 0;
          numEdges = 0;
          break;
      }
    }
    
    int countCellDim(int vtkGeomType)
    {
      assert(geometryTypeIsSupported(vtkGeomType));
      switch(vtkGeomType)
      {
        case vtkLine:
          return 1;
        case vtkTriangle:
          return 2;
        case vtkquadTriangle:
          return 2;
        case vtkTetrahedron:
          return 3;
        case vtkquadTetrahedron:
          return 3;
        default:
          return 0;
      }
      return 0;
    }

    std::pair<std::size_t, std::size_t> getVTKLocalCorners(int vtkGeomType, std::size_t edge)
    {
      assert(geometryTypeIsSupported(vtkGeomType));
      
      // map local edge-indices to pairs of local corner indices (according to vtk file-format)
      static std::pair<std::size_t, std::size_t>
      quadTriangleEdges[] = {{0, 1}, {1, 2}, {0, 2}},
      quadTetrahedronEdges[] = {{0, 1}, {1, 2}, {0, 2}, {0, 3}, {1, 3}, {2, 3}};
      
      switch(vtkGeomType)
      {
        case vtkLine:
          break;
        case vtkquadTriangle:
          return quadTriangleEdges[edge];
        case vtkquadTetrahedron:
          return quadTetrahedronEdges[edge];
        
        default:
          break;
      }
      return std::pair<std::size_t, std::size_t>(0, 0);
    }
    
    std::vector<std::size_t> getIndicesPerCell(std::vector<std::size_t> const& corners, std::vector<std::size_t> const& edges, int vtkGeomType)
    {
      assert(geometryTypeIsSupported(vtkGeomType));
      
      std::vector<std::size_t> const& c = corners;
      std::vector<std::size_t> const& e = edges;
      
      if (vtkGeomType == vtkTriangle)
      {
        std::vector<std::size_t> indices = {c[0], c[2], c[1]};
        return indices;
      }
      else if (vtkGeomType == vtkquadTriangle)
      {
        std::vector<std::size_t> indices = {c[0], e[1], c[2], e[0], e[2], c[1]};
        return indices;
      }
      else if (vtkGeomType == vtkTetrahedron)
      {
        std::vector<std::size_t> indices = {c[0], c[3], c[2], c[1]};
        return indices;
      }
      else if (vtkGeomType == vtkquadTetrahedron)
      {
        std::vector<std::size_t> indices = {c[0], e[3], c[3], e[1], e[5], c[2], e[0], e[4], e[2], c[1]};
        return indices;
      }
      
      return std::vector<std::size_t>();
    }
    
    std::pair<std::size_t, std::size_t> sortIndexPair(std::pair<std::size_t, std::size_t> const& p)
    {
      return p.first < p.second ? p : std::make_pair(p.second, p.first);
    }
    
    bool CompareTwoSet::operator ()(std::pair<std::size_t, std::size_t> const& p1, std::pair<std::size_t, std::size_t> const& p2) const
    {
      std::pair<std::size_t, std::size_t> twoSet1 = sortIndexPair(p1), twoSet2 = sortIndexPair(p2);
      return twoSet1.first == twoSet2.first ? twoSet1.second < twoSet2.second : twoSet1.first < twoSet2.first;
    }
    
    std::size_t findEdgeIndex(EdgeTable const& edgeTable, std::size_t a, std::size_t b)
    {
      auto it = edgeTable.find(std::make_pair(a, b));
      assert(it != edgeTable.end());
      return it->second;
    }
    
    // --------------------------------------------------------------------------------------------
    
    // This is the main reader, that actually interprets the XML file structure.
    void loadVTKFromPTree(VTKData& vtk, boost::property_tree::ptree const& root, bool onlyData)
    {
      using namespace boost::property_tree;
      using namespace VTKReaderDetail;
      using OptTree = boost::optional<ptree const&>;
      
      ptree const& vtkFile = root.get_child("VTKFile");
      std::string vtkFileType = vtkFile.get<std::string>("<xmlattr>.type");
      std::string vtkByteOrder = vtkFile.get<std::string>("<xmlattr>.byte_order");
      
      if (!(littleEndian() || bigEndian()))
        throw Kaskade::FileIOException("Could not determine own endianness.", vtk.filename, __FILE__, __LINE__);
      
      assert(littleEndian() != bigEndian());
      
      if (vtkByteOrder == "BigEndian")
        vtk.reverseByteOrder = littleEndian();
      else if (vtkByteOrder == "LittleEndian")
        vtk.reverseByteOrder = bigEndian();
      else
      {
        std::ostringstream out;
        out << "Unknown byte order: " << vtkByteOrder;
        throw Kaskade::FileIOException(out.str(), vtk.filename, __FILE__, __LINE__);
      }
      
      // File format 1.0 provides an additional attribute specifying the binary header type. This
      // defaults to UInt32. See https://www.paraview.org/Wiki/VTK_XML_Formats
      vtk.headerType = vtkFile.get("<xmlattr>.header_type","UInt32");
      
      
      if (vtkFileType == "UnstructuredGrid")
      {
        ptree const& unstructuredGrid = vtkFile.get_child("UnstructuredGrid");
        std::size_t numPieces = unstructuredGrid.count("Piece");
        if (numPieces > 1)
        {
          // warning
          std::cout << formatWarning(vtk.filename, "Only one Piece will be loaded.") << std::endl;
        }
        
        ptree const& piece = unstructuredGrid.get_child("Piece");
        
        if (!onlyData)
        {
          vtk.numPoints = piece.get<std::size_t>("<xmlattr>.NumberOfPoints");
          vtk.numCells = piece.get<std::size_t>("<xmlattr>.NumberOfCells");
          loadPointsAndCells(vtk, piece.get_child("Points"), piece.get_child("Cells"));
        }
        
        // To points an cells there may be data, i.e. coefficients, be provided.
        // Read them if available.
        OptTree pData = piece.get_child_optional("PointData");
        if (pData)
        {
          loadCoefficientData(vtk, *pData, pointData);
        }
        
        OptTree cData = piece.get_child_optional("CellData");
        if (cData)
        {
          loadCoefficientData(vtk, *cData, cellData);
        }
      }
      else
      {
        std::ostringstream out;
        out << "Cannot load VTK-file of type " << vtkFileType << ".";
        throw Kaskade::FileIOException(out.str(), vtk.filename, __FILE__, __LINE__);
      }
    }
    
    // --------------------------------------------------------------------------------------------

    void loadPointsAndCells(VTKData& vtk, boost::property_tree::ptree const& points, 
                            boost::property_tree::ptree const& cells)
    {
      using namespace boost::property_tree;
      
      // read coordinates
      ptree const& coordinates = points.get_child("DataArray");
      vtk.dimVertex = coordinates.get<std::size_t>("<xmlattr>.NumberOfComponents");

      try
      {
        readDataArray(coordinates, vtk.coordinates, vtk.reverseByteOrder, vtk.headerType);
      }
      catch(DataArrayError const& e)
      {
        std::ostringstream out;
        out << "In coordinates: " << e.what();
        throw Kaskade::FileIOException(out.str(), vtk.filename, __FILE__, __LINE__);
      }
      
      // read connectivity, offsets, types
      for(auto const& child: cells)
      {
        if (child.first == "DataArray")
        {        
          std::string name = child.second.get<std::string>("<xmlattr>.Name");
          try
          {
            if (name == "connectivity")
              readDataArray(child.second, vtk.connectivity, vtk.reverseByteOrder, vtk.headerType, &vtk.numPoints);
            else if (name == "offsets")
            {
              std::size_t upperBound = vtk.connectivity.size() + 1;
              readDataArray(child.second, vtk.offsets, vtk.reverseByteOrder, vtk.headerType, &upperBound);
            }
            else if (name == "types")
              readDataArray(child.second, vtk.types, vtk.reverseByteOrder, vtk.headerType);
          }
          catch (DataArrayError const& e)
          {
            std::ostringstream out;
            out << "In " << name << ": " << e.what();
            throw Kaskade::FileIOException(out.str(), vtk.filename, __FILE__, __LINE__);
          }
        }
      }
    }
    
    // --------------------------------------------------------------------------------------------

    void loadCoefficientData(VTKData& vtk, boost::property_tree::ptree const& data, CoefficientDataType type)
    {
      using namespace VTKReaderDetail;
      
      for (auto const& child: data)
      {
        if (child.first == "DataArray")
        {
          std::string name = child.second.get<std::string>("<xmlattr>.Name");

          // Check if coefficients with same name have already been read.
          auto it = vtk.coefficientData.find(name);
          if (it != vtk.coefficientData.end())
          {
            std::cout << formatWarning(vtk.filename, "More than one function is named \"" + name + "\".") << std::endl;
          }
          else
          {
            
            CoefficientData& cd = vtk.coefficientData[name];
            cd.type = type;
            // For reading coefficients we default to scalars, if no number of components is given.
            cd.numComponents = child.second.get("<xmlattr>.NumberOfComponents",1);
            try
            {
              readDataArray(child.second, cd.data, vtk.reverseByteOrder, vtk.headerType);
            }
            catch(DataArrayError const& e)
            {
              std::ostringstream out;
              out << "In " << name << ": " << e.what();
              throw Kaskade::FileIOException(out.str(), vtk.filename, __FILE__, __LINE__);
            }
          }
        }
      }
    }
    
    // --------------------------------------------------------------------------------------------

    bool vtkIsConsistent(VTKData const& vtk, std::ostream& problems)
    {
      bool vtkOK = true;
      
      std::size_t coordinatesExpected = vtk.dimVertex * vtk.numPoints;
      if (vtk.coordinates.size() != coordinatesExpected)
      {
        problems << "Expected " << coordinatesExpected << " coordinates, got " << vtk.coordinates.size() << ".\n";
        vtkOK = false;
      }
      if (vtk.offsets.size() != vtk.numCells)
      {
        problems << "Expected " << vtk.numCells << " offsets, got " << vtk.offsets.size() << ".\n";
        vtkOK = false;
      }
      if (vtk.types.size() != vtk.numCells)
      {
        problems << "Expected " << vtk.numCells << " types, got " << vtk.types.size() << ".\n";
        vtkOK = false;
      }
      // connectivity not handled here
      
      for(auto const& p : vtk.coefficientData)
      {
        CoefficientData const& cd = p.second;
        std::size_t valuesExpected = 0;
        if (cd.type == pointData)
        {
          valuesExpected = cd.numComponents * vtk.numPoints;
        }
        else if (cd.type == cellData)
        {
          valuesExpected = cd.numComponents * vtk.numCells;
        }
        
        if (cd.data.size() != valuesExpected)
        {
          problems << "\n" << (cd.type == pointData? "PointData" : "CellData")  << " \"" << p.first << "\": Expected " << valuesExpected << " values, got " << cd.data.size() << ".\n";
          vtkOK = false;
        }
      }
      return vtkOK;
    }
    
  } // namespace VTKReaderDetail
  
  std::ostream& operator<<(std::ostream& stream, VTKReader::FunctionInfo const& info)
  {
    if (info.dataProvided)
    {
      stream << (info.continuous? "continuous" : "discontinuous") << ", order = " << info.order << ", number of components = " << info.numComponents;
    }
    else
    {
      stream << "no data provided";
    }
    return stream;
  }
  
  //            VTKReader
    
  VTKReader::VTKReader(std::string const& filename)
  {
    load(filename);
  }
  
  void VTKReader::load(std::string const& filename)
  {
    try
    {
      loadVTKInternal(filename, loadOptionEverything);
    }
    catch(...)
    {
      reset();
      throw;
    }
    gridWasLoaded = true;
  }
  
  void VTKReader::loadOnlyData(std::string const& filename)
  {
    try
    {
      loadVTKInternal(filename, loadOptionOnlyData);
    }
    catch(...)
    {
      clearCoefficientData();
      throw;
    }
    gridWasLoaded = false;
  }
  
  VTKReader::FunctionInfo VTKReader::getFunctionInfo(std::string const& functionName) const
  {
    using namespace VTKReaderDetail;
    
    FunctionInfo info;
    
    auto it = vtk.coefficientData.find(functionName);
    if (it != vtk.coefficientData.end())
    {
      CoefficientData const& cd = it->second;
      info.dataProvided = true;
      info.numComponents = cd.numComponents;
      if (cd.type == cellData)
      {
        info.order = 0;
        info.continuous = false;
      }
      else if (cd.type == pointData)
      {
        info.order = polynomialOrder;
        info.continuous = (coefficientMapping == continuous);
      }
    }
    return info;
  }
  
  void VTKReader::loadVTKInternal(std::string const& filename, LoadOption option)
  {
    using namespace boost::property_tree;
    
    ptree root;
    int flags = xml_parser::trim_whitespace;
    
    try
    {
      read_xml(filename, root, flags);
    }
    catch(xml_parser::xml_parser_error const& e)
    {
      std::ostringstream out;
      if (e.line() > 0)
      {
        out << "XML error (line " << e.line() << "): ";
      }
      out << e.message();
      throw Kaskade::FileIOException(out.str(), filename, __FILE__, __LINE__);
    }
    
    if (option == loadOptionEverything)
      reset();
    else if (option == loadOptionOnlyData)
      clearCoefficientData();
    
    vtk.filename = filename;
    try
    {
      loadVTKFromPTree(vtk, root, option == loadOptionOnlyData? true : false);
    }
    catch(ptree_bad_data const& e)
    {
      std::ostringstream out;
      out << "Conversion of data failed for \"" << e.data<std::string>() << "\".";
      throw Kaskade::FileIOException(out.str(), vtk.filename, __FILE__, __LINE__);
    }
    catch(ptree_error const& e)
    {
      throw Kaskade::FileIOException(e.what(), vtk.filename, __FILE__, __LINE__);
    }
    
    std::ostringstream problems;
    if (!vtkIsConsistent(vtk, problems))
    {
      throw Kaskade::FileIOException(problems.str(), vtk.filename, __FILE__, __LINE__);
    }
  }
  
  void VTKReader::prepareGridCreation(std::vector<VTKReaderDetail::VTKCell>& vtkCells, std::vector<char>& isCorner)
  {
    using namespace VTKReaderDetail;
    
    pointDataIndices.clear();
    cellDataIndices.clear();
    gridWasCreated = false;
    polynomialOrder = 0;
    
    std::size_t cellStart = 0, numCorners = 0, numEdges = 0;
    for(std::size_t i = 0; i < vtk.numCells; ++i)
    {
      int vtkGeomType = vtk.types[i];
      
      if (!geometryTypeIsSupported(vtkGeomType))
      {
        std::ostringstream out;
        out << "Cell type " << vtkGeomType << " is not supported.";
        throw Kaskade::FileIOException(out.str(), vtk.filename, __FILE__, __LINE__);
      }
      
      countCornersAndEdges(vtkGeomType, numCorners, numEdges);
      
      std::string noMixOrders = "Mixing polynomial orders is not supported.";
      
      if (numEdges == 0)
      {
        if (polynomialOrder == 2)
        {
          throw Kaskade::FileIOException(noMixOrders, vtk.filename, __FILE__, __LINE__);
        }
        polynomialOrder = 1;
      }
      else if (numEdges > 0)
      {
        if (polynomialOrder == 1)
        {
          throw Kaskade::FileIOException(noMixOrders, vtk.filename, __FILE__, __LINE__);
        }
        polynomialOrder = 2;
      }
      
      std::size_t offsCorners = cellStart, offsEdges = cellStart + numCorners;
      
      // out-of-bounds check
      if (cellStart + numCorners + numEdges > vtk.connectivity.size())
      {
        std::ostringstream out;
        out << "Connectivity error in cell " << i << ".";
        throw Kaskade::FileIOException(out.str(), vtk.filename, __FILE__, __LINE__);
      }
      
      VTKCell& vtkCell = vtkCells[i];
      vtkCell.vtkGeomType = vtkGeomType;
      vtkCell.corners.resize(numCorners, 0);
      vtkCell.edges.resize(numEdges, 0);
      
      for(std::size_t k = 0; k < numCorners; ++k)
      {
        std::size_t cornerIndex = vtk.connectivity[offsCorners + k];
        vtkCell.corners[k] = cornerIndex;
        // mark as corner
        isCorner[cornerIndex] = 1;
      }
      
      for(std::size_t k = 0; k < numEdges; ++k)
      {
        vtkCell.edges[k] = vtk.connectivity[offsEdges + k];
      }
      cellStart = vtk.offsets[i];
    }
  }
  
  void VTKReader::reset()
  {
    gridWasLoaded = gridWasCreated = false;
    polynomialOrder = 0;
    vtk = VTKReaderDetail::VTKData();
    pointDataIndices.clear();
    cellDataIndices.clear();
  }
  
  void VTKReader::clearCoefficientData()
  {
    vtk.coefficientData.clear();
  }
} // namespace Kaskade

#ifdef UNITTEST

#include <numeric> // std::iota
#include <iostream>

#include <dune/grid/uggrid.hh>
#include <dune/alugrid/grid.hh>

#include <fem/assemble.hh>
#include <fem/gridmanager.hh>
#include <fem/lagrangespace.hh>
#include <fem/variables.hh>
#include <io/vtk.hh>
#include <utilities/gridGeneration.hh> //  createUnitSquare



// switch between square and cube grid depending on dimension
template<class Grid, int dim>
struct CreateTestGrid
{};
template<class Grid>
struct CreateTestGrid<Grid, 2>
{
  std::unique_ptr<Grid> operator()() const {return Kaskade::createUnitSquare<Grid>();}
};
template<class Grid>
struct CreateTestGrid<Grid, 3>
{
  std::unique_ptr<Grid> operator()() const {return Kaskade::createUnitCube<Grid>();}
};

// convenience typedef
template<int dim>
using AluGrid = Dune::ALUGrid<dim, dim, Dune::ALUGridElementType::simplex, Dune::ALUGridRefinementType::conforming>;

template<class Space>
void testVTKReader(int order)
{
  // Writes Variables to a file, then loads the file and compares loaded version and original.
  
  using namespace Kaskade;
  using namespace boost;
  
  using Spaces = fusion::vector<Space const*>;
  // VariableDescriptions: Note that writeVTK() will always output 1 or 3 data-components.
  using VariableDescriptions = fusion::vector<Variable<SpaceIndex<0>, Components<1>>,
                                              Variable<SpaceIndex<0>, Components<3>>>;
  using VariableSetDesc = VariableSetDescription<Spaces, VariableDescriptions>;
  using Grid = typename Space::Grid;
  constexpr int dim = Grid::dimension;
  constexpr int codimVertex = dim, codimEdge = dim - 1;
  
  // Create grid, spaces and a VariableSet varSet1 containing the FSEs u and v for output.
  CreateTestGrid<Grid, dim> create;
  GridManager<Grid> gridMan1(create()); // create unit-square or unit-cube
  Space space1(gridMan1, gridMan1.grid().leafGridView(), order);
  Spaces spaces1(&space1);
  std::vector<std::string> varNames = {"u", "v"};
  
  VariableSetDesc varSetDesc1(spaces1, varNames);
  typename VariableSetDesc::VariableSet varSet1(varSetDesc1);
  
  // Fill the coefficient vectors of varSet1 with unique values.
  std::vector<double> values;
  fusion::for_each(varSet1.data, [&](auto& fse)
  {
    values.resize(fse.size() * fse.dim());
    std::iota(values.begin(), values.end(), 0);
    vectorFromSequence(fse, values.begin());
  });
  
  // Save varSet1.
  std::ostringstream out;
  out << "unittest_dim" << dim << "_order" << order << "_" << (Space::continuity >= 0? "continuous" : "discontinuous");
  std::string filename = out.str();
  // Need to set order and data-mode for output.
  writeVTK(varSet1, std::string("graph/") + filename, IoOptions().setOrder(order).setDataMode(IoOptions::inferred).setOutputType(IoOptions::binary));
  
  // Load the file.
  VTKReader vtk(std::string("graph/") + filename + ".vtu");
  
  // Create a new grid from the vtk-file.
  GridManager<Grid> gridMan2(vtk.createGrid<Grid>());
  
  // Use the new grid to create a new space and a new VariableSet.
  Space space2(gridMan2, gridMan2.grid().leafGridView(), order);
  Spaces spaces2(&space2);
  VariableSetDesc varSetDesc2(spaces2, varNames); // (reusing varNames)
  typename VariableSetDesc::VariableSet varSet2(varSetDesc2);
  
  // Fill varSet2 with the loaded data. 
  vtk.getCoefficients(varSet2);
  
  // Compare old and new VariableSets.
  using Vertex = Dune::FieldVector<typename Grid::ctype, dim>;
  std::vector<Vertex> positions;
  // Going through the cells of one grid. Functions will be evaluated at global positions (independent of grid).
  for(auto const& cell : elements(gridMan2.grid().leafGridView()))
  {
    Vertex centerPos = cell.geometry().center();
    int numCorners = cell.subEntities(codimVertex), numEdges = cell.subEntities(codimEdge);
    // For each corner of a cell
    for(int i = 0; i < numCorners; ++i)
    {
      // add a point halfway between the cell's center and this corner.
      auto vertexEntity = cell.template subEntity<codimVertex>(i);
      Vertex cornerPos = vertexEntity.geometry().center();
      Vertex position = 0.5 * (centerPos + cornerPos);
      positions.push_back(position);
    }
    // The same for edge-midpoints.
    for(int i = 0; i < numEdges; ++i)
    {
      auto edgeEntity = cell.template subEntity<codimEdge>(i);
      Vertex edgePos = edgeEntity.geometry().center();
      Vertex position = 0.5 * (centerPos + edgePos);
      positions.push_back(position);
    }
  }
  
  // Store the norm of the difference between old and new functions at these positions.
  std::vector<double> errors;
  for(Vertex const& position : positions)
  {
    fusion::for_each(fusion::zip(varSet1.data, varSet2.data), [&](auto const& fsePair)
    {
      auto const &fse1 = fusion::at_c<0>(fsePair), &fse2 = fusion::at_c<1>(fsePair);
      auto error = fse1.value(position) - fse2.value(position);
      errors.push_back(error.two_norm());
    });
  }
  
  // Print the maximal error (should be 0).
  double maxError = *std::max_element(errors.begin(), errors.end());
  std::cout << filename << ": max error = " << maxError << std::endl;
}

int main()
{
  using namespace Kaskade;
  
  std::cout << "Testing UGGrid:" << std::endl;
  
  {
    using Grid2D = Dune::UGGrid<2>;
    using Grid3D = Dune::UGGrid<3>;
    using H1Space2D = FEFunctionSpace<ContinuousLagrangeMapper<double, Grid2D::LeafGridView>>;
    using H1Space3D = FEFunctionSpace<ContinuousLagrangeMapper<double, Grid3D::LeafGridView>>;
    using L2Space2D = FEFunctionSpace<DiscontinuousLagrangeMapper<double, Grid2D::LeafGridView>>;
    using L2Space3D = FEFunctionSpace<DiscontinuousLagrangeMapper<double, Grid3D::LeafGridView>>;
    
    // Test all combinations of 2D/3D, order 1/2, cont./disc.
    testVTKReader<H1Space2D>(1);
    testVTKReader<H1Space2D>(2);
    testVTKReader<H1Space3D>(1);
    testVTKReader<H1Space3D>(2);
    testVTKReader<L2Space2D>(1);
    testVTKReader<L2Space2D>(2);
    testVTKReader<L2Space3D>(1);
    testVTKReader<L2Space3D>(2);
    
    // Test discontinuous, order 0.
    testVTKReader<L2Space2D>(0);
    testVTKReader<L2Space3D>(0);
  }
  
  std::cout << "\nTesting ALUGrid:" << std::endl;
  
  {
    using Grid2D = AluGrid<2>;
    using Grid3D = AluGrid<3>;
    using H1Space2D = FEFunctionSpace<ContinuousLagrangeMapper<double, Grid2D::LeafGridView>>;
    using H1Space3D = FEFunctionSpace<ContinuousLagrangeMapper<double, Grid3D::LeafGridView>>;
    using L2Space2D = FEFunctionSpace<DiscontinuousLagrangeMapper<double, Grid2D::LeafGridView>>;
    using L2Space3D = FEFunctionSpace<DiscontinuousLagrangeMapper<double, Grid3D::LeafGridView>>;
    
    // Test all combinations of 2D/3D, order 1/2, cont./disc.
    testVTKReader<H1Space2D>(1);
    testVTKReader<H1Space2D>(2);
    testVTKReader<H1Space3D>(1);
    testVTKReader<H1Space3D>(2);
    testVTKReader<L2Space2D>(1);
    testVTKReader<L2Space2D>(2);
    testVTKReader<L2Space3D>(1);
    testVTKReader<L2Space3D>(2);
    
    // Test discontinuous, order 0.
    testVTKReader<L2Space2D>(0);
    testVTKReader<L2Space3D>(0);
  }
  
  return 0;
}
  
#endif





