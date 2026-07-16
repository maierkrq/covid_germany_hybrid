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

#ifndef AMIRA_HH
#define AMIRA_HH

#include <string>
#include <iostream>
#include <fstream>

#define HAVE_AMIRAMESH 1
#undef max
#include <dune/grid/io/file/amirameshwriter.hh>
#include "io/iobase.hh"

#include <boost/fusion/include/for_each.hpp>
#include <boost/fusion/include/zip_view.hpp>
#include <boost/fusion/include/vector.hpp>

#include "utilities/enums.hh"
#include "fem/assemble.hh"
#include "fem/spaces.hh"

#include "io/vtk.hh"
#include "io/matlab.hh"
#include "io/amirameshreader.hh"
#include "io/amira.hh"
#include "io/amiramesh.hh"


/**
 * @file
 * @brief  Output of mesh and solution for visualization software
 * @author Martin Weiser, Bodo Erdmann
 *
 * This file contains functionality needed to write mesh data and
 * solution data for Amira.
 */

//---------------------------------------------------------------------

namespace Kaskade
{
  template<class GridType>
  class MaterialAmiraMeshWriter
    : public Dune::AmiraMeshWriter<typename GridType::LeafGridView>
  {
    
  public:
    typedef typename GridType::LeafGridView GridView;
    typedef typename Dune::BlockVector<Dune::FieldVector<double, 1> > DataContainer;
  
    /** \brief Default constructor */
    MaterialAmiraMeshWriter() {}

    /** \brief Constructor which initializes the AmiraMesh object with a given leaf grid */
    MaterialAmiraMeshWriter(const GridType& grid) {
      this->addGrid(grid.leafGridView());
    }
    
    /** 
     *  This method writes the data of the material of each cell to the undelying AmiraMesh. 
     * 
     *  @param gridView     gridView, i.e. LeafGridView
     *  @param materialData Dune::BlockVector containing data on the material of the elements in the grid.
     * 
     *  */
    void addMaterialData(const GridView& gridView, DataContainer& materialData) {
      
      const typename GridView::IndexSet& indexSet = gridView.indexSet();
      int noOfElements  = indexSet.size(0);
      
      AmiraMesh::Data* data = Dune::AmiraMeshWriter<GridView>::amiramesh_.findData("Tetrahedra",HxDOUBLE,1,"Materials");
      
      for(int i=0; i<noOfElements; i++)
        ((unsigned char*)data->dataPtr())[i] = materialData[i];
    }
    
    
  };
  
  
  /**
   * \cond internals
   */
  namespace IoDetail {

    template <class GridView, class Filter, class RAIter>
    struct AddDataTo<GridView,Kaskade::MaterialAmiraMeshWriter<typename GridView::Grid>,Filter,RAIter>
    {
      typedef typename GridView::Grid Grid;
      typedef typename Kaskade::MaterialAmiraMeshWriter<typename GridView::Grid> AmiraWriter;

      AddDataTo(GridView const& gridView_, AmiraWriter& writer_, Filter const& filter_, RAIter const& names_):
        gridView(gridView_),
        writer(writer_),
        filter(filter_),
        names(names_)
      {}


      template <class Pair>
      void operator()(Pair const& pair) const
      {
        using namespace boost::fusion;

        int const id = boost::remove_reference<typename result_of::value_at_c<Pair,0>::type>::type::id;

        if (!filter(names[id]))
          return;

        typedef typename boost::remove_reference<typename result_of::value_at_c<Pair,1>::type>::type Function;

        Function const& f = at_c<1>(pair);
        //    std::string const& name = names[id];

        // Get the index set. We hope this is actually the leaf index set
        // (probably yes, but not guaranteed). Unfortunately there seems
        // to be no possibility to check index sets for equality.
        typedef typename Function::Space::IndexSet IndexSet;
        IndexSet const& indexSet = f.space().indexSet();


        // Create evaluator. Todo: Create and use a vertex-based
        // quadrature rule, create and use a shape function cache based on
        // this quadrature rule.
        typename Function::Space::Evaluator eval(f.space());


        if (Function::Space::continuity >= 0) 
        {
          // Continuous function. Written as values at grid vertices.
          Dune::BlockVector<typename Function::ValueType> data(indexSet.size(GridView::dimension));

          // Step through all cells.
          for (auto const& cell: Dune::elements(f.space().gridView()))
          {
            // Obtain reference element.
            auto const& refElem = Dune::ReferenceElements<typename Grid::ctype, GridView::dimension>::general(cell.type());
            // Move evaluator to current cell.
            eval.moveTo(cell);

            // Step through all vertices of the current cell.
            for (int i=0; i<cell.subEntities(GridView::dimension); ++i) 
            {
              // Move evaluator to current vertex.
              eval.evaluateAt(refElem.position(i,GridView::dimension));
              // Store value of function in data vector (inefficient, occurs multiple times...).
              data[indexSet.subIndex(cell,i,GridView::dimension)] = f.value(eval);
            }
          }

          writer.addVertexData(data,gridView);
        } 
        else 
        {
          // Discontinuous function. Written as values at cells.
          Dune::BlockVector<typename Function::ValueType> data(indexSet.size(0));

          // Step through all cells.
          for (auto const& cell: Dune::elements(f.space().grid().leafGridView()))
          {
            // Obtain reference element.
            auto const& refElem = Dune::ReferenceElements<typename Grid::ctype, GridView::dimension>::general(cell.type());

            data[indexSet.index(cell)] = f.value(cell,refElem.position(0,0));
          }

          writer.addCellData(data,gridView);
        }
      }

    private:
      GridView const& gridView;
      AmiraWriter&  writer;
      Filter const& filter;
      RAIter        names;

    };

  } // End of namespace IoDetail

  /**
   * \endcond
   */

  //---------------------------------------------------------------------
   
   
     /**
   *  This method copies/adapts the header of an Amira-file
   *  containing additional information about the materials. (e.g which material 
   *  id belongs to which material)
   *  @param gridView gridView, i.e. LeafGridView
   *  @param file_in file, from which the additional material data is read (header in Amira-file)
   *  @param file_out file to be written in. Should be the same file as in the call of writeAMIRAFile
   * 
   */
   template <class GridView>
  void copyHeader(GridView const& gridView, std::string const& file_in, std::string const& file_out){
    
    int nTet   = gridView.indexSet().size(0);
    int nNodes = gridView.indexSet().size(int(GridView::dimensionworld));
    
    std::string find1("nNodes ");
    std::string find2("nTetrahedra ");
    std::ofstream tmpFile("tmp.txt");
    std::ifstream inFile(file_in);
    std::ifstream outFile(file_out);
    std::string   line;
    bool start = false;
    
    if (!file_in.empty())
    {                    //get amira-file header from input file_in
      while( getline(inFile, line) )
      {
        if( line.substr(0,7).compare("Nodes {") == 0 )
          break;
        else if( line.substr(0,7).compare(find1) == 0)
          tmpFile << find1 << nNodes << std::endl;
        else if( line.substr(0,12).compare(find2) == 0)
          tmpFile << find2 << nTet << std::endl;
        else
          tmpFile << line << std::endl;
      }
    }
    else
      start = true;

    inFile.close();

    //copying data section of file_out
    while( getline(outFile,line) )
    {
      if (line.substr(0,7).compare("Nodes {") == 0)
        start = true;        //Start of the data section

      if (start)
        tmpFile << line << std::endl;
    }
    outFile.close();
    tmpFile.close();
    
    //rename
    if (std::rename("tmp.txt" , file_out.c_str()) == 0)
      ;
    else
      perror( "Error renaming file in call of copyHeader() in writeAMIRAFile" );
    
    std::remove("tmp.txt");
  }


  /**
   *   This procedure writes data in AmiraMesh with cell material data to ascii/binary format, e.g. for visualization in AMIRA software.
   *   @param gridView    gridView, i.e. LeafGridView
   *   @param description variable set description
   *   @param vars        data of variables
   *   @param file_out    name of output file
   *   @param options     defines the format for the output,e.g., ascii/binary
   *   @param material    material data of underlying grid. Typically called as material.coefficients(),
   *                      where material is a is a FunctionSpaceElement.
   *   @param file_in     the Amira-input file the material data was read from
   */
  template <class GridView, class VariableSet>
  void writeAMIRAFile(GridView const& gridView,
                      VariableSet const& vars,
                      std::string file_out,
                      IoOptions options=ioOptions_default,
                      Dune::BlockVector<Dune::FieldVector<double, 1> > materialData = Dune::BlockVector<Dune::FieldVector<double, 1>>(),
                      std::string file_in = "")
  {
    //check for .am extension
    if( (!file_in.empty()) && (file_in.substr( file_in.length()-3) != ".am") ) file_in.append(".am");
    if( (file_out.substr( file_out.length()-3) != ".am") ) file_out.append(".am");
    
    typedef Kaskade::MaterialAmiraMeshWriter<typename GridView::Grid> AmiraWriter;
    //  AmiraWriter writer(gridView.grid());
    AmiraWriter writer(vars.descriptions.gridView.grid());
    
    if( materialData.size() != 0) writer.addMaterialData(vars.descriptions.gridView, materialData);

    std::string fullname = file_out;

    if (options.outputType == IoOptions::ascii)
      writePartialFile(vars.descriptions.gridView,writer,vars,fullname,UnaryTrue(),1);
    else if (options.outputType == IoOptions::binary)
      writePartialFile(vars.descriptions.gridView,writer,vars,fullname,UnaryTrue(),0);
    else std::cout << "Amira:  graphic output for nonvalid driver" << std::endl;
    
    if(!file_in.empty()) copyHeader(vars.descriptions.gridView, file_in, file_out);
    
  }
   
   
} // end namespace Kaskade

#endif
