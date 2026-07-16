/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2022 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef VTK_HH
#define VTK_HH

#include "io/check_endianness.hh"
#include "io/iobase.hh"
#include "io/vtkwriter.hh"
#include "fem/variables.hh"
#include "utilities/deprecation.hh"

/**
 * @file
 * @brief  Output of mesh and solution for visualization software
 * @author Martin Weiser, Bodo Erdmann
 *
 * This file contains functionality needed to export mesh data and
 * solution data.
 *
 * - \a writeVTKFile: write data to VTK XML file, e.g. for visualization in Paraview software
 * - \a TimeSeriesVTKWriter: A class for writing a (time-) sequence of variable sets into VTK
 * XML files and creating a Paraview description file for playing the
 * animation. The Paraview description file is completed and closed on
 * destruction of the writer object.

 *
 */

//---------------------------------------------------------------------

namespace Kaskade
{

  /**
   * \ingroup IO
   * \brief Writes a single finite element function to a VTK file.
   * 
   * \tparam Function the finite element function type
   *
   * \param[in] f the function to plot
   * \param[in] filename the name of the output file (.vtu is appended automatically)
   * \param[in] options
   * \param[in] fname name of the function in the VTK output file. If empty, the filename is used.
   */
  template <class Function>
  void writeVTK(Function const& f, std::string const& filename, IoOptions options, std::string fname)
  {
    using namespace boost::fusion;
    
    if (fname.empty())
      fname = filename;
    
    // Construct a variable set (trivially consisting of only one function) to be handed to the
    // VTK writer.
    auto vsDescription = makeVariableSetDescription(
                                makeSpaceList(&f.space()),
                                make_vector(Variable<SpaceIndex<0>,
                                                     Components<Function::components>>(fname)));
    auto fs = vsDescription.variableSet();
    component<0>(fs) = f;
    
    // write the VTK file
    writeVTK(fs,filename,options);
  }


  //---------------------------------------------------------------------

  /**
   * \ingroup IO
   *
   * A class for writing a (time-) sequence of variable sets into VTK
   * XML files and creating a Paraview description file for playing the
   * animation. The Paraview description file is completed and closed on
   * destruction of the writer object.
   */
  class TimeSeriesVTKWriter
  {
  public:
    TimeSeriesVTKWriter(std::string const& basename, IoOptions options_= ioOptions_default)
    : base(basename), pvd((basename+".pvd").c_str()), count(0), options(options_)
    {
      pvd << "<?xml version=\"1.0\"?>\n"
          << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"" <<  VTKWriterDetail::vtkEndianness() << "\">\n"
          << "  <Collection>\n";
    }

    /**
     * \brief Adds frame with values given by \a vars at time \a time to the series.
     */
    template <class VariableSet>
    void add(double time, VariableSet const& vars)
    {
      auto fname = base+"-"+paddedString(count,4);
      writeVTK(vars,fname,options);

      pvd << "    <DataSet timestep=\"" << time << "\" file=\"" << fname << ".vtu\" />\n";
      ++count;
    }

    ~TimeSeriesVTKWriter()
    {
      pvd << "  </Collection>\n"
          << "</VTKFile>\n";
    }

  private:
    std::string base;
    std::ofstream pvd;
    int count;
    IoOptions options;
  };
  
  
  //---------------------------------------------------------------------
  //---------------------------------------------------------------------
  


  
} // end namespace Kaskade
//---------------------------------------------------------------------

#endif
