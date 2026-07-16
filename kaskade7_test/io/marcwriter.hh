/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2020-2021 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef MARCWRITER_HH
#define MARCWRITER_HH

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <fstream>
#include <map>

#include <utilities/detailed_exception.hh>

namespace Kaskade
{

  /// \cond internals
  namespace MarcWriterDetail
  {

    int marcElementIdFromGeometryType(Dune::GeometryType gt)
    {
      // These are listed in Marc Volume B as "heat conduction elements"
      if (gt.isTriangle())     return 37;
      if (gt.isTetrahedron())  return 135;
      throw LookupException("Unhandled cell type encountered in Marc output.",__FILE__,__LINE__);
      return 0;
    }

    // MARC wants floating point numbers as -1.234+5, i.e. without the usual 'e' denoting the
    // exponent.
    std::string marcFloatNumber(double x)
    {
      int const n = 50;
      char buf[n];
      std::snprintf(buf,n,"%E",x);
//       std::remove(buf,buf+n,'e');
      return std::string(buf);
    }
  }
  /// \endcond

  /**
   * \ingroup io
   * \brief Writes grid and material IDs given as FE function to a Marc input file.
   *
   * This creates an input file for the Marc finite element solver by MSC Software. The created
   * file is intended for heat simulations (not mechanics), and provides dummy thermal material
   * parameters which need to be fixed before submitting the file to MARC.
   *
   * The underlying grid must contain just one type of cell.
   *
   * \param material a scalar integer-valued piecewise constant FE function defining the
   *                 material id
   * \param filename file name (without trailing ".marc" suffix)
   * \param scale a scaling factor for the geometry. All lengths and positions are scaled by this factor.
   */
  template <class Material>
  void writeMarc(Material const& material, std::string const& filename, double const scale=1.0, IoOptions options=IoOptions())
  {
    using namespace MarcWriterDetail;

    // open file
    std::ofstream out(filename+".marc");
    if (!out)
      throw FileIOException("File " + filename + ".marc could not be opened for writing Marc output.",
                            filename+".marc",__FILE__,__LINE__);


    using GridView = typename Material::GridView;
    auto const& gridView = material.space().gridView();
    auto const& is = gridView.indexSet();
    constexpr int dim = GridView::dimension;
    constexpr int dimworld = GridView::dimensionworld;

    int const nCells = gridView.size(0);
    int const nVertices = gridView.size(dim);

    // Make sure the material function is piecewise constant and scalar
    assert(material.space().mapper().maxOrder()==0);
    static_assert(Material::components==1);
    assert(material.dim()==nCells);


    // Get the geometry (of first cell, which must coincide with all other cells)
    Dune::GeometryType gt = gridView.template begin<0>()->geometry().type();
    if (gridView.size(gt) != nCells)
      throw GridException("Marc output implemented for grids with a single cell geometry type only.",
                          __FILE__,__LINE__);

    int elementId = marcElementIdFromGeometryType(gt);

    // Marc requires floating point data to be formatted as float.
    out << std::scientific << std::setprecision(10);

    // Marc input is organized in three groups: parameter data, model definition data,
    // and history definition data. First write parameter data.

    out << "TITLE " << "Kaskade 7 output file " << filename << "\n";
    out << "EXTENDED\n";                          // necessary for long float mantissa
    out << "ELEMENTS," << elementId << ",\n";     // number of grid cells
    out << "SIZING\n";
    out << "VERSION, 15, 1, 0, 1\n";              // necessary for correct interpretation by MARC
    out << "HEAT, 0, 0, 1, 0, 1, 1\n";            // necessary for interpretation as heat problem (not mechanics)
    out << "END\n";


    // Next gather cells by their material id.
    std::map<int,std::vector<std::vector<int>>> id;
    for (auto const& cell: Dune::elements(gridView))
    {
      int cellIdx = is.index(cell);
      std::vector<int> c;
      c.push_back(cellIdx+1);
      for (int i=0; i<cell.subEntities(dim); ++i)
        c.push_back(is.index(cell.template subEntity<dim>(i))+1);

      int idi = static_cast<int>(std::round(material.coefficients()[cellIdx]));
      id[idi].push_back(c);
    }


    // Next the model definition data. First the grid topology, with one "CONNECTIVITY" section
    // per material id. Marc appears not to accept arbitrary material ids, only contiguous numbers
    // starting at 1. We therefore skip the material id and count the number up.
    int materialNumber = 1;
    for (auto const& [_,cs]: id)
    {
      out << "CONNECTIVITY\n";
      out << "0, 0, 1, 0, 0, " << materialNumber << ", \n";  // this encodes the material id of subsequent cells
      for (auto const& c: cs)
      {
        out << c[0] << ", " << elementId;
        for (int i=1; i<c.size(); ++i)
          out << ", " << c[i];
        out << "\n";
      }
      ++materialNumber;
    }

    // Write out vertex coordinates
    out << "COORDINATES\n";
    out << dimworld << ", " << nVertices << "\n";
    for (auto const& v: Dune::entities(gridView,Dune::Codim<dim>()))
    {
      out << is.index(v)+1;
      auto x = v.geometry().center();
      for (int i=0; i<dimworld; ++i)
        out << ", " << marcFloatNumber(scale*x[i]);
      out << "\n";
    }

    // Then write out dummy (thermal) material properties for each material.
    // The odd blank lines in the output appear to be necessary for MARC not to choke.
    for (int i=1; i<=size(id); ++i)
      out << "ISOTROPIC\n\n"
          << i << ", 0, 0, 0, 0, mati" << i << "\n"
          << "8.00000-2, 1.00000+3, 4.70000+2, 0.00000+0, 9.00000-1, 0.00000+0, 0.00000+0\n\n";


    out << "END OPTION\n";

  }

}

#endif
