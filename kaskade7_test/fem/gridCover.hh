/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2018-2018 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef GRIDCTOOLS_HH
#define GRIDCTOOLS_HH

#include <algorithm>
#include <utility>

#include "fem/fixdune.hh"
#include "utilities/gridGeneration.hh"

namespace Kaskade
{
  /**
   * \ingroup grid 
   * \brief Computes a the bounding box of the given grid view.
   * 
   * \return a pair (low,high) of vectors denoting the smallest and largest corner of
   *         the bounding box
   * Keep in mind that level grid views need not cover the whole domain.
   */
  template <class GridView>
  std::pair<Dune::FieldVector<typename GridView::ctype,GridView::dimensionworld>,
            Dune::FieldVector<typename GridView::ctype,GridView::dimensionworld>>
  boundingBox(GridView const& gv)
  {
    using ctype = typename GridView::ctype;
    int const dim = GridView::dimensionworld;
    
    Dune::FieldVector<ctype,GridView::dimensionworld> low(std::numeric_limits<ctype>::max()), 
                                                      high(std::numeric_limits<ctype>::lowest());
    
    for (auto const& cell: Dune::elements(gv))
    {
      auto const& geo = cell.geometry();
      for (int i=0; i<geo.corners(); ++i)
      {
        auto x = geo.corner(i);
        low = Dune::min(low,x);
        high = Dune::max(high,x);
      }
    }
    
    return std::make_pair(low,high);
  }
  
  // ----------------------------------------------------------------------------------------------
  
  /**
   * \ingroup grid
   * \brief Creates a cartesian structured grid occupying the bounding box of the given grid view.
   * 
   * Starting from a uniform grid, the cover grid is locally refined until its cells are not larger
   * than the provided grid view cells by a certain factor.
   * 
   * \param gv the provided grid view
   * \param volumeRatio the maximum local ratio of cell volumes of the cover grid and the provided grid
   */
  template <class GridView>
  std::unique_ptr<typename GridView::Grid> createCoverGrid(GridView const& gv, double volumeRatio=10)
  {
    using Grid = typename GridView::Grid;
    
    // Create a uniform cartesian structured coarse grid that covers the bounding box of the given grid.
    // The edge length of the cubes is chosen maximal, i.e. the smallest width of the bounding box.
    auto box = boundingBox(gv);
//     std::cout << "Grid bounding box: " << box.first << " --- " << box.second << "\n";
    auto d = box.second - box.first;
    auto h = d;
    while (true)
    {
      auto& maxh = *std::max_element(h.begin(),h.end());
      auto& minh = *std::min_element(h.begin(),h.end());
      if (maxh > 2*minh)
        maxh /= 2;
      else
        break;
    }
//     std::cout << "aiming at bounding box: " << box.first << " --- " << box.first+d << "\n";
    
    std::unique_ptr<Grid> grid = createCuboid<Grid>(box.first,d,h);
    auto coverBox = boundingBox(grid->leafGridView());
//     std::cout << "Grid bounding box: " << coverBox.first << " --- " << coverBox.second << "\n";
    
    // Refine the cover grid until its cells have a volume larger than the given volumeRatio of
    // any source grid cells they cover.
    while (true)
    {
      bool refine = false;
      for (auto const& cell: Dune::elements(gv))
      {
        auto x = cell.geometry().center();
        auto const covercell = findCell(grid->leafGridView(),x);
        if (covercell.geometry().volume() > volumeRatio*cell.geometry().volume())
        {
          grid->mark(1,covercell);
          refine = true;
        }
      }
      
      // stop refining if there is no cell to be refined
      if (!refine)
        break;
      
      grid->preAdapt();
      grid->adapt();
      grid->postAdapt();
    }
    
    return grid;
  }

} /* end of namespace Kaskade */
#endif
