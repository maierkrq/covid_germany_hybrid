/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2016-2018 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef FEMSUPPORTS_HH
#define FEMSUPPORTS_HH

#include <algorithm>
#include <vector>

#include "scalarspace.hh"

namespace Kaskade
{

  // ----------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------


  /**
   * \ingroup fem
   * \brief Computes the global indices of all ansatz functions whose support is a subset of the union of given cells.
   *
   * These ansatz functions span the space of FE functions on the patch with zero trace.
   *
   * \tparam Space a FEFunctionSpace instatiation
   * \tparam Range an STL range type with grid cells as value type
   */
  template <class Space, class Range>
  std::vector<size_t> algebraicPatch(Space const& space, Range const& cells)
  {
    using std::begin; using std::end;
    std::vector<size_t> dofs;

    // find all ansatz functions with support intersecting the patch (these are too many, of course)
    for (auto const& cell: cells)
    {
      auto const& range = space.mapper().globalIndices(cell);
      dofs.insert(end(dofs),begin(range),end(range));
    }
    assert(dofs.size()>0);

    // remove duplicates
    std::sort(begin(dofs),end(dofs));
    dofs.erase(std::unique(begin(dofs),end(dofs)),end(dofs));

    // We still have too many ansatz functions, as there are some whose support extends to cells outside
    // the given patch. Thus, for each cell in the patch, consider every neighbor that is not contained
    // in the patch, and remove each ansatz function that lives on that patch.
    // TODO: consider using std::set_difference, but this needs two sorted ranges (could be done efficiently by using sortedIndices),
    //       and the output range may not overlap the input range as of C++14 (why? appears to be overspecified)
    auto const& gridview = space.gridView();
    for (auto const& cell: cells)                                               // step through all cells
      for (auto fi=gridview.ibegin(cell); fi!=gridview.iend(cell); ++fi)        // step through all its faces
        if (fi->neighbor())                                                     // if there's a neighbor
        {
          auto neighbor = fi->outside();
          if (std::find(begin(cells),end(cells),neighbor)==end(cells))          // and this neigbor is not in the patch
            for (size_t i: space.mapper().globalIndices(neighbor))              // consider all ansatz functions with support on the neighbor
            {
              auto j = std::lower_bound(begin(dofs),end(dofs),i);               // and if they had been entered into our dofs
              if (j!=end(dofs) && *j==i)
              {
                dofs.erase(j);                                                  // remove them - they have support outside the patch
                assert(dofs.size()>0);                                          // (erasure retains the sorting of indices)
              }
            }
        }

    return dofs;
  }





  /**
   * \ingroup fem
   * \brief Computes the global indices of all ansatz functions whose support is a subset of the union of given cells.
   *
   * These ansatz functions span the space of FE functions on the patch with zero trace.
   * This version also treats patches that touch the domain boundary correctly (the above version takes too many indices on these patches).
   *
   * \tparam Space a FEFunctionSpace instatiation
   * \tparam Range an STL range type with grid cells as value type
   * \param vertexIdx the global index of patch middle vertex
   * \param indicesOnFacets contains for each facet (by global facet index) the indices of all ansatz function that do not vanish on that facet
   */
  template <class Space, class Range>
  std::vector<size_t> algebraicPatch(Space const& space, Range const& cells,
                                     size_t vertexIdx, std::vector<std::vector<size_t>> const& indicesOnFacets)
  {
    using std::begin; using std::end;
    std::vector<size_t> dofs;
    int const dim = Space::Grid::dimension;

    // find all ansatz functions with support intersecting the patch (these are too many, of course)
    for (auto const& cell: cells)
    {
      auto const& range = space.mapper().globalIndices(cell);
      dofs.insert(end(dofs),begin(range),end(range));
    }
    assert(dofs.size()>0);

    // remove duplicates
    std::sort(begin(dofs),end(dofs));
    dofs.erase(std::unique(begin(dofs),end(dofs)),end(dofs));

    // We still have too many ansatz functions, as there are some whose support extends to cells outside
    // the given patch. Thus, for each cell in the patch, consider the facet that lies on the patch boundary
    // and remove each ansatz function that does not vanish on that face.
    for (auto const& cell: cells)                                               // step through all cells
    {
      int localVertexIndex = 0; // determine which cell-local index lies in the middle of the patch
      for (int ci=0; ci<dim+1; ++ci)
        if(space.indexSet().subIndex(cell,ci,dim) == vertexIdx)
          localVertexIndex = ci;

      for (auto const& intersection : intersections(space.gridView(), cell))        // step through all the faces
      {
        int indexInInside = intersection.indexInInside();
        if (!ScalarSpaceDetail::onSimplexFace(dim,dim, localVertexIndex, indexInInside)) // if face lies on boundary of patch
        {
          for (size_t i: indicesOnFacets[space.indexSet().subIndex(cell,indexInInside,1)]) // consider all ansatz functions which are nonzero on this face
          {
            auto j = std::lower_bound(begin(dofs),end(dofs),i);               // and if they have not already been removed from dofs
            if (j!=end(dofs) && *j==i)
            {
              dofs.erase(j);                                                  // remove them - they have support outside the patch
              assert(dofs.size()>0);                                          // (erasure retains the sorting of indices)
            }
          }
        }
      }
    }

    return dofs;
  }


  /**
   * \ingroup fem
   * \brief Computes the patches around vertices of the grid view.
   *
   * \return a random access sequence (indexed by vertex index) containing a sequence of cells forming the patch
   *         around that cell
   */
  template <class GridView>
  auto computePatches(GridView const& gridview)
  {
    using Cell = typename GridView::template Codim<0>::Entity;
    int const dim = GridView::dimension;
    std::vector<std::vector<Cell>> patches(gridview.size(dim));                 // initialize with size of number of vertices

    for (auto const& cell: elements(gridview))                                  // step through the cells
      for (int i=0; i<cell.subEntities(dim); ++i)                               // and their vertices
        patches[gridview.indexSet().subIndex(cell,i,dim)].push_back(cell);      // add the cell to the patch of that vertex

    return patches;
  }

}

#endif
