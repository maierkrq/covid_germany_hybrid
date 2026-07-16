/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2019-2022 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef GRIDBASICS_HH_
#define GRIDBASICS_HH_

#include <mutex>

// forward declarations of grids
#include "dune/common/fvector.hh"
#include "dune/grid/geometrygrid/declaration.hh"

#include "utilities/threading.hh"



// ------------------------------------------------------------------------------------------------

namespace Kaskade
{
  /**
   * \ingroup grid
   * \brief The type of cells (entities of codimension 0) in the grid view.
   * \tparam GridView the type of grid view used (a grid type works as well)
   */
  template <class GridView>
  using Cell = typename GridView::template Codim<0>::Entity;

 /**
   * \ingroup grid
   * \brief The type of faces in the grid view.
   */
  template <class GridView>
  using Face = typename GridView::Intersection;

  /**
   * \ingroup grid
   * \brief The type of global positions within the grid view.
   * \tparam GridView a grid view type (a grid type works as well)
   */
  template <class GridView>
  using GlobalPosition = Dune::FieldVector<typename GridView::ctype,GridView::dimensionworld>;


  /**
   * \ingroup grid
   * \brief The type of local positions within the grid view.
   * \tparam GridView a grid view type (a grid type works as well)
   */
  template <class GridView>
  using LocalPosition = Dune::FieldVector<typename GridView::ctype,GridView::dimension>;
  

  // ----------------------------------------------------------------------------------------------

  /**
   * \brief Grid locking information based on grid type.
   *
   * Some Dune grids appear not to be thread-safe even for read operations. Based on the
   * grid type, we provide this information and provide a mutex which can be used to lock
   * individual calls to the grid interface.
   */
  template <class Grid>
  struct GridLocking: public std::false_type
  {
    /**
     * \brief Returns a globally unique lockable object that shall be used for locking access to
     *        operations on the grid.
     *
     * For grids which are (supposed to be) thread-safe, this returns a dummy lockable, which does
     * exactly nothing (as it should be). For non-safe grids this returns a recursive mutex. Thus,
     * the grid access can be made mutually exclusive without caring whether a calling routine has
     * already acquired the mutex.
     */
    static GridLocking<Grid>& mutex()
    {
      static GridLocking<Grid> m;
      return m;
    }

    void lock()
    {}

    bool try_lock()
    {
      return true;
    }

    void unlock()
    {}
  };



  template <class HostGrid, class Deformation, class Allocator>
  struct GridLocking<Dune::GeometryGrid<HostGrid,Deformation,Allocator>>: public std::true_type
  {
    static std::recursive_mutex& mutex()
    {
      static std::recursive_mutex m;
      return m;
    }
  };

  // ----------------------------------------------------------------------------------------------


}

#endif
