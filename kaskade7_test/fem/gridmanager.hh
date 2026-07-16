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
#ifndef GRIDMANAGER_HH_
#define GRIDMANAGER_HH_

#include <ctime>
#include <iostream>
#include <mutex>
#include <numeric>
#include <utility>

#include <boost/signals2.hpp>
#include <boost/range/iterator_range.hpp>

#ifdef MULTITHREAD
#include <boost/thread.hpp>
#endif

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/gridview.hh>
#include <dune/grid/io/file/dgfparser/gridptr.hh>

#ifdef MULTITHREAD
#include "utilities/threadpool/threadpool-0_2_5-src/threadpool/boost/threadpool.hpp"
#endif

#include "fem/gridBasics.hh"
#include "utilities/threading.hh"
#include "utilities/timing.hh"



namespace Kaskade
{

  // forward declaration
  template<class G, class T> class CellData;

  // ----------------------------------------------------------------------------------------------

  /**
   * \ingroup grid
   * \brief A class that provides access to signals that are emitted from the grid manager on various occasions.
   */
  struct GridSignals
  {
    /**
     * \ingroup grid
     * Signal group names for use with the informAboutRefinement
     * signal. Freeing resources should occur before allocating new
     * resources.
     */
    enum { freeResources=0, allocResources=1 };

    /**
     * \ingroup grid
     * \brief The argument type of the signal that is emitted before and after grid adaptation.
     */
    enum Status { BeforeRefinement, AfterRefinement, TransferCompleted };


    /**
     * \ingroup grid
     * \brief A signal that is emitted thrice on grid modifications, once before
     * adaptation takes place and twice after it is completed.
     * 
     * The current state of grid refinement and FE coefficient transfer is signaled by the status argument.
     */
    boost::signals2::signal<void (Status)> informAboutRefinement;
  };
  
  // ----------------------------------------------------------------------------------------------

  /**
   * \ingroup grid
   * \brief A class storing cell ranges of roughly equal size for multithreaded mesh traversal.
   */
  template <class GridView>
  class CellRanges
  {
  public:
    /**
     * \brief The type of iterator to step through the cells of the given view.
     */
    typedef typename GridView::template Codim<0>::Iterator CellIterator;
    
    /**
     * \brief Default constructor.
     * This creates a pretty useless empty cell ranges object.
     */
    CellRanges() = default;
    
    /**
     * \brief Move constructor.
     */
    CellRanges(CellRanges&&) = default;
    
    /**
     * \brief Constructs cell ranges.
     * \param nmax Compute and store cell ranges for up to nmax groups
     * \param gridView the grid view to use (usually a leaf view). The iterators returned by the grid view
     *                 have to be valid for the lifetime of the CellRanges object - usually this means the 
     *                 grid shall not be modified/adapted.
     */
    CellRanges(int nmax, GridView const& gridView): cellRanges(nmax+1)
    {
      // We know how many iterators we will insert - allocate space beforehand. 
      for (int n=1; n<=nmax; ++n)
        cellRanges[n].reserve(n+1);
      
      // number of cells
      size_t const size = gridView.size(0);
      
      // Step through all cells, storing the iterator in case it is at a boundary of cell ranges.
      // The iteration through the cells is monotone, hence the ranges we construct cannot overlap. If few
      // cells exist, some ranges may be empty, however. 
      // The following postconditions hold: 
      // (i) cellRanges[n].size() <= n. 
      // Proof: Assume cellRanges[n].size()==n. Then there is no cell with j >= size (it is at most size-1), where a further push_back could happen.
      // (ii) cellRanges[n][0] == begin. 
      // Proof: In the first iterations j==0 and cellRanges[n].size()==0 holds, such that 0>=0 evaluates true and triggers the push_back(begin).
      // (iii) any range contains less than size/n + 1 cells.
      // Proof: Assume a range [r,s[. Then there is k with r>=k*size/n and t<(k+1)*size/n for all t<s. Thus s-1-r < ((k+1)-k)*size/n = size/n,
      // such that the length is s-r < size/n + 1.
      size_t j = 0; // cell counter
      for (CellIterator i=gridView.template begin<0>(); i!=gridView.template end<0>(); ++i, ++j)
        for (int n=1; n<=nmax; ++n)
          if (j >= cellRanges[n].size()*size/n)
            cellRanges[n].push_back(i);

      // Now we have to make sure that every range vector contains the number of ranges it's supposed to hold. Add empty ranges as needed.
      // As up to here the size of cellRanges[n] is at most n (not n+1), we automatically take care of terminating the last range with 
      // the end iterator.
      for (int n=1; n<=nmax; ++n)
        while (cellRanges[n].size()<n+1)
          cellRanges[n].push_back(gridView.template end<0>());
    }
    
    /**
     * \brief Move assignment.
     */
    CellRanges<GridView>& operator=(CellRanges<GridView>&&) = default;
    
    /**
     * \brief Obtain cell range.
     * \param n the number of cell ranges (1<=n<=maxRanges())
     * \param k the k-th range will be returned (0<=k<n)
     */
    boost::iterator_range<CellIterator> range(int n, int k) const
    {
      assert(1<=n && n<cellRanges.size());
      assert(0<=k && k<n);
      return boost::iterator_range<CellIterator>(cellRanges[n][k],cellRanges[n][k+1]);
    }
    
    /**
     * \brief Obtain the maximum number of ranges.
     * Default constructed "empty" cell ranges return a negative value.
     */
    int maxRanges() const
    {
      return static_cast<int>(cellRanges.size())-1;
    }
    
  private:
    std::vector<std::vector<CellIterator>> cellRanges;
  };
  
  //----------------------------------------------------------------------------------------

  // Start of implementation for better numbering of entities. Very much in progress, still.
  template <class GridView>
  class EntityNumbering
  {
    using IndexSet = typename GridView::IndexSet;
    constexpr static int dim = GridView::dimension;
    
  public:
    using Index = typename IndexSet::IndexType;
    
    /**
     * \brief Defines the available orderings
     * \todo shall this be extendable by deriving from the class, or by policies, or ...?
     */
    enum SortOrder { 
      /**
       * \brief sorts cells first, then faces, edges, and vertices last
       */
      SORTBYCODIM, 
      /**
       * \brief sorts cells according to the cell iterator of the grid view, and places lower 
       *        codim entities close to incident cells
       * This provides very good locality of access for all algorithms that step through the cells of
       * the grid view. Depending on the grid implementation, it does *not* mean that stiffness matrices
       * have a small band width.
       * 
       * NOT YET IMPLEMENTED
       */
      SEQUENTIAL, 
      /**
       * \brief sorts all entities according to their adjacency by Cuthill-McKee
       * This provides small bandwidth of stiffness matrices and hence good locality of access in
       * matrix-vector products.
       * 
       * NOT YET IMPLEMENTED
       */
      CUTHILLMCKEE };

    EntityNumbering(GridView const& gv, SortOrder order)
    {
      IndexSet const& is = gv.indexSet();
      for (int codim=0; codim<=dim; ++codim)
        indices[codim].resize(is.size(codim));
      
      switch (order)
      {
        case SORTBYCODIM:
        {
          Index n = 0;
          for (int codim=0; codim<=dim; ++codim)
          {
            std::iota(begin(indices[codim]),end(indices[codim]),n);
            n += indices[codim].size();
          }
          break;
        }
        case SEQUENTIAL:
          abort();
          break;
        case CUTHILLMCKEE:
          abort();
          break;
      }
    }
    
    Index operator()(int codim, Index i) const 
    {
      assert(0<=codim && codim<=dim);
      assert(0<=i && i<indices[codim].size());
      return indices[codim][i];
    }
    
  private:
    std::array<std::vector<Index>,dim+1> indices;
  };
  
  //----------------------------------------------------------------------------------------
  
  /**
   * \ingroup grid
   * \brief Basic functionality for managing grids and their refinement.
   */
  template<class Grd>
  class GridManagerBase
  {
  public:    
    /**
     * \brief the type of managed grid
     */
    using Grid = Grd;
    
    typedef GridManagerBase<Grid> Self;
    
    GridManagerBase(Dune::GridPtr<Grid> grid_, bool verbose_, bool enforceConcurrentReads = false)
    : GridManagerBase(grid_.release(),verbose_,enforceConcurrentReads)
    {}

    GridManagerBase(std::unique_ptr<Grid>&& grid_, bool verbose_, bool enforceConcurrentReads = false)
    : GridManagerBase(grid_.release(),verbose_,enforceConcurrentReads)
    {}

    GridManagerBase(Grid*&& grid_, bool verbose_, bool enforceConcurrentReads = false):
      markedAny(false), gridptr(grid_), 
      gridIsThreadSafe_(Dune::Capabilities::viewThreadSafe<Grid>::v), enforceConcurrentReads_(enforceConcurrentReads),
      levelRanges(gridptr->maxLevel()+1), verbose(verbose_)
    {}

    virtual ~GridManagerBase() {};
    
    /**
     * \name Access to the managed grid.
     * @{
     */

    /// Returns a const reference on the owned grid
    Grid const& grid() const
    {
      return *gridptr;
    }
    
    /// Returns a non-const reference on the owned grid
    Grid& grid_non_const()
    {
      return *gridptr;
    }

    /**
     * \brief Provides shared ownership management. 
     * 
     * Note that only an immutable grid is provided, such that the responsibility for 
     * grid modifications remains with the GridManager object. This is necessary for the 
     * GridManager to trigger prolongation/interpolation on grid modification.
     */
    std::shared_ptr<Grid const> gridShared() const
    {
      return gridptr;
    }

    /**
     * @}
     */
    
    /**
     * \name Grid adaptation
     * @{
     */
    
    /**
     * \brief Marks the given cell for refinement or coarsening
     * \param refCount if positive, marks for (multiple) refinement; if negative, for (multiple) coarsening
     * \param cell the affected cell
     */
    bool mark(int refCount, typename Grid::template Codim<0>::Entity const& cell) 
    {
      if(refCount != 0)
      {
        markedAny = true;
        return gridptr->mark(refCount,cell);
      }
      return false;
    }

    /// Marks all Entities of a grid according to the corresponding entries in a CellData
    template<class S>
    void mark(CellData<Grid,S> const& cellData)
    {
      this->flushMarks();
      typename CellData<Grid,S>::CellDataVector::const_iterator imark;
      typename CellData<Grid,S>::CellDataVector::const_iterator imark_end = cellData.end();
      for(imark=cellData.begin(); imark != imark_end; ++imark)
        if(!this->mark(imark->first,imark->second) && imark->first!=0) 
          std::cout << "GRIDMANAGER: Warning: Cannot mark element!" << std::endl;
      this->markedAny = true;
    }

    /**
     * \brief Reports the number of marked elements.
     * 
     * Elements may be marked for multiple refinement. This is reflected in the sum. The exact 
     * number returned is \f[ \sum_i r_i, \f] where \f$ r_i \f$ is the refinement count for 
     * cell \f$ i \f$.
     */
    size_t countMarked() const
    {
      size_t marks(0);
      for (auto const& cell: elements(gridptr->leafGridView()))
        marks += abs(gridptr->getMark(cell));
      return marks;
    }

    /// Calls grid.preAdapt().
    bool preAdapt() { 
      if(markedAny) 
        return gridptr->preAdapt(); 
      else 
        return false; 
    };

    /// Calls grid.postAdapt().
    void postAdapt()
    {
      gridptr->postAdapt();

      if(verbose)
        std::cout << "after refinement: " << gridptr->size(0) << " cells, " << gridptr->size(Grid::dimension) << " vertices." << std::endl;
    };


    /**
     * \brief DEPRECATED Performs grid refinement without prolongating registered FE functions
     * 
     * As this can easily lead to inconsistent states of FE functions, usage is not recommended. 
     * Probably this method will be removed soon.
     */
    bool adaptAtOnce()
    {
      this->preAdapt();
      bool result = adapt();
      #if !defined(NDEBUG)
      if(!result) std::cout << "GRIDMANAGER: Nothing has been refined!" << std::endl;
      #endif
      this->postAdapt();
      return result;
    }
    
    /// Remove all cell marks from the grid
    void flushMarks()
    {
      for(auto const& element: elements(gridptr->leafGridView()))
        mark(0,element);

      markedAny = false;
    }

    /**
     * \brief Refines each marked leaf cell of the grid and prolongates registered FE functions to the new leaf grid view
     */
    bool adapt()
    {
      bool result;
      if(this->markedAny) 
      {
        Timings& timer = Timings::instance();
        ScopedTimingSection refineSection("grid refinement");
        
        timer.start("before refinement");
        this->beforeRefinement();
        timer.stop("before refinement");
        
        if (verbose)  std::cout << std::endl << "GRIDMANAGER: mesh is refined from " << this->gridptr->size(0) << " cells " << std::flush;
        
        timer.start("adapt");
        result = adaptGrid();
        update();
        timer.stop("adapt");
        
        if(verbose) 
          std::cout << "to " << this->gridptr->size(0) << " cells" << std::endl;
        
        timer.start("after refinement");
        this->afterRefinement();
        timer.stop("after refinement");
        
        this->markedAny = false;
      } 
      else
        result = false;
      return result;
    }
    
    /**
     * \brief  Refines each leaf cell of the grid and prolongates registered FE functions to the new leaf grid view
     */
    void globalRefine(int refCount) 
    {
      if(refCount<=0) 
        return;
      
      Timings& timer = Timings::instance();
      ScopedTimingSection refineSection("grid refinement");
      timer.start("before refinement");
      this->beforeRefinement();
      timer.stop("before refinement");
      
      timer.start("adapt");
      refineGrid(refCount);
      update();
      timer.stop("adapt");
      
      timer.start("after refinement");
      this->afterRefinement();
      timer.stop("after refinement");
      
      this->markedAny=false;
    };
    
    /**
     * \brief Must be overloaded by derived classes to adjust internal state immediately after mesh refinement.
     */
    virtual void update() = 0;
    
    /**
     * @}
     */
    
    /**
     * \brief Tells the grid manager that concurrent reads of the grid are safe.
     */
    Self& enforceConcurrentReads(bool enforceConcurrentReads) 
    { 
      enforceConcurrentReads_ = enforceConcurrentReads; 
      return *this; 
    }


    /**
     * \brief sets the verbosity level of the grid manager
     *
     * The grid manager prints status and progress information to std::cout. The amount of
     * information depends on the verbosity:
     * - 0: no messages are printed
     * - 1: mesh size change on refinement is reported
     */
    Self& setVerbosity(const bool verbosity) 
    { 
      verbose = verbosity; return *this; 
    }

    /**
     * \brief Returns true if concurrent read accesses to the grid do not lead to data races.
     * 
     * This is either the case if the grid implementation itself reports that it's thread-safe, or if the 
     * user has explicitly stated that concurrent reads are safe (using enforceConcurrentReads) - even if they
     * actually are not.
     * 
     * The thread-safety reported by this method is used e.g. by the assembler to decide how many threads 
     * to use during assembly.
     * 
     * If the grid access may not be thread-safe, consider using GridLocking.
     */
    bool gridIsThreadSafe() const
    {
#ifdef NDEBUG
      return gridIsThreadSafe_ || enforceConcurrentReads_;
#else
      return false;
#endif
    }

    template<class TM, class FSE>
    struct TransFSE
    {
      TransFSE(TM const& mat_, FSE& fse_) : mat(mat_), fse(fse_) {}

      void operator()()
      {
        fse.transfer(mat);
      }

      TM const& mat;
      FSE& fse;
    };


    template<class TM, class FSE>
    void transferFSE(TM const& mat, FSE& fse) const
    {
#ifdef MULTITHREAD
      pool->schedule(TransFSE<TM, FSE>(mat,fse));
#else
      fse.transfer(mat);
#endif
    }

    template<class TD, class Space>
    struct TransBefore
    {
      TransBefore(TD* td_, Space const& space_) 
      : td(td_), space(space_) 
      {}

      void operator()()
      {
        td = new TD(space);
      }

      TD* td;
      Space const& space;
    };




    template<class TransferD, class Space>
    void getTransferData(TransferD*& result, Space const& space)
    {
#ifdef MULTITHREAD
      pool->schedule(TransBefore<TransferD, Space>(result,space));
#else
      result = (new TransferD(space));
#endif
    }


    template<class TD, class TM, class Space>
    struct TransAfter
    {
      TransAfter(TD* td_, TM* result_, Space const& space_, GridManagerBase<Grd>const & gm_) : td(td_), result(result_), space(space_), gm(gm_) {}

      void operator()()
      {
        result=td->transferMatrix().release();
        space.requestToRefine(*result,gm);
      }

      TD* td;
      TM* result;
      Space const& space;
      GridManagerBase<Grd> const& gm;
    };

    template<class TransferD, class Space>
    void getTransferMatrix(TransferD* transferData, typename TransferD::TransferMatrix*& result, Space const& space)
    {
#ifdef MULTITHREAD
      pool->schedule(TransAfter<TransferD,typename TransferD::TransferMatrix, Space>(transferData,result,space,*this));
#else
      result = transferData->transferMatrix().release();
      space.requestToRefine(*result,*this);
#endif
    }

    /**
     * \brief Returns a CellRanges object for the given grid view.
     * The cell ranges are created on demand. The reference is valid up to the next mesh modification. 
     * This method is thread-safe.
     */
    CellRanges<typename Grid::LeafGridView> const& cellRanges(typename Grid::LeafGridView const&) const
    {
      std::lock_guard<std::mutex> lock(cellRangeMutex);
      if (!leafRanges)
        leafRanges.reset(new CellRanges<typename Grid::LeafGridView>(2*NumaThreadPool::instance().cpus(),grid().leafGridView()));
      
      return *leafRanges;
    }
    
    /**
     * \brief Returns a CellRanges object for the given grid view.
     * \param gridView a level grid view of the grid managed by this grid manager
     * 
     * The cell ranges are created on demand. The reference is valid up to the next mesh modification. 
     * This method is thread-safe.
     */
    CellRanges<typename Grid::LevelGridView> const& cellRanges(typename Grid::LevelGridView const& gridView) const
    {
      int const level = gridView.template begin<0>()->level(); // cells have the level of their level view (inferred from the Dune tutorial)
      assert(level<=grid().maxLevel());
      
      std::lock_guard<std::mutex> lock(cellRangeMutex);
      if (levelRanges[level].maxRanges()<0)
        levelRanges[level] = CellRanges<typename Grid::LevelGridView>(2*NumaThreadPool::instance().cpus(),grid().levelGridView(level));
        // we could have used gridView directly here. But can we make sure this is a level view of the grid we keep?
      return levelRanges[level];
    }
    
    /**
     * \brief Returns an EntityNumbering object for the given grid view.
     * \todo allow to specify the sorting
     */
    EntityNumbering<typename Grid::LeafGridView> const& entityNumbering(typename Grid::LeafGridView const&) const
    {
      using Numbering = EntityNumbering<typename Grid::LeafGridView>;
      
      std::lock_guard<std::mutex> lock(cellRangeMutex); // maybe we should use our own mutex here?
      if (!leafEntityNumbers)
        leafEntityNumbers = std::make_unique<Numbering>(grid().leafGridView(),Numbering::SORTBYCODIM);
      
      return *leafEntityNumbers;
    }
    
    
  private:
    // TODO: docme!
    virtual void refineGrid(int refcount) = 0;
    virtual bool adaptGrid() = 0;
    
    bool markedAny; // true if any cells have been marked for refinement/coarsening, otherwise false

  protected:
    void beforeRefinement()
    {
      // trigger creation of transfer data structures
#ifdef MULTITHREAD
      pool.reset(new boost::threadpool::pool(16));
#endif
      
      signals.informAboutRefinement(GridSignals::BeforeRefinement);
      
#ifdef MULTITHREAD
      pool->wait(0);
#endif
      
      // invalidate the cell range data structures
      leafRanges.reset();
      levelRanges.clear();
      leafEntityNumbers.reset();
    }

    void afterRefinement()
    {
      levelRanges.resize(grid().maxLevel()+1);
      signals.informAboutRefinement(GridSignals::AfterRefinement);
#ifdef MULTITHREAD
      pool->wait(0);
#endif
      signals.informAboutRefinement(GridSignals::TransferCompleted);
#ifdef MULTITHREAD
      pool->wait(0);
      pool.reset();
#endif
    }

    template<class LocalToGlobalMapper> friend class FEFunctionSpace;
    template<class GridMan, class ErrEst> friend class AdaptiveGrid;

    /**
     * \brief The grid itself.
     */
    std::shared_ptr<Grid> gridptr;
    
#ifdef MULTITHREAD
    std::unique_ptr<boost::threadpool::pool> pool;
#endif

  public:
    GridSignals signals;

  private:
    const bool gridIsThreadSafe_;
    
    /// enforce concurrent reads for grids that possibly do not support concurrency
    bool enforceConcurrentReads_;
    
    // cell ranges and entity numberings are constructed on demand
    mutable std::mutex cellRangeMutex;
    mutable std::unique_ptr<CellRanges<typename Grid::LeafGridView>> leafRanges;
    mutable std::vector<CellRanges<typename Grid::LevelGridView>> levelRanges;
    mutable std::unique_ptr<EntityNumbering<typename Grid::LeafGridView>> leafEntityNumbers;

  protected:
    bool verbose;
  };

  // -------------------------------------------------------------------
  
  /** 
   * \ingroup grid
   * \brief Class that takes care of the transfer of data (of each FunctionSpaceElement) after modification of the grid
   *
   * Class that owns a grid. It takes care that all grid functions stay
   * consistent, when the grid is refined. A Gridmanager<Grid> is constructed with a
   * GridPointer<Grid>,  a std::unique_ptr<Grid>&& or an r-value reference on a Grid*. After this, all
   * modifications on the grid have to be done via the GridManager, in
   * particular refinement. If the grid is refined, the GridManager sends a
   * boost::signal to each FEFunctionSpace, which themselves send a
   * boost::signal to each FunctionSpaceElement, such that their transfer
   * from one grid to the other is organized.
   * 
   * The grid manager can print status and progress information to std::cout, see \ref setVerbosity.
   *
   * <b> Examples for usage </b>
   *
   * - Construction
   *
   * \verbatim
   * Dune::GridPtr<Grid> gridptr("MyDomain.dgf");
   * GridManager<Grid> gm(gridptr);
   * \endverbatim
   *
   * - Alternative construction
   *
   * \verbatim
   * std::unique_ptr<Grid> gridptr(new ...);
   * GridManager<Grid> gm(gridptr);
   * \endverbatim
   *
   * - Construction of FEFunctionSpaces with GridManager:
   *
   * \verbatim
   * MySpace mySpace(gm,gm.grid().leafIndexSet(),order) //instead of: MySpace mySpace(grid,grid().leafIndexSet(),order)
   * \endverbatim
   *
   * - Uniform Refinement:
   *
   * \verbatim
   * gm.globalRefine(nTimes) //instead of: grid.globalRefine(nTimes)
   * \endverbatim
   *
   * - Adaptivity
   *
   * Cf. detailed description of averaging_errorest.hh
   */
  template <class Grd>
  class GridManager : public GridManagerBase<Grd>
  {
  public:
    typedef Grd Grid;
    /// Construction from a GridPtr<Grid>
    /** Construct a GridManager<Grid> from a GridPointer<Grid>. The
     * ownership of the grid is transferred to the GridManager. Afterwards
     * the grid can only be modified via the GridManager.
     */
    explicit GridManager(Dune::GridPtr<Grid> grid_, bool verbose=false) 
    : GridManagerBase<Grid>(grid_,verbose)
    {}

    /// Construction from a unique_ptr<Grid>
    /** Construct a GridManager<Grid> from a std::unique_ptr<Grid>. The
     * ownership of the grid is transferred to the GridManager. Afterwards
     * the grid can only be modified via the GridManager.
     */
    explicit GridManager(std::unique_ptr<Grid>&& grid_, bool verbose=false) 
    : GridManagerBase<Grid>( std::move(grid_),verbose)
    {}

    /// Construction from a Grid*
    /** Construct a GridManager<Grid> from a std::unique_ptr<Grid>. The
     * ownership of the grid is transferred to the GridManager. Afterwards
     * the grid can only be modified via the GridManager.
     */
    explicit GridManager(Grid*&& grid_, bool verbose=false) 
    : GridManagerBase<Grid>( std::move(grid_),verbose)
    {}


    virtual void update() {}

  private:
    virtual void refineGrid(int refCount)
    {
      this->gridptr->globalRefine(refCount);
    };


    virtual bool adaptGrid()
    {
      return this->gridptr->adapt();
    };
  };


  //---------------------------------------------------------------------

  /**
   * \ingroup grid
   * \brief A convenience routine that refines all cells incident to a boundary face. 
   * 
   * This is useful to create a priori boundary-refined meshes.
   */
  template <class Grid>
  void refineAtBoundary(GridManagerBase<Grid>& gridManager)
  {
    for (auto const& cell: elements(gridManager.grid().leafGridView()))
      if (cell.hasBoundaryIntersections())
        gridManager.mark(1,cell);
    gridManager.adaptAtOnce();
  }
}

#endif /* GRIDMANAGER_HH_ */
