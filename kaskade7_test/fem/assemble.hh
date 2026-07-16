/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2021 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef MULTIVARIATIONAL_ASSEMBLER_HH
#define MULTIVARIATIONAL_ASSEMBLER_HH

#include <iostream>
#include <numeric>
#include <memory>
#include <mutex>
#include <tuple>
#include <new>
#include <type_traits>

#include <boost/version.hpp>
#include <boost/mpl/accumulate.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/if.hpp>
#include <boost/timer/timer.hpp>

#include <boost/exception/diagnostic_information.hpp>
#include <boost/fusion/algorithm.hpp>
#include <boost/fusion/sequence.hpp>
#include <boost/fusion/include/find_if.hpp>

#include "dune/istl/bcrsmatrix.hh"
#include "dune/istl/bdmatrix.hh"
#include "dune/common/fvector.hh"
#include "dune/common/fmatrix.hh"
#include "dune/geometry/type.hh"
#include "dune/grid/common/capabilities.hh"

#include "fem/firstless.hh"
#include "fem/fixfusion.hh"
#include "fem/functional_aux.hh"
#include "fem/gridmanager.hh"
#include "fem/quadrature.hh"
#include "fem/variables.hh"

#include "linalg/dynamicMatrix.hh"
#include "linalg/localMatrices.hh"
#include "linalg/threadedMatrix.hh"

#include "utilities/detailed_exception.hh"
#include "utilities/threading.hh"
#include "utilities/timing.hh"
#include "utilities/typeTraits.hh"


// this is omitted in boost/fusion/sequence/container/vector/detail/vector_n.hpp
#undef N

namespace Kaskade
{

  /**
   * @file
   * @brief  Compute local and global stiffness matrices for variational problems with multiple variables.
   * @author Martin Weiser
   */

  static int const timerStatistics = 0x1;       // gather global execution time statistics
  static int const scatterStatistics = 0x2;
  static int const lockStatistics = 0x4;
  static int const localTimerStatistics = 0x0;  // gather local (per element) time statistics
  /**
   * \brief flags for gathering performance and timing statistics
   * 
   * This is intended as a debugging aid. Gathering timing statistics for different assembler parts 
   * is difficult, as the granularity is small: The outer loop runs over all cells, and thus 
   * the different tasks occupy very small time spans on each cell. The timing statistics here 
   * is tailored for minimal overhead.
   * 
   * For general timing statistics, use Kaskade::Timings.
   */
  static int const statistics = 0;

  //---------------------------------------------------------------------

  /**
   * \cond internals
   *
   * This namespace contains template metaprogramming details used in
   * the VariationalFunctionalAssembler.
   */
  namespace AssemblyDetail {

    using namespace boost::fusion;


    /**
     * \brief A class that supports exclusive access to several row groups of a matrix
     *        or column vector by providing appropriate locking.
     *
     * Essentially, this is a structured group of mutexes, and therefore noncopyable.
     */
    class RowGroupManager
    {
    public:
      /**
       * \brief Default constructor creates a (pretty useless) empty row group manager
       */
      RowGroupManager();

      /**
       * \brief Destructor
       */
      virtual ~RowGroupManager();

      // Noncopyable since we contain mutexes.
      RowGroupManager(RowGroupManager const& m) = delete;
      RowGroupManager& operator=(RowGroupManager const& m) = delete;

      /**
       * \brief initializes the row group manager to row groups of approximately the same number of rows
       * \param[in] nrg the desired number of row groups (typically the number of hardware threads)
       * \param[in] rows the total number of rows
       */
      void init(int nrg, size_t rows);

      /**
       * \brief initializes the row group manager to row groups of approximately 
       *        the same number of matrix entries
       * \param[in] nrg the desired number of row groups (typically the number of hardware threads)
       * \param[in] colIndex a container of column indices for matrix entries
       */
      void init(int nrg, std::vector<std::vector<size_t>> const& colIndex);

      /**
       * \brief initializes the row group manager to row groups of approximately the same weight
       * \param[in] nrg the desired number of row groups (typically the number of hardware threads)
       * \param[in] weight a container of row weights
       *
       * In a sparse matrix, the row weight might just be the number of entries in the row.
       */
      void init(int nrg, std::vector<size_t> const& weight);

      /**
       * \brief Queries the row range in group \a n.
       * \return the half open range of row indices belonging to this group.
       */
      std::pair<size_t,size_t> range(int n) const;

      /**
       * \brief Provides access to the mutex of row group \a n.
       */
      std::mutex& mutex(int n) const;

      /**
       * \brief Returns the number of managed row groups.
       */
      int size() const;


    private:
      // start indices of row groups
      std::vector<size_t> rowGroupStart;

      // Support for simultaneous write access to different row groups.
      std::mutex* mutexes;

      // deletes the mutexes
      void clear();
    };

    // ----------------------------------------------------------------------------------------------------

    /**
     * \brief A class that describes a matrix block and holds the associated global matrix.
     *
     * \tparam Policy \todo docme
     * \tparam RowVar \todo docme
     * \tparam ColVar \todo docme
     *
     * Despite holding the global stiffness matrix by a shared pointer,
     * the copy and assignment are deep (i.e. expensive).
     */
    template <class Policy, class RowVar, class ColVar, class SparseIdx>
    struct MatrixBlock
    {
      static int const rowId = RowVar::id;
      static int const colId = ColVar::id;
      static int const rowSpaceIndex = RowVar::spaceIndex;
      static int const colSpaceIndex = ColVar::spaceIndex;
      static bool const symmetric    = Policy::template BlockInfo<rowId,colId>::symmetric;
      static bool const mirror       = Policy::template BlockInfo<rowId,colId>::mirror;
      static bool const lumped       = Policy::template BlockInfo<rowId,colId>::lumped;
      static bool const makePositive = Policy::template BlockInfo<rowId,colId>::makePositive;

      using SparseIndex = SparseIdx;
      typedef MatrixBlock<Policy,RowVar,ColVar,SparseIndex> Self;
      typedef typename Policy::Scalar Scalar;
      typedef typename Policy::Spaces Spaces;
      typedef Dune::FieldMatrix<Scalar,RowVar::m,ColVar::m> BlockType;
      typedef NumaBCRSMatrix<BlockType, SparseIndex> Matrix;

      typedef typename SpaceType<Spaces,RowVar::spaceIndex>::type RowSpace;
      typedef typename SpaceType<Spaces,ColVar::spaceIndex>::type ColSpace;
      typedef typename RowSpace::GridView GridView;
      typedef typename GridView::template Codim<0>::Iterator CellIterator;

      static int const dim = ColSpace::dim;

      /**
       * \brief Default constructor creates a (pretty useless) empty matrix block.
       */
      MatrixBlock() {}

      /**
       * \brief Deep copy constructor.
       */
      MatrixBlock(MatrixBlock const& mb)
      : matrix(new Matrix(*mb.matrix))
      , isDense(mb.isDense)
      {}

      /**
       * \brief Deep copy.
       */
      MatrixBlock& operator=(MatrixBlock const& mb)
      {
        *matrix = *mb.matrix;
        isDense = mb.isDense;
      }

      MatrixBlock& operator=(Scalar x)
      {
        assert(matrix.get());
        *matrix = x;
        return *this;
      }

      MatrixBlock& operator+=(MatrixBlock const& mb)
      {
        *matrix += *mb.matrix;
        return *this;
      }


      /**
       * \brief Provides read access to the global matrix.
       */
      Matrix& globalMatrix()
      {
        assert(matrix.get());
        return *matrix;
      }

      /**
       * \brief Provides write access to the global matrix.
       */
      Matrix const& globalMatrix() const
      {
        assert(matrix.get());
        return *matrix;
      }

      /**
       * \brief Provides persistent read access to the matrix.
       *
       * The matrix remains
       * valid even on destruction of the assembler or calling flush(),
       * but may be overwritten by the next assembly call.
       */
      std::shared_ptr<Matrix const> globalMatrixPointer() const { return matrix; }


      /**
       * \brief Computes the sparsity pattern of the matrix block. For each cell
       * (codim 0 entity) of the mesh, all shape functions on the cell are
       * assumed to interact, hence their associated global indices induce
       * a scattering of nonzero entries into the matrix.
       *
       * \param considerFace a callable that, for any face in the grid view, returns
       *                     true if and only if integration shall be performed on
       *                     that face. If it returns false, no matrix entries for
       *                     coupling the dofs on the incident cells are created.
       *
       * For symmetric matrix blocks, only the lower triangular part is stored.
       */
      template <class FaceOracle>
      void init(Spaces const& spaces, CellIterator first, CellIterator last,
                FaceOracle const& considerFace)
      {
        RowSpace const& rowSpace = *at_c<RowVar::spaceIndex>(spaces);
        ColSpace const& colSpace = *at_c<ColVar::spaceIndex>(spaces);

        // @todo: respect *dynamic* presence flag. static one is respected
        // on construction of sequence of matrix blocks. Respect symmetric flag.
        if (!Policy::template BlockInfo<rowId,colId>::present)
        {
          assert("Aieee! Nonpresent matrix block detected!\n"==0);
          abort();
        }

        boost::timer::cpu_timer timer;

        size_t nnz = 0;
        size_t const rows = rowSpace.degreesOfFreedom();
        size_t const cols = colSpace.degreesOfFreedom();

        // If mass lumping is desired, we just create a diagonal matrix.
        // Otherwise, we need to establish the connectivity pattern of
        // the ansatz functions.
        if (Policy::template BlockInfo<rowId,colId>::lumped)
        {
          assert(RowVar::spaceIndex==ColVar::spaceIndex && rows==cols);
          NumaCRSPatternCreator<SparseIndex> creator(rows,cols,false,1);
          for (size_t i=0; i<rows; ++i)
            creator.addElements(&i,&i+1,&i,&i+1);
          matrix.reset(new Matrix(creator)); // zero init
          nnz = rows;
        }
        else
        {
          // In case of global support of shape functions of one of the involved spaces (usually ConstantSpace),
          // the resulting matrix is dense. Use a shortcut to define the sparsity structure.
          bool const dense = RowSpace::Mapper::globalSupport || ColSpace::Mapper::globalSupport;


          // Allocate sparse matrix and define the sparsity structure. Estimate the required number of entries per row
          // for efficient preallocation.
          // Ws 2016-01-22: Zero nnz init is actually faster on P1 pattern, probably because otherwise the
          // chunks put a high burden on the memory system by requesting memory at the same time. This appears
          // to parallelize quite badly and hence is turned off in parallel mode.
          int preallocateEntriesPerRow = 0;
          if (NumaThreadPool::instance().nodes()==1)                       // sequential allocation in NumaCRSPatternCreator
          {
            int localSize = colSpace.mapper().globalIndices(0).size();     // # entries per row of local matrices (rough estimate)
            preallocateEntriesPerRow = dense? 0:                           // for dense matrices we'll use a specialized version
                                       dim<=1? 2*localSize:                // make room for (roughly) all elemental matrices appended...
                                       dim==2? 8*localSize:                // ... this is overprovisioning, but faster due to less ...
                                               20*localSize;               // ... allocations and data movement.
          }
          NumaCRSPatternCreator<SparseIndex> creator(rows,cols,symmetric,preallocateEntriesPerRow);

          if (statistics & timerStatistics)
            std::cout << "initial creator construction (" << rowId << "," << colId << "):  " << timer.format();
          timer.start();


          if (dense)
          {
            creator.addAllElements();
            if (statistics & timerStatistics)
              std::cout << "entering indices for (" << rowId << "," << colId << "):        " << timer.format();
            timer.start();
          }
          else
          {
            // Traditional local FE spaces - that's the hard case. We have to traverse the grid and pick up all
            // interactions of (local) basis functions.

            // Iteration over all cells. Note that the grid views of row and
            // column spaces need not cover the whole grid (e.g., in case of
            // local support). It would be best to iterate over the grid
            // view with smaller support, but for simplicity we choose the
            // row space.

            size_t const nCells = rowSpace.indexSet().size(0);

            std::vector<typename RowSpace::Mapper::GlobalIndexRange> rIndices(nCells,rowSpace.mapper().initGlobalIndexRange());
            std::vector<typename ColSpace::Mapper::GlobalIndexRange> cIndices(nCells,colSpace.mapper().initGlobalIndexRange());
            for (size_t i=0; i<nCells; ++i)
            {
              rIndices[i] = rowSpace.mapper().globalIndices(i);
              cIndices[i] = colSpace.mapper().globalIndices(i);
            }
            if (statistics & timerStatistics)
              std::cout << "gathering indices for (" << rowId << "," << colId << "):         " << timer.format();
            timer.start();

            creator.addElements(rIndices,cIndices); // enter all entries in one go

            // We may need to assemble integrals on inner faces - then there are matrix entries coupling degrees of freedom 
            // associated to different cells (e.g., in discontinuous Galerkin methods). Add these entries, too.
            if(Policy::considerInnerFaces)
            {
              GridView const& gridView = rowSpace.gridView();
              
              for (CellIterator ci=first; ci!=last; ++ci)
              {
                auto const& rIndices = rowSpace.mapper().globalIndices(*ci);

                for (auto const& face: intersections(gridView,*ci))
                {
                  // only true inner faces are considered, and only if requested
                  if (face.boundary() || !face.neighbor() || considerFace(face)==false)
                    continue;
                  
                  auto const& cIndices = colSpace.mapper().globalIndices(face.outside());

                  creator.addElements(begin(rIndices),end(rIndices),
                                      begin(cIndices),end(cIndices));
                } // end iteration over faces
              } // end iteration over cells
            } // end Policy::considerInnerFaces

            if (statistics & timerStatistics)
              std::cout << "entering indices for (" << rowId << "," << colId << "):          " << timer.format();
            timer.start();

            creator.balance();
            if (statistics & timerStatistics)
              std::cout << "balancing for (" << rowId << "," << colId << "):                 " << timer.format();
            timer.start();
          }

          // Now that the pattern creator is filled, create pattern and matrix.

          // Count number of nonzeros.
          nnz = creator.nonzeroes();

          // Create the matrix (zero initialization of entries).
          matrix.reset(new Matrix(creator));
          if (statistics & timerStatistics)
            std::cout << "creating pattern & matrix for (" << rowId << "," << colId << "): " << timer.format();
          timer.start();
        }

        // Check density (nonzeros>(rows*cols/2). Be careful not to get
        // overflow for very large and sparse matrices!
        if (std::min(rows,cols)>0)
          isDense = nnz/std::min(rows,cols) > std::max(rows,cols)/2;
        else
          isDense = false; // well, empty matrices are as sparse as one can hope for ;)

        if (statistics & timerStatistics)
            std::cout << "init cleanup for (" << rowId << "," << colId << "):              " << timer.format() << "\n";
      }

      // Returns true if the matrix is essentially dense (more than 50% nonzeros).
      bool dense() const { return isDense; }

    private:
      std::shared_ptr<Matrix> matrix;
      bool                    isDense;
    };

    // ---------------------------------------------------------------------------------------

    // a boost::fusion predicate telling us whether the current matrix block has the given row and column indices.
    template <int row, int col>
    struct IsBlock
    {
      template <class MatrixBlock>
      struct apply
      {
        using type = std::conditional_t<MatrixBlock::rowId==row && MatrixBlock::colId==col,boost::mpl::true_,boost::mpl::false_>;
      };
    };

    // ---------------------------------------------------------------------------------------

    /**
     * \brief defines a flat sequence of matrix blocks
     */
    template <class Policy, class AnsatzVariables, class TestVariables, class SparseIndex>
    class BlockArray
    {
      // Retain only those blocks which are present in the functional
      struct PresenceFilter
      {
        template <class Block>
        struct apply
        {
          typedef boost::mpl::bool_<Policy::template BlockInfo<Block::rowId,Block::colId>::present> type;
        };
      };

      static auto create(AnsatzVariables a, TestVariables t)
      {
        auto joiner = [=](auto blocks, auto rowVar)
        {
          using RowVar = decltype(rowVar);
          auto thisRow = transform(a,[](auto avar){ return MatrixBlock<Policy,RowVar,decltype(avar),SparseIndex>(); });
          return join(blocks,thisRow);
        };
        auto allBlocks = accumulate(t,vector<>(),joiner);
        return as_vector(filter_if<PresenceFilter>(allBlocks));
      }

    public:
      typedef decltype(create(AnsatzVariables(),TestVariables())) type;
    };


    //---------------------------------------------------------------------

    /**
     * \brief A class that defines one block of the right hand side
     */
    template <class Policy, class RV>
    struct RhsBlock
    {
      typedef typename std::remove_reference<RV>::type RowVar;

      static int const rowId = RowVar::id;

      typedef typename Policy::Scalar Scalar;
      typedef typename Policy::Spaces Spaces;
      typedef Dune::FieldVector<Scalar,RowVar::m> BlockType;
      typedef typename Dune::BlockVector<BlockType> Rhs;
      typedef RhsBlock<Policy,RowVar> Self;

      static int const rowSpaceIndex = RowVar::spaceIndex;
      typedef typename SpaceType<Spaces,RowVar::spaceIndex>::type RowSpace;
      typedef typename RowSpace::GridView::template Codim<0>::Iterator CellIterator;


      /**
       * \brief Initializes the rhs block info.
       */
      void init(Spaces const& spaces, CellIterator /*first*/, CellIterator /*last*/, int nrg)
      {
        RowSpace const& rowSpace = *at_c<RowVar::spaceIndex>(spaces);
        size_t rows = rowSpace.degreesOfFreedom();

        rowGroupManager = std::make_unique<RowGroupManager>();
        rowGroupManager->init(nrg,rows);
      }

      RowGroupManager& rowGroup() { return *rowGroupManager; }

      // A class storing a local element rhs vector and the corresponding indices
      struct LocalVector
      {
        // just because the ranges are not default constructible
        LocalVector(): ridx(RowSpace::Mapper::initSortedIndexRange()) {}

        std::vector<BlockType> vector;

        typedef typename RowSpace::Mapper::SortedIndexRange SortedRowIdx;

        /**
         * \brief A container of (global row, local row) indices.
         */
        SortedRowIdx ridx;
      };

      // A structure for holding a sequence of several local vectors to be filled sequentially
      // and to be scattered together.
      struct LocalVectors
      {
        typedef Self RhsBlock;
        typedef std::vector<BlockType> LocalVectorType;

        LocalVectors(int n, RhsBlock& rb_): localVectors(n), currentLocalVector(0), rb(&rb_) {}

        // Storage for local rhs vectors.
        std::vector<LocalVector> localVectors;

        // Which is the current local vector to be filled?
        int currentLocalVector;

        // Pointing back to its rhs block
        RhsBlock* rb;
      };

    private:
      std::unique_ptr<RowGroupManager> rowGroupManager;
    };

    /**
     * \brief Defines a flat sequence of present rhs blocks
     */
    template <class Policy, class TestVariables>
    class RhsBlockArray
    {
      // Define a mapping from Variable to RhsBlock for Variable as row variable.
      // Retain only those blocks which are present in the functional
      struct PresenceFilter
      {
        template <class Block>
        struct apply
        {
          typedef boost::mpl::bool_<Policy::template RhsBlockInfo<Block::rowId>::present> type;
        };
      };

      static auto create(TestVariables t)
      {
        auto allBlocks = transform(t,[](auto tvar) { return RhsBlock<Policy,decltype(tvar)>(); });
        return as_vector(filter_if<PresenceFilter>(allBlocks));
      }

    public:
      typedef decltype(create(TestVariables())) type;
    };

    //---------------------------------------------------------------------

    /**
     * \brief A boost::fusion functor that for a given test variable defines the type of local rhs blocks.
     * \TODO rewrite this to use multiple local rhs vectors to be scattered simultaneously (as already done for matrix blocks)
     */
    template <class Policy>
    struct RhsLocalData
    {
      /**
       * \param[in] n_ desired number of local vectors to complete before scattering
       */
      RhsLocalData(int n_): n(n_) {}

      template <class RhsBlock>
      typename RhsBlock::LocalVectors operator()(RhsBlock const& rb) const
      {
        // due to boost::fusion::transform providing only const references to the source sequence elements, we
        // need a const cast here :(
        return typename RhsBlock::LocalVectors(n,removeConst(rb));
      }

    private:
      int n;
    };
    
    /**
     * \brief A boost::fusion functor that defines for each matrix block the local matrices storage for assembly.
     */
    struct MatrixLocalData
    {
      /**
       * \brief Constructor.
       * \param s The memory size in bytes the local matrices should not exceed.
       */
      MatrixLocalData(size_t s_)
      : s(s_) {}

      template <class MatrixBlock>
      auto operator()(MatrixBlock const& mb) const
      {
        // When switching from boost 1.51 to 1.57, it seems that the filter_view only presents
        // const references to the elements (why?), but we need non-const ones. Taking non-const
        // reference as arguments leads to failure of overload resolution, omitting this
        // functor, and ultimately to error. As a workaround, we take const reference and cast
        // it away.
        return LocalMatrices<typename MatrixBlock::BlockType,
                             MatrixBlock::lumped,
                             typename MatrixBlock::RowSpace::Mapper::SortedIndexRange,
                             typename MatrixBlock::ColSpace::Mapper::SortedIndexRange,
                             MatrixBlock>(removeConst(mb).globalMatrix(),s);
      }

    private:
      size_t s;
    };

    //---------------------------------------------------------------------

    template <class Spaces>
    class GetMaxDerivativeBase
    {
    public:
      GetMaxDerivativeBase(Spaces const& spaces_, std::map<void const*,int>& deriv_): spaces(spaces_), deriv(deriv_) {}

    private:
      Spaces const&        spaces;
      std::map<void const*,int>& deriv;

    protected:
      template <int spaceIndex>
      void enterSpace(int d) const
      {
        void const* space = at_c<spaceIndex>(spaces);
        auto i = deriv.find(space);
        if (i==deriv.end())
          deriv[space] = d;
        else
          i->second = std::max(i->second,d);
      }
    };

    template <class Functional, class Spaces>
    class GetMaxDerivativeMatrix: public GetMaxDerivativeBase<Spaces>
    {
    public:
      GetMaxDerivativeMatrix(Spaces const& spaces, std::map<void const*,int>& deriv): GetMaxDerivativeBase<Spaces>(spaces,deriv)  {}

      template <class LocalMatrices>
      void operator()(LocalMatrices const&) const
      {
        using MatrixBlock = typename LocalMatrices::InfoType;
        int d = Functional::template D2<MatrixBlock::rowId,MatrixBlock::colId>::derivatives;
        this->template enterSpace<MatrixBlock::rowSpaceIndex>(d);
        this->template enterSpace<MatrixBlock::colSpaceIndex>(d);
      }
    };


    template <class Functional, class Spaces>
    class GetMaxDerivativeRhs: public GetMaxDerivativeBase<Spaces>
    {
    public:
      GetMaxDerivativeRhs(Spaces const& spaces, std::map<void const*,int>& deriv): GetMaxDerivativeBase<Spaces>(spaces,deriv)  {}

      template <class RBlock>
      void operator()(RBlock const&) const
      {
        int d = Functional::template D1<RBlock::RhsBlock::rowId>::derivatives;
        this->template enterSpace<RBlock::RhsBlock::rowSpaceIndex>(d);
      }
    };


    /**
     * \brief Computes the maximum occuring derivative for each space.
     *
     * \returns a map from space address to the highest occuring derivative order
     */
    template <class Functional, class Spaces, class LocalMatrices, class RBlocks>
    std::map<void const*,int> derivatives(Functional const& f, Spaces const& spaces, LocalMatrices const& localMatrices, RBlocks const& rblocks)
    {
      std::map<void const*,int> deriv; // map from space address to maximum occuring derivative

      for_each(localMatrices,GetMaxDerivativeMatrix<Functional,Spaces>(spaces,deriv));
      for_each(rblocks,GetMaxDerivativeRhs<Functional,Spaces>(spaces,deriv));
      return deriv;
    }

    //---------------------------------------------------------------------

    /**
     * \brief A boost::fusion functor that, given a test variable, initializes the current local rhs vector tor correct size and 0 values.
     */
    template <class Evaluators, class TestVariables>
    struct ClearLocalRhs
    {
      ClearLocalRhs(Evaluators const& evaluators_):
        eval(evaluators_)
      {}

      template <class LocalVectors>
      void operator()(LocalVectors& lv) const
      {
        // Deletes all entries from vector and inserts the needed number
        // of zero elements. Thus, finite elements with different number
        // of shape functions are supported. Memory reallocations do
        // rarely occur.
        assert(lv.currentLocalVector <  lv.localVectors.size());

        int const rSpaceIdx = result_of::value_at_c<TestVariables,LocalVectors::RhsBlock::rowId>::type::spaceIndex;

        auto& reval = at_c<rSpaceIdx>(eval);                    // evaluator

        auto& rc = lv.localVectors[lv.currentLocalVector];        // current local vector structure

        // set current local vector to correct size and value 0
        rc.vector.resize(reval.size());
        std::fill(rc.vector.begin(),rc.vector.end(),0);

        // retrieve sorted indices
        rc.ridx = reval.sortedIndices();
      }

      Evaluators const&  eval;
    };

    /**
     * \brief A boost::fusion functor that creates a new local stiffness matrix initialized to zero.
     */
    template <class Evaluators, class AnsatzVariableSetDescription, class TestVariableSetDescription>
    struct NewLocalMatrix
    {
      NewLocalMatrix(Evaluators const& reval_, Evaluators const& ceval_)
        : reval(reval_), ceval(ceval_)
      {}

      template <class LocalMatrices>
      void operator()(LocalMatrices& lm) const
      {
        using MatrixBlock = typename LocalMatrices::InfoType;

        int const cSpaceIdx = spaceIndex<AnsatzVariableSetDescription,MatrixBlock::colId>;
        int const rSpaceIdx = spaceIndex<TestVariableSetDescription  ,MatrixBlock::rowId>;

        lm.push_back(at_c<rSpaceIdx>(reval).sortedIndices(),at_c<cSpaceIdx>(ceval).sortedIndices());
      }

      Evaluators const& reval;
      Evaluators const& ceval;
    };


    //---------------------------------------------------------------------

    /**
     * \brief A boost::fusion functor that evaluates local right hand side contributions
     */
    template <class TestVariableSetDescription, class Evaluators, class Real, class Cache>
    struct UpdateLocalRhs
    {
      UpdateLocalRhs(Evaluators const& evaluators_, Real integrationFactor_, Cache const& cache_)
      : evaluators(evaluators_)
      , integrationFactor(integrationFactor_)
      , cache(cache_)
      {}

      template <class LocalVectors>
      void operator()(LocalVectors& lv) const
      {
        typedef typename LocalVectors::RhsBlock RhsBlock;

        int const rSpaceIndex = spaceIndex<TestVariableSetDescription,RhsBlock::rowId>;
        auto& rc = lv.localVectors[lv.currentLocalVector].vector; // current local vector
        size_t const rows = rc.size();
        for (size_t i=0; i<rows; ++i)
        {
          // Check that the cache returns the correct type. The typedefs just serve for
          // creating readable error messages.
          using rhsEntryType = std::remove_reference_t<decltype(rc[i])>;
          using d1ResultType = decltype(cache.template d1<RhsBlock::rowId>(at_c<rSpaceIndex>(evaluators).evalData[i]));
          static_assert(std::is_same<rhsEntryType,d1ResultType>::value,
                        "Types don't match. Check return type of d1() in problem formulation.");

          rc[i] += integrationFactor * cache.template d1<RhsBlock::rowId>(at_c<rSpaceIndex>(evaluators).evalData[i]);
          assert(!std::isnan(rc[i][0])); // check that first entry in rhs vector is valid
        }
      }

      Evaluators const&    evaluators;
      Real                 integrationFactor;
      Cache const&         cache;
    };


    /**
     * \brief A boost::fusion functor which adds up contributions to the local stiffness matrices during the integration
     */
    template <class AnsatzVariableSetDescription, class TestVariableSetDescription, class Evaluators, class Real, class Cache>
    struct UpdateLocalMatrix
    {
      UpdateLocalMatrix(Evaluators const& evaluators_, Real integrationFactor_, Cache const& cache_):
        evaluators(evaluators_), integrationFactor(integrationFactor_), cache(cache_)
      {}

      template <class LocalMatrices>
      void operator()(LocalMatrices& lm) const
      {
        using MatrixBlock = typename LocalMatrices::InfoType;

        // Get the evaluated ansatz (column) and test (row) functions
        int const rowId =    MatrixBlock::rowId;
        int const colId =    MatrixBlock::colId;
        int const rSpaceIndex = spaceIndex<TestVariableSetDescription,rowId>;
        int const cSpaceIndex = spaceIndex<AnsatzVariableSetDescription,colId>;
        auto const& rEval = at_c<rSpaceIndex>(evaluators).evalData;
        auto const& cEval = at_c<cSpaceIndex>(evaluators).evalData;

        auto& A = lm.back();  // the current local Galerkin matrix to contribute to

        // retrieve number of localized ansatz functions
        int const rows = A.ridx().size();
        int const cols = A.cidx().size();

        // loop over all combinations of localized ansatz and test functions
        for (int i=0; i<rows; ++i)
        {
          int begin = MatrixBlock::lumped? i: 0;
          int end = (MatrixBlock::symmetric||MatrixBlock::lumped)? i+1: cols;
          for (int j=begin; j<end; ++j)
          {
            static_assert(std::is_same<std::remove_reference_t<decltype(A(i,j))>,
                                       decltype(cache.template d2<rowId,colId>(rEval[i],cEval[j]))
                          >::value, "Types don't match. Check return type of d2() in problem formulation.");
            A(i,j) += integrationFactor * cache.template d2<rowId,colId>(rEval[i], cEval[j]);
            assert(!std::isnan(A(i,j)[0][0])); // check that top left entry is valid
          }
        }
      }

      Evaluators const& evaluators;
      Real              integrationFactor;
      Cache const&      cache;
    };


    // TODO: factorize common functionality out of both UpdateLocalMatrix classes
    template <class AnsatzVariableSetDescription, class TestVariableSetDescription,
              class Evaluators, class Real, class Cache>
    struct UpdateLocalMatrixFromInnerBoundaryCache
    {
      UpdateLocalMatrixFromInnerBoundaryCache(Evaluators const& insideEval, Evaluators const& outsideEval, 
                                              Real integrationFactor_, Cache const& cache_, bool sameCell_)
      : insideEvaluator(insideEval)
      , outsideEvaluator(outsideEval)
      , integrationFactor(integrationFactor_)
      , cache(cache_)
      , sameCell(sameCell_)
      {}

      template <class LocalMatrices>
      void operator()(LocalMatrices& lm) const
      {
        using MatrixBlock = typename LocalMatrices::InfoType;

        int const rSpaceIndex = spaceIndex<TestVariableSetDescription,  MatrixBlock::rowId>;
        int const cSpaceIndex = spaceIndex<AnsatzVariableSetDescription,MatrixBlock::colId>;

        auto& A = lm.localMatrices.back();

        // retrieve number of localized ansatz functions
        int const rows = A.ridx().size();
        int const cols = A.cidx().size();

        // loop over all combinations of localized ansatz and test functions
        for (int i=0; i<rows; ++i)
        {
          int begin = MatrixBlock::lumped? i: 0;
          int end = ((MatrixBlock::symmetric && sameCell)||MatrixBlock::lumped)? i+1: cols;
          for (int j=begin; j<end; ++j)
          {
            auto tmp = cache.template d2<MatrixBlock::rowId,MatrixBlock::colId>(
                at_c<rSpaceIndex>(insideEvaluator).evalData[i],
                at_c<cSpaceIndex>(outsideEvaluator).evalData[j], sameCell);
            A(i,j) += integrationFactor * tmp;
            assert(!std::isnan(A(i,j)[0][0])); // check that top left entry is valid
          }
        }
      }

      Evaluators const& insideEvaluator, outsideEvaluator;
      Real              integrationFactor;
      Cache const&      cache;
      bool              sameCell; 
    };
    //---------------------------------------------------------------------

    /**
     * \brief A boost::fusion functor that scatters a set of local rhs vectors into the global ones
     */
    template <class GlobalRhs>
    struct ScatterLocalRhs
    {
      /**
       * \brief Constructor.
       * \param globalRhs
       * \param immediate if true, scatter all available local rhs immediately. Otherwise,
       *                  scatter only if the local rhs buffer is full.
       */
      ScatterLocalRhs(GlobalRhs& globalRhs_, bool immediate_ = false):
        globalRhs(globalRhs_), immediate(immediate_)
      {}

      template <class LocalVectors>
      void operator()(LocalVectors& lv) const
      {
        auto& gRhs = at_c<LocalVectors::RhsBlock::rowId>(globalRhs.data);

        auto const& rowGroups = lv.rb->rowGroup();

        if (immediate || lv.currentLocalVector+1==lv.localVectors.size())
        {
          int const last = immediate ? lv.currentLocalVector : lv.localVectors.size();

          // step through all row groups and scatter all the included local vectors' components
          for (int rg=0; rg<rowGroups.size(); ++rg)
          {
            // obtain and lock the row group
            auto const [firstRow,lastRow] = rowGroups.range(rg);
            std::lock_guard lock(rowGroups.mutex(rg));

            // step through all local vectors
            for (int i=0; i<last; ++i)
            {
              auto& ri = lv.localVectors[i].vector;
              auto& ridx = lv.localVectors[i].ridx;

              assert(ri.size()==ridx.size());

              for (int r=0; r<ridx.size(); ++r)
              {
                size_t const rGlobalIndex = ridx[r].first;

                // check if this row is inside the row group
                if (firstRow<=rGlobalIndex && rGlobalIndex<lastRow)
                {
                  size_t const rLocalIndex = ridx[r].second;
                  gRhs[rGlobalIndex] += ri[rLocalIndex];
                }
              }
            }
          }

          // everything's scattered - clean up
          lv.currentLocalVector = 0;
        }
        else
          ++lv.currentLocalVector; // prepare next local vector for updating
      }

      GlobalRhs&        globalRhs;
      bool              immediate;
    };

    //---------------------------------------------------------------------------------------------

    /**
     * \brief A boost::fusion functor for symmetrizing local Galerkin matrices.
     *
     * This symmetrizes the local Galerkin matrix in case it should be
     * symmetric (in which case only the lower triangular part is
     * assembled).
     *
     * This functor does only work on local Galerkin
     * matrices and is therefore thread-safe.
     */
    template <class Functional>
    struct SymmetrizeLocalMatrix
    {
      template <class LocalMatrices>
      void operator()(LocalMatrices& lm) const
      {
        using MatrixBlock = typename LocalMatrices::InfoType;

        // If a block is symmetric, only the lower part of the local
        // matrix is filled. Mirror to the upper part for easy access in
        // the following loop. For lumped (purely diagonal) blocks, of course,
        // nothing is to do.
        if constexpr (!MatrixBlock::lumped && MatrixBlock::symmetric)
        {
          auto& m = lm.back();
          int const N = m.ridx().size();
          assert(N==m.cidx().size());

          // copy lower triangular part to upper triangular part, transposing
          // each entry
          for (int i=0; i<N-1; ++i)
            for (int j=i+1; j<N; ++j)
              m(i,j) = transpose(m(j,i));
            
          // If we are asked to enforce positive semidefiniteness of these local
          // matrices (symmetric and non-lumped), we do this using an eigendecomposition,
          // and shifting all nonpositive eigenvalues to zero.
          if constexpr (MatrixBlock::makePositive)
          {
            assert(f);
            constexpr int row = MatrixBlock::rowId;
            constexpr int col = MatrixBlock::colId;

            // First we need to extract the local matrix.
            using Entry = typename MatrixBlock::BlockType;
            DynamicMatrix<Entry> Aloc(N,N);
            for (int j=0; j<N; ++j)
              for (int i=0; i<N; ++i)
                Aloc[i][j] = m(i,j);

            // now we set all negative eigenvalues to zero in the eigendecomposition
            bool changed = makePositiveSemiDefinite(Aloc,f->template makePositiveThreshold<row,col>());

            // now we copy the matrix back
            if (changed)
              for (int j=0; j<N; ++j)
                for (int i=0; i<N; ++i)
                  m(i,j) = Aloc[i][j];
          }
        }
      }

      Functional const* f;
    };

    //---------------------------------------------------------------------
    //---------------------------------------------------------------------

    /**
     * \brief provides some meta information for assembly problems (here for variational functionals)
     *
     * In contrast to weak formulations, variational functionals provide a method d0() to evaluate the
     * scalar value of the functional. The matrix is always symmetric (it's a second derivative).
     */
    template <class Functional>
    struct VariationalFunctionalPolicy
    {
      typedef typename Functional::Scalar                Scalar;
      typedef typename Functional::AnsatzVars::Spaces    Spaces;
      static constexpr bool considerInnerFaces = hasInnerBoundaryCache<Functional>;

      template <int rowId>
      struct RhsBlockInfo
      {
        static bool const present = Functional::template D1<rowId>::present;
      };


      template <int rowId, int colId>
      struct BlockInfo
      {
        static bool const present = Functional::template D2<rowId,colId>::present && rowId>=colId;
        static bool const symmetric = (rowId==colId || Functional::template D2<rowId,colId>::symmetric)
                      && (result_of::value_at_c<typename Functional::TestVars::Variables,rowId>::type::spaceIndex
                          == result_of::value_at_c<typename Functional::AnsatzVars::Variables,colId>::type::spaceIndex);
        static bool const mirror = rowId>colId;
        static bool const lumped = symmetric && Functional::template D2<rowId,colId>::lumped;
        static bool const makePositive = symmetric && Functional::template D2<rowId,colId>::makePositive;
      };

      /**
       * Returns the functional's value that has been assembled in a
       * previous call to assemble() with flags|VALUE==true. If no such
       * call to assemble() occured before, the value is undefined.
       */
      Scalar functional() const { return fValue; }

      static constexpr bool hasValue = true;

      Scalar fValue;

      template <class Cache>
      static Scalar d0(Cache const& cache) { return cache.d0(); }
    };

    /**
     * \brief provides some meta information for assembly problems (here for weak formulations)
     *
     * In contrast to variational functional problems, weak formulations are not induced by real-valued
     * functionals. As a consequence, there is no functional value d0() to be evaluated, and the matrix
     * need not be symmetric.
     */
    template <class Functional>
    struct WeakFormulationPolicy
    {
      typedef typename Functional::Scalar                Scalar;
      typedef typename Functional::AnsatzVars::Spaces    Spaces;
      static constexpr bool considerInnerFaces = hasInnerBoundaryCache<Functional>;

      template <int rowId>
      struct RhsBlockInfo
      {
        static bool const present = Functional::template D1<rowId>::present;
      };


      template <int rowId, int colId>
      struct BlockInfo
      {
        static bool const present = Functional::template D2<rowId,colId>::present;
        static bool const symmetric = Functional::template D2<rowId,colId>::symmetric
            && (result_of::value_at_c<typename Functional::TestVars::Variables,rowId>::type::spaceIndex
                == result_of::value_at_c<typename Functional::AnsatzVars::Variables,colId>::type::spaceIndex);
        static bool const mirror = false;
        static bool const lumped = symmetric && Functional::template D2<rowId,colId>::lumped;
        static bool const makePositive = symmetric && Functional::template D2<rowId,colId>::makePositive;
      };

    public:
      static constexpr bool hasValue = false;

      Scalar fValue = 0;

      template <class Cache>
      static Scalar d0(Cache const& cache) { return 0; }
    };


    /**
     * \ingroup assembly 
     * \brief A base class for the assembler adapted to either variational functionals or weak formulations.
     */
    template <class Problem>
    using FormulationPolicy = std::conditional_t<Problem::type==VariationalFunctional,
                                                 VariationalFunctionalPolicy<Problem>,
                                                 WeakFormulationPolicy<Problem>>;

    // -------------------------------------------------------------------------------------------------------

    /**
     * \ingroup assembly
     * \brief A trivial block filter: take all blocks
     */
    struct TakeAllBlocks
    {
      template <int row>          struct D1 { static bool const assemble = true; };
      template <int row, int col> struct D2 { static bool const assemble = true; };
    };

    // Mapping a block filter to an MPL predicate
    template <class BlockFilter>
    struct MatrixBlockFilter
    {
      template <class X>
      struct apply
      {
        typedef boost::mpl::bool_<BlockFilter::template D2<X::rowId,X::colId>::assemble> type;
      };
    };

    template <class> struct Fill;
    template <class> struct FillPointer;
    template <class,int,int> struct FillPointerWithBlock;

    template <class Block, bool, bool>
    struct CopyBlock
    {
      template <class OtherBlock>
      static void apply(OtherBlock const&, std::unique_ptr<Block>&, bool, size_t){}
    };

    template <class Scalar, int n, class Allocator, bool symmetric>
    struct CopyBlock<Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,n,n>,Allocator>,symmetric,true>
    {
      typedef Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,n,n>,Allocator> Block;

      static void apply(Block const& from, std::unique_ptr<Block>& to, bool extractOnlyLowerTriangle, size_t nnz)
      {
        // copy
        if(extractOnlyLowerTriangle == symmetric)
        {
          to.reset(new Block(from));
          return;
        }

        // only lower triangle is stored in the assembler -> restore full matrix
        if( (!extractOnlyLowerTriangle && symmetric) || (extractOnlyLowerTriangle && !symmetric) )
        {
          std::vector<std::vector<size_t> > bcrsPattern(from.N());
          for(size_t i=0; i<from.N(); ++i)
            for(size_t j=0; j<=i; ++j)
            {
              if(from.exists(i,j))
              {
                bcrsPattern[i].push_back(j);
                if(i!=j && !extractOnlyLowerTriangle)
                  bcrsPattern[j].push_back(i);
              }
            }
          for(std::vector<size_t>& row : bcrsPattern) std::sort(row.begin(),row.end());

          to.reset(new Block(from.N(),from.M(),nnz,Block::row_wise));

          // init sparsity pattern
          for (typename Block::CreateIterator row=to->createbegin(); row!=to->createend(); ++row)
            for(size_t col : bcrsPattern[row.index()]) row.insert(col);

          // read data
          for(size_t row=0; row<bcrsPattern.size(); ++row)
            for(size_t col : bcrsPattern[row])
            {
              if(row >= col) (*to)[row][col] = from[row][col];
              else if(!extractOnlyLowerTriangle)
              {
                for(size_t i=0; i<n; ++i)
                  for(size_t j=0; j<n; ++j)
                    (*to)[row][col][i][j] = from[col][row][j][i];
              }
            }
        }
      }
    };


    template <class BCRSMat, int rowIndex, int columnIndex>
    struct BlockToBCRS
    {
      explicit BlockToBCRS(std::unique_ptr<BCRSMat>& m, bool extractOnlyLowerTriangle_, size_t nnz_) : mat(m), extractOnlyLowerTriangle(extractOnlyLowerTriangle_), nnz(nnz_)
      {}

      template <class MatrixBlock>
      void operator()(MatrixBlock const& mb) const
      {
        CopyBlock<BCRSMat,MatrixBlock::symmetric,MatrixBlock::rowId==rowIndex && MatrixBlock::colId==columnIndex>::apply(mb.globalMatrix(),mat,extractOnlyLowerTriangle,nnz);
      }

    private:
      std::unique_ptr<BCRSMat>& mat;
      bool extractOnlyLowerTriangle;
      size_t nnz;
    };


    template <class Scalar, int n, int m, class Allocator, int row, int column>
    struct FillPointerWithBlock<Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,n,m>,Allocator>,row,column>
    {
      typedef Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,n,m>,Allocator> Matrix;

      template <class MatrixBlockArray>
      static std::unique_ptr<Matrix> apply(MatrixBlockArray const& block, bool extractOnlyLowerTriangle, size_t nnz)
        {
          std::unique_ptr<Matrix> result;

          for_each(block,BlockToBCRS<Matrix,row,column>(result,extractOnlyLowerTriangle,nnz));
          return result;
        }
    };
  } // End of namespace AssemblyDetail
  /// \endcond
  //---------------------------------------------------------------------
  //---------------------------------------------------------------------


  /**
   * \ingroup assembly
   * A boundary detector to be used with the
   * VariationalFunctionalAssembler. This detector checks the existence
   * of boundary intersections on construction and simply looks it up
   * when questioned.
   *
   * If multiple assembly passes are performed on the same mesh, this
   * tends to be faster than the ForwardingBoundaryDetector or even the
   * TrivialBoundaryDetector.
   */
  template <class GridView>
  class CachingBoundaryDetector
  {
    using Cell = Kaskade::Cell<GridView>;
    typedef CachingBoundaryDetector<GridView> Self;

  public:
    CachingBoundaryDetector(GridSignals& signals, GridView const& gridView_) :
      gridView(gridView_), indexSet(gridView.indexSet())
    {
      refConnection = signals.informAboutRefinement.connect(int(GridSignals::allocResources),
                                                            [=](GridSignals::Status when){this->update(when);});
      update(GridSignals::AfterRefinement);
    }

    CachingBoundaryDetector(GridView const& gridView_) :
      gridView(gridView_), indexSet(gridView.indexSet())
    {
      update(GridSignals::AfterRefinement);
    }

    bool hasBoundaryIntersections(Cell const& cell) const
    {
      assert(indexSet.size(0)==boundaryFlag.size()); // check for grid adaptation...
      return boundaryFlag[indexSet.index(cell)];
    }

  private:
    GridView const&                    gridView;
    typename GridView::IndexSet const& indexSet;
    std::vector<bool>                  boundaryFlag;
    boost::signals2::scoped_connection refConnection;

    void update(int const status)
    {
      using namespace Dune;

      if (status == GridSignals::AfterRefinement)
      {
        boundaryFlag.resize(indexSet.size(0));
        // TODO: run this loop in parallel (no so easy because of concurrent access to the bool vector...)
        for (auto const& cell: elements(gridView))
          boundaryFlag[indexSet.index(cell)] = cell.hasBoundaryIntersections();
      }
    }
  };

  /**
   * \ingroup assembly
   * A boundary detector to be used with the
   * VariationalFunctionalAssembler. This detector simply forwards the
   * questioning to the entity.
   */
  class ForwardingBoundaryDetector
  {
  public:
    template <class GridView>
    ForwardingBoundaryDetector(GridSignals const&, GridView const&)
    {}

    template <class GridView>
    ForwardingBoundaryDetector(GridView const&)
    {}

    template <class Entity>
    bool hasBoundaryIntersections(Entity const& cell) const
    {
      return cell.hasBoundaryIntersections();
    }
  };

  /**
   * \ingroup assembly
   * A boundary detector to be used with the
   * VariationalFunctionalAssembler. This detector simply says every
   * cell might have boundary intersections.
   */
  class TrivialBoundaryDetector
  {
  public:
    template <class GridView>
    TrivialBoundaryDetector(GridSignals const&, GridView const&)
    {}

    template <class GridView>
    TrivialBoundaryDetector(GridView const&)
    {}

    template <class Entity>
    bool hasBoundaryIntersections(Entity const& cell) const
    {
      return true;
    }
  };


  //---------------------------------------------------------------------

  //---------------------------------------------------------------------
  //---------------------------------------------------------------------

  namespace Assembler
  {
    /**
     * \ingroup assembly
     * \brief What to assemble.
     *
     * Binary flags to be or'ed for supplying assembly requests. An integer bitfield is supplied when calling the assemble() method of any 
     * \ref VariationalFunctionalAssembler, to tell the assembler which parts (right hand side, matrix, functional value) are to be 
     * computed. This information is passed on to the DomainCache and BoundaryCache objects of the VariationalFunctional, such that
     * they can restrict the computation to actually required values.
     */
    enum { VALUE=0x1, RHS=0x2, MATRIX=0x4, EVERYTHING=0x7 };
  }

  // ----------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------

  /**
   * \ingroup assembly
   * \brief A class for assembling Galerkin representation matrices and right hand sides for
   *        variational functionals depending on multiple variables.
   *
   * Template parameters:
   * \tparam F: The variational functional. This has to be a model of the VariationalFunctional concept.
   * \tparam BoundaryDetector: A policy class type for detecting cells
   *   incident to the boundary. Currently CachingBoundaryDetector,
   *   ForwardingBoundaryDetector, and TrivialBoundaryDetector may be used.
   * \tparam QuadRule: allows choice of quadrature formula. Default: provided by DUNE, some special formulas: special_quadrature.hh
   *
   * After a grid refinement an assembler is not valid anymore and has to be
   * constructed again. This is, because the assembler copies the index set
   * at construction, which changes, if the grid changes
   *
   */
  template <class F,
            class SparseIndex = size_t,
            class BoundaryDetector = CachingBoundaryDetector<typename F::AnsatzVars::GridView>,
            class QuadRule =  Dune::QuadratureRule<typename F::AnsatzVars::Grid::ctype, F::AnsatzVars::Grid::dimension>>
  class VariationalFunctionalAssembler
  : public AssemblyDetail::FormulationPolicy<F>
  {
  protected:
    typedef typename AssemblyDetail::FormulationPolicy<F>               Policy;

  public:
    typedef VariationalFunctionalAssembler<F,SparseIndex,BoundaryDetector,QuadRule> Self;

    /// functional
    typedef          F                            Functional;
    /// ansatz variables 
    typedef typename Functional::AnsatzVars       AnsatzVariableSetDescription;

    /// test variables 
    typedef typename Functional::TestVars         TestVariableSetDescription;

    /// grid
    typedef typename AnsatzVariableSetDescription::Grid      Grid;

    /**
     * \brief The grid view on which the variables live.
     */
    typedef typename AnsatzVariableSetDescription::GridView  GridView;

    /// spaces
    typedef typename AnsatzVariableSetDescription::Spaces    Spaces;
    /// index set
    typedef typename AnsatzVariableSetDescription::IndexSet  IndexSet;

    /**
     * \brief Underlying field type.
     */
    typedef typename Functional::Scalar           field_type;
    
    using Scalar = field_type;

  protected:
    typedef typename AnsatzVariableSetDescription::Variables AnsatzVariables;
    typedef typename TestVariableSetDescription::Variables   TestVariables;

    static_assert(std::is_same<typename AnsatzVariableSetDescription::Spaces, typename TestVariableSetDescription::Spaces>::value,
                  "VariableSetDescriptions for ansatz spaces and test spaces must provide the same space list.");
    static_assert(Functional::type==WeakFormulation || std::is_same<AnsatzVariableSetDescription,TestVariableSetDescription>::value,
                  "For problem type VariationalFunctional, ansatz and test variables must be the same.");

    // Whether we need to assemble on inner faces.
    constexpr static bool innerBoundaries = hasInnerBoundaryCache<Functional>;

  public:

    /**
     * \brief A boost::fusion sequence of AssemblyDetail::MatrixBlock elements for present matrix blocks.
     *
     * The elements of this sequence type include both block info and the global matrices.
     */
    using MatrixBlockArray = typename AssemblyDetail::BlockArray<Policy,AnsatzVariables,TestVariables,SparseIndex>::type;

    /**
     * \brief A boost::fusion sequence of AssemblyDetail::RhsBock elements for present rhs blocks.
     *
     * The elements of this sequence type include the rhs block infos.
     */
    typedef typename AssemblyDetail::RhsBlockArray<Policy,TestVariables>::type RhsBlockArray;

    /**
     * \brief A LinearProductSpace type of right hand side coefficient vectors.
     */
    typedef typename TestVariableSetDescription::template CoefficientVectorRepresentation<>::type RhsArray;

  private:
    VariationalFunctionalAssembler(GridManagerBase<Grid>& gridManager_, Spaces const& spaces)
    : spaces_(spaces),
      gridManager(gridManager_),
      gridView(boost::fusion::at_c<0>(spaces)->gridView()),
      indexSet(gridView.indexSet()),
      boundaryDetector(gridManager_.signals,gridView),
      nSimultaneousBlocks(80),
      rowBlockFactor(2.0),
      localStorageSize(128000) // roughly 128 kB
    {
      static_assert(symmetryCheck<Functional,std::min(AnsatzVariableSetDescription::noOfVariables,
                                                      TestVariableSetDescription::noOfVariables)>,
                    "If a problem is of type VariationalFunctional, the present diagonal blocks need to be marked as symmetric in D2.");
      refConnection = gridManager_.signals.informAboutRefinement.connect(int(GridSignals::freeResources),
                                                                         [this](GridSignals::Status when){this->reactToRefinement(when);});
    }
    
    
  public:
    /// Construct an empty assembler, gridManager is implicitly obtained from the first space
    explicit VariationalFunctionalAssembler(Spaces const& spaces)
    : VariationalFunctionalAssembler(boost::fusion::deref(boost::fusion::begin(spaces))->gridManager(),spaces)
    {
    }

    // injected here for backward compatibility and convenience. TODO: Remove in the long run?
    /// DEPRECATED, enums in the Assembler namespace
    static unsigned int const VALUE  = Assembler::VALUE;
    static unsigned int const RHS    = Assembler::RHS;
    static unsigned int const MATRIX = Assembler::MATRIX;
    static unsigned int const EVERYTHING = Assembler::EVERYTHING;

    /**
     * \name Tuning parameters
     * @{
     */
    
    /**
     * \brief Defines how many cells are assembled locally before scattering them together into the global data structures.
     *
     * Higher numbers reduce the overhead of locking for mutually exclusive write access to global matrix/rhs and improve
     * memory access locality, but increase the amount of thread-local computation and storage. The default value is 42, but the
     * performance appears to be not very sensitive to this parameter in the range [20,60].
     */
    Self& setNSimultaneousBlocks(int n)
    {
      nSimultaneousBlocks = n;
      return *this;
    }

    /**
     * \brief Defines how many more row blocks in each matrix are used compared to the number of threads.
     *
     * When scattering the local matrices/rhs into the global data structures, blocks of sequential rows are
     * processed in turn. Each block of rows is equipped with a lock in order to guarantee exclusive write access.
     *
     * A high number of row blocks allows fine granularity and improves cache locality and simultaneous scattering into different
     * parts of the global data structures, but incurs some inefficiency due to more frequent locking and more
     * thread-local processing. The default value is to use two times the number of threads.
     */
    Self& setRowBlockFactor(double a)
    {
      rowBlockFactor = a;
      return *this;
    }

    /**
     * \brief Defines how many memory the locally assembled matrices may occupy before they are scattered
     *
     * The amount of memory (in bytes) should scale with CPU cache, roughly a quarter of the second level cache size should be
     * a reasonable value. Note that this is merely a hint, not a hard limit. The actual local matrix storage
     * can be slightly larger.
     */
    Self& setLocalStorageSize(size_t s)
    {
      localStorageSize = s;
      return *this;
    }

    /**
     * \brief Whether to gather timing statistics using Kaskade::Timings.
     *
     * Gathering timings can be switched on and off. Note that if several assembler
     * calls are done in parallel, it *must* be switched off, since the timings
     * functionality is not thread-safe. Defaults to off.
     */
    Self& setGatherTimings(bool gt)
    {
      gatherTimings = gt;
      return *this;
    }
    
    /**
     * @}
     */

    /// Assembly without block filter or cell filter
    void assemble(F const& f, unsigned int flags=Assembler::EVERYTHING, int nThreads=0, bool verbose=false)
    {
      if (flags)
        assemble<AssemblyDetail::TakeAllBlocks>(f,[](auto const& cell) { return true; },
                                                flags,nThreads,verbose);
    }

  protected:
    // Define iterator and entity type for stepping through all cells of the grid.
    typedef typename GridView::template Codim<0>::Iterator CellIterator;
    typedef typename CellIterator::Entity Entity;


  public:

    /**
     * \brief Create data in assembler
     *
     * Assembles the value and/or derivative and/or Hessian (according
     * to the given flags, default is to assemble all three) of the
     * FE-discretized variational functional given by f. The block
     * filter can be used for partial assembly of specified blocks in
     * the matrix and rhs.
     *
     * The data that is not concerned by the given flags and block
     * filter is left unmodified. Be careful: This allows to drive the
     * Galerkin operator representation into a semantically inconsistent
     * state. This need not be a problem (instead it can be useful and more efficient in
     * certain situations), but one has to be aware of that fact.
     * The data that is not concerned by the cell filter is cleared before assembly
     * (must be, because new values are partially added up to existing data).
     *
     * \tparam BlockFilter a traits class implementing the BlockFilter concept
     * \tparam CellFilter  a callable type for deciding on which cells to assemble
     *
     * The assembled data can afterwards be accessed via the functional,
     * rhs, and matrix methods.
     *
     * \param f the variational functional or weak formulation to be assembled
     * \param cellFilter a callable mapping a grid cell to a boolean. Assembly is performed only on those cells
     *                   for which the cell filter returns true. This can be used, e.g., for subassembly in
     *                   non-overlapping domain decomposition such as BDDC, or for partial update of stiffness
     *                   matrices/rhs in case of spatially local changes in Newton's method.
     * \param flags a bitset determining what to assemble
     * \param nTasks number of tasks to submit to the thread pool for parallel execution. 0 (default) means to use the
     *                 default hardware concurrency reported by the system. Sequential computation is performed if the grid is reported
     *                 not to be thread-safe. You can call \ref GridManagerBase::enforceConcurrentReads
     *                 to lie about thread-safety in case you know (or believe) the implementation to be
     *                 thread safe.
     * \param verbose whether to output status messages to standard output
     *
     * \see NumaThreadPool
     */
    template <class BlockFilter, class CellFilter>
    void assemble(F const& f, CellFilter const& cellFilter, unsigned int flags=Assembler::EVERYTHING, int nTasks=0, bool verbose=false)
    {
      using namespace Dune;
      using namespace boost::fusion;
      using namespace AssemblyDetail;

      auto& timings = Timings::instance();
      ScopedTimingSection assemblyTiming("assembly",gatherTimings);

      if (gatherTimings) timings.start("get cell ranges");
      auto const& cellRanges = gridManager.cellRanges(gridView);
      if (gatherTimings) timings.stop("get cell ranges");

      // Choose the number of tasks
      if (nTasks < 1)
        nTasks = NumaThreadPool::instance().cpus();

      if (!gridManager.gridIsThreadSafe())
      {
        if (verbose) 
          std::cout << "Grid is not thread safe. Reducing number of used threads to 1." << std::endl;
        nTasks = 1;
      }
      else 
        if(verbose) std::cout << "#tasks:" << nTasks << " " << std::endl;
      assert(nTasks>=1);

      nTasks = std::min(nTasks,cellRanges.maxRanges());

      if (gatherTimings) timings.start("creating rhs/matrix data structures");
      // set entries to zero
      auto clearGlobalData = [](auto& x) { x = 0; };

      if (flags & RHS)
        for_each(getRhs().first.data,clearGlobalData);

      if (flags & MATRIX)
      {
        typedef filter_view<MatrixBlockArray,MatrixBlockFilter<BlockFilter>> MatrixBlockFilter;
        for_each(MatrixBlockFilter(getMatrix(f)),clearGlobalData);
      }
      if (gatherTimings) timings.stop("creating rhs/matrix data structures");

      // If there is only one task, things are somewhat simpler and more efficient...
      if (nTasks == 1)
      {
        ScopedTimingSection("perform sequential assembly",gatherTimings);
        Scalar retval =
            assembleCellRange<BlockFilter>(f,cellFilter,flags,gridView.template begin<0>(),gridView.template end<0>(),0);
        if (flags & VALUE)
          this->fValue = retval;
      }
      else
      {
        ScopedTimingSection parAssTime("perform parallel assembly",gatherTimings);
        // Partition the cells in ranges of almost equal size and do the assembly in parallel.
        std::vector<Scalar> values(nTasks,0.0);

        parallelFor([&cellRanges,this,&f,flags,&values,&cellFilter](int i, int n) {
                      auto range = cellRanges.range(n,i);
                      values[i] = this->assembleCellRange<BlockFilter>(f,cellFilter,flags,range.begin(),range.end(),i);
                    },
                    nTasks);

        // At the end, add all the functional values up
        if (flags & VALUE)
          this->fValue = std::accumulate(values.begin(),values.end(),0.0);
      }

      // announce validity of assembled data
      validparts |= flags;
    }

  private:


    /**
     * \brief Assembles stiffness matrix and/or right hand side on a given part of the cells.
     *
     * \param[in] f the variational functional (or weak formulation) to be used
     * \param[in] flags what to assemble
     * \param[in] first start of cell range
     * \param[in] last one behind end of cell range
     * \param[in] threadNo thread identifier (for debugging/timing output)
     *
     * This method is intended to be thread-safe.
     */
    template <class BlockFilter, class CellFilter>
    Scalar assembleCellRange(F const& f, CellFilter const& cellFilter,
                             unsigned int flags, CellIterator cell, CellIterator last, int threadNo)
    {
      using namespace Dune;
      using namespace boost::fusion;
      using namespace AssemblyDetail;

      // timer for gathering performance statistics
      boost::timer::cpu_timer scatterTimer, evalTimer, localTimer, setupTimer, totalTimer;
      if (statistics & localTimerStatistics)
      {
        evalTimer.stop();
        scatterTimer.stop();
        localTimer.stop();
      }
      // Obtain global data structures. Note that here is no race condition even if the
      // data is constructed on demand, as the required data structures have been created
      // and initialized before.
      MatrixBlockArray* globalMatrix    = flags & MATRIX? & getMatrix()    : nullptr;
      RhsArray*         globalRhs       = flags & RHS   ? & getRhs().first : nullptr;
      RhsBlockArray*    globalRhsBlocks = flags & RHS   ? & getRhs().second: nullptr;
      

      using MatrixBlockFilter = filter_view<MatrixBlockArray,AssemblyDetail::MatrixBlockFilter<BlockFilter>>;

      // Local contribution to global value.
      Scalar fValue = 0;

      int const dim = Grid::dimension;
      typedef typename Grid::ctype CoordType;

      // Allocate local rhs vectors that can be reused on each
      // element. Since all variables may have different element types
      // (i.e. different number of components), this has to be one local
      // rhs per variable.
      auto localRhs = as_vector(transform(*globalRhsBlocks,RhsLocalData<Policy>(nSimultaneousBlocks))); // assemble several cells en bloc before scattering

      // We assemble local matrices into a buffer from where they are scattered simultaneously. The buffers
      // scatter automatically when they are full, such that access to previous local matrices is not possible.
      // For DG methods, there are several local matrices that are assembled simultaneously in a particular 
      // access pattern, such that we use two buffers. The total buffer size is localStorageSize.
      using LocalMatrices = decltype(as_vector(transform(MatrixBlockFilter(*globalMatrix),MatrixLocalData(0))));
      LocalMatrices localMatrices, localNeighborMatrices;
      if (flags & MATRIX)
      {
        localMatrices = as_vector(transform(MatrixBlockFilter(*globalMatrix),
                                            MatrixLocalData(innerBoundaries? localStorageSize/2: localStorageSize)));
        localNeighborMatrices = as_vector(transform(MatrixBlockFilter(*globalMatrix),
                                                    MatrixLocalData(innerBoundaries? localStorageSize/2: 0)));
      }

      // element-local functional value
      Scalar floc = 0;

      // Shape function cache. Remember that every thread has to use its own cache, thus we create our own here.
      typedef ShapeFunctionCache<Grid,Scalar> SfCache;
      SfCache sfCache;

      // Quadrature rule caches. Remember that every thread has to use its own quadrature caches, thus we create our own here.
      QuadratureTraits<QuadRule>                          quadratureCache;
      QuadratureTraits<QuadratureRule<CoordType,dim-1>>   faceQuadratureCache;

      // Evaluators for shape functions of all FE spaces. Remember that
      // every thread has to use its own collection of evaluators (for both current cell and neighbouring cells in case
      // we assemble on inner boundaries)!
      auto evaluators = getEvaluators(spaces(),&sfCache,AssemblyDetail::derivatives(f,spaces(),localMatrices,localRhs));
      using Evaluators = decltype(evaluators);
      auto neighbourEvaluators = evaluators;
      // Create the "caches" responsible for evaluating the variational functional / weak formulation.
      // The inner boundary cache need not exist, thus we evaluate it as 0 if not present. Lambda function
      // is a trick to allow using constexpr if.
      auto domainCache = f.createDomainCache(flags);
      auto boundaryCache = f.createBoundaryCache(flags);
      auto innerBoundaryCache = [&](){ if constexpr (innerBoundaries) return f.createInnerBoundaryCache(flags);
                                       else                           return 0; }();
      // The symmetrizer functor. We pass the heterogeneous array of global matrices
      // such that run-time meta information, e.g., on enforcing positive definiteness,
      // can be used.
      SymmetrizeLocalMatrix<Functional> symmetrizeLocalMatrix{&f};

      // Dummy variable sets to be fed to boost::fusion algorithms.
      typename AnsatzVariableSetDescription::Variables ansatzVars;
      typename TestVariableSetDescription::Variables testVars;

      if (statistics & localTimerStatistics)
        setupTimer.stop();

      // Assemble representation: step through all elements, assemble
      // the local matrix and rhs, and scatter their entries into the
      // global data structures.
      if (statistics & localTimerStatistics)
        evalTimer.resume();

      while(cell!=last)
      {
        // If the cell is not of interest, skip
        if (!cellFilter(*cell))
        {
          ++cell;
          continue;
        }

        // tell the problem on which cell we are
        domainCache.moveTo(*cell);

        // for all spaces involved, compute the shape functions and
        // their global indices, which are needed for evaluating the
        // functional's derivative.
        moveEvaluatorsToCell(evaluators,*cell,indexSet.index(*cell));
        int const shapeFunctionMaxOrder = maxOrder(evaluators);

        if (statistics & localTimerStatistics) localTimer.resume();
        if (statistics & localTimerStatistics) evalTimer.stop();

        // create new local matrices and right hand sides.
        floc = 0;
        if (flags & RHS)
          for_each(localRhs,ClearLocalRhs<Evaluators,TestVariables>(evaluators));

        if (flags & MATRIX)
          for_each(localMatrices,NewLocalMatrix<Evaluators,AnsatzVariableSetDescription,
                                                TestVariableSetDescription>(evaluators,evaluators));
          
        // loop over all quadrature points and integrate local stiffness matrix and rhs
        int const integrationOrderOnCell = f.integrationOrder(*cell,shapeFunctionMaxOrder,false);
        GeometryType gt = cell->type();

        QuadRule const qr = quadratureCache.rule(gt,integrationOrderOnCell);
        useQuadratureRuleInEvaluators(evaluators,qr,0);

        size_t nQuadPos = qr.size();
        for (size_t g=0; g<nQuadPos; ++g)
        {
          // pos of integration point
          FieldVector<CoordType,dim> const& quadPos = qr[g].position();

          // for all spaces involved, update the evaluators associated to this quadrature point
          moveEvaluatorsToIntegrationPoint(evaluators,quadPos,qr,g,0);

          // prepare evaluation of functional
          domainCache.evaluateAt(quadPos,evaluators);

          CoordType weightTimesDetJac(cell->geometry().integrationElement(quadPos)); // determinant of jacobian
          assert(weightTimesDetJac > 0);
          // weight of quadrature point
          weightTimesDetJac *= qr[g].weight();

          // step through all local matrix and rhs blocks and update their data
          if (this->hasValue && (flags & VALUE))
            floc += weightTimesDetJac*Policy::d0(domainCache);
          if (flags & RHS)
            for_each(localRhs,UpdateLocalRhs<TestVariableSetDescription,Evaluators,Scalar,
                                             typename F::DomainCache>(evaluators,weightTimesDetJac,domainCache));
          if (flags & MATRIX)
            for_each(localMatrices,UpdateLocalMatrix<AnsatzVariableSetDescription,TestVariableSetDescription,
                                                     Evaluators,Scalar,typename F::DomainCache>(
                                                  evaluators,weightTimesDetJac,domainCache));
        }


        // Handle Robin (Cauchy) type boundary conditions and interior face terms. Loop over all
        // faces (codim 1) that are exterior boundaries or are relevant due to interior face terms and integrate
        // the contribution. We assume that intersections on the domain boundary always consist of the whole face (codim 1
        // subentity).
        if(innerBoundaries || boundaryDetector.hasBoundaryIntersections(*cell))
        {
          int const integrationOrderOnFace = f.integrationOrder(*cell,shapeFunctionMaxOrder,true);

          // This would be nice, but since there is no way to obtain an IntersectionIterator from an
          // Intersection, we'd need to pass the Intersection directly to boundaryCache.moveTo()
          // - which would be nicer, but is a breaking change. Consider this for the next major revision.
          //
          // for (auto const& face: intersections(gridView,*cell))
          //
          // Instead, we use the old style.
          auto faceEnd = gridView.iend(*cell);
          for (auto face=gridView.ibegin(*cell); face!=faceEnd; ++face)
          {
            // @warning: make SURE the following does always hold: (i)
            // On boundary faces, there is exactly one intersection, and
            // this covers the whole face. (ii) When an intersection
            // covers the whole codim 1 subentity, its geometry in the
            // cell is the same as specified in the reference
            // element. Currently this is not guaranteed by the Dune
            // interface. Maybe there is a relation to the ominous
            // TwistUtility from Freiburg?

            bool const onBoundary = face->boundary();

            // Assemble only if needed.
            if (innerBoundaries || onBoundary)
            {
              GeometryType gt = face->type();
              QuadratureRule<CoordType,dim-1> const qr = faceQuadratureCache.rule(gt,integrationOrderOnFace);
              size_t const nQuadPos = qr.size();

              // Move the caches to the current face (and neighboring evaluators to the adjacent cell).
              if (onBoundary)
                boundaryCache.moveTo(face);
              else if constexpr (innerBoundaries) // the constexpr if just protects the otherwise ill-formed body
              {
                assert(face->neighbor());
                assert(gt == face->geometryInOutside().type());

                // If the problem claims the face is not relevant, skip assembly
                if (! f.considerFace(*face))
                  continue;

                innerBoundaryCache.moveTo(face);
                moveEvaluatorsToCell(neighbourEvaluators,face->outside());

                // for each neighbor (i.e. face), there is one local matrix that is computed by quadrature. Let's
                // create a new local matrix if needed. The row (test variable) is given by the current cell, and the
                // column (ansatz variable) determined by the neighbor.
                if (flags & MATRIX)
                  for_each(localNeighborMatrices,
                           NewLocalMatrix<Evaluators,AnsatzVariableSetDescription,
                                                     TestVariableSetDescription>(evaluators,neighbourEvaluators));
              }


              for (size_t g=0; g<nQuadPos; ++g)
              {
                // Evaluation by quadrature point index cannot be used on faces due to inconsistent numbering
                // of quadrature points relative to the cell (see below). Hence we don't need to announce the quadrature
                // rule to the evaluators here.
                // useQuadratureRuleInEvaluators(evaluators,qr,face->indexInInside());

                // position of integration point
                Dune::FieldVector<CoordType,dim-1> quadPos = qr[g].position();
                CoordType weightTimesDetJac = qr[g].weight() * face->geometry().integrationElement(quadPos);

                // evaluate values at integration points for all spaces
                // involved, update the evaluators associated to this
                // quadrature point
                //
                // TODO: check whether this call is
                // thread-safe. (currently as of 2010-02-11 this is not
                // guaranteed by the Dune interface for UG). Copying the
                // quadrature position instead of obtaining a const
                // reference minimizes the time window for concurrent
                // access.
                // Note that global() here returns the *cell-local* position.
                FieldVector<CoordType,dim> const quadPosInCell = face->geometryInInside().global(quadPos);

                // Move evaluators to integration point
                // quadrature points on faces are not consistently indexed (if even symmetric...) - hence evaluating
                // shape functions based on quadrature point index is impossible
                // moveEvaluatorsToIntegrationPoint(evaluators,quadPosInCell,qr,g,face->indexInInside());
                moveEvaluatorsToIntegrationPoint(evaluators,quadPosInCell);

                if (onBoundary)     // assemble on domain boundary
                {
                  boundaryCache.evaluateAt(quadPos,evaluators);

                  // step through all rhs blocks and update their local rhs
                  if (this->hasValue && (flags & VALUE))
                    floc += weightTimesDetJac*Policy::d0(boundaryCache);
                  if (flags & RHS)
                    for_each(localRhs,UpdateLocalRhs<TestVariableSetDescription,Evaluators,Scalar,
                      typename F::BoundaryCache>(evaluators,weightTimesDetJac,boundaryCache));
                  if (flags & MATRIX)
                    for_each(localMatrices,UpdateLocalMatrix<AnsatzVariableSetDescription,TestVariableSetDescription,
                      Evaluators,Scalar,typename F::BoundaryCache>(evaluators,weightTimesDetJac,boundaryCache));

                }  // done boundary faces
                else if constexpr (innerBoundaries)   // assemble on inner boundary ("constexpr if"
                {                                     // only because the body may be ill-formed)
                  // Move in addition the neighboring evaluators to the correct evaluation point
                  moveEvaluatorsToIntegrationPoint(neighbourEvaluators,face->geometryInOutside().global(quadPos));

                  innerBoundaryCache.evaluateAt(quadPos,evaluators,neighbourEvaluators);

                  // Update value. Remember that each face is visited twice, so we add here only our half.
                  if (this->hasValue && (flags & VALUE))
                    floc += 0.5 * weightTimesDetJac*FormulationPolicy<Functional>::d0(innerBoundaryCache);

                  using InnerBoundaryCache = typename Functional::InnerBoundaryCache;

                  // Step through all rhs blocks and update their local rhs. Remember that even if the face
                  // is visited twice, these visits correspond to different rows (the row being associated
                  // to the current center cell).
                  // Hence we do not sum up the samve values twice and have to omit the factor 0.5.
                  if (flags & RHS)
                    boost::fusion::for_each(localRhs,UpdateLocalRhs<TestVariableSetDescription,Evaluators,
                                            Scalar,InnerBoundaryCache>
                                            (evaluators,weightTimesDetJac,innerBoundaryCache));

                  // Update local matrices. As with the rhs, we have to omit the factor 0.5.
                  if (flags & MATRIX)
                  {
                    // update diagonal block, contribution of coupling terms on cell (boolean true)
                    boost::fusion::for_each(localMatrices,UpdateLocalMatrixFromInnerBoundaryCache<AnsatzVariableSetDescription,
                                            TestVariableSetDescription,Evaluators,Scalar,InnerBoundaryCache>
                                            (evaluators,evaluators,weightTimesDetJac,innerBoundaryCache,true));
                    // update off-diagonal block, contribution of coupling terms with neighbouring cell (boolean false)
                    boost::fusion::for_each(localNeighborMatrices,UpdateLocalMatrixFromInnerBoundaryCache<AnsatzVariableSetDescription,
                                            TestVariableSetDescription,Evaluators,Scalar,InnerBoundaryCache>
                                            (evaluators,neighbourEvaluators,weightTimesDetJac,innerBoundaryCache,false));
                  }
                }
              }
            }
          } // done face loop
        } // done boundary conditions

        // symmetrize local matrix (only the center cell matrices, not the off-diagonal neighbor matrices which are 
        // not symmetric even for symmetric problems).
        if (flags & MATRIX)
          for_each(localMatrices,symmetrizeLocalMatrix);

        if (statistics & localTimerStatistics) scatterTimer.resume();
        if (statistics & localTimerStatistics) localTimer.stop();

        // Move to next cell.
        ++cell;
        // Scatter local data into global data (note that this occurs blockwise, i.e.
        // the actual scatter may be delayed). Local matrix buffers scatter automatically.
        fValue += floc;
        if (flags & RHS)
          for_each(localRhs,ScatterLocalRhs<RhsArray>(*globalRhs));

        if (statistics & localTimerStatistics) evalTimer.resume(); // count cell iterator increment as well
        if (statistics & localTimerStatistics) scatterTimer.stop();
      } // end iteration over cells

      if (statistics & localTimerStatistics) scatterTimer.resume();
      if (statistics & localTimerStatistics) evalTimer.stop();

      // scattering of local vectors may not be complete (there may be local data left over in the buffers). Again, local matrices
      // scatter automatically on destruction.
      if (flags & RHS)
        for_each(localRhs,ScatterLocalRhs<RhsArray>(*globalRhs,true));

      if (statistics & localTimerStatistics) scatterTimer.stop();

      // report statistics
      if (statistics & localTimerStatistics)
      {
        std::cerr << "--------------\nThread " << threadNo << ":\n";
        double totalTime = setupTimer.elapsed().wall + evalTimer.elapsed().wall + scatterTimer.elapsed().wall + localTimer.elapsed().wall;
        std::cerr.precision(1);
        std::cerr.setf(std::ios_base::fixed,std::ios_base::floatfield);
        std::cerr << "Total time:    " << 100* totalTimer.elapsed().wall / totalTime << "% -- " << totalTimer.format()
        << "Setup time:    " << 100* setupTimer.elapsed().wall / totalTime << "% -- " << setupTimer.format()
        << "Eval times:    " << 100* evalTimer.elapsed().wall / totalTime << "% -- " << evalTimer.format()
        << "Scatter times: " << 100* scatterTimer.elapsed().wall / totalTime << "% -- " << scatterTimer.format()
        << "Local times:   " << 100 * localTimer.elapsed().wall / totalTime << "% -- " << localTimer.format();
        std::cerr << "--------------\n";
      }

      return fValue;
    }

  public:
    /**
     * \brief Returns a bitfield specifyign which of the parts have been assembled since
     * construction or flush() according to the format (VALUE|RHS|MATRX).
     */
    int valid() { return validparts; }

    /**
     * \brief Destructs parts of the assembled quatities according to the format (VALUE|RHS|MATRX)
     *
     * The default flag leads to destruction of all three values.
     */
    void flush(int flags=(VALUE|RHS|MATRIX))
    {
      validparts &= !flags;
      if (flags & VALUE)  this->fValue = 0;
      if (flags & RHS)
      {
        rhss.reset();
        rhsBlocks.reset();
      }
      if (flags & MATRIX) matrixBlocks.reset();
    }

    /**
     * \brief the size of a matrix block (in terms of scalar rows/columns)
     * \return a pair (nrows,ncols)
     */
    std::pair<size_t,size_t> size(int row, int col) const
    {
      return std::make_pair(TestVariableSetDescription::degreesOfFreedom(spaces(),row,row+1),
                            AnsatzVariableSetDescription::degreesOfFreedom(spaces(),col,col+1));
    }

    /**
     * \name Submatrix access
     * @{
     */

    /**
     * \brief Writes the submatrix defined by the half-open block ranges
     * [rbegin,rend), [cbegin,cend) to a newly created std::unique_ptr<MatrixType>.
     *
     * Use this method only if your matrix type does not support move-semantics.
     *
     * If extractOnlyLowerTriangle is true, the Galerkin operator is treated as
     * if its upper triangle was structurally zero.
     *
     * In order to be able to read the assembled matrix into the smart pointer
     * a specialization of AssemblyDetail::FillPointer, i.e. AssemblyDetail::FillPointer<MatrixType>
     * must be available!!!
     *
     * Kaskade provides spezializations of AssemblyDetail::Fill for
     *  - MatrixType = MatrixAsTriplet<Scalar>
     *  - MatrixType = Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,1,1> >
     *
     */
    template <class MatrixType, int rbegin=0, int rend=TestVariableSetDescription::noOfVariables,
              int cbegin=0, int cend=AnsatzVariableSetDescription::noOfVariables>
    std::unique_ptr<MatrixType> getPointer(bool extractOnlyLowerTriangle) const __attribute__((deprecated)) /*definition at end of file*/;

    /**
     * \brief Extracts the submatrix defined by the half-open block ranges
     * [rbegin,rend), [cbegin,cend).
     *
     * If extractOnlyLowerTriangle is true, the Galerkin operator is treated as
     * if its upper triangle was structurally zero.
     *
     * In order to be able to read the assembled matrix into the smart pointer
     * a specialization of AssemblyDetail::Fill, i.e. AssemblyDetail::Fill<MatrixType>
     * must be available!!!
     *
     * Kaskade provides spezializations of AssemblyDetail::Fill for
     *  - MatrixType = MatrixAsTriplet<Scalar>
     *  - MatrixType = Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,1,1> >
     *
     * \param extractOnlyLowerTriangle if true, and the matrix is symmetric, only the lower triangular part is extracted
     * \param rbegin number first block row to store
     * \param rend   number of one past the last block row to store
     * \param cbegin number first block column to store
     * \param cend   number of one past the last block column to store
     */
    template <class Matrix>
    Matrix get(bool extractOnlyLowerTriangle=false, int rbegin=0, int rend=TestVariableSetDescription::noOfVariables,
                                                    int cbegin=0, int cend=AnsatzVariableSetDescription::noOfVariables) const
    {
      // block sizes
      auto cSizes = AnsatzVariableSetDescription::variableDimensions(spaces());
      auto rSizes = TestVariableSetDescription::variableDimensions(spaces());

      // row and column offsets
      std::vector<size_t> rowOff(rend-rbegin,0), colOff(cend-cbegin,0);
      std::partial_sum(rSizes.begin()+rbegin,rSizes.begin()+rend-1,rowOff.begin()+1);
      std::partial_sum(cSizes.begin()+cbegin,cSizes.begin()+cend-1,colOff.begin()+1);

      return AssemblyDetail::Fill<Matrix>::apply(getMatrix(), rbegin, cbegin, rowOff, colOff, extractOnlyLowerTriangle,
                                                 nnz(rbegin, rend, cbegin, cend, extractOnlyLowerTriangle), nrows(rbegin,rend), ncols(cbegin,cend));
    }

    /**
     * \brief Extracts a raw single submatrix block indexed by row and column. 
     * 
     * Note that for symmetric matrices, only the lower triangular part is stored.
     * 
     * Compilation will fail if the referenced matrix block is not present.
     */
    template <int row, int col>
    auto const& get() const
    {
      return (*boost::fusion::find_if<AssemblyDetail::IsBlock<row,col>>(getMatrix())).globalMatrix();
    }

    /**
     * \brief Extracts the submatrix defined by the half-open block ranges given as template parameter.
     *
     * Provide the block information via IstlInterface_Detail::BlockInfo.
     */
    template <class MatrixType, class BlockInformation>
    MatrixType get(bool extractOnlyLowerTriangle) const
    {
      using BI = BlockInformation;
      return get<MatrixType,BI::firstRow,BI::lastRow,BI::firstCol,BI::lastCol>(extractOnlyLowerTriangle);
    }

    /**
     * @}
     */

    /**
     * \brief Writes the subvector defined by the half-open block range [rbegin,rend) to the output iterator xi.
     */
    template <class DataOutIter>
    void toSequence(int rbegin, int rend, DataOutIter xi) const
    {
      // WARNING: this assumes that for_each processes the rhsArray in
      // correct order from front to back!
      for_each(getRhs().first.data,BlockToSequence<DataOutIter>(rbegin,rend,xi));
    }


    /**
     * \brief Returns the number of structurally nonzero entries in the
     * submatrix defined by the half-open block ranges [rbegin,rend),
     * [cbegin,cend).
     *
     * This is the number of elements written by the
     * corresponding call of toTriplet().
     */
    size_t nnz(size_t rbegin=0, size_t rend=TestVariableSetDescription::noOfVariables,
               size_t cbegin=0, size_t cend=AnsatzVariableSetDescription::noOfVariables, bool extractOnlyLowerTriangle=false) const
    {
      size_t n = 0;
      boost::fusion::for_each(getMatrix(),CountNonzeros(n,rbegin,rend,cbegin,cend,extractOnlyLowerTriangle));
      return n;
    }

    /**
     * \brief Returns the number of scalar rows in the row block range [firstBlock,lastBlock[.
     */
    size_t nrows(int firstBlock=0, int lastBlock=TestVariableSetDescription::noOfVariables) const
    {
      return TestVariableSetDescription::degreesOfFreedom(spaces(),firstBlock,lastBlock);
    }

    /**
     * \brief Returns the number of scalar cols in the column block range [firstBlock,lastBlock[.
     */
    size_t ncols(int firstBlock=0, int lastBlock=AnsatzVariableSetDescription::noOfVariables) const
    {
      return AnsatzVariableSetDescription::degreesOfFreedom(spaces(),firstBlock,lastBlock);
    }

    /**
     * \brief Returns the list of spaces used.
     */
    Spaces const& spaces() const { return spaces_; }


    /**
     * \brief Returns a mutable reference to the sequence of matrix blocks.
     *
     *
     *
     * The data structure for the matrix is created on the fly if needed.
     *
     * \todo check const correctness - why is this method const and returns a mutable matrix?
     */
    MatrixBlockArray& getMatrix() const 
    {
      if (!matrixBlocks)
        createMatrix([](auto const& face) { return true; });
      return *matrixBlocks;
    }

  private:
    // Returns the matrix, creating it on the fly if necessary.
    MatrixBlockArray& getMatrix(F const& f)
    {
      if (!matrixBlocks)
        createMatrix([&](auto const& face)
        {
          if constexpr (innerBoundaries)
            return f.considerFace(face);
          else
            return true;
        });

      return *matrixBlocks;
    }


    template <class FaceOracle>
    void createMatrix(FaceOracle const& considerFace) const
    {
      matrixBlocks.reset(new MatrixBlockArray());
      boost::fusion::for_each(*matrixBlocks,[&](auto& block)
      {
        block.init(this->spaces(),gridView.template begin<0>(),gridView.template end<0>(),
                   considerFace);
      });
    }

  public:
    /**
     * \brief Returns a contiguous subrange of the rhs coefficient vectors.
     *
     * \return an element of a linear product space with components of coefficient vectors.
     */
    template <int first=0, int last=TestVariableSetDescription::noOfVariables>
    typename TestVariableSetDescription::template CoefficientVectorRepresentation<first,last>::type rhs() const
    {
      using Rhs = typename TestVariableSetDescription::template CoefficientVectorRepresentation<first,last>::type;
      return Rhs(make_range<first,last>(getRhs().first.data));
    }

  protected:
    std::pair<RhsArray&,RhsBlockArray&> getRhs() const
    {
      if(rhss.get()==0)
      {
        // construct global vectors
// old version, should be equivalent to the following line.
//        rhss.reset(new RhsArray(Variables_Detail::VariableRangeCreator<TestVariables,ConstructCoefficientVector<Spaces>>
//            ::apply(ConstructCoefficientVector<Spaces>(spaces()))));
        rhss.reset(new RhsArray(TestVariableSetDescription::template CoefficientVectorRepresentation<>::init(spaces())));
        // construct block infos
        rhsBlocks.reset(new RhsBlockArray());
        // Choose the number of row groups in the matrix blocks. The obvious choice would be 1 in sequential mode, but it
        // appears that a somewhat larger number is beneficial, probably due to better memory access locality. For multithreaded
        // mode, we choose the number of hardware threads as the default, since then all threads can (in principle) be active
        // simultaneously.
        int nrg = std::max((int)(rowBlockFactor*boost::thread::hardware_concurrency()),4);
        boost::fusion::for_each(*rhsBlocks,[&](auto& block)
        {
          block.init(this->spaces(),gridView.template begin<0>(),gridView.template end<0>(),nrg);
        });
      }

      return std::pair<RhsArray&,RhsBlockArray&>(*rhss,*rhsBlocks);
    }


    // Resource management: delete matrix and rhs data structures as well as cell ranges on
    // refinement since they are no longer meaningful on the new grid.
    void reactToRefinement(GridSignals::Status const status)
    {
      if (status == GridSignals::BeforeRefinement)
        flush();
    }


    // resource management
    boost::signals2::scoped_connection refConnection;

    Spaces const&                      spaces_;
    GridManagerBase<Grid> const&       gridManager;
    GridView const&                    gridView;
    IndexSet const&                    indexSet;

    mutable std::shared_ptr<MatrixBlockArray> matrixBlocks; // matrix block info including global matrices
    mutable std::shared_ptr<RhsArray>         rhss;         // global rhs vectors
    mutable std::shared_ptr<RhsBlockArray>    rhsBlocks;    // rhs block info
    BoundaryDetector                          boundaryDetector;
    int                                       validparts;


    int nSimultaneousBlocks;
    double rowBlockFactor;
    size_t localStorageSize;
    bool gatherTimings = false;                 // if true, use default Kaskade timings


    template <class IdxOutIter, class DataOutIter>
    struct BlockToTriplet
    {
      BlockToTriplet(size_t rbegin_, size_t cbegin_, std::vector<size_t> const& rowOff_, std::vector<size_t> const& colOff_,
          IdxOutIter& ri_, IdxOutIter& ci_, DataOutIter& xi_, bool extractOnlyLowerTriangle_):
            rbegin(rbegin_), cbegin(cbegin_), rowOff(rowOff_), colOff(colOff_),
            ri(ri_), ci(ci_), xi(xi_), extractOnlyLowerTriangle(extractOnlyLowerTriangle_)
      {}

      template <class MatrixBlock>
      void operator()(MatrixBlock const& mb) const
      {
        // Check if block is in requested range
        if (inRange(mb.rowId,mb.colId))
          Matrix_to_Triplet<typename MatrixBlock::Matrix>::call(mb.globalMatrix(),ri,ci,xi,
              rowOff[MatrixBlock::rowId-rbegin],colOff[MatrixBlock::colId-cbegin],
              mb.rowId==mb.colId && extractOnlyLowerTriangle,
              mb.symmetric);
        // For mirror blocks, the transposed block needs to be written
        // if the transposed block is in the requested
        // range. Transposition is implicitly achieved by swapping
        // column and row index output iterators.
        if (MatrixBlock::mirror && inRange(mb.colId,mb.rowId) && extractOnlyLowerTriangle==false)
          Matrix_to_Triplet<typename MatrixBlock::Matrix>::call(mb.globalMatrix(),ci,ri,xi,
              colOff[MatrixBlock::rowId-cbegin],rowOff[MatrixBlock::colId-rbegin],
              false,
              mb.symmetric);
      }

      bool inRange(size_t r, size_t c) const
      {
        return r>=rbegin && r<rbegin+rowOff.size() && c>=cbegin && c<cbegin+colOff.size();
      }


      size_t rbegin, cbegin;
      std::vector<size_t> const& rowOff;
      std::vector<size_t> const& colOff;
      IdxOutIter& ri;
      IdxOutIter& ci;
      DataOutIter& xi;
      bool extractOnlyLowerTriangle;
    };

    struct CountNonzeros
    {
      CountNonzeros(size_t& n_, size_t rbegin_, size_t rend_, size_t cbegin_, size_t cend_, bool onlyLowerTriangle_):
        n(n_), rbegin(rbegin_), rend(rend_), cbegin(cbegin_), cend(cend_), onlyLowerTriangle(onlyLowerTriangle_)
      {}

      template <class MatrixBlock>
      void operator()(MatrixBlock const& mb) const
      {
        size_t myN = Matrix_to_Triplet<typename MatrixBlock::Matrix>::nnz(mb.globalMatrix(),
                                                                          (mb.rowId==mb.colId) && onlyLowerTriangle,
                                                                          mb.symmetric);
        // Check if block is in requested range
        if (inRange(MatrixBlock::rowId,MatrixBlock::colId))
          n += myN;

        // For subdiagonal blocks, the transposed block needs to be
        // counted if onlyLowerTriangle is false and the transposed
        // block is in the requested range and we need to mirror this
        // block.
        if (MatrixBlock::mirror && inRange(mb.colId,mb.rowId) && onlyLowerTriangle==false)
          n += myN;
      }

      bool inRange(size_t r, size_t c) const
      {
        return r>=rbegin && r<rend && c>=cbegin && c<cend;
      }



      size_t& n;
      size_t rbegin, rend, cbegin, cend;
      bool onlyLowerTriangle;
    };

    template <class DataOutIter>
    struct BlockToSequence
    {
      BlockToSequence(int& rbegin_, int& rend_, DataOutIter& out_): rbegin(rbegin_), rend(rend_), out(out_) {}
      template <class VectorBlock> void operator()(VectorBlock const& v) const
      {
        if (rbegin<=0 && rend>0)
          out = vectorToSequence(v,out);
        --rbegin;
        --rend;
      }
    private:
      int& rbegin;
      int& rend;
      DataOutIter& out;
    };

  };

  // ----------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------

  template <class F, class SparseIndex, class BoundaryDetector,class QuadRule>
  template <class MatrixType, int rbegin, int rend, int cbegin, int cend>
  std::unique_ptr<MatrixType> VariationalFunctionalAssembler<F,SparseIndex,BoundaryDetector,QuadRule>::getPointer(bool extractOnlyLowerTriangle) const
  {
    using namespace boost::fusion;

    // block sizes
    auto cSizes = AnsatzVariableSetDescription::variableDimensions(spaces());
    auto rSizes = TestVariableSetDescription::variableDimensions(spaces());

    // row and column offsets
    std::vector<size_t> rowOff(rend-rbegin,0), colOff(cend-cbegin,0);
    std::partial_sum(rSizes.begin()+rbegin,rSizes.begin()+rend-1,rowOff.begin()+1);
    std::partial_sum(cSizes.begin()+cbegin,cSizes.begin()+cend-1,colOff.begin()+1);

    return AssemblyDetail::FillPointer<MatrixType>::apply(getMatrix(), rbegin, cbegin, rowOff, colOff, extractOnlyLowerTriangle,
        nnz(rbegin, rend, cbegin, cend, extractOnlyLowerTriangle), nrows(rbegin,rend), ncols(cbegin,cend));
  }
  // End of VariationalFunctionalAssembler
  //---------------------------------------------------------------------
  //---------------------------------------------------------------------
} /* end of namespace Kaskade */

#endif
