/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2012-2022 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef THREADED_MATRIX
#define THREADED_MATRIX

#include <iostream>
#include <boost/timer/timer.hpp>


#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/operators.hh>

#include "fem/firstless.hh"
#include "fem/fixdune.hh"
#include "linalg/dynamicMatrix.hh"
#include "utilities/duneInterface.hh"
#include "utilities/threading.hh"
#include "utilities/detailed_exception.hh"



namespace Kaskade {
  /**
   * \internal
   */
  // forward declarations
  template <class Entry, class Index>
  class NumaBCRSMatrix;
  
  template <class Scalar, class Index>
  class MatrixAsTriplet;

  template <class Target, class Source, class RowIndices, class ColIndices>
  Target submatrix(Source const& A, RowIndices const& ri, ColIndices const& ci);

  
    
  namespace ThreadedMatrixDetail 
  {
    // Computes for each row the number of entries in that row. The count is stored in the 
    // vector rowCount, which has to be big enough.
    template <class Entry, class Allocator, class Index>
    void getRowCount(Dune::BCRSMatrix<Entry,Allocator> const& matrix, bool isTransposed, std::vector<Index>& rowCount)
    {
      for (auto ri=matrix.begin(); ri!=matrix.end(); ++ri)
        if (isTransposed)
          for (auto ci=ri->begin(); ci!=ri->end(); ++ci)
            ++rowCount[ci.index()];
        else
          rowCount[ri.index()] = ri->size();
    }
    
    template <class Entry, class Index2, class Index>
    void getRowCount(NumaBCRSMatrix<Entry,Index2> const& matrix, bool isTransposed, std::vector<Index>& rowCount)
    {
      for (auto ri=matrix.begin(); ri!=matrix.end(); ++ri)
        if (isTransposed)
          for (auto ci=ri->begin(); ci!=ri->end(); ++ci)
            ++rowCount[ci.index()];
        else
          rowCount[ri.index()] = ri->size();
    }
    
    template <class Entry, class Index, class Index2>
    void getRowCount(MatrixAsTriplet<Entry,Index2> const& matrix, bool isTransposed, std::vector<Index>& rowCount)
    {
      auto const& indices = isTransposed? matrix.cidx: matrix.ridx;
      for (auto r: indices)
        ++rowCount[r];
    }
    
    
    // Class with partial specialization for copying from a Dune ISTL sparse matrix into ThreadedMatrix Chunks. 
    // Supported modes of copying are (i) one-to-one direct copy (ii) transposed copy (iii) source is symmetric
    // (superdiagonal part is never accessed).
    template <class Entry, class Matrix, bool symmetric, bool transposed, class Index>
    struct CopyMatrixToChunk 
    {
      // The simple case: not symmetric, not transposed
      static void init(Index first, Index last, Matrix const& matrix, std::vector<size_t,NumaAllocator<size_t>>& colStart, 
                       std::vector<Index,NumaAllocator<Index>>& cols, std::vector<Entry,NumaAllocator<Entry>>& values,
                       std::vector<Index> const& /* nRowEntries */)
      {
        for (Index i=0; i<last-first; ++i)
        {
          Index j=0;
          for (auto ci=matrix[first+i].begin(); ci!=matrix[first+i].end(); ++ci, ++j)
          {
            values[colStart[i]+j] = *ci; 
            cols[colStart[i]+j] = ci.index();
          }
        }
      }
    };
    
    template <class Entry, class Matrix, bool transposed, class Index>
    struct CopyMatrixToChunk<Entry, Matrix,true,transposed,Index>
    {
      // The hard case: symmetric with only lower triangular part to be accessed
      static void init(Index first, Index last, Matrix const& matrix, std::vector<size_t,NumaAllocator<size_t>>& colStart, 
		       std::vector<Index,NumaAllocator<Index>>& cols, std::vector<Entry,NumaAllocator<Entry>>& values,
		       std::vector<Index> const& nRowEntries) 
      {
        // Symmetric case is a mixture of regular and transposed. First we insert the regular entries,
        // subsequently the transposed ones.
        
        // Start with regular entries (including the diagonal)
        std::vector<Index> entriesInRow(last-first,0);
        for (Index i=0; i<last-first; ++i)
        {
          Index j=0;
          for (auto ci=matrix[first+i].begin(); ci!=matrix[first+i].end(); ++ci, ++j)
          {
            values[colStart[i]+j] = *ci; 
            cols[colStart[i]+j] = ci.index();
          }
          entriesInRow[i] = matrix[first+i].size();
          assert(entriesInRow[i]<=nRowEntries[first+i]);
        }
        
        // Do the transposed entries (excluding the diagonal). Note that as we only want entries with
        // column index >= first+1, we can start (for the transposed part) at row first+1
        for (Index r=first+1; r!=matrix.N(); ++r)
          for (auto c=matrix[r].begin(); c!=matrix[r].end(); ++c)
          {
            // get transposed indices 
            Index ridx = c.index();
            if (ridx>=first && ridx<last)
            {
              Index cidx = r;
              
              if (cidx > ridx) // take only superdiagonal entries 
                {
                  assert(ridx-first<colStart.size() && ridx-first < entriesInRow.size());
                  Index pos = colStart[ridx-first] + entriesInRow[ridx-first];
                  assert(pos<colStart[ridx-first+1]);
                  cols[pos] = cidx;
                  values[pos] = transpose(*c);
                  ++entriesInRow[ridx-first];
                  assert(entriesInRow[ridx-first]<=nRowEntries[ridx]);
                }
            }
          }
      }
    };
    
    
    template <class Entry, class Matrix, class Index>
    struct CopyMatrixToChunk<Entry, Matrix,false,true,Index>
    {
      // The case of transposed, nonsymmetric supplied matrices.
      static void init(Index first, Index last, Matrix const& matrix, std::vector<size_t,NumaAllocator<size_t>>& colStart, 
		       std::vector<Index,NumaAllocator<Index>>& cols, std::vector<Entry,NumaAllocator<Entry>>& values,
		       std::vector<Index> const& nRowEntries) 
      {
        // Essentially we have to exchange row and column indices. If we go through the matrix
        // row by row, the entries of the transposed matrix to be stored in our CRS chunk come
        // with random row indices but increasing column indices. Hence we can push the entries
        // at the back of each row. We only have to remember where the current end of each row is.
        std::vector<Index> entriesInRow(last-first,0);
        
        for (auto r=matrix.begin(); r!=matrix.end(); ++r)
          for (auto c=r->begin(); c!=r->end(); ++c)
          {
            // get transposed indices
            Index ridx = c.index();
            if (ridx>=first && ridx<last)
            {
              Index cidx = r.index();
              
              assert(ridx-first<colStart.size() && ridx-first < entriesInRow.size());
              Index pos = colStart[ridx-first] + entriesInRow[ridx-first];
              assert(pos<colStart[ridx-first+1]);
              cols[pos] = cidx;
              values[pos] = transpose(*c);
              ++entriesInRow[ridx-first];
              assert(entriesInRow[ridx-first]<=nRowEntries[ridx]);
            }
          }
      }
    };
    
    //-------------------------------------------------------------------------

    // Class for copying matrix entries. If the entry is to be transposed during
    // copying (e.g. because only the lower triangular part of the source block matrix 
    // has been stored), the matrix dimensions may change due to transposition 
    // (if they are not equal). This case is covered here - we know the entry is
    // to be transposed.
    template <class To, class From, class scalarsMatch = std::false_type>
    struct MatrixEntry
    {
      static void copy(From const& from, To& to, bool isTransposed)
      {
        assert(isTransposed);
        to = transpose(from);
      }
    };

    // This specialization is for the case of different floating point types as scalar entries in the matrix blocks.
    // Currently this is only implemented for blocks of quadratic shape and non transposed blocks.
    template <class To, class From>
    struct MatrixEntry<To, From, typename std::is_same<typename To::field_type, typename From::field_type>::type>
    {
      static void copy(From const& from, To& to, bool isTransposed)
      {
        assert(!isTransposed);
        static_assert(To::rows == From::rows && To::cols == From::cols, "Mismatch of dimensions.");
        for(int i=0;i<To::rows;++i)
          for(int j=0;j<To::cols;++j)
            to[i][j] = static_cast<typename To::field_type>(from[i][j]);
      }
    };
    
    // The case of matching dimensions of source and target is covered here.
    // Both cases can occur, transposition or not. Partial specialization of
    // this handler class is necessary due to the runtime switch between
    // transposed or not - this dynamic switch is only legal in C++ if the 
    // dimensions are the same in both cases.
    template <class Entry>
    struct MatrixEntry<Entry,Entry>
    {
      static void copy(Entry const& from, Entry& to, bool isTransposed)
      {
        if (isTransposed && from.rows > 1 && from.cols > 1 && from.rows != from.cols)
        {
          std::cerr << "Transposing rectangular Dune::FieldMatrix not implemented!" << std::endl
                    << "#rows=" << from.rows << ", #cols=" << from.cols << std::endl;
          exit(1);
        }
         if (isTransposed && from.rows == from.cols)
//           to = transpose(from);
//       the following transposing code also covers the case that 
//       the memory locations of from and to are the same
         {
            for (int i=0;i<from.rows;i++)
            {
              to[i][i]=from[i][i];
              for (int j=i+1;j<from.cols;j++)
              {
                typename Entry::field_type s=from[j][i];
                to[j][i]=from[i][j];
                to[i][j]=s;
              }
            }
         }
         else
           to = from;
      }
    };
    
    //-------------------------------------------------------------------------

    /**
     * \ingroup linalgbasic
     * \brief A base class representing basic meta information about sparsity patterns of NUMA matrix chunks.
     */
    template <class Index=size_t>
    class CRSChunkPatternInfo
    {
    public:
      /**
       * \brief Constructor.
       * \param firstRow index of first row 
       * \param lastRow index one behind the las row (half-open STL style range)
       * \param cols the number of columns
       * \param symmetric if true, store only lower triangular part
       * \param node NUMA node on which to hold the memory
       */
      CRSChunkPatternInfo(Index first_, Index last_, Index cols_, bool symmetric, int node)
      : firstRow(first_), lastRow(last_), cols(cols_), numaNode(node), symm(symmetric)
      {}
      
      /**
       * \brief start of the covered row range
       */
      Index first() const { return firstRow; }
      
      /**
       * \brief end of the half-open row range
       */
      Index last() const { return lastRow; }
      
      /**
       * \brief number of columns in the matrix
       */
      Index columns() const { return cols; }
      
      /**
       * \brief the NUMA node on which to allocate the memory
       */
      int node() const { return numaNode; }
      
      /**
       * \brief if true, only the lower triangular part is stored
       */
      bool symmetric() const { return symm; } 

    protected:
      Index firstRow, lastRow; // first and one behind last row
      Index cols;              // number of columns
      
    private:      
      int numaNode;
      bool symm;
    };
    
    //-------------------------------------------------------------------------

    /**
     * \ingroup linalgbasic
     * \brief A class supporting two-stage construction of sparsity patterns of NUMA matrix chunks.
     * \tparam Index integral type for representing row/column indices (usually int, long, or size_t).
     * 
     * This class allows to insert nonzero entry positions efficiently one by one. Lateron, a CRSChunkPattern
     * can be created (also efficiently) from this pattern creator.
     * 
     * Design rationale: Maintaining and updating a CRS sparsity pattern data structure during
     * insertion of elements is quite expensive. The alternative would be a stateful interface for
     * two-stage construction as in Dune::BCRSMatrix (but stateful interfaces are bad).
     */
    template <class Index=size_t>
    class CRSChunkPatternCreator: public CRSChunkPatternInfo<Index>
    {
    public:
      /**
       * \brief An STL container holding a sequence of indices.
       */
      typedef std::vector<Index> IndexArray; // column indices are indexed by user-provided Index type          
      // We could use a Numa allocator here, but testing as of 2014-03-12 revealed that the performance is 
      // *much* worse with Numa allocators, even though the Numa allocator is actually faster in
      // synthetic tests of a structure similar to here. This reason for this disparity is up to now unknown.
      
      
      /**
       * \brief Constructor.
       * \param firstRow index of first row 
       * \param lastRow index one behind the las row (half-open STL style range)
       * \param ncols number of columns
       * \param symmetric if true, store only lower triangular part
       * \param node NUMA node on which to hold the memory
       */
      CRSChunkPatternCreator(Index firstRow, Index lastRow, Index ncols, bool symmetric, int node)
      : CRSChunkPatternInfo<Index>(firstRow,lastRow,ncols,symmetric,node), 
        cols(lastRow-firstRow)
      {
      }
      
      /**
       * \brief Reserves a certain amount of entries per row.
       * \param nnzPerRow hint for the expected number of nonzeros per row
       * 
       * Use this to prevent frequent reallocations.
       */
      void reserve(Index nnzPerRow=8)
      {
        // perform allocation
        for (auto& r: cols)
          r.reserve(nnzPerRow);
      }
      
      /**
       * \brief Clears the complete data, handing back the memory.
       */
      void clear()
      {
        for (auto& r: cols)
        {
          r.clear();
          r.shrink_to_fit();
        }
      }
      
      /**
       * \brief Enters elements into the sparsity pattern.
       * \tparam IterRow forward iterator for a range of row indices
       * \tparam IterCol random access iterator for a range of column indices
       * \param fromRow start of row indices
       * \param toRow one behind of last row index
       * \param fromCol start of sorted column indices
       * \param toCol one behind of last column index
       * \param colIsSorted true if the provided range of column indices is sorted ascendingly
       * 
       * All entries \f$ (i,j) \f$ with \f$ i \f$ in the given row indices and \f$ j \f$ in the 
       * given column indices are entered into the sparsity structure. 
       * 
       * For performance reasons, no index in either range should occur twice.
       */
      template <class IterRow, class IterCol>
      void addElements(IterRow fromRow, IterRow const toRow,
                       IterCol const fromCol, IterCol const toCol, bool colIsSorted=false)
      {
        // sort column indices if needed, and remove duplicates (set_union below does NOT remove duplicates...)
        sortedCols.clear();
        sortedCols.insert(end(sortedCols),fromCol,toCol);
        if (!colIsSorted)
          std::sort(begin(sortedCols),end(sortedCols));
        sortedCols.erase(std::unique(begin(sortedCols),end(sortedCols)),end(sortedCols));
        
        // For each row, merge the specified column indices into the column indices that are already present.
        for ( ; fromRow != toRow; ++fromRow)                         // step through all affected rows
          if (this->first() <= *fromRow && *fromRow < this->last())  // but treat only those that actually lie in our chunk
          {
            IndexArray& c = cols[*fromRow-this->first()];
            auto top = this->symmetric()?
                              std::upper_bound(begin(sortedCols),end(sortedCols),*fromRow):
                              end(sortedCols); // omit superdiagonal elements if symmetric
            tmp.resize(c.size() + std::distance(begin(sortedCols),top));         // reserve space for the union
            tmp.erase(std::set_union(begin(c),end(c),begin(sortedCols),top,      // compute the union
                                     begin(tmp)),end(tmp));                      // and cut the container at the proper end
            std::swap(tmp,c);                                                    // write the sorted range into its target
          }
        // Note that instead of a temporary buffer tmp we could append the new column indices and use 
        // std::inplace_merge. This, however, uses a temporary buffer of its own behind the scenes, one that is held
        // by the STL library implementation somewhere. It is highly probable that this memory is *not* located
        // on our NUMA node, but frequently accessed from all nodes. Hence, false sharing is bound to occur with 
        // inplace_merge even if it performs proper locking.
      }
      
      /**
       * \brief Enters all possible elements (defining a dense chunk).
       */
      void addAllElements(Index columns) 
      {
        for (Index i=0; i<cols.size(); ++i)
        {
          cols[i].resize(this->symmetric()? i+this->first()+1: columns);
          std::iota(begin(cols[i]),end(cols[i]),static_cast<Index>(0));
        }
      }
      
      /**
       * \brief Returns the sorted and unique column indices of elements in row \f$ i \f$.
       * \param i row index in the range [first,last[
       */
      IndexArray const& row(Index i) const { return cols[i-this->first()]; }
      
      /**
       * \brief Returns the number of stored entries (structurally nonzero elements).
       */
      size_t nonzeroes() const 
      {  
        size_t nnz = 0;
        for (auto const& r: cols)
          nnz += r.size();
        return nnz;
      }
      
      /**
       * \brief Moves rows in and out of the chunk in order to equilibrate the number of nonzeroes.
       * \param[in] covered      the number of nonzeroes contained in previous (lower node) chunks without the moved rows
       * \param[in] nnz          the total number of nonzeroes
       * \param[in] chunks       the number of chunks
       * \param[in,out] moveRows on entry, contains the rows to prepend. On exit, contains the rows removed from the end 
       * \return                 the number of nonzeroes contained in this and previous chunks without the moved rows
       * 
       * Explicit instantiations for int and size_t are defined in are defined in threadedMatrix.cpp.
       */
      size_t balanceForward(size_t const covered, size_t const nnz, int chunks, std::vector<IndexArray>& moveRows);

      /**
       * \brief Moves rows in and out of the chunk in order to equilibrate the number of nonzeroes.
       * \param[in] covered      the number of nonzeroes contained in this and previous (lower node) chunks including the ones in moveRows
       * \param[in] nnz          the total number of nonzeroes
       * \param[in] chunks       the number of chunks
       * \param[in,out] moveRows on entry, contains the rows to append. On exit, contains the rows removed from the beginning
       * \return                 the number of nonzeroes contained in this and previous chunks
       * 
       * Explicit instantiations for int and size_t are defined in are defined in threadedMatrix.cpp.
       */
      size_t balanceBackward(size_t covered, size_t const nnz, int chunks, std::vector<IndexArray>& moveRows);

    private:
      std::vector<IndexArray> cols;      // For each row a std::vector of column indices  
      IndexArray tmp;                    // some temporary buffer
      IndexArray sortedCols;             // some temporary buffer for sorted column indices
    };
    
    //-------------------------------------------------------------------------
    
    /**
     * \ingroup linalgbasic
     * \brief This class maintains the sparsity structure of a couple of matrix rows (a NUMA matrix chunk).
     * \tparam Index integral type for representing row/column indices (usually int, long, or size_t).
     */
    template <class Index=size_t>
    class CRSChunkPattern: public CRSChunkPatternInfo<Index>
    {
    public: 
      
      /**
       * \brief Constructor
       * 
       * Constructs the pattern from a given pattern creator.
       */
      CRSChunkPattern(CRSChunkPatternCreator<Index> const& creator);
      
      /**
       * \brief Constructor.
       * \tparam Expanded an array type with value type convertible to Index
       * \tparam Condensed an array type with value type convertible to Index
       * 
       * 
       * \param first    the first row/col of the condensed indices to include
       * \param last     one behind the last row/col of the condensed indices to include in this chunk
       * \param eIndices sorted global array of expanded indices
       * \param cIndices sorted global array of condensed indices
       * \param mat      a matrix with BCRSMatrix interface
       * 
       * The expanded indices are the ones in the matrix mat, the condensed refer to the matrix to
       * be constructed. It must hold that cIndices[eIndices[i]] == i for all firstRow <= i < lastRow.
       */
      template <class Expanded, class Condensed, class Matrix>
      CRSChunkPattern(Index first, Index last, Expanded const& eIndices, Condensed const& cIndices, 
	              Matrix const& mat, int node)
      : CRSChunkPatternInfo<Index>(first,last,eIndices.size(),false,node)
      , colStarts(last-first+1,0,NumaAllocator<size_t>(node)), cols(NumaAllocator<Index>(node))
      {
        cols.reserve(last-first);                  // prevent frequent reallocation (heuristic)

        for (Index i=first; i<last; ++i)           // scan all our rows
        {
          colStarts[i-first] = cols.size();        // note where this row starts

          auto const& row = mat[eIndices[i]];      // get a handle to the column indices in that row

          auto cend = row.end();
          for (auto c = row.begin(); c!=cend; ++c) // step through all column entries in that row
          {
            Index ci = cIndices[c.index()];         // extract the condensed index
            if (ci<eIndices.size())                 // if that is in the condensed matrix range...
              cols.push_back(ci);                   // ...include it
          }
        }

        colStarts.back() = cols.size();            // sentinel

        // investigate sparsity
        nnzPerRow = this->last()==this->first()? 0 : nonzeroes() / (this->last()-this->first());
      }
      
      /**
       * \brief Constructor
       * \tparam Matrix the supplied Dune::ISTL matrix
       * 
       * \param matrix the matrix to be copied
       * \param symmetric whether the supplied matrix is symmetric and only its lower triangular part is stored
       * \param transposed whether the supplied matrix is transposed
       * \param firstRow first row index in our chunk
       * \param lastRow one behind last row index in our chunk
       * \param symmetric whether to store only the lower triangular part
       * \param node the NUMA node on which to store the data
       */
      template <class Matrix>
      CRSChunkPattern(Matrix const& matrix, bool isSymmetric, bool isTransposed, 
                      Index firstRow, Index lastRow, bool symmetric, int node)
      : CRSChunkPatternInfo<Index>(firstRow,lastRow,matrix.M(),symmetric,node)
      , colStarts(NumaAllocator<size_t>(node)), cols(NumaAllocator<Index>(node))
      {
        if (isSymmetric)         // Make sure that transposed is only flagged if
          isTransposed = false;  // it's not symmetric.
        
        // Compute start of row indices. Step through the supplied matrix and count all elements that 
        // fall in one of our rows.
        std::vector<size_t,NumaAllocator<size_t>> rowCount(this->last()-this->first(),0,NumaAllocator<size_t>(node));
        
        // Compute the row range of the supplied matrix we have to scan.
        Index fromRow = this->first(), toRow = this->last();
        if (isSymmetric && !symmetric) // If the supplied matrix is symmetrically stored but we are not, we have to
          toRow = matrix.N();          // scan the later rows for the elements that appear in later columns here.
        if (isTransposed)                   // If the supplied matrix is transposed, we have to scan all the rows
          fromRow = 0, toRow = matrix.N();  // in order to cover all the columns in our range.
        
        for (Index ri=fromRow; ri<toRow; ++ri)
        {
          auto cend = matrix[ri].end();
          for (auto ci=matrix[ri].begin(); ci!=cend; ++ci)
          {
            Index r = ri, c = ci.index();
            if (isTransposed) // supplied matrix is transposed, swap row and column indices
              std::swap(r,c); // to get the row/col indices in our world
            if (this->first()<=r && r<this->last() && (!symmetric || c<=r)) // yep, this element is in our chunk
              ++rowCount[r-this->first()];                                  // 
            if (isSymmetric && !symmetric && c<r && this->first()<=c && c<this->last()) // if supplied matrix is stored symmetrically but we are not
              ++rowCount[c-this->first()];                                              // copy truely subdiagonal entries to the upper triangle
          }
        }
        
        // compute the partial sums of number of elements per row.
        colStarts.resize(rowCount.size()+1,0);
        std::partial_sum(rowCount.begin(),rowCount.end(),colStarts.begin()+1);
        
        // Get space for all nonzeros
        cols.resize(colStarts.back());
        
        // Copy the column indices of elements.
        for (Index ri=fromRow; ri<toRow; ++ri)
        {
          auto cend = matrix[ri].end();
          for (auto ci=matrix[ri].begin(); ci!=cend; ++ci)
          {
            Index r = ri, c = ci.index();
            if (isTransposed) // supplied matrix is transposed, swap row and column indices
              std::swap(r,c); // to get the row/col indices in our world
            if (this->first()<=r && r<this->last() && (!symmetric || c<=r))    // yep, this element is in our chunk
            {
              cols[colStarts[r-this->first()]+rowCount[r-this->first()]-1] = c; // enter the column index
              --rowCount[r-this->first()];                                      // one less element to come
            }
            if (isSymmetric && !symmetric && c<r && this->first()<=c && c<this->last()) // if supplied matrix is stored symmetrically but we are not
            {
              cols[colStarts[c-this->first()]+rowCount[c-this->first()]-1] = r;         // copy truely subdiagonal entries to the upper triangle
              --rowCount[c-this->first()];                                              // one less element to come
            }
          }
        }
        
        // sort the column indices in each row 
        for (Index ri=0; ri<rowCount.size(); ++ri) {
          assert(rowCount[ri] == 0);
          std::sort(cols.begin()+colStarts[ri],cols.begin()+colStarts[ri+1]);
        }
                
        // investigate sparsity
        nnzPerRow = this->last()==this->first()? 0 : nonzeroes() / (this->last()-this->first());
      }
      
      /**
       * \brief Constructor
       * \param matrix the matrix to be copied (in triplet format)
       * \param symmetric whether the supplied matrix is symmetric and only its lower triangular part is stored
       * \param transposed whether the supplied matrix is transposed
       * \param firstRow first row index in our chunk
       * \param lastRow one behind last row index in our chunk
       * \param symmetric whether to store only the lower triangular part
       * \param node the NUMA node on which to store the data
       */
      template <class Scalar, class Index2>
      CRSChunkPattern(MatrixAsTriplet<Scalar,Index2> const& matrix, bool isSymmetric, bool isTransposed,  // TODO: fuse with constructor above
                      Index firstRow, Index lastRow, bool symmetric, int node)
      : CRSChunkPatternInfo<Index>(firstRow,lastRow,matrix.M(),symmetric,node)
      , colStarts(NumaAllocator<size_t>(node)), cols(NumaAllocator<Index>(node))
      {
        if (isSymmetric)         // Make sure that transposed is only flagged if
          isTransposed = false;  // it's not symmetric.
        
        // Compute start of row indices. Step through the supplied matrix and count all elements that 
        // fall in one of our rows.
        std::vector<size_t,NumaAllocator<size_t>> rowCount(this->last()-this->first(),0,NumaAllocator<size_t>(node));
        
        for (size_t i=0; i<matrix.ridx.size(); ++i)
        {
          Index r = matrix.ridx[i], c = matrix.cidx[i];
          if (isTransposed) // supplied matrix is transposed, swap row and column indices
            std::swap(r,c); // to get the row/col indices in our world
          if (this->first()<=r && r<this->last() && (!symmetric || c<=r)) // yep, this element is in our chunk
            ++rowCount[r-this->first()];                                  // 
          if (isSymmetric && !symmetric && c<r && this->first()<=c && c<this->last()) // if supplied matrix is stored symmetrically but we are not
            ++rowCount[c-this->first()];                                              // copy truely subdiagonal entries to the upper triangle
        }
        
        // compute the partial sums of number of elements per row.
        colStarts.resize(rowCount.size()+1,0);
        std::partial_sum(rowCount.begin(),rowCount.end(),colStarts.begin()+1);
        
        // Get space for all nonzeros
        cols.resize(colStarts.back());
        
        // Copy the column indices of elements.
        for (size_t i=0; i<matrix.ridx.size(); ++i)
        {
          Index r = matrix.ridx[i], c = matrix.cidx[i];
          if (isTransposed) // supplied matrix is transposed, swap row and column indices
            std::swap(r,c); // to get the row/col indices in our world
          if (this->first()<=r && r<this->last() && (!symmetric || c<=r))    // yep, this element is in our chunk
          {
            cols[colStarts[r-this->first()]+rowCount[r-this->first()]-1] = c; // enter the column index
            --rowCount[r-this->first()];                                      // one less element to come
          }
          if (isSymmetric && !symmetric && c<r && this->first()<=c && c<this->last()) // if supplied matrix is stored symmetrically but we are not
          {
            cols[colStarts[c-this->first()]+rowCount[c-this->first()]-1] = r;         // copy truely subdiagonal entries to the upper triangle
            --rowCount[c-this->first()];                                              // one less element to come
          }
        }
        
        // sort the column indices in each row 
        for (Index ri=0; ri<rowCount.size(); ++ri)    // TODO: do this in parallel on our node?
        {
          assert(rowCount[ri] == 0);
          auto first = cols.begin()+colStarts[ri];
          auto last = cols.begin()+colStarts[ri+1];
          std::sort(first,last);                      // sort in ascending order
          last = std::unique(first,last);             // duplicate entries may come from triplet matrices, remove
          rowCount[ri] = last-first;                  // now we have the correct number of entries.
        }
        
        // Merging duplicate column indices in a row may have reduced the number of entries in that row,
        // hence there may be a gap to the start of the next row. Compactify the column index array now.
        for (Index ri=1; ri<rowCount.size(); ++ri)
        {
          Index newStart = colStarts[ri-1] + rowCount[ri-1];
          if (newStart != colStarts[ri])                // entries need to be shifted
          {
            std::copy(cols.begin()+colStarts[ri],cols.begin()+colStarts[ri]+rowCount[ri],cols.begin()+newStart);
            colStarts[ri] = newStart;
          }
        }
        if (!rowCount.empty())                                  // set the sentinel to point just behind the last entry...
          colStarts.back() = colStarts[rowCount.size()-1]+rowCount.back();
        cols.erase(cols.begin()+colStarts.back(),cols.end());   // ...and drop the rubbish left over in the tail
        
                
        // investigate sparsity
        nnzPerRow = this->last()==this->first()? 0 : nonzeroes() / (this->last()-this->first());
      }
      
      // returns the column index of entry at position idx
      Index col(size_t idx) const { return cols[idx]; }
      
      /**
       * \brief returns the index from which on the entries in given local row are stored 
       * \param row row index in the range [0, last-first[
       */
      size_t colStart(Index row) const { return colStarts[row]; }
      
      /**
       * \brief an iterator pointing to the start of the column index for the given local row 
       * \param row row index in the range [0, last-first[
       */
      typename std::vector<Index,NumaAllocator<Index>>::const_iterator colStartIterator(Index row) const 
      { 
        return cols.begin()+colStart(row); 
      }
      
      /**
       * \brief returns the position of an element with given global row and column indices
       * \param r row index in the range [first,last[
       * \param c column index in the range 
       * 
       * \return the position of the entry within the column and value arrays of the chunk pattern if the entry exists,
       *         std::numeric_limits<size_t>::max() otherwise.
       */
      size_t position(Index r, Index c) const
      {
        auto b = cols.begin();
        r -= this->first();
        assert(0 <= r && r+1 < colStarts.size());
        auto p = std::lower_bound(b+colStarts[r],b+colStarts[r+1],c); // col indices are sorted in each row
        if (p==b+colStarts[r+1] || *p!=c)                             // return sentinel in case the entry is not here
          return std::numeric_limits<size_t>::max();
        return p-b;                                                   // return offset
      }

      /**
       * \brief returns true if the entry with global row and column indices is present in the sparsity pattern
       */
      bool exists(Index r, Index c) const
      {
        return position(r,c) != std::numeric_limits<size_t>::max();
      }
        
      /**
       * \brief Returns the number of stored entries.
       * 
       * This differs from nonzeroes() for symmetric storage.
       */
      size_t storage() const { return cols.size(); }
      
      /**
       * \brief Returns the number of structurally nonzero entries.
       * 
       * In case of symmetric storage, subdiagonal entries are counted twice.
       */
      size_t nonzeroes() const;
      
      /**
       * \brief Returns the average number of nonzeroes per row.
       */
      Index nonzeroesPerRow() const { return nnzPerRow; }
      
    private:
      
      // Raw data storage. The allocators used are NUMA allocators, guaranteeing local memory access.
      std::vector<size_t,NumaAllocator<size_t>> colStarts; // memory positions are indexed by size_t (since there may be MANY)
      std::vector<Index,NumaAllocator<Index>> cols;        // column indices are indexed by user-provided Index type      
      Index nnzPerRow;                                     // average number of nonzeroes per row (rounded)

    };
    
    
    //-------------------------------------------------------------------------

    /**
     * \brief An iterator stepping through all entries in a row.
     * \todo Implement a proper iterator, e.g., using boost::iterator_adaptor
     */
    template <class Entry, class Index>
    class NumaBCRSMatrixConstRowIterator
    {
    public:
      typedef typename std::vector<Index,NumaAllocator<Index>>::const_iterator ColIterator;
      typedef typename std::vector<Entry,NumaAllocator<Entry>>::const_iterator ValueIterator;
      using Self = NumaBCRSMatrixConstRowIterator<Entry,Index>;
      
      // typedefs required by std::iterator_traits
      using value_type = Entry;
      using difference_type = std::ptrdiff_t;
      using pointer = Entry const*;
      using reference = Entry const&;
      using iterator_category = std::random_access_iterator_tag;
      
      NumaBCRSMatrixConstRowIterator(ColIterator col_, ValueIterator val_): col(col_), val(val_) {}
      
      Index index() const { return *col; }
      
      /**
       * \name Access
       * \{
       */
      Entry const& operator*() const { return *val; }
      Entry const* operator->() const { return &*val; }
      /**
       * \}
       */
      
      /**
       * \name Iterator motion
       * \{
       */
      void operator++() { ++col; ++val; }
      void operator--() { --col; --val; }
      void operator+=(Index i) { col += i; val += i; }
      void operator-=(Index i) { col -= i; val -= i; }
      /**
       * \}
       */
      
      bool operator==(Self const& it) const { return val==it.val; }
      bool operator!=(Self const& it) const { return ! (*this==it); }
      difference_type operator-(Self const& it) const { return val-it.val; }
      
    protected:
      ColIterator col;
      ValueIterator val;
    };

    template <class Entry, class Index>
    class NumaBCRSMatrixRowIterator: public NumaBCRSMatrixConstRowIterator<Entry,Index>
    {
    public:
      typedef typename std::vector<Index,NumaAllocator<Index>>::const_iterator ColIterator;
      typedef typename std::vector<Entry,NumaAllocator<Entry>>::iterator ValueIterator;
      
      NumaBCRSMatrixRowIterator(ColIterator col_, ValueIterator val_): NumaBCRSMatrixConstRowIterator<Entry,Index>(col_,val_) {}
      
      /**
       * \name Mutable access
       * \{
       */
      Entry& operator*() const { return const_cast<Entry&>(*this->val); }
      Entry* operator->() const { return &(**this); }
      /**
       * \}
       */
    };
    
    
 
    //-------------------------------------------------------------------------
    
    template <class Arguments, class Operation>
    class NumaBCRSMatrixExpressionChunk
    {
    public:
      class iterator
      {
      };
      
      iterator begin(size_t) const
      {
      }
      
      iterator end(size_t) const
      {
      }
    };
    
    template <class Arguments, class Operation>
    class NumaBCRSMatrixExpression
    {
    public:
      NumaBCRSMatrixExpressionChunk<Arguments,Operation> const& operator[](int i) const
      {
      }
    };
    
    //-------------------------------------------------------------------------

    /**
     * \ingroup linalgbasic
     * \brief This class stores a couple of compressed row storage rows in memory allocated locally on a NUMA node.
     * \tparam Entry type of matrix elements (usually Dune::FieldMatrix<...>)
     * \tparam Index type of row/column indices (usually int, long, or size_t)
     * 
     * The sparsity structure (in form of a CRSChunkPattern) is hold using std::shared_ptr and
     * hence can be shared between multiple CRSChunk objects. The chunk itself owns just the actual 
     * matrix entries, not their indices.
     */
    template <class Entry, class Index=size_t>
    class CRSChunk
    {
      typedef CRSChunk<Entry,Index> Self;
      
    public:
      using Scalar = typename EntryTraits<Entry>::field_type;
//       typedef typename GetScalar<Entry>::type Scalar;
      
      /**
       * \brief Constructor initializing an empty Chunk. 
       * 
       * Note that this chunk is then pretty useless, and in fact no methods other than (move) assignment 
       * are allowed. Do only use if you need to (almost) default construct a chunk, and fill it directly 
       * afterwards (e.g., in parallel).
       */
      CRSChunk(int node)
      : values(NumaAllocator<Entry>(node))
      {
      }
      
      /**
       * \brief Constructor.
       */
      CRSChunk(std::shared_ptr<CRSChunkPattern<Index>> const& pattern_, Entry const& init)
      : pat(pattern_), scatterMutex(4), values(pat->nonzeroes(),init,NumaAllocator<Entry>(pat->node()))
      {
      }
      
      /**
       * \brief Constructor.
       * \tparam Matrix the type of the supplied matrix to copy
       * 
       * \param pattern the sparsity pattern of our row chunk
       * \param matrix the matrix to be copied
       * \param isSymmetric whether the supplied matrix is symmetric
       * \param isTransposed whether the supplied matrix is transposed
       * 
       * The sparsity pattern to use can be provided here, such that sharing of patterns between multiple
       * chunk objects is possible.
       */
      template <class Matrix>
      CRSChunk(std::shared_ptr<CRSChunkPattern<Index>> const& pattern_, Matrix const& matrix, bool isSymmetric, bool isTransposed)
      : CRSChunk(pattern_,Entry(0))
      {
        typedef typename Matrix::block_type SuppliedEntry;
        
        Index first = pat->first(), last = pat->last();
        
        if (isSymmetric)         // Make sure that transposed is only flagged if
          isTransposed = false;  // it's not symmetric.

        // Compute the row range of the supplied matrix we have to scan.
        Index fromRow = first, toRow = last;
        if (isSymmetric && !pat->symmetric()) // If the supplied matrix is symmetrically stored but we are not, we have to
          toRow = matrix.N();                 // scan the later rows for the elements that appear in later columns here.
        if (isTransposed)                     // If the supplied matrix is transposed, we have to scan all the rows
          fromRow = 0, toRow = matrix.N();    // in order to cover all the columns in our range.

        // Copy the column indices of elements.
        for (Index ri=fromRow; ri<toRow; ++ri)
        {
          auto cend = matrix[ri].end();
          for (auto ci=matrix[ri].begin(); ci!=cend; ++ci)
          {
            Index r = ri, c = ci.index();
            if (isTransposed) // supplied matrix is transposed, swap row and column indices
              std::swap(r,c); // to get the row/col indices in our world
            if (first<=r && r<last && (!pat->symmetric() || c<=r))                                   // yep, this element is in our chunk
            {
              assert(pat->position(r,c) < values.size());
              MatrixEntry<Entry,SuppliedEntry>::copy(*ci,values[pat->position(r,c)],isTransposed);   // copy, transpose if needed
            }
            if (isSymmetric && !pat->symmetric() && c<r && first<=c && c<last)                       // if supplied matrix is stored symmetrically but we are not..
              MatrixEntry<Entry,SuppliedEntry>::copy(*ci,values[pat->position(c,r)],true);           // copy, transpose 
          }
        }
      }
      
      /**
       * \brief Constructor.
       * 
       * \param pattern the sparsity pattern of our row chunk
       * \param matrix the matrix to be copied
       * \param isSymmetric whether the supplied matrix is symmetric
       * \param isTransposed whether the supplied matrix is transposed
       * 
       * The sparsity pattern to use can be provided here, such that sharing of patterns between multiple
       * chunk objects is possible.
       */
      template <class Scalar, class Index2>
      CRSChunk(std::shared_ptr<CRSChunkPattern<Index>> const& pattern_, MatrixAsTriplet<Scalar,Index2> const& matrix, bool isSymmetric, bool isTransposed)
      : CRSChunk(pattern_,Entry(0))
      {
        Index first = pat->first(), last = pat->last();
        
        if (isSymmetric)         // Make sure that transposed is only flagged if
          isTransposed = false;  // it's not symmetric.


        // Copy the column indices of elements.
        for (size_t i=0; i<matrix.ridx.size(); ++i)
        {
          Index r = matrix.ridx[i], c = matrix.cidx[i];
          if (isTransposed) // supplied matrix is transposed, swap row and column indices
            std::swap(r,c); // to get the row/col indices in our world
          if (first<=r && r<last && (!pat->symmetric() || c<=r))                                   // yep, this element is in our chunk
          {
            assert(pat->position(r,c) < values.size());
            values[pat->position(r,c)] = matrix.data[i];                                           // copy (no transpose as triplet entries are scalar)
          }
          if (isSymmetric && !pat->symmetric() && c<r && first<=c && c<last)                       // if supplied matrix is stored symmetrically but we are not..
            values[pat->position(c,r)] = matrix.data[i];                                           // copy (no transpose as triplet entries are scalar)
        }
      }
      
      /**
       * \brief Constructor.
       * \tparam Matrix the type of the supplied matrix to copy
       * 
       * \param matrix the matrix to be copied
       * \param isSymmetric whether the supplied matrix is symmetric
       * \param isTransposed whether the supplied matrix is transposed
       * \param firstRow
       * \param lastRow
       * \param node
       * 
       * This is a convenience constructor that creates a new sparsity pattern from scratch.
       */
      template <class Matrix>
      CRSChunk(Matrix const& matrix, bool isSymmetric, bool isTransposed, Index firstRow, Index lastRow, bool symmetric, int node)
      : CRSChunk(std::make_shared<CRSChunkPattern<Index>>(matrix,isSymmetric,isTransposed,firstRow,lastRow,symmetric,node),
                 matrix,isSymmetric,isTransposed)
      {}
      
      /**
       * \brief Constructor.
       * \tparam Matrix the type of the supplied matrix to copy
       * 
       * \param pattern the sparsity pattern of our row chunk
       * \param matrix the matrix to be copied
       * \param isSymmetric whether the supplied matrix is symmetric
       * \param isTransposed whether the supplied matrix is transposed
       * 
       * The sparsity pattern to use can be provided here, such that sharing of patterns between multiple
       * chunk objects is possible.
       */
      template <class Expanded, class Condensed,  class Matrix>
      CRSChunk(std::shared_ptr<CRSChunkPattern<Index>> const& pattern_, 
	             Expanded const& eIndices, Condensed const& cIndices, Matrix const& matrix)
      : CRSChunk(pattern_,Entry(0))
      {
        Index first = pat->first(), last = pat->last();
        for (Index r=first; r<last; ++r)
        {
          auto row = matrix[eIndices[r]];
          auto p = begin(values) + pat->colStart(r-first);
          for (auto c=row.begin(); c!=row.end(); ++c)
            if (cIndices[c.index()]<eIndices.size())
            {
              *p = *c;
              assert(!isnan(*p));
              ++p;
            }

          assert(p==begin(values)+pat->colStart(r-first+1));
        }
      }
      
      /**
       * \brief Copy assignment.
       */
      Self& operator=(Self const& c) = default;

      /**
       * \brief Assigns the given value to each entry.
       */
      template <class Value>
      Self&  operator=(Value const& a)
      {
        std::fill(begin(values),end(values),a);
        return *this;
      }
      
      
      /**
       * \brief Assigns the given matrix expression componentwise.
       */
      template <class Arguments, class Operation>
      Self& operator=(ThreadedMatrixDetail::NumaBCRSMatrixExpressionChunk<Arguments,Operation> const& e)
      {
        Index first = pat->first();
        for (Index r=0; r<pat->last()-first; ++r)
        {
          auto eend = e.end(r);
          auto ci = pat->colStartIterator(r);
          auto cend = pat->colStartIterator(r+1);
          auto vi = begin(values)+pat->colStart(r);
          for (auto ei = e.begin(r); ei != eend; ++ei) // step through the expression row
          {
            auto pos = std::find(ci,cend,ei.index());  // find our entry
            assert(pos != cend);                       // our pattern must be a superset of the expression pattern
            vi += pos-ci;                              // advance to found entry
            ci = pos;
            *vi = *ei;                                 // assign
          }
        }
      }
      
      /**
       * \brief Matrix-vector multiplication \f$ y \leftarrow y + a Ax \f$ or \f$ y \leftarrow aAx \f$ and \f$ x^T y \f$.
       * \param initialize if nonzero, the result y is set to zero before updating, realizing \f$ y \leftarrow aAx \f$.
       */
      template <class Domain, class Range>
      Scalar apply(Scalar a, Domain const& x, Range& y, bool initialize) const
      {
        Scalar dp = 0;
        auto const& p = *pat;
        Index first = p.first();
        Index last = p.last();
        
        // Perform matrix vector multiplication. Step through all rows in our chunk...
        // On Numa nodes with several cores, one may think that using multiple threads
        // on this node could improve the performance. Alas, this appears not to be the 
        // case as of 2014-05 (both aged AMD Opteron 32cores/8nodes and a newer 32/4 Intel
        // Xeon). It seems that a single thread saturates the memory channel(s). 
        // Hence we stay with the sequential version here.
        for (Index row=first; row<last; ++row)
        {
          // and compute the scalar product of (sparse) row and vector
          typename Range::block_type z(0);
          
          size_t jEnd = p.colStart(row+1-first);
          for (size_t j=p.colStart(row-first); j<jEnd; ++j)
            if (Entry::rows==1 && Entry::cols==1)  
              z.axpy(values[j][0][0],x[p.col(j)]);              // this is necessary as otherwise scalar matrix entries and 
            else                                                // vecotr-valued rhs entries do not work together
              values[j].umv(x[p.col(j)],z);
          
          // scale result as requested
          z *= a;
          if (initialize)   y[row] = z;
          else              y[row] += z;
          
          // accumulate duality product
          if (Entry::rows==Entry::cols && row<x.size())
            dp += y[row]*x[row];
        }

        return dp;
      }
      
      /**
       * \brief Transpose matrix-vector multiplication \f$ y \leftarrow y + a A^Tx \f$ or \f$ y \leftarrow aA^Tx \f$ 
       * \warning CURRENTLY THIS IS NOT THREAD-SAFE! DON'T CALL THIS FOR CHUNKS IN PARALLEL!
       */
      template <class Domain, class Range>
      void applyTransposed(Scalar a, Domain const& x, Range& y) const
      {
        auto const& p = *pat;
        Index first = p.first(); // first row in our chunk
        Index last = p.last();   // last row in our chunk
        
        assert(x.size()>=last);
        
        // We compute y_i = a sum_j A_ji^T x_j, with j in the outer loop
        for (Index j=first; j<last; ++j)                        // looping over the rows of A, i.e. the columns of A^T
        {
          size_t jEnd = p.colStart(j+1-first);
          for (size_t k=p.colStart(j-first); k<jEnd; ++k)       // looping over the columns in the row
          {
            auto i = p.col(k);                                  // the column index in A, i.e. the row index in A^T
            if (Entry::rows==1 && Entry::cols==1)               
              y[i] += a*values[k][0][0] * x[j];                 // this is necessary as otherwise scalar matrix entries and 
            else                                                // vecotr-valued rhs entries do not work together
              values[k].usmtv(a,x[j],y[i]);                     // y_i += a * (A^T)_ij x_j = a * (A_ji)^T x_j
          }
        }
      }
      
      /**
       * \brief Matrix-vector multiplication of transposed parts
       * 
       * This is required for subdiagonal parts in symmetric (lower triangular) storage .
       * 
       * \param block the chunk of which the transpose has to be considered
       * \param subdiagonal consider only true subdiagonal elements
       */
      template <class Domain, class Range>
      void gatherMirrored(Scalar a, Domain const& x, Range& y, CRSChunk<Entry,Index> const& block, bool subdiagonal, bool initToZero) const
      {
        auto const& p = *block.pat;
        Index firstRow = p.first();
        Index lastRow = p.last();
        Index firstCol = pat->first();
        Index lastCol = pat->last();

        assert(lastRow<=y.size());
        
        if (initToZero)
          for (Index row=firstRow; row<lastRow; ++row)
            y[row] = 0;
        
        for (Index row=firstRow; row<lastRow; ++row)
        {
          auto const& xrow = x[row];
          
          auto cend = p.colStartIterator(row+1-firstRow);
          auto ci = firstCol==0?                                                       // We consider all elements with column index in our
                    p.colStartIterator(row-firstRow):                                  // row chunk. Compute the starting point by binary search
                    std::lower_bound(p.colStartIterator(row-firstRow),cend,firstCol);  // (or used the shortcut if we know we start from the beginning).
          auto vi = block.values.begin() + (ci-p.colStartIterator(0));
          if (subdiagonal)
            for ( ; ci!=cend && *ci<lastCol && *ci<row; ++ci, ++vi)                    // Then step through the row until we hit the upper bound,
              vi->usmtv(a,xrow,y[*ci]);                                                // which is either the upper row range end or the diagonal,
          else                                                                         // depending on the matrix symmetry (encoded in subdiagonal).
            for ( ; ci!=cend && *ci<lastCol; ++ci, ++vi)
              vi->usmtv(a,xrow,y[*ci]);
        }
      }      
      
      /**
       * \brief Returns the sparsity pattern of this chunk.
       */
      CRSChunkPattern<Index> const& pattern() const { return *pat; }
      
      /**
       * \brief Returns an iterator to the start of the values for the given local row index
       * \param row the row index in range [0,last-first[
       */
      auto valStartIterator(Index row) 
      { 
        return values.begin()+pat->colStart(row); 
      }      
      
      /**
       * \brief Scatters given sub-matrices into the chunk by adding up their entries.
       */
      template <class LMIterator>
      void scatter(LMIterator first, LMIterator last)
      {
        Index firstRow = pat->first();
        Index lastRow = pat->last();
        
        if (first==last || firstRow==lastRow)    // empty data or empty chunk
          return;                                // -> nothing to do
        
        Index nnzPerRow = pat->nonzeroesPerRow();

        // scatter into uniform ranges
        Index n = scatterMutex.size();
        for (Index i=0; i<n; ++i)
        {
          #ifndef KASKADE_SEQUENTIAL
          boost::lock_guard<boost::mutex> lock(scatterMutex[i].get());
          #endif
          
          Index rowRangeStart = uniformWeightRangeStart(i,  n,lastRow-firstRow)+firstRow;
          Index rowRangeEnd   = uniformWeightRangeStart(i+1,n,lastRow-firstRow)+firstRow;  
          
          if (nnzPerRow == pat->columns())                     // completely dense
            scatterRowRange<0>(first,last,firstRow,rowRangeStart,rowRangeEnd);
          else if (nnzPerRow > 10*(first->cidx().size()))      // many more entries than we scatter right now
            scatterRowRange<1>(first,last,firstRow,rowRangeStart,rowRangeEnd);
          else                                                 // rather few entries per row
            scatterRowRange<2>(first,last,firstRow,rowRangeStart,rowRangeEnd);
        }
      }
      
      /**
       * \brief Performs entry-wise operations on our own matrix and some other.
       * 
       * Only entries present in both matrices are touched: a_ij = f(a_ij,b_ij)
       */
      template <class F, class EntryB, class IndexB>
      void entrywiseOp(F const& f, NumaBCRSMatrix<EntryB,IndexB> const& B) 
      {
        auto const& p = *pat;
        Index first = p.first();  // first row index 
        Index last = p.last();    // end row index
        
        for (Index row=first; row<last; ++row)
        {
          auto const& Br = B[row];                              // get the other matrix row
          size_t jEnd = p.colStart(row+1-first);                // step through all entries in this row
          auto bri = Br.begin();
          auto brend = Br.end();
          for (size_t j=p.colStart(row-first); j<jEnd; ++j)
          {
            while (bri!= brend && bri.index()<p.col(j))         // advance other row pointer until we find our entry (but do not overshoot)
              ++bri;
            if (bri!=brend && bri.index()==p.col(j))            // if the column indices coincide, apply the functor
              values[j] = f(values[j],*bri);
          }
        }
      }
      
    private:
      std::shared_ptr<CRSChunkPattern<Index>> pat;
      std::vector<Mutex> scatterMutex;
      
      // Raw data storage. The allocators used are NUMA allocators to guarantee local memory access.
      std::vector<Entry,NumaAllocator<Entry>> values;     // matrix entries
      
      // scatters given local matrices into the provided row range. The expected sparsity is
      // encoded as static template parameter in order to allow an efficient if-statement in 
      // the innermost loop
      template <int sparsity, class LMIterator>
      void scatterRowRange(LMIterator first, LMIterator last, Index startRow, Index firstRow, Index lastRow)
      {
        // simple-minded scatter implementation
        for ( ; first!=last; ++first)    // step through all local matrices
        {           
          for (auto rgl: first->ridx())  // step through all global rows affected by this matrix
          {
            Index rg = rgl.first;
            if (rg >= lastRow)           // global row indices in rgl.first are sorted ascendingly,
              break;                     // if we pass our last row, we can stop with this local matrix
            if (firstRow<=rg)
              scatterRow<sparsity>(*first,rg-startRow,rgl.second);
          }
        }
      }

      template <int sparsity, class LM, class LIndex>
      void scatterRow(LM const& a, Index rg, LIndex rl)
      {
        auto cbegin = pat->colStartIterator(rg);
        auto cend = pat->colStartIterator(rg+1);
        auto vbegin = valStartIterator(rg);
        
        if (LM::lumped) // only diagonal, i.e. column index is row index both in global and local matrix
        {
          auto pos = std::find(cbegin,cend,rg+pat->first());
          assert(pos != cend);
          vbegin[pos-cbegin] += a(rl,rl);
        }
        else 
        {
          for (auto cgl: a.cidx())
          {
            auto pos = sparsity==0? std::lower_bound(cbegin,cend,cgl.first) // find entry in row with global column number cgl.first 
                     : sparsity==1? std::lower_bound(cbegin,cend,cgl.first) // TODO: specialized find for dense matrices.
                     :              std::find(cbegin,cend,cgl.first); 
                     
            if (pos==cend)
              return;                                    // reached end of row - skip the rest (there may be more, e.g. if we are symmetric)
              
            vbegin += pos-cbegin;                        // move on to the found entry
            cbegin =  pos;
            
            *vbegin += a(rl,cgl.second);                 // add up appropriate entry from local matrix
          }
        }
      }
    };
    
    //---------------------------------------------------------------------------------------------
    
    template <class Entry, class Index>
    class NumaBCRSRow
    {
    public:
      typedef Index size_type;
      typedef Entry block_type;
      
      typedef NumaBCRSMatrixRowIterator<Entry,Index>      Iterator;
      typedef NumaBCRSMatrixConstRowIterator<Entry,Index> ConstIterator;
      
      /**
       * \brief Returns the number of entries in the current row.
       */
      size_type size() const { return std::distance(vBegin,vEnd); }
      
      /**
       * \brief Low-level access to array of (sorted) column indices in this row.
       */
      Index const* getindexptr() const { return &*cBegin; }
      
      /**
       * \brief Low-level access to array of values in this row.
       */
      block_type const* getptr() const { return &*vBegin; }
      block_type*       getptr()       { return &*vBegin; }
      
      /**
       * \brief Start of row entries.
       */
      ConstIterator begin() const { return ConstIterator(cBegin,vBegin); }
      ConstIterator end() const { return ConstIterator(cEnd,vEnd); }
      
      Iterator begin() { return Iterator(cBegin,vBegin); }
      Iterator end() { return Iterator(cEnd,vEnd); }
      
      /**
       * \brief Looks up a particular entry specified by column index.
       * \param c the column index
       * \return an iterator pointing to the found entry, or end() if none exists
       */
      ConstIterator find(Index c) const
      {
        auto ci = std::lower_bound(cBegin,cEnd,c);
        if (ci==cEnd || *ci!=c)
          return end();
        return ConstIterator(ci,vBegin+(ci-cBegin));
      }
      
      /**
       * \brief Random read access to row entries by global column index.
       * 
       * Note that this is
       * - of logarithmic complexity in the row length
       * - not the same as advancing an iterator c times.
       * 
       * In case there is no entry with column index c, a reference to 
       * a zero entry is returned.
       */
      Entry const& operator[](Index c) const
      {
        static Entry const zero(0);
        
        auto ci = std::lower_bound(cBegin,cEnd,c);
        if (ci==cEnd || *ci!=c)
          return zero;

        return vBegin[ci-cBegin];
      }
      
      
      /**
       * \brief Random write access to row entries by global column index.
       * 
       * This is of logarithmic complexity in the row length. In case there is no entry 
       * with column index c, a LookupException is thrown.
       */
      Entry& operator[](Index c) 
      {
        auto ci = std::lower_bound(cBegin,cEnd,c);
        if (ci==cEnd || *ci!=c)
          throw LookupException("write access to nonexistent matrix entry",__FILE__,__LINE__);
        
        return vBegin[ci-cBegin];
      }


    protected:
      typedef typename std::vector<Index,NumaAllocator<Index>>::const_iterator ColIterator;
      typedef typename std::vector<Entry,NumaAllocator<Entry>>::iterator       ValueIterator;
      
      ColIterator cBegin, cEnd;                               // iterators giving the range of column indices
      ValueIterator vBegin, vEnd;                             // iterators giving the range of values
      
      // Extracts the row begin and end iterators at the current position (if it is a valid position)
      void update(NumaBCRSMatrix<Entry,Index>& matrix, Index row, int chunk)
      {
        if (0<=row && row<matrix.N())                         // we have a valid row in a valid chunk in front of us 
        {
          auto& c = matrix.chunk(chunk);
          auto const& p = c.pattern();
          assert(p.first()<=row && row<p.last());             // check that the row is contained in the chunk
          Index first = p.first();
          cBegin = p.colStartIterator(row-first);
          cEnd = p.colStartIterator(row-first+1);
          vBegin = c.valStartIterator(row-first);
          vEnd = c.valStartIterator(row-first+1);
        }
      }
    };
    
    // -----------------------------------------------------------------------------------
    
    template <class Entry, class Index>
    class NumaBCRSMatrixConstIterator: public NumaBCRSRow<Entry,Index> // TODO: derive from const row type!
    {
      typedef NumaBCRSMatrixConstIterator<Entry,Index> Self;
      typedef NumaBCRSRow<Entry,Index> Row;
      using value_type = Row;
      
    public:
      NumaBCRSMatrixConstIterator(NumaBCRSMatrix<Entry,Index> const& matrix_, Index row_)
      : matrix(const_cast<NumaBCRSMatrix<Entry,Index>*>(&matrix_)), row(row_), chunk(getChunk(row_))
      {
        assert(row>=0 && row<=matrix->N());
        this->update(*matrix,row,chunk);
      }
      
      /// \{
      void operator++() 
      {
        ++row;                                                   // next row
        if (row==matrix->N())                                    // we've reached the end 
          return;                                                // -> skip anything else
        while (matrix->getPattern()->rowStart()[chunk+1]==row)   // oops, we tripped into the next chunk...
          ++chunk;                                               // which may be empty (zero rows), then go on
        this->update(*matrix,row,chunk);
      }
      
      void operator--()
      {
        --row;
        while (matrix->getPattern()->rowStart()[chunk] > row) // oops, we stepped before our chunk
          --chunk;                                            // previous one may be empty, then go on
        this->update(*matrix,row,chunk);
      }
      
      void operator+=(Index inc)
      {
        row += inc;
        if (row > matrix->N())
        {
          row = matrix->N();                  // we're behind the end - move to end 
          return;                             // such that we compare equal to end.
        }
        chunk = getChunk(row);
        this->update(*matrix,row,chunk);
      }
      
      bool operator==(Self const& it) const { return matrix==it.matrix && row==it.row; }
      bool operator!=(Self const& it) const { return !(*this==it); }
      /// \}
      
      Index index() const { return row; }
      
      ///\{
      /**
       * \brief Dereferentiation yields the row.
       * 
       * The row object is just the iterator, which provides the row interface as well. As no row objects 
       * exist independently of the iterators, the row objects are volatile. Don't store any references
       * to rows!
       */
      Row const& operator*() const { return *this; }  
      Row const* operator->() const { return this; }
      ///\}
      
    private:
      NumaBCRSMatrix<Entry,Index>* matrix;
      Index row;
      int chunk;
      
      int getChunk(Index r) const { return matrix->getPattern()->chunk(r); }
    };
    
    
    
    template <class Entry, class Index>
    class NumaBCRSMatrixIterator: public NumaBCRSMatrixConstIterator<Entry,Index>
    {
      typedef NumaBCRSMatrixIterator<Entry,Index> Self;
      typedef NumaBCRSRow<Entry,Index> Row;
      
    public:
      NumaBCRSMatrixIterator(NumaBCRSMatrix<Entry,Index>& matrix_, Index row_)
      : NumaBCRSMatrixConstIterator<Entry,Index>(matrix_,row_) {}
      
     
      Row& operator*()  { return *this; }  
      Row* operator->()  { return this; }
    };
    
    

 
  } // end of namespace ThreadedMatrixDetail
  /**
   * \endinternal
   */
  
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  /**
   * \ingroup linalgbasic
   * \brief A NUMA-aware creator for matrix sparsity patterns.
   * 
   * This supports a two-phase creation of sparsity pattern. Nonzero entry positions can be
   * added to the creator at will. At any time, a \ref NumaCRSPattern sparsity pattern can be constructed 
   * efficiently by supplying this creator.
   * 
   * \tparam Index the integral type used for row/column indices (defaults to size_t)
   */
  template <class Index=size_t>
  class NumaCRSPatternCreator
  {
    typedef ThreadedMatrixDetail::CRSChunkPatternCreator<Index> ChunkCreator;
    
  public:
    /**
     * \brief Constructs a rows times cols matrix sparsity structure creator.
     * 
     * \param rows the number of rows.
     * \param cols the number of columns.
     * \param symmetric if true, the matrix is symmetric and only the lower triangular part is actually stored (only for square matrices).
     * \param nnzPerRow hint for the expected number of nonzeroes per row
     * 
     * The rows are distributed evenly among the NUMA nodes in case of unsymmetric patterns, and 
     * roughly proportional to \f$ i^{-1/2} \f$ for symmetric ones. This balances the number of 
     * elements in each chunk if there is the same number of elements in each row.
     * 
     * Providing a sufficiently large nnzPerRow hint can reduce the number of memory reallocations performed
     * during element insertion and hence speed up the insertion significantly (a factor of 40 has been observed!)
     * at the cost of potentially increased memory footprint. A guideline to what is "sufficiently large": The 
     * memory buffer of a row needs to temporarily store the number of entries in this row plus the entries 
     * that are inserted by the \ref addElements operation to the row. E.g., for 2D linear finite elements
     * on a triangular grid, this would be 7 (number of entries per row in a regular grid, say 8 or 9 as a 
     * precaution for unstructured grids) + 3 (elemental matrices are 3x3, so 3 elements are added 
     * by an addElements operation), in total 10 to 12. If memory consumption is not the prime concern, err on 
     * the larger side.
     * 
     * The effect need not be present, and in fact, specifying a nonzero nnzPerRow may actually be slower since
     * all chunks request memory at the same time and put a high load on the memory management system. This has
     * been observed to parallelize rather badly. 
     * 
     * In the end, choose a value based on profiling *your* application.
     *
     * \warning Memory usage may be exhaustive if there are few essentially dense rows into which
     *          entries are added alternatingly with entries into many sparse rows!
     */
    NumaCRSPatternCreator(Index rows_, Index cols, bool symmetric=false, int nnzPerRow=0)
    : rows(rows_), columns(cols)
    , sym(symmetric && rows==columns)
    {
      int nodes = static_cast<Index>(NumaThreadPool::instance().nodes());
      creators.reserve(nodes);
      
      for (int i=0; i<nodes; ++i)
        if (sym)
          creators.push_back(ChunkCreator(i*rows/nodes,(i+1)*rows/nodes,cols,sym,i)); // TODO: more balanced 
        else
          creators.push_back(ChunkCreator(i*rows/nodes,(i+1)*rows/nodes,cols,sym,i));
        
      if (nnzPerRow > 0)
        parallelForNodes([this,nnzPerRow](int i, int n)
        {
          this->creators[i].reserve(nnzPerRow);
        },nodes);
    }
    
    // The default destructor does not parallelize the release of memory on the NUMA chunks.Thus we 
    // define our own destructor.
    ~NumaCRSPatternCreator()
    {
      int nodes = NumaThreadPool::instance().nodes();
      parallelForNodes([this](int i, int n)
      {
        this->creators[i].clear();
      },nodes);
    }
    
    /**
     * \brief Enters a single entry into the sparsity pattern.
     * 
     * Note that entering individual entries has a significant overhead. Use one of the other methods if possible.
     */
    void addElement(Index row, Index col)
    {
      addElements(&row,&row+1,&col,&col+1,true);
    }

    /**
     * \brief Enters a contiguous dense block into the sparsity pattern.
     *
     * \param fromRow first row in the block
     * \param toRow   one behind last row in the block
     * \param fromCol first column in the block
     * \param toCol   one behind last column in the block
     */
    void addDenseBlock(Index fromRow, Index toRow, Index fromCol, Index toCol)
    {
      assert(fromRow>=0 && toRow<=rows && fromRow<=toRow);
      assert(fromCol>=0 && toCol<=columns && fromCol<=toCol);

      for (int r=fromRow; r<toRow; ++r)           // TODO: this is probably inefficient and not
        for (int c=fromCol; c<toCol; ++c)         //       NUMA-parallelized. Consider a more
          addElement(r,c);                        //       efficient implementation
    }
    
    /**
     * \brief Enters entries into the sparsity pattern.
     * \tparam IterRow random access iterator for a range of row indices
     * \tparam IterCol random access iterator for a range of column indices
     * \param fromRow start of row indices
     * \param toRow one behind of last row index
     * \param fromCol start of sorted column indices
     * \param toCol one behind of last column index
     * \param colIsSorted true if the column index range is sorted ascendingly
     * 
     * All entries \f$ (i,j) \f$ with \f$ i \f$ in the given row indices and \f$ j \f$ in the 
     * given column indices are entered into the sparsity structure. Note that the column indices
     * have to be sorted. No index in the row range shall occur twice.
     */
    template <class IterRow, class IterCol>
    void addElements(IterRow const fromRow, IterRow const toRow, IterCol const fromCol, IterCol const toCol, bool colIsSorted=false)
    {
      assert(fromCol==toCol || *std::max_element(fromCol,toCol) < columns);
      assert(fromRow==toRow || *std::max_element(fromRow,toRow) < rows);
      for (auto& c: creators)
        c.addElements(fromRow,toRow,fromCol,toCol,colIsSorted);    
    }
    
    /**
     * \brief Enters elements into the sparsity pattern.
     * \tparam RowRangeSequence a sequence type with value type representing a row index range
     * \tparam ColRangeSequence a sequence type with value type representing a column index range
     * \param rrs the sequence of row ranges. 
     * \param crs the sequence of column ranges. The size has to be the same as that of rrs.
     * \param colsAreSorted true if all the column index ranges are sorted ascendingly
     * 
     * For each corresponding pair of row and column ranges in the sequences \arg rrs and \arg crs, the elements
     * of the cartesian product of the row and column ranges are entered into the sparsity structure. 
     * 
     * The insertion of indices is done in parallel on the NUMA chunks. In order to maintain efficiency, the number
     * or size of the ranges provided in the sequences should not be too small (a couple of ten ranges should be ok).
     */
    template <class RowRangeSequence, class ColRangeSequence>
    void addElements(RowRangeSequence const& rrs, ColRangeSequence const& crs, bool colsAreSorted=false)
    {
      assert(rrs.size() == crs.size());

      for (auto c=std::begin(crs); c!=std::end(crs); ++c)
        assert((*c).empty() || (*std::max_element(std::begin(*c),std::end(*c)) < columns)); // check that all column indices are ok
      
      parallelForNodes([&rrs,&crs,colsAreSorted,this](int i, int nodes) { 
          auto c = std::begin(crs);
          for (auto r=std::begin(rrs); r!=std::end(rrs); ++r, ++c)
            this->creators[i].addElements(std::begin(*r),std::end(*r),std::begin(*c),std::end(*c),colsAreSorted);
        },nodes());
    }
    
    /**
     * \brief Enters all possible elements (defining a dense matrix).
     * 
     * Sometimes (e.g. with spatially constant variables), FE matrices are actually dense. In this case we
     * need not create a sparsity pattern, but can simplify the process.
     * \todo This should be statically encoded in a different matrix type.
     */
    void addAllElements()
    {
      parallelForNodes([this](int i, int nodes) {
        this->creators[i].addAllElements(this->columns);
      }, nodes());
    }
    
    /**
     * \brief Enters the diagonal elements.
     */
    void addDiagonal()
    {
      for (Index i=0; i<std::min(rows,columns); ++i)
        addElement(i,i);
    }
    
    /**
     * \brief The number of structurally nonzero elements.
     * 
     * For symmetric storage without superdiagonal elements stored, this counts subdiagonal elements twice.
     */
    size_t nonzeroes() const
    {
      size_t nnz = 0;
      for (auto const& c: creators)
        nnz += c.nonzeroes();
      return nnz;
    }
    
    /**
     * \brief The number of columns.
     */
    Index cols() const { return columns; }
    
    /**
     * \brief Redistributes the rows to the NUMA chunks in order to have the same number of entries in each chunk.
     * 
     * The association of matrix rows to chunks is recalculated in order to have approximately the same number of 
     * nonzeroes in every chunk. While this is not an extremely expensive operation, it does several reallocations.
     * Call this once before creating a sparsity pattern.
     */
    void balance();
    
    /**
     * \brief Returns the number of NUMA nodes/chunks used.
     * 
     * This is guaranteed not to exceed the number of NUMA nodes, but can be less.
     */
    int nodes() const { return creators.size(); }
    
    /**
     * \brief Returns the chunk creator.
     * \param node the number of the NUMA node / chunk (0<=node<nodes()).
     */
    ChunkCreator const& creator(int node) const { return creators[node]; }
    
    /**
     * \brief Returns the symmetry status of the pattern.
     * 
     * If true, the matrix is symmetric, and only the lower triangular entries are actually stored.
     */
    bool isSymmetric() const { return sym; }
    
  private:
    Index rows, columns;
    std::vector<ChunkCreator> creators;
    bool sym;
  };

  //---------------------------------------------------------------------------

  /**
   * \ingroup linalgbasic
   * \brief A NUMA-aware compressed row storage sparsity pattern.
   * \tparam Index integral type for row/column indices
   */
  template <class Index=size_t>
  class NumaCRSPattern
  {
    typedef ThreadedMatrixDetail::CRSChunkPattern<Index> ChunkPattern;
    
  public:
    
    /**
     * \brief Constructs an empty 0x0 pattern.
     */
    NumaCRSPattern()
    : NumaCRSPattern(NumaCRSPatternCreator<Index>(0,0))
    {}
    
    /**
     * \brief Constructor creating a sparsity pattern from the given creator.
     */
    NumaCRSPattern(NumaCRSPatternCreator<Index> const& creator)
    : patterns(creator.nodes()), sym(creator.isSymmetric()), rowSt(creator.nodes()+1), cols(creator.cols())
    {
      parallelForNodes([this,&creator](int i, int n)
      {
        this->patterns[i] = std::make_shared<ChunkPattern>(creator.creator(i)); 
        this->rowSt[i] = this->patterns[i]->first();
      },creator.nodes());
      rowSt[patterns.size()] = patterns.back()->last();
    }
    
    /**
     * \brief Constructor.
     * 
     * This works like Matlab A(idx,idx), where idx is given by eIndices.
     * 
     * \tparam Expanded an array type with value type convertible to Index
     * \tparam Condensed an array type with value type convertible to Index
     * 
     * 
     * \param eIndices sorted global array of expanded indices
     * \param cIndices sorted global array of condensed indices
     * \param mat      a matrix with BCRSMatrix interface
     */
    template <class Expanded, class Condensed, class Matrix>
    NumaCRSPattern(Expanded const& eIndices, Condensed const& cIndices, Matrix const& mat)
    : patterns(NumaThreadPool::instance().nodes()), sym(false), rowSt(eIndices.size(),0), cols(eIndices.size())
    {
      assert(mat.M()==mat.N()); // works only for quadratic matrices
      assert(eIndices.size() <= cIndices.size());
      
      for (Index r=0; r<eIndices.size(); ++r)          // extract the number of entries in each 
      {                                                // row of the condensed matrix
        auto row = mat[eIndices[r]];
        for (auto ci=row.begin(); ci!=row.end(); ++ci)
          if (cIndices[ci.index()] < eIndices.size())  // only count those entries that fall in our col range
            ++rowSt[r];
      }
      
      equalWeightRanges(rowSt,patterns.size());        // balance the number of entries in each chunk
      
      for (int i=0; i<patterns.size(); ++i)            // create patterns
        patterns[i] = std::make_shared<ChunkPattern>(rowSt[i],rowSt[i+1],eIndices,cIndices,mat,i);
    }
    
    /**
     * \brief Constructor extracting the sparsity pattern of a given matrix (usually a Dune::BCRSMatrix).
     * \tparam Matrix the type of the supplied matrix to copy
     * 
     * \param matrix the matrix to be copied
     * \param isSymmetric whether the supplied matrix is symmetric
     * \param isTransposed whether the supplied matrix is transposed
     * \param symmetric whether only the lower triangular part should be stored
     * 
     * This copies the sparsity pattern of the provided Dune ISTL matrix. The number of chunks created
     * is the number of NUMA nodes as reported by the \ref NumaThreadPool.
     * 
     * Note that !isSymmetric && symmetric is highly questionable and hence not allowed.
     */
    template <class Matrix>
    NumaCRSPattern(Matrix const& matrix, bool isSymmetric, bool isTransposed, bool symmetric)
    : sym(symmetric), cols(isTransposed? matrix.N(): matrix.M())
    {
      assert(isSymmetric || !symmetric);
      
      // Compute an equilibrated distribution of rows to chunks (attempt as many chunks as NUMA nodes).
      Index rows = isTransposed? matrix.M(): matrix.N();
      std::vector<size_t> rowCount(rows,0);
      ThreadedMatrixDetail::getRowCount(matrix,isTransposed,rowCount);
      equalWeightRanges(rowCount,NumaThreadPool::instance().nodes());  
      
      // Create chunk patterns
      patterns.reserve(rowCount.size()-1);
      for (int i=0; i<rowCount.size()-1; ++i)
        patterns.push_back(std::make_shared<ChunkPattern>(matrix,isSymmetric,isTransposed, 
                                                          rowCount[i],rowCount[i+1],symmetric,i));
        
      rowSt.resize(patterns.size()+1);
      for (int i=0; i<patterns.size(); ++i)
        rowSt[i] = patterns[i]->first();
      rowSt[patterns.size()] = patterns.back()->last();
    }
        
    /**
     * \brief The number of rows.
     */
    Index N() const { return rowSt.back(); }
    
    /**
     * \brief The number of columns.
     */
    Index M() const { return cols; }
    
    /**
     * \brief Returns the number of stored entries.
     * 
     * This differs from nonzeroes() for symmetric storage.
     */
    size_t storage() const 
    {  
      size_t nnz = 0;
      for (auto const& p: patterns)
        nnz += p->storage();
      return nnz;
    }
    
    /**
     * \brief Returns the number of structurally nonzero entries.
     */
    size_t nonzeroes() const 
    {  
      size_t nnz = 0;
      for (auto const& p: patterns)
        nnz += p->nonzeroes();
      return nnz;
    }

    /**
     * \brief Returns the number of NUMA nodes/chunks used.
     */
    int nodes() const { return patterns.size(); }
    
    /**
     * \brief Returns the individual patterns.
     */
    std::shared_ptr<ChunkPattern> pattern(int i) const { return patterns[i]; }
    
    /**
     * \brief Returns the symmetry status of the pattern.
     * 
     * If true, the matrix is symmetric, and only the lower triangular entries are actually stored.
     */
    bool isSymmetric() const { return sym; }
    
    /**
     * \brief Returns the number of the chunk containing the given row.
     * \param row the index of the row in the range [0,N]
     * 
     * For row==N, the total number of chunks is returned (i.e. one behind the last chunk).
     */
    int chunk(Index row) const 
    { 
      assert(0<=row && row<=N());
      auto it = std::upper_bound(rowSt.begin(),rowSt.end(),row);  // use upper bound here, with lower_bound rows 0 and 1 would end up
      int c = it - rowSt.begin() - 1;                             // in different chunks...
      assert(0<=c && c<=patterns.size());
      assert(c==patterns.size() || (patterns[c]->first() <= row && row < patterns[c]->last()));
      return c;
    }
    
    /**
     * \brief Returns the limiting row indices between the chunks.
     * 
     * For \f$ 0\le i < n \f$, the chunk \f$ i \f$ contains the rows in the half-open range
     * [ rowStart()[i], rowStart()[i+1] [.
     */
    std::vector<Index> const& rowStart() const { return rowSt; }

    /**
     * \brief queries whether an entry is present in the sparsity pattern
     */
    bool exists(Index r, Index c) const
    {
      return patterns[chunk(r)]->exists(r,c);
    }

  private:
    std::vector<std::shared_ptr<ChunkPattern>> patterns;
    bool sym;
    std::vector<Index> rowSt;   // row start indices
    Index cols;                 // how many columns the pattern has
  };
  
  template <class Index, class Index2>
  NumaCRSPattern<Index> operator+(NumaCRSPattern<Index> const& pa, NumaCRSPattern<Index2> const& pb)
  {
    assert(pa.N()==pb.N() && pa.M()==pb.M());
    Index nnzPerRow = (pa.nonzeroes()+pb.nonzeroes())/pa.N();
    NumaCRSPatternCreator<Index> creator(pa.N(),pa.M(),pa.isSymmetric()&&pb.isSymmetric(),nnzPerRow);
    
    for (int i=0; i<pa.nodes(); ++i)
    {
      bool symmetrize = pa.isSymmetric() && !pb.isSymmetric();
      auto const& chunkPattern = *pa.pattern(i);
      Index first = chunkPattern.first();
      for (Index row=first; row<chunkPattern.last(); ++row)
        for (auto it=chunkPattern.colStartIterator(row-first); it!=chunkPattern.colStartIterator(row+1-first); ++it)
        {
          creator.addElement(row,*it);
          if (symmetrize)
            creator.addElement(*it,row);
        }
    }
    
    for (int i=0; i<pb.nodes(); ++i)
    {
      bool symmetrize = pb.isSymmetric() && !pa.isSymmetric();
      auto const& chunkPattern = *pb.pattern(i);
      Index2 first = chunkPattern.first();
      for (Index2 row=first; row<chunkPattern.last(); ++row)
        for (auto it=chunkPattern.colStartIterator(row-first); it!=chunkPattern.colStartIterator(row+1-first); ++it)
        {
          creator.addElement(static_cast<Index>(row),*it);
          if (symmetrize)
            creator.addElement(*it,static_cast<Index>(row));
        }
    }
    
    return NumaCRSPattern<Index>(creator);
  }

  //---------------------------------------------------------------------------

  /**
   * \ingroup linalgbasic
   * \brief A NUMA-aware compressed row storage matrix adhering mostly to the Dune ISTL interface (to complete...)
   * 
   * This class distributes the matrix rows block-wise across available NUMA nodes to exploit the higher
   * memory bandwidth in parallelized matrix-vector products and similar operations. 
   * 
   * As the distribution of matrix entries to the NUMA nodes is optimized for equidistribution of storage size 
   * and computational effort in matrix vector products, transpose times vector multiplications are much slower 
   * than plain matrix-vector multiplications (roughly factor 6).
   * 
   * If the matrix is a priori specified as being symmetric, only the lower triangular elements are stored and 
   * can be accessed. Nevertheless, all linear algebra operations (i.e., matrix vector products) act on the 
   * complete symmetric matrix. Note that symmetric storage incurs a significant performance penalty for matrix-vector
   * products (roughly factor 3), and should only be used if those operations are rare. In contrast, transpose 
   * times vector multiplications do not suffer.
   * 
   * For creating a NumaBCRSMatrix, a sparsity pattern (\ref NumaCRSPattern) has to be provided. This in turn is
   * created using a creator object (\ref NumaCRSPatternCreator). This multi-stage creation allows efficient 
   * creation while avoiding stateful interfaces (incompletely created objects). The matrix creation can be done
   * as follows:
   * \code
   * // define the pattern creator for a 9x9 matrix
   * NumaCRSPatternCreator<> creator(9,9);
   * // add a diagonal entry
   * creator.addElement(2,2);
   * // create the matrix
   * NumaBCRSMatrix<Dune::FieldMatrix<double,1,1>> matrix(creator);
   * // fill the entry
   * matrix[2][2] = 1.0;
   * \endcode
   * 
   * Concurrent calls to one NumaBCRSMatrix object are not safe.
   * 
   * \tparam Entry the type of matrix elements (in general Dune::FieldMatrix<double,n,m>)
   * \tparam Index integral type for row/column indices (usually int, long, or size_t), defaults to size_t
   */
  template <class Entry, class Index=size_t>
  class NumaBCRSMatrix
  {
    typedef ThreadedMatrixDetail::CRSChunk<Entry,Index> Chunk;
    typedef NumaBCRSMatrix<Entry,Index> Self;
    
  public:
    typedef ScalarType<Entry> Scalar;
    typedef Scalar field_type;        // for compatibility with Dune
    typedef Entry block_type;
    using value_type = Entry;         // or should this be row_type?
    using size_type = Index;          // compatibility with Dune::BCRSMatrix
    
    /**
     * \brief iterator type stepping through the rows
     */
    typedef ThreadedMatrixDetail::NumaBCRSMatrixIterator<Entry,Index>      iterator;
    typedef ThreadedMatrixDetail::NumaBCRSMatrixConstIterator<Entry,Index> const_iterator;
    typedef iterator       row_type;
    typedef const_iterator const_row_type;
    
    typedef iterator ConstRowIterator;
    typedef typename iterator::ConstIterator ConstColIterator;
    

    /**
     * \name Constructors
     * @{
     */
    
    /**
     * \brief Constructs an empty 0x0 matrix.
     */
    NumaBCRSMatrix()
    : NumaBCRSMatrix(std::make_shared<NumaCRSPattern<Index>>())
    {}
    
    /**
     * \brief Copy constructor.
     */
    NumaBCRSMatrix(Self const& A) = default; // TODO: do this in parallel
    
    /**
     * \brief Move constructor.
     */
    NumaBCRSMatrix(Self&& A) = default;  // move should be efficient without parallelization
        
    /**
     * \brief Constructor creating a matrix from a given sparsity pattern creator.
     * 
     * This is a convenience constructor that creates a new pattern under the hood.
     * 
     * \param creator the sparsity pattern creatorof the matrix to be constructed
     * \param init    initialization value for the entries
     */
    NumaBCRSMatrix(NumaCRSPatternCreator<Index> const& creator, Entry const& init=Entry(0))
    : NumaBCRSMatrix(std::make_shared<NumaCRSPattern<Index>>(creator),init)
    {}
    
    /**
     * \brief Constructor creating a matrix from a given sparsity pattern.
     * 
     * \param pattern the sparsity pattern of the matrix to be constructed
     * \param init    initialization value for the entries
     */
    NumaBCRSMatrix(std::shared_ptr<NumaCRSPattern<Index>> const& pattern_, Entry const& init=Entry(0))
    : pattern(pattern_)
    {
      for (int i=0; i<pattern->nodes(); ++i)
        chunks.push_back(Chunk(i));
      
      parallelForNodes([this,&init] (int i, int n) 
      { 
          this->chunks[i] = std::move(Chunk(this->pattern->pattern(i),init));
      },chunks.size());    

      dps.resize(chunks.size());      
    }
    
    /**
     * \brief Deprecated. Use simpler version without cIndices instead.
     *
     * This constructor will be made private after 2022-12.
     */
    template <class Expanded, class Condensed, class Matrix>
    NumaBCRSMatrix(Expanded const& eIndices, Condensed const& cIndices, Matrix const& mat)
    : pattern(std::make_shared<NumaCRSPattern<Index>>(eIndices,cIndices,mat))   // TODO: move this to private, and have a convenience constructor creating cIndices under the hood
    {
      chunks.reserve(pattern->nodes());
      for (int i=0; i<pattern->nodes(); ++i)
        chunks.push_back(Chunk(pattern->pattern(i),eIndices,cIndices,mat)); // TODO: do this in parallel
      dps.resize(chunks.size());
    }

  private:

    // Creates a vector cIndices that for each index in the source matrix tells
    // the index in the extracted submatrix (or a too large sentinel value if the
    // index is not contained in the submatrix). I.e. result[eIndices[i]]==i
    // for 0 <= i < eIndices.size(), and result[j]>=eIndices.size() otherwise.
    template <class Expanded>
    static std::vector<Index> createCondensed(Expanded const& eIndices, Index n)
    {
      std::vector<Index> cIndices(n,eIndices.size());
      for (Index i=0; i<eIndices.size(); ++i)
        cIndices[eIndices[i]] = i;
      return cIndices;
    }

  public:

    /**
     * \brief Indexed submatrix constructor.
     *
     * This works like Matlab A(idx,idx), where idx is given by eIndices.
     * It creates a \f$ n \times n \f$ matrix \f$ a \f$  from the \f$ N \times \f$ N matrix \f$ A \f$, where
     * \f$ a_{ij} = A_{e_ie_j} \f$. Of course, \f$ n \le N \f$ has to hold.
     *
     * \tparam Expanded an array type with value type convertible to Index
     *
     * \param eIndices sorted global array of size \f$ n \f$ of expanded indices
     * \param mat      a square matrix with BCRSMatrix interface
     */
    template <class Expanded, class Matrix>
    NumaBCRSMatrix(Expanded const& eIndices, Matrix const& mat)
    : NumaBCRSMatrix(eIndices,createCondensed(eIndices,mat.N()),mat)
    { }

    
    /**
     * \brief Constructor copying a given matrix.
     * \tparam Matrix the type of the supplied matrix to copy (usually a Dune::BCRSMatrix 
     *                or a NumaBCRSMatrix with different scalar type).
     * 
     * \param pattern the sparsity pattern to use
     * \param matrix the matrix to be copied (shall have the given sparsity pattern)
     * \param isSymmetric whether the supplied matrix is symmetric (then only its lower triangular part is accessed)
     * \param isTransposed whether the supplied matrix is transposed
     */
    template <class Matrix>
    NumaBCRSMatrix(std::shared_ptr<NumaCRSPattern<Index>> const& pattern_, Matrix const& matrix, 
                   bool isSymmetric, bool isTransposed)
    : pattern(pattern_)
    {
      chunks.reserve(pattern->nodes());
      for (int i=0; i<pattern->nodes(); ++i)
        chunks.push_back(Chunk(pattern->pattern(i),matrix,isSymmetric,isTransposed)); // TODO: do this in parallel
      dps.resize(chunks.size());
    }
    
    /**
     * \brief Constructor copying a given matrix.
     * 
     * \tparam OtherEntry The entry type of the matrix to be copied. The entries shall have 
     *                    the same dimension (but possibly different scalar type).
     * \param matrix the matrix to be copied (shall conform to the Dune::BCRSMatrix interface)
     * 
     * Use this for switching between different entry types (e.g. float vs double). 
     * The sparsity pattern is shared between both matrices.
     * 
     * The number of chunks created is the number of NUMA nodes as reported by \ref NumaThreadPool.
     */
    template <class OtherEntry>
    NumaBCRSMatrix(NumaBCRSMatrix<OtherEntry,Index> const& matrix)
    : NumaBCRSMatrix(matrix.getPattern(),matrix,matrix.getPattern()->isSymmetric(),false)
    { }

    /**
     * \brief Constructor.
     * \tparam Matrix the type of the supplied matrix to copy (shall conform to the Dune::BCRSMatrix interface)
     * 
     * \param matrix the matrix to be copied
     * \param isSymmetric whether the supplied matrix is symmetric (then only its lower triangular part 
     *                    is accessed, requires square matrix and entries)
     * \param isTransposed whether the supplied matrix is transposed
     * \param symmetric whether only the lower triangular part should be stored (requires a square matrix and entries)
     * 
     * This is a convenience constructor that creates a new sparsity pattern from scratch. The number of chunks
     * created is the number of NUMA nodes as reported by \ref NumaThreadPool.
     * 
     * Note that !isSymmetric && symmetric is highly questionable and hence not allowed.
     */
    template <class Matrix>
    NumaBCRSMatrix(Matrix const& matrix, bool isSymmetric, bool isTransposed=false, bool symmetric=false)
    : NumaBCRSMatrix(std::make_shared<NumaCRSPattern<Index>>(matrix,isSymmetric,isTransposed,symmetric),
                     matrix,isSymmetric,isTransposed)
    {
      assert(matrix.N()==matrix.M() || !(symmetric||isSymmetric));   // only go into symmetric mode if the matrix is square
      assert(Entry::rows==Entry::cols || !(symmetric||isSymmetric)); // symmetry only if the entries are square

      // Make sure symmetry is handled correctly (even if not in debug mode). TODO: shall we take this as a feature instead of complaining via assert?
      symmetric &= Entry::rows==Entry::cols && matrix.N()==matrix.M();
      isSymmetric &= Entry::rows==Entry::cols && matrix.N()==matrix.M();
    }
    
    /**
     * \brief Constructor.
     * 
     * \param matrix the matrix to be copied (shall conform to the Dune::BCRSMatrix interface)
     * \param isSymmetric whether the supplied matrix is symmetric (then only its lower triangular part is accessed, requires square matrix and entries)
     * \param isTransposed whether the supplied matrix is transposed
     * \param symmetric whether only the lower triangular part should be stored (requires a square matrix and entries)
     * 
     * This is a convenience constructor that creates a new sparsity pattern from scratch. Use it for switching between symmetric/normal storage patterns
     * and for creating transposed matrices.
     * 
     * The number of chunks created is the number of NUMA nodes as reported by \ref NumaThreadPool.
     * 
     * Note that !isSymmetric && symmetric is highly questionable and hence not allowed.
     */
    NumaBCRSMatrix(NumaBCRSMatrix<Entry,Index> const& matrix, bool isSymmetric, bool isTransposed=false, bool symmetric=false)
    : NumaBCRSMatrix(isTransposed==false && isSymmetric==symmetric? matrix.pattern            // just copy
                                                                  : std::make_shared<NumaCRSPattern<Index>>(matrix,isSymmetric,isTransposed,symmetric),
                     matrix,isSymmetric,isTransposed)
    {
      assert(matrix.N()==matrix.M() || !(symmetric||isSymmetric));   // only go into symmetric mode if the matrix is square
      assert(Entry::rows==Entry::cols || !(symmetric||isSymmetric)); // symmetry only if the entries are square
      assert(isTransposed? (matrix.N()==M() && matrix.M()==N()): (matrix.N()==N() && matrix.M()==M()));
      assert(!isTransposed || Entry::rows==Entry::cols);
    }
    
    
    /// @}
    
    /**
     * \name Assignment and simple in-place modifications
     * @{
     */
    
    /**
     * \brief Copy assignment.
     * 
     * The matrix becomes a copy of the given one, no matter what its previous structure or size was. 
     * The sparsity pattern is shared between both matrices.
     */
    Self& operator=(Self const& mat) = default;
    
    /**
     * \brief Move assignment.
     * 
     * The matrices need not be of same size or sparsity pattern.
     */
    Self& operator=(Self&& mat) = default;
    
    /**
     * \brief Assigns the given value to each entry.
     */
    Self& operator=(Entry const& a)
    {
      if (chunks.size()>0) 
        parallelForNodes([this,&a] (int i, int n) 
        { 
          this->chunks[i] = a;  
        },chunks.size());
      
      return *this;
    }
    
    /**
     * \brief Assigns the given scalar value to each entry.
     */
    Self& operator=(typename Entry::field_type const& a)
    {
      return *this = Entry(a);
    }
    
    /**
     * \brief NOT YET IMPLEMENTED Assigns the given Numa matrix expression.
     */
    template <class Arguments, class Operation>
    Self& operator=(ThreadedMatrixDetail::NumaBCRSMatrixExpression<Arguments,Operation> const& e)
    {
      if (chunks.size()>0) 
        parallelForNodes([this,&e] (int i, int n) 
        { 
          this->chunks[i] = e[i];  
        },chunks.size());
      
      return *this;
    }
    
    /**
     * \brief Adds a sparse matrix to this one.
     * 
     * The sparsity pattern of B must be a subset of our sparsity pattern.
     */
    template <class EntryB, class IndexB>
    Self& operator+=(NumaBCRSMatrix<EntryB,IndexB> const& B)
    {
      entrywiseOp([](Entry const& a, Entry const& b) { return a+b; }, B);
      return *this;
    }
    
    /**
     * \brief Subtracts a sparse matrix from this one.
     * 
     * The sparsity pattern of B must be a subset of our sparsity pattern.
     */
    template <class EntryB, class IndexB>
    Self& operator-=(NumaBCRSMatrix<EntryB,IndexB> const& B)
    {
      entrywiseOp([](Entry const& a, Entry const& b) { return a-b; }, B);
      return *this;
    }
    
    /**
     * \brief Multiplication with a "scalar".
     * 
     * The factor a can be a scalar or a quadratic matrix of compatible static size.
     */
    template <class Factor>
    Self& operator*=(Factor a)
    {
      entrywiseOp([=](Entry const& e, Entry const&)  { return a*e; }, *this);
      return *this;
    }
    
    /// @}
    
    /**
     * \name Element access
     *
     * Matrix entries can be accessed via iterators or by direct access. Iterators
     * allow row-major scanning of existing matrix entries, with an iterator iterating
     * over the rows and for each row an iterator iterating over the entries in that
     * row.
     *
     * Direct access is provided by the subscript operator [], which returns a row
     * (reference), which in turn provides a subscript operator for selecting the
     * column, e.g.,
     * \code
     * A[i][j] = 42;
     * \endcode
     * Note that the second (column) subscript has logarithmic complexity in the number
     * of entries in the row.
     * @{
     */
    
    /**
     * \brief returns an iterator to the first row
     * 
     * In symmetric storage, the iterators iterate only through the lower triangular part
     * that is actually stored.
     */
    iterator       begin()       { return iterator(*this,0); }
    const_iterator begin() const { return const_iterator(*this,0); }
    
    /**
     * \brief returns an iterator to the first row
     */
    iterator       end()       { return iterator(*this,N()); }
    const_iterator end() const { return const_iterator(*this,N()); }
    
    /**
     * \brief Subscript operator allowing random access to rows.
     *
     * This has constant complexity.
     */
    row_type       operator[](Index r)       { return iterator(*this,r); }
    const_row_type operator[](Index r) const { return const_iterator(*this,r); }
    
    /// @}

    /**
     * \name Submatrix extraction
     * @{
     */
    template <class RowIndices, class ColIndices>
    Self operator()(RowIndices const& ri, ColIndices const& ci) const
    {
      return submatrix<Self>(*this,ri,ci);
    };

    /**
     * @}
     */
    
    /**
     * \name General matrix information.
     * @{
     */
    
    /**
     * \brief The number of rows.
     */
    Index N() const { return pattern->N(); }
    
    /**
     * \brief The number of columns.
     */
    Index M() const { return pattern->M(); }
    
    /**
     * \brief Returns the number of structurally nonzero elements.
     * 
     * In symmetric storage, the number of actually stored entries is smaller,
     * usually by a factor of almost two.
     */
    size_t nonzeroes() const { return pattern->nonzeroes(); }
    
    /**
     * \brief Returns a pointer to the sparsity pattern.
     */
    std::shared_ptr<NumaCRSPattern<Index>> getPattern() const { return pattern; }
    
    /**
     * \brief Obtains a reference to the given chunk.
     */
    Chunk& chunk(int i) { return chunks[i]; }

    /**
     * \brief returns true if (r,c) is structurally nonzero
     */
    bool exists(Index r, Index c) const
    {
      return pattern->exists(r,c);
    }

    /**
     * \brief Computes the Frobenius norm \f$ \|A\|_F = \sqrt{\mathrm{tr}(A^TA)} \f$.
     * \todo parallelize
     */
    field_type frobenius_norm() const
    {
      return sqrt(frobenius_norm2());
    }

    /**
     * \brief Computes the square of the Frobenius norm \f$ \|A\|_F^2 = \mathrm{tr}(A^TA) \f$.
     * \todo parallelize
     */
    field_type frobenius_norm2() const
    {
      field_type sum = 0;
      for (auto ri=begin(); ri!=end(); ++ri)
        for (auto ci=ri->begin(); ci!=ri->end(); ++ci)
          sum += ci->frobenius_norm2();
      return sqrt(sum);
    }

    /// @}
    
    /**
     * \name Matrix-vector operations
     * @{
     */
    
    /**
     * \brief Matrix-vector multiplication \f$ y \leftarrow Ax \f$ with computation of \f$ y^T x \f$ if A is square.
     * \param x a vector of size M()
     * \param y a vector of size N()
     */
    template <class X, class Y>
    field_type mv(X const& x, Y& y) const { return doMv(1.0,x,y,true); }
    
    /**
     * \brief Matrix-vector multiplication \f$ y \leftarrow A^Tx \f$.
     * \param x a vector of size M()
     * \param y a vector of size N()
     * 
     * Note that transpose times vector operations are quite expensive (a factor of 6 slower than matrix-vector operations)
     * and should be avoided. If you need transpose times vector products, consider storing a transposed matrix directly.
     */
    template <class X, class Y>
    void mtv(X const& x, Y& y) const { doMtv(1.0,x,y,true); }
    
    /**
     * \brief Matrix-vector multiplication \f$ y \leftarrow -Ax \f$.
     */
    template <class X, class Y>
    void mmv(X const& x, Y& y) const { doMv(-1.0,x,y,true); }
    
    /**
     * \brief Matrix-vector multiplication \f$ y \leftarrow aAx \f$.
     */
    template <class X, class Y>
    field_type smv(field_type const& a, X const& x, Y& y) const { return doMv(a,x,y,true); }
    
    /**
     * \brief Matrix-vector multiplication \f$ y \leftarrow aA^Tx \f$.
     * 
     * Note that transpose times vector operations are quite expensive (a factor of 6 slower than matrix-vector operations)
     * and should be avoided. If you need transpose times vector products, consider storing a transposed matrix directly.
     */
    template <class X, class Y>
    void smtv(field_type const& a, X const& x, Y& y) const { doMtv(a,x,y,true); }
    
    /**
     * \brief Matrix-vector multiplication \f$ y \leftarrow y + Ax \f$ and subsequent computation of \f$ y^T x \f$ if A is square.
     */
    template <class X, class Y>
    field_type umv(X const& x, Y& y) const { return doMv(1.0,x,y,false); }
    
    /**
     * \brief Matrix-vector multiplication \f$ y \leftarrow y + A^Tx \f$.
     * 
     * Note that transpose times vector operations are quite expensive (a factor of 6 slower than matrix-vector operations)
     * and should be avoided. If you need transpose times vector products, consider storing a transposed matrix directly.
     */
    template <class X, class Y>
    void umtv(X const& x, Y& y) const { doMtv(1.0,x,y,false); }
    
    /**
     * \brief Matrix-vector multiplication \f$ y \leftarrow y + aAx \f$ and subsequent computation of \f$ y^T x \f$ if A is square.
     */
    template <class X, class Y>
    field_type usmv(field_type const& a, X const& x, Y& y) const { return doMv(a,x,y,false); }
    
    /**
     * \brief Matrix-vector multiplication \f$ y \leftarrow y + aA^Tx \f$.
     * 
     * Note that transpose times vector operations are quite expensive (a factor of 6 slower than matrix-vector operations)
     * and should be avoided. If you need transpose times vector products, consider storing a transposed matrix directly.
     */
    template <class X, class Y>
    void usmtv(field_type const& a, X const& x, Y& y) const { doMtv(a,x,y,false); }
    
    /// @}
    
    
    /**
     * \brief Scatters given sub-matrices into the matrix by adding up their entries.
     * 
     * \tparam LMIterator a forward iterator with value type of LocalMatrix (requirements see below).
     * 
     * Local matrices \f$ a \f$ in the range [first,last[ are scattered into the global matrix as follows:
     * \f[ A_{I_k,J_l} \leftarrow A_{I_k,J_l} + a_{i_k,j_l} \quad 0\le k < r, 0\le l < c \f]
     * The possibly different ordering of entries in the local and global matrices is taken into account by
     * the indirect indexing via \f$ (I_k,J_k) \f$ for the global indices and \f$ (i_k,j_l) \f$ for the local
     * indices.
     * 
     * Objects a of LocalMatrix type have to provide the following access:
     * - a(ik,jl) is the value of the local matrix entry \f$ a_{i_k,j_l} \f$
     * - ridx() returns a range of \f$ (I_k,i_k) \f$ pairs (of type std::pair), sorted ascendingly according to \f$ I_k \f$
     * - cidx() returns a range of \f$ (J_l,j_l) \f$ pairs (of type std::pair), sorted ascendingly according to \f$ J_l \f$
     * 
     * In case LocalMatrix::lumped is true, then ridx()==cidx() shall hold and only the diagonal of \f$ a \f$ is to be scattered.
     * 
     * In contrast to the other methods, this method can (and should) be called concurrently from different threads.
     */
    template <class LMIterator>
    void scatter(LMIterator first, LMIterator last)
    {
      for (auto& c: chunks)
        c.scatter(first,last);
    }
    
  private:
    std::shared_ptr<NumaCRSPattern<Index>> pattern;
    std::vector<Chunk>                     chunks;
    mutable std::vector<field_type>        dps;     // could be local to doMv, but we don't like frequent reallocations
    
    
    // Perform matrix-vector multiplication y = (initToZero? 0: y) + a A x and return  a x^T Ax
    template <class X, class Y>
    Scalar doMv(Scalar const& a, X const& x, Y& y, bool initToZero) const
    {
      assert(y.size()==N());
      assert(x.size()==M());
      
      if (chunks.size()>1) 
      {
        parallelForNodes([this,a,&x,&y,initToZero] (int i, int n) 
        { 
          if (pattern->isSymmetric()) // if it's symmetric, we have to add the transpose of the other blocks
            for (int k=i; k<n; ++k)
              this->chunks[i].gatherMirrored(a,x,y,this->chunks[k],true,initToZero&&k==i);
          this->dps[i] = this->chunks[i].apply(a,x,y,initToZero&&!pattern->isSymmetric());  
        },chunks.size());
        
        return std::accumulate(dps.begin(),dps.end(),0.0);
      }
      else 
      {
        if (pattern->isSymmetric())
          chunks[0].gatherMirrored(a,x,y,chunks[0],true,initToZero);
        return chunks[0].apply(a,x,y,initToZero&&!pattern->isSymmetric());
      }
    }
    
    // Perform transpose matrix-vector multiplication y = a A^T x and return y^T x
    template <class X, class Y>
    void doMtv(Scalar const& a, X const& x, Y& y, bool initToZero) const
    {
      assert(y.size()==M());
      assert(x.size()==N());
      
      if (pattern->isSymmetric())  // matrix is symmetric, hence we can multiply with A instead of A^T
        doMv(a,x,y,initToZero);
      else
      {
        if (initToZero)
          y = 0;
        for (auto const& chunk: chunks)  // WARNING: applyTransposed is not thread-safe due to global unstructured scattering
          chunk.applyTransposed(a,x,y);  // into y. It is possible to parallelize this, in a window-shingled fashion, but is
      }                                  // it worth the implementation effort?
    }
    
    // performs elementwise operations on two matrices (most useful for axpy-like operations):
    // a_ij = f(a_ij,b_ij)
    // This affects exactly the entries present in *both* matrices. Others are unchanged. Matrices have to have the same shape.
    template <class F, class EntryB, class IndexB>
    void entrywiseOp(F const& f, NumaBCRSMatrix<EntryB,IndexB> const& B)
    {
      assert(N()==B.N() && M()==B.M());
      parallelForNodes([this,&f,&B](int i, int n)
      {
        this->chunks[i].entrywiseOp(f,B);
      },chunks.size());
    }
  };
  
  // ------------------------------------------------------------------------------------------
  
  /**
   * \ingroup linalgbasic
   * \brief Writes a NumaBCRSMatrix to an output stream.
   * 
   * The matrix is written as text in triplet format, one entry per row.
   * \relates NumaBCRSMatrix
   */
  template <class Entry, class Index>
  std::ostream& operator <<(std::ostream& out, NumaBCRSMatrix<Entry,Index> const& A)
  {
    for (auto ri=A.begin(); ri!=A.end(); ++ri)
      for (auto ci=ri->begin(); ci!=ri->end(); ++ci)
        out << ri.index() << ' ' << ci.index() << ' ' << *ci << std::endl;
    return out;
  }

  // ----------------------------------------------------------------------------------------------
  
  /**
   * \ingroup linalgbasic
   * \relates NumaBCRSMatrix
   * \brief Converts a sparse NumaBCRSMatrix to a dense matrix.
   * 
   * Be aware that, depending on the number of nonzeros in the sparse matrix, the memory requirement of
   * the dense matrix may be *much* higher. Use this only for matrices up to 1000x1000, and rather less.
   */
  template <class SparseMatrix>
  DynamicMatrix<typename SparseMatrix::block_type> full(SparseMatrix const& A)
  {
    DynamicMatrix<typename SparseMatrix::block_type> B(A.N(),A.M());
    for (auto ri=A.begin(); ri!=A.end(); ++ri)
      for (auto ci=ri->begin(); ci!=ri->end(); ++ci)
        B[ri.index()][ci.index()] = *ci;
    return B;
  }

  /**
   * \ingroup linalgbasic
   * \relates NumaBCRSMatrix
   * \brief Converts a subrange of a sparse matrix to a dense matrix.
   * 
   * This is equivalent to Matlab's full(A(rows,cols)) (except that we index 0-based).
   *
   * \tparam SparseMatrix a sparse matrix type adhering to the Dune::BCRSMatrix interface (e.g., NumaBCRSMatrix)
   * 
   * \param rows an STL range of row indices
   * \param cols an STL range of column indices
   *
   * The row and column index ranges need not be sorted, but cols must not contain duplicates.
   * 
   * Be aware that, depending on the number of nonzeros in the sparse matrix, the memory requirement of
   * the dense matrix may be *much* higher. Use this only for matrices up to 1000x1000, and rather less.
   */
  template <class SparseMatrix, class RowRange, class ColRange>
  DynamicMatrix<typename SparseMatrix::block_type> full(SparseMatrix const& A, RowRange const& rows, ColRange const& cols)
  {
    using std::begin;
    using std::end;
    using Index = typename ColRange::value_type;

    // First we prepare a lookup structure for efficient determination where
    // a sparse matrix entry should end up.
    std::vector<std::pair<Index,Index>> cidx;
    cidx.reserve(cols.size());

    Index i=0;
    for (auto c: cols)
    {
      cidx.push_back(std::make_pair(c,i));
      ++i;
    }
    std::sort(begin(cidx),end(cidx),FirstLess());

    // Now create the matrix and extract the entries.
    DynamicMatrix<typename SparseMatrix::block_type> B(rows.size(),cols.size());

    for (size_t r=0; r<rows.size(); ++r)
      for (auto ci=A[rows[r]].begin(); ci!=A[rows[r]].end(); ++ci)
      {
        auto it = std::lower_bound(begin(cidx),end(cidx),
                                   std::make_pair(ci.index(),ci.index()),FirstLess());
        if (it!=end(cidx) && it->first==ci.index())       // this is an entry from our column range
          B[r][it->second] = *ci;
      }

    return B;
  }
  
  // ----------------------------------------------------------------------------------------------
  
  /**
   * \ingroup linalgbasic
   * \relates NumaBCRSMatrix
   * \brief Matrix addition \f$ (A,B) \mapsto A+B \f$.
   * The sparsity patterns of both matrices can be different. The size of the matrices have to be the same. 
   * Both must not have symmetric storage.
   */
  template <class Entry, class Index, class Index2>
  NumaBCRSMatrix<Entry,Index> operator+(NumaBCRSMatrix<Entry,Index> const& A, NumaBCRSMatrix<Entry,Index2> const& B)
  {
    auto pat = std::make_shared<NumaCRSPattern<Index>>(*A.getPattern()+*B.getPattern());
    NumaBCRSMatrix<Entry,Index> AB(pat);
    assert(!A.getPattern()->isSymmetric() && !B.getPattern()->isSymmetric());
    
    // TODO: parallelize
    for (Index row=0; row<AB.N(); ++row)
    {
      auto ABrow = AB[row];
      for (auto ci=A[row].begin(); ci!=A[row].end(); ++ci)
        ABrow[ci.index()] += *ci;
      for (auto ci=B[row].begin(); ci!=B[row].end(); ++ci)
        ABrow[static_cast<Index>(ci.index())] += *ci;
    }
    return AB;
  }
  
  // ----------------------------------------------------------------------------------------------

  
  /**
   * \ingroup linalgbasic
   * \brief creates a unit matrix in NUMA block compressed row storage format
   * 
   * The entries are fixed-dimensional unit matrices.
   * 
   * \tparam Scalar the field type of entries
   * \tparam n dimension of the entries
   * 
   * \param N the size of the matrix to create (N blocks of size n each)
   */
  template <class Scalar, int n, class Index=size_t>
  NumaBCRSMatrix<Dune::FieldMatrix<Scalar,n,n>,Index> sparseUnitMatrix(Index N)
  {
    NumaCRSPatternCreator<Index> creator(N,N,false,1);
    for (Index i=0; i<N; ++i)
      creator.addElement(i,i);
    return NumaBCRSMatrix<Dune::FieldMatrix<Scalar,n,n>,Index>(creator,unitMatrix<Scalar,n>());
  }
  
  
  // ----------------------------------------------------------------------------------------------

  
    /**
   * \ingroup linalgbasic
   * \brief creates a zero matrix in NUMA block compressed row storage format
   * 
   * The entries are fixed-dimensional zero matrices.
   * 
   * \tparam Scalar the field type of entries
   * \tparam n dimension of the entries
   * 
   * \param N,M the size of the matrix to create (NxM blocks of size n each)
   */
  template <class Scalar, int n, int m, class Index=size_t>
  NumaBCRSMatrix<Dune::FieldMatrix<Scalar,n,m>,Index> sparseZeroMatrix(Index N, Index M)
  {
    NumaCRSPatternCreator<Index> creator(N,M,false,1);
    return NumaBCRSMatrix<Dune::FieldMatrix<Scalar,n,m>,Index>(creator);
  }
  
  
  // ----------------------------------------------------------------------------------------------

  /**
   * \ingroup linalgbasic
   * \relates NumaBCRSMatrix
   * \brief Creates the transposed sparse matrix \f$ A^T \f$.
   */
  template <class Entry, class Index>
  auto transpose(NumaBCRSMatrix<Entry,Index> const& A)
  {
    NumaBCRSMatrix<typename EntryTraits<Entry>::transpose_type,Index> At(A,false,true,false);
    return At;
  }

  // ----------------------------------------------------------------------------------------------

  namespace ThreadedMatrixDetail
  {
    template <class Target>
    class Submatrix;
    
    template <class Entry, class Allocator>
    class Submatrix<Dune::BCRSMatrix<Entry,Allocator>>
    {
    public:
      template <class Source, class RowIndices, class ColIndices>
      static Dune::BCRSMatrix<Entry,Allocator> apply(Source const& A,
                                                     RowIndices const& ri, ColIndices const& ci)
      {
        using Target = Dune::BCRSMatrix<Entry,Allocator>;
        using Index = typename Source::size_type;

        // Extract nonzero entries in the requested row and column ranges and store them
        // in sorted triplet format.
        std::vector<std::tuple<Index,Index>> a;
        std::vector<Entry> av;
        for (Index r=0; r<ri.size(); ++r)
        {
          auto row = A[ri[r]];
          for (Index c=0; c<ci.size(); ++c)
          {
            auto coli = row.find(ci[c]);
            if (coli != row.end())
            {
              a.push_back(std::make_tuple(r,c));
              av.push_back(*coli);
            }
          }
        }

        // Define the sparsity pattern.
        Target S(ri.size(),ci.size(),a.size(),Target::row_wise);
        using CreateIter = typename Target::CreateIterator;
        auto ai = begin(a);
        for(CreateIter row=S.createbegin(); row!=S.createend(); ++row)
          for ( ; ai!=end(a) && std::get<0>(*ai)==row.index(); ++ai)
            row.insert(std::get<1>(*ai));

        // Enter the matrix entries
        auto avi = begin(av);
        for (auto r=S.begin(); r!=S.end(); ++r)
          for (auto c=r->begin(); c!=r->end(); ++c, ++avi)
            *c = *avi;

        return S;
      }
    };
    
    
    template <class Entry, class Index>
    class Submatrix<NumaBCRSMatrix<Entry,Index>>
    {
    public:
      template <class Source, class RowIndices, class ColIndices>
      static NumaBCRSMatrix<Entry,Index> apply(Source const& A, 
                                               RowIndices const& ri, ColIndices const& ci)
      {
        NumaCRSPatternCreator<Index> creator(ri.size(), ci.size());

        // define sparsity pattern
        // WARNING: This has O(n*m) run time if the target matrix is n x m.
        // TODO: Devise a more sophisticated algorithm with O(nnz+n) run time.
        for(Index r=0; r<ri.size(); ++r)                                // iterate over rows in ri
        {
          auto row = A[ri[r]];
          for (Index c=0; c<ci.size(); ++c)                             // iterate over columns in ci
            if (auto coli = row.find(ci[c]); coli != row.end())
              creator.addElement(r,c);
        }

        NumaBCRSMatrix<Entry,Index> S(creator);

        // Enter the matrix entries
        for (auto r=S.begin(); r!=S.end(); ++r)
          for (auto c=r->begin(); c!=r->end(); ++c)
            *c = A[ri[r.index()]][ci[c.index()]];

        return S;
      }
    };
  }
  
  /**
   * \ingroup linalgbasic
   * \brief extracts a sparse submatrix
   *
   * This works as Matlab's A(ri,ci) submatrix extraction..
   *
   * \tparam Target the target matrix format , usually Dune::BCRSMatrix.
   * \tparam Source the source matrix, usually Dune::BCRSMatrix or NumaBCRSMatrix.
   * \tparam RowIndices an STL sequence type
   * \tparam ColIndices an STL sequence type
   *
   * # Complexity
   * For target Dune::BCRSMatrix, the complexity is \f$ n_r n_c \log n_s \f$, where
   * \f$ n_r \f$ is the size of ri, \f$ n_c \f$ the size of ci, and \f$ n_s \f$ the
   * number of nonzeros per row in A.
   *
   * \relates NumaBCRSMatrix
   */
  template <class Target, class Source, class RowIndices, class ColIndices>
  Target submatrix(Source const& A, RowIndices const& ri, ColIndices const& ci)
  {
    return ThreadedMatrixDetail::Submatrix<Target>::apply(A,ri,ci);
  }
  
  /**
   * \ingroup linalgbasic
   * \brief "deletes" rows/columns by extracting a copy of the matrix keeping the non-deleted rows/columns.
   *
   *
   * \tparam Target the target matrix format , usually Dune::BCRSMatrix.
   * \tparam Source the source matrix, usually Dune::BCRSMatrix or NumaBCRSMatrix.
   * \tparam RowIndices an STL sequence type
   * \tparam ColIndices an STL sequence type
   *
   *
   * \relates NumaBCRSMatrix
   */
  template <class Target, class Source, class RowIndices, class ColIndices>
  Target eraseRowsNCols(Source const& A, RowIndices const& ri, ColIndices const& ci)
  {
    // rows and columns to keep
    std::vector<size_t> rows(A.N()), cols(A.M());
    std::iota(rows.begin(),rows.end(),0);           // vector of rows to keep           
    std::iota(cols.begin(),cols.end(),0);           // vector of columns to keep    

    if(!ri.empty())
      for(int i=0; i<ri.size();++i)
      {
        rows.erase(rows.begin + ri[i]);               // erase the corresponding row indices 
      } 

    if(!ci.empty())
      for(int i=0; i<ci.size();++i)
      {
        cols.erase(cols.begin + ci[i]);               // erase the corresponding columns indices
      }

    return ThreadedMatrixDetail::Submatrix<Target>::apply(A,rows,cols);
  }
  
  /**
   * \ingroup linalgbasic
   * \brief "deletes" rows by extracting a copy of the matrix keeping only the non-deleted rows.
   *
   *
   * \tparam Target the target matrix format , usually Dune::BCRSMatrix.
   * \tparam Source the source matrix, usually Dune::BCRSMatrix or NumaBCRSMatrix.
   * \tparam RowIndices an STL sequence type
   *
   *
   * \relates NumaBCRSMatrix
   */
  template <class Target, class Source, class RowIndices>
  Target eraseRows(Source const& A, RowIndices const& ri)
  {
    // empty column vector (no columns to delete)
    std::vector<size_t> cols;
    return eraseRowsNCols(A,ri,cols);
  }
  
  /**
   * \ingroup linalgbasic
   * \brief "deletes" columns by extracting a copy of the matrix keeping only the non-deleted columns.
   *
   *
   * \tparam Target the target matrix format , usually Dune::BCRSMatrix.
   * \tparam Source the source matrix, usually Dune::BCRSMatrix or NumaBCRSMatrix.
   * \tparam ColIndices an STL sequence type
   *
   *
   * \relates NumaBCRSMatrix
   */
  template <class Target, class Source, class ColIndices>
  Target eraseCols(Source const& A, ColIndices const& ci)
  {
    // empty row vector (no rows to delete)
    std::vector<size_t> rows;
    return eraseRowsNCols(A,rows,ci);
  }
  
  
  /**
   * \ingroup linalgbasic
   * \brief returns all indices of nonzero columns.
   *
   *
   * \tparam Source the source matrix, usually Dune::BCRSMatrix or NumaBCRSMatrix.
   *
   *
   * \relates NumaBCRSMatrix
   */
  template <class Source>
  std::vector<size_t> nonZeroColumns(Source const& A)
  {
    // flag vector:
    // 1 = non zero column
    // 0 = zero column
    std::vector<size_t> nzFlags(A.M(),0);
    for(auto r=A.begin(); r!=A.end(); ++r)
      for(auto c=r->begin(); c!= r->end(); ++c)
        nzFlags.at(c.index()) = 1;
    
    // build column index vector
    std::vector<size_t> nzCols;
    for(int i=0; i<nzFlags.size(); ++i)
      if(nzFlags.at(i) == 1) nzCols.push_back(i);
    
    return nzCols;
  }
  
  
    /**
   * \ingroup linalgbasic
   * \brief   reshapes NumaBRCSMAtrix entry block structure
   * 
   * Converts a NumaBCRSMAtrix with entry-type Dune::FieldMatrix<Scalar,row1,col1> to entry type
   * Dune::FieldMatrix<Scalar,row2,col2>, i.e. the size of the contained blocks are reshaped.
   * Row size row1 must be divisible by row2 and column size col1 must be divisible by col2.
   * 
   * \tparam row1 number of rows in input block
   * \tparam col1 number of columns in input block
   * \tparam row2 desired number of block-rows
   * \tparam col2 desired number of column-rows
   * 
   * \param A NumaBRCSMatrix whose entry-type is to be reshaped
   * 
   * \relates NumaBCRSMatrix
   */
  template<int row2, int col2, int row1, int col1, class Scalar, class Index>
  NumaBCRSMatrix<Dune::FieldMatrix<Scalar,row2,col2>,Index>
  reshapeBlocks(NumaBCRSMatrix<Dune::FieldMatrix<Scalar,row1,col1>,Index> const& A)
  {
    static_assert(row1 % row2 == 0);
    static_assert(col1 % col2 == 0);
    int const qN = row1/row2;
    int const qM = col1/col2;
    
    Index N = A.N()*qN;
    Index M = A.M()*qM;
    NumaCRSPatternCreator<Index> creator(N, M, A.getPattern()->isSymmetric());
    
    for(auto row=A.begin(); row!=A.end(); ++row)
      for(auto col=row->begin(); col!=row->end(); ++col)                                        //loop over non-zero elements in A
        for(int i=0; i<qN; ++i)
          for(int j=0; j<qM; ++j)
            creator.addElement(row.index()*qN + i, col.index()*qM + j);                         //create qN*qM non-zero entries in B
            
    
    NumaBCRSMatrix<Dune::FieldMatrix<Scalar,row2,col2>, Index> B(creator);                      //empty matrix with sparsity pattern
    
    for(auto row=A.begin(); row!=A.end(); ++row)
      for(auto col=row->begin(); col!=row->end(); ++col)                                        // loop over non-zero elements in A
      {
        Dune::FieldMatrix<Scalar,row1,col1> mat = *col;                                         // get block entry of A
        for(int k=0; k<row1; ++k)
          for(int l=0; l<col1; ++l)
            B[row.index()*qN + k/row2][col.index()*qM + l/col2][k%row2][l%col2] = mat[k][l];    // add data at corresponding entries of B
      }

    return B;
  }
  

  
  /**
   * \ingroup linalgbasic
   * \brief   concatenate two matrices horizontally
   * 
   * Requires the matrices to have same number of rows and entry types
   * of Dune::FieldMatrix<Scalar,blockrows,blockcols>
   * 
   * \param A NumaBRCSMatrix, first submatrix (left)
   * \param B NumaBRCSMatrix, second matrix (right)
   * 
   * \related NumaBCRSMatrix
   */
  template<int blockrows, int blockcols, class Scalar, class Index>
  NumaBCRSMatrix<Dune::FieldMatrix<Scalar,blockrows,blockcols>, Index>
  horzcat(NumaBCRSMatrix<Dune::FieldMatrix<Scalar,blockrows,blockcols>, Index> const& A,
          NumaBCRSMatrix<Dune::FieldMatrix<Scalar,blockrows,blockcols>, Index> const& B)
  {
    assert(A.N() == B.N());   
    Index cols = (A.M() + B.M());
    Index rows = A.N();
    
    NumaCRSPatternCreator<Index> creator(rows, cols, false);
    
    //loop over non-zero elements in A        
    for(auto row=A.begin(); row!=A.end(); ++row)
      for(auto col=row->begin(); col!=row->end(); ++col)        
        creator.addElement(row.index(), col.index());
    //loop over non-zero elements in B        
    for(auto row=B.begin(); row!=B.end(); ++row)
      for(auto col=row->begin(); col!=row->end(); ++col)     
        creator.addElement(row.index(), col.index()+A.M()); // append columns on the right
            
    
    NumaBCRSMatrix<Dune::FieldMatrix<Scalar,blockrows,blockcols>,Index> result(creator);            //empty matrix with sparsity pattern
    
    for(auto row=A.begin(); row!=A.end(); ++row)
      for(auto col=row->begin(); col!=row->end(); ++col)                                            // loop over non-zero elements in A
        result[row.index()][col.index()] = *col;
      
    for(auto row=B.begin(); row!=B.end(); ++row)
      for(auto col=row->begin(); col!=row->end(); ++col)   
        result[row.index()][col.index()+A.M()] = *col;

          
    return result;
  }
  
  /**
   * \ingroup linalgbasic
   * \brief   concatenate two matrices vertically
   * 
   * Requires the matrices to have same number of rows and entry types
   * of Dune::FieldMatrix<Scalar,blockrows,blockcols>
   * 
   * \param A NumaBRCSMatrix, first submatrix (top)
   * \param B NumaBRCSMatrix, second matrix (bottom)
   * 
   */
  template<int blockrows, int blockcols, class Scalar, class Index>
  NumaBCRSMatrix<Dune::FieldMatrix<Scalar,blockrows,blockcols>, Index>
  vertcat(NumaBCRSMatrix<Dune::FieldMatrix<Scalar,blockrows,blockcols>, Index> const& A,
          NumaBCRSMatrix<Dune::FieldMatrix<Scalar,blockrows,blockcols>, Index> const& B)
  {
    assert(A.M() == B.M());   
    Index cols = A.M();
    Index rows = (A.N()+B.N());
    NumaCRSPatternCreator<Index> creator(rows, cols, false);
    
    //loop over non-zero elements in A        
    for(auto row=A.begin(); row!=A.end(); ++row)
      for(auto col=row->begin(); col!=row->end(); ++col)                                        
          creator.addElement(row.index(), col.index());
    //loop over non-zero elements in B        
    for(auto row=B.begin(); row!=B.end(); ++row)
      for(auto col=row->begin(); col!=row->end(); ++col) 
            creator.addElement(row.index()+A.N(), col.index()); // append columns on the right
            
    
    NumaBCRSMatrix<Dune::FieldMatrix<Scalar,blockrows,blockcols>, Index> result(creator);                      //empty matrix with sparsity pattern
    
    for(auto row=A.begin(); row!=A.end(); ++row)
      for(auto col=row->begin(); col!=row->end(); ++col)                                            // loop over non-zero elements in A
        result[row.index()][col.index()] = *col;

    for(auto row=B.begin(); row!=B.end(); ++row)
      for(auto col=row->begin(); col!=row->end(); ++col)                                            // loop over non-zero elements in B
        result[row.index()+A.N()][col.index()] = *col;
          
    return result;
  }
  
  
  /**
   * \ingroup linalgbasic
   * \brief   concatenate two matrices diagonally,
   *          resulting matrix is zero on the off block-diagonal
   * 
   * Requires the matrices to have same number of entry types
   * of Dune::FieldMatrix<double,blockrows,blockcols>
   * 
   * \param A NumaBRCSMatrix, first submatrix (top left)
   * \param B NumaBRCSMatrix, second matrix (bottom right)
   * 
   */
  template<int blockrows, int blockcols, class Index=size_t>
  NumaBCRSMatrix<Dune::FieldMatrix<double,blockrows,blockcols>, Index> diagcat(NumaBCRSMatrix<Dune::FieldMatrix<double,blockrows,blockcols>, Index> const& A,
          NumaBCRSMatrix<Dune::FieldMatrix<double,blockrows,blockcols>, Index> const& B)
  {   
    Index cols = (A.M()+B.M());
    Index rows = (A.N()+B.N());
    NumaCRSPatternCreator<Index> creator(rows, cols, false);
    
    //loop over non-zero elements in A        
    for(auto row=A.begin(); row!=A.end(); ++row)
      for(auto col=row->begin(); col!=row->end(); ++col)                                        
            creator.addElement(row.index(), col.index());
    //loop over non-zero elements in B        
    for(auto row=B.begin(); row!=B.end(); ++row)
      for(auto col=row->begin(); col!=row->end(); ++col) 
            creator.addElement(row.index()+A.N(), col.index()+A.M()); // append columns on the right
            
    
    NumaBCRSMatrix<Dune::FieldMatrix<double,blockrows,blockcols>, Index> result(creator);                      //empty matrix with sparsity pattern
    
    for(auto row=A.begin(); row!=A.end(); ++row)
      for(auto col=row->begin(); col!=row->end(); ++col)                                            // loop over non-zero elements in A
      {
        auto mat = *col;
        for(int k=0; k < blockrows; ++k)
          for(int l= 0; l < blockcols; ++l)
            result[row.index()][col.index()][k][l] =mat[k][l];
      }

    for(auto row=B.begin(); row!=B.end(); ++row)
      for(auto col=row->begin(); col!=row->end(); ++col)                                            // loop over non-zero elements in B
      {
        auto mat = *col;
        for(int k=0; k < blockrows; ++k)
          for(int l= 0; l < blockcols; ++l)
            result[row.index()+A.N()][col.index()+A.M()][k][l] = mat[k][l];
      }

          
    return result;
  }

}   // end namespace Kaskade

#endif
