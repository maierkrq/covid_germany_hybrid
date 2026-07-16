/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2017-2017 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef MATRIX_BLOCKS_HH
#define MATRIX_BLOCKS_HH

#include <vector>
#include <cassert>
#include <sstream>
#include <algorithm>

#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>

#include <linalg/threadedMatrix.hh>
#include <utilities/detailed_exception.hh>

namespace Kaskade
{
/**
 * @brief The GlobalLocalIndexRelationship class provides relationships between local indices of (e.g. matrix) blocks and global indices (of the corresponding block matrix)
 */
template <class Index>
class GlobalLocalIndexRelationship
{
public:
  GlobalLocalIndexRelationship() {}

  /**
   * @brief GlobalLocalIndexRelationship
   * @param localSizes are sizes of involved blocks
   */
  GlobalLocalIndexRelationship(std::vector<Index> const& localSizes) : offset_(localSizes.size()+1){
    offset_[0] = 0;
    for (Index i=0; i<localSizes.size(); ++i) {
      assert(localSizes[i]>=0);
      offset_[i+1] = offset_[i] + localSizes[i];
    }
  }

  /**
   * @brief offset
   * @param block
   * @return global start index of some block
   */
  Index offset(Index block) const {
    assert(block>=0 && block<numberBlocks());
    return offset_[block];
  }

  /**
   * @brief globalIndex
   * @param block
   * @param localIndex
   * @return global index by block number and block local index
   */
  Index globalIndex(Index block, Index localIndex) const {
    return offset(block) + localIndex;
  }

  /**
   * @brief localIndex
   * @param globalIndex
   * @return pair containing block number and local index by global index
   */
  std::pair<Index, Index> localIndex(Index globalIndex) const {
    std::pair<Index, Index> res = findOffset(globalIndex);
    res.second = globalIndex - res.second;
    return res;
  }

  /**
   * @brief findOffset
   * @param globalIndex
   * @return pair containing block number and global start index of block by global index
   */
  std::pair<Index, Index> findOffset(Index globalIndex) const {
    assert(globalIndex>=0);
    Index block = 0;
    for(; block<numberBlocks(); ++block) {
      if(offset_[block+1] > globalIndex) break;
      assert(block < numberBlocks()-1);
    }
    return std::make_pair(block, offset(block));
  }

  /**
   * @brief globalSize
   * @return number of global indices
   */
  Index globalSize() const {
    return offset_.back();
  }

  /**
   * @brief localSize
   * @param block
   * @return size of some block
   */
  Index localSize(Index block) const {
    return -offset(block) + offset_[block+1]; // this way assert in offset() is reused
  }

  /**
   * @brief numberBlocks
   * @return number of involved blocks
   */
  Index numberBlocks() const {
    return offset_.size()-1;
  }

private:
  // contains the global start indices of the blocks, last entry is number of global indices
  std::vector<Index> offset_ {0};
};

/**
 * Add indices of a sparse matrix block to a global sparse matrix (creator).
 */
template <class Entry, class Index>
void insertMatrixBlockIndices(NumaCRSPatternCreator<Index> & globalCreator, NumaBCRSMatrix<Entry, Index> const& matrixBlock,
                              typename NumaBCRSMatrix<Entry, Index>::size_type startRow, typename NumaBCRSMatrix<Entry, Index>::size_type startCol) {
//  if(globalMatrix.N() < startRow+matrixBlock.N()) {
//    std::ostringstream message;
//    message << "Block (" << matrixBlock.N() << " rows and starting at row " << startRow << ") does not fit in global matrix (" << globalMatrix.N() << " rows). ";
//    throw LinearAlgebraException(message.str(),__FILE__,__LINE__);
//  }
  if(globalCreator.cols() < startCol+matrixBlock.M()) {
    std::ostringstream message;
    message << "Block (" << matrixBlock.M() << " columns and starting at column " << startCol << ") does not fit in global matrix (creator) (" << globalCreator.cols() << " columns). ";
    throw LinearAlgebraException(message.str(),__FILE__,__LINE__);
  }
  if(matrixBlock.getPattern()->isSymmetric() && (startCol != startRow)) {
    throw LinearAlgebraException("Symmetric stored blocks can only be inserted on diagonal.",__FILE__,__LINE__);
  }
  if(matrixBlock.getPattern()->isSymmetric() && !globalCreator.isSymmetric()) {
    throw LinearAlgebraException("Symmetric stored blocks can only be inserted into symmetric stored Matrices.",__FILE__,__LINE__);
  }
  for(Index i=0; i<matrixBlock.N(); ++i) {
    auto const& rowMb = matrixBlock[i];
    std::vector<Index> colIdxs(rowMb.size());
    std::transform(rowMb.getindexptr(), rowMb.getindexptr()+rowMb.size(), colIdxs.begin(), [&](Index const& colIdx) {
      return startCol+colIdx;
    });
    Index rowIdx = startRow+i;
    globalCreator.addElements(&rowIdx, (&rowIdx)+1, colIdxs.begin(), colIdxs.end());
  }
}

/**
 * Add indices of a sparse matrix block to a global sparse matrix (creator).
 */
template <class Entry>
void insertMatrixBlockIndices(Dune::MatrixIndexSet & globalIndexSet, Dune::BCRSMatrix<Entry> const& matrixBlock, Dune::MatrixIndexSet::size_type startRow, Dune::MatrixIndexSet::size_type startCol) {
  if(globalIndexSet.rows() < startRow+matrixBlock.N()) {
    std::ostringstream message;
    message << "Block (" << matrixBlock.N() << " rows and starting at row " << startRow << ") does not fit in global matrix (index set)(" << globalIndexSet.rows() << " rows). ";
    throw LinearAlgebraException(message.str(),__FILE__,__LINE__);
  }
//  if(globalIndexSet.cols() < startCol+matrixBlock.M()) {
//    std::ostringstream message;
//    message << "Block (" << matrixBlock.M() << " columns and starting at column " << startCol << ") does not fit in global matrix (creator) (" << globalCreator.cols() << " columns). ";
//    throw LinearAlgebraException(message.str(),__FILE__,__LINE__);
//  }
  globalIndexSet.import(matrixBlock, startRow, startCol);
}

/**
 * Add entries of a sparse matrix block to a global sparse matrix, which already has the suitable sparsity pattern assigned.
 */
template <class BMatrix>
void insertMatrixBlockEntries(BMatrix & globalMatrix, BMatrix const& matrixBlock, typename BMatrix::size_type startRow, typename BMatrix::size_type startCol) {
  if(globalMatrix.N() < startRow+matrixBlock.N()) {
    std::ostringstream message;
    message << "Block (" << matrixBlock.N() << " rows and starting at row " << startRow << ") does not fit in global matrix (" << globalMatrix.N() << " rows). ";
    throw LinearAlgebraException(message.str(),__FILE__,__LINE__);
  }
  if(globalMatrix.M() < startCol+matrixBlock.M()) {
    std::ostringstream message;
    message << "Block (" << matrixBlock.M() << " columns and starting at column " << startCol << ") does not fit in global matrix (" << globalMatrix.M() << " columns). ";
    throw LinearAlgebraException(message.str(),__FILE__,__LINE__);
  }
  for(typename BMatrix::size_type i=0; i<matrixBlock.N(); ++i) {
    auto const& rowMb = matrixBlock[i];
    if(rowMb.size() > 0) std::copy(rowMb.getptr(), rowMb.getptr()+rowMb.size(), &(globalMatrix[startRow+i][startCol+rowMb.begin().index()]));
  }
}

/**
 * Add entries of a vector block to a global vector.
 */
template <class BVector>
void insertVectorBlock(BVector & globalVector, BVector const& vectorBlock, typename BVector::size_type startIdx) {
  if(globalVector.N() < startIdx+vectorBlock.N()) {
    std::ostringstream message;
    message << "Block (of size " << vectorBlock.N() << " and starting at index " << startIdx << ") does not fit in global vector (of size " << globalVector.N() << "). ";
    throw LinearAlgebraException(message.str(),__FILE__,__LINE__);
  }
  for(typename BVector::size_type i=0; i<vectorBlock.N(); ++i) {
    globalVector[startIdx+i] = vectorBlock[i];
  }
}
}

#endif
