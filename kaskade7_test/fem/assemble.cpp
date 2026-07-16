/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2012-2020 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


#include "dune/grid/config.h"

#include "fem/assemble.hh"

namespace Kaskade
{
  namespace AssemblyDetail {
    
    RowGroupManager::RowGroupManager()
    : rowGroupStart(1,0)
    , mutexes(nullptr) 
    {}
    
    RowGroupManager::~RowGroupManager() 
    {
      clear();
    }

    void RowGroupManager::clear()
    {
      if (mutexes)
        delete[] mutexes;
    }
    
    
    void RowGroupManager::init(int nrg, size_t rows)
    {
      // no more row groups than rows...
      nrg = std::min((size_t)nrg,rows);        
      rowGroupStart.resize(nrg+1);
      for (size_t i=0; i<nrg; ++i)
        rowGroupStart[i] = i*rows/nrg;
      rowGroupStart[nrg] = rows;
      
      clear();
      mutexes = new std::mutex[size()];
    }
    
    void RowGroupManager::init(int nrg, std::vector<size_t> const& rowSize)
    {
      size_t const rows = rowSize.size();
      size_t nnz = 0;
      for (size_t i=0; i<rows; ++i)
        nnz += rowSize[i];
      
      
      // no more row groups than rows...
      nrg = std::min((size_t)nrg,rows);        
      rowGroupStart.resize(nrg+1);
      
      // first and last entry are clear
      rowGroupStart[0] = 0;
      rowGroupStart[nrg] = rows;
      size_t count = 0; // count the number of entries in group i-1
      
      int i = 1;
      while (count<nnz)
      {
        assert(i<=nrg);
        
        rowGroupStart[i] = rowGroupStart[i-1]+1;    // at least one row per group
        count += rowSize[rowGroupStart[i-1]];       // and its elements are done

        // add more rows to the current group until we're at the desired position
        while (count < (i*nnz)/nrg)
        {
          count += rowSize[rowGroupStart[i]];
          ++rowGroupStart[i];
          assert(rowGroupStart[i]<=rows);
        }
        
        ++i;
      }
      
      // maybe we ended up with fewer row groups than intended.
      rowGroupStart.erase(rowGroupStart.begin()+i,rowGroupStart.end());
      
      clear();
      mutexes = new std::mutex[size()];
    }
    
    void RowGroupManager::init(int nrg, std::vector<std::vector<size_t> > const& colIndex)
    {
      std::vector<size_t> rowSize(colIndex.size());
      std::transform(begin(colIndex),end(colIndex),begin(rowSize),
                     [](auto const& cidx) { return cidx.size(); });
      init(nrg,rowSize);
    }
    
    std::pair<size_t,size_t> RowGroupManager::range(int n) const
    {
      assert(0<=n && n<size());
      return std::pair(rowGroupStart[n],rowGroupStart[n+1]);
    }
    
    std::mutex& RowGroupManager::mutex(int n) const
    {
      assert(0<=n && n<size());
      return mutexes[n];
    }
    
    int RowGroupManager::size() const 
    {
      return rowGroupStart.size()-1;
    }
  }
}
