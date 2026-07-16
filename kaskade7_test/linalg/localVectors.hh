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

#ifndef LOCALVECTORS_HH
#define LOCALVECTORS_HH

#include <vector>
#include <mutex>

namespace Kaskade
{
  
  /**
   * \ingroup linalgbasic
   * \brief Providing a vector interface for contiguously stored entries.
   *
   * \tparam Entry the type of vector entries
   * \tparam SortedIdx a range of indices (global idx, local idx) sorted ascendingly by global idx
   */
  template <class Entry, class SortedIdx>
  class LocalVector
  {
  public:

    LocalVector(SortedIdx const& idx, Entry* data_)
    : idx_(idx), data(data_)
    {
    }
    
    /**
     * \brief A sequence of indices, sorted ascendingly by global idx.
     */
    SortedIdx idx() const { return idx_; }

    /**
     * \brief The number of entries stored.
     */
    size_t size() const
    {
      return idx().size();
    }
    
    /**
     * \brief Resets the data pointer.
     */
    void relocate(Entry* newData)
    {
      data = newData;
    }
    
    /**
     * \brief Access the matrix entries
     */
    Entry& operator[](int idx)
    {
      return data[idx];
    }
    
    Entry const& operator[](int idx) const
    {
      return data[idx];
    }
    
  private:
    SortedIdx idx_;
    
    // a pointer to the raw storage for the entries
    Entry* data;
  };
  
  /**
   * \ingroup linalgbasic
   * \brief A structure for holding a sequence of several local vectors to be filled sequentially
   *        and to be scattered jointly.
   * 
   * This realizes a container of vector,
   * with the critical feature that the entries for all the vectors are stored in a contiguous memory block. This
   * should improve locality of memory access during assembly and scattering of local vectors.
   * 
   * \tparam Entry the type of matrix entries
   * \tparam SortedIdx a range of indices (global idx, local idx) sorted ascendingly by global idx
   * \tparam Vector the type of the global Vector
   */
  template <class Entry, class SortedIdx, class Vector>
  class LocalVectors
  {
  public:
    using value_type = LocalVector<Entry,SortedIdx>;
    
    /**
     * \param globalVector global vector where computed local vectors are to be scattered to.
     * \param scatterMutexes_ several mutexes (associated with ranges in the global vector), required to lock ranges in the global vector during scattering.
     *                        By the number of given mutexes it can be controlled into how many ranges the global vector is to be subdivided during scattering.
     * \param maxStorage the desired size of the buffer for data of local vectors, in bytes.
     */
    LocalVectors(Vector& globalVector, std::vector<std::mutex> & scatterMutexes_, size_t maxStorage=256*1024)
    : globalVector(&globalVector), scatterMutexes(&scatterMutexes_)
    {
      // Allocate memory as desired. The vector will only grow if there is a local matrix with more entries than defined here.
      localData.reserve(maxStorage/sizeof(Entry));
      if(scatterMutexes_.empty())
      {
        throw DetailedException("At least one mutex must be given.",__FILE__,__LINE__);
      }
    }
    
    /**
     * \brief Default constructor. 
     * 
     * This creates an essentially useless empty local matrix buffer. 
     */
    LocalVectors()
    {}
    
    /**
     * \brief Destructor.
     * 
     * The destructor scatters all remaining local matrices into the global matrix.
     */
    ~LocalVectors()
    {
      scatter();
    }

    /**
     * \brief Appends another (zero-initialized) local vector.
     * 
     * In case that appending the new local vector exceeds the maximum storage size, the existing local vectors
     * are scattered into the global vector before the new vector is created.
     * 
     * \param idx a range of (global idx, local idx) pairs sorted by global index
     */
    void push_back(SortedIdx const& idx)
    {
      size_t size = idx.size();
      
      // Appending the new vector exceeds the storage capacity. Scatter the existing vectors to their target and
      // remove them.
      if (localData.capacity() < localData.size()+size)
        scatter();
      
      // Append the new local vector. Note that either the capacity of localData is sufficient, and
      // the localVectors data pointers remain valid, or the local vectors have been scattered,
      // such that localVectors is empty and no pointer exists that could be invalidated.
      Entry* pos = &*localData.insert(localData.end(),size,Entry(0));
      localVectors.emplace_back(idx,pos);
    }

    /**
     * \brief A reference to the n-th local matrix.
     */
    value_type      & operator[](int n)       { return localVectors[n]; }
    value_type const& operator[](int n) const { return localVectors[n]; }
    
    /**
     * \brief A reference to the last pushed local matrix.
     */
    value_type      & back()       { return localVectors.back(); }
    value_type const& back() const { return localVectors.back(); }
    
    /**
     * \brief The number of local vectors that are stored.
     */
    size_t size() const { return localVectors.size(); }
    
    /**
     * \brief Scatters the local vectors into the global one and removes them from the local container.
     * 
     * This is generally not necessary, as the LocalVectors class scatters automatically as soon as its buffer
     * is filled, and also on destruction. Use this only if you need to control the specific time point at which
     * scattering occurs.
     */
    void scatter()
    {
      if (localData.empty()) return;

      if (!globalVector)
      {
        throw DetailedException("No global vector given to scatter to",__FILE__,__LINE__);
      }

      // step through all global vector ranges and scatter all the included local vectors' components
      unsigned n = scatterMutexes->size();
      for (unsigned ri=0; ri<n; ++ri)
      {
        // lock global vector range
        std::lock_guard<std::mutex> lock((*scatterMutexes)[ri]);

        // step through all local vectors
        for (auto const& vec : localVectors)
        {
          size_t globalStartIdx = uniformWeightRangeStart(ri,   n, globalVector->size());
          size_t globalEndIdx   = uniformWeightRangeStart(ri+1, n, globalVector->size());

          for (auto const& glIdx : vec.idx())
          {
            // check which rows are inside the global vector range
            if (glIdx.first >= globalEndIdx) break; // vec.idx() are sorted by global index
            if (globalStartIdx <= glIdx.first)
            {
              (*globalVector)[glIdx.first] += vec[glIdx.second];
            }
          }
        }
      }

      // everything's scattered - clean up
      clear();
    }
    
    /**
     * \brief clears all the data, leaving an empty state
     */
    void clear() 
    {
      localData.clear();
      localVectors.clear();
    }
    
    /**
     * \brief reports the size of the local vectors storage in bytes
     */
    size_t storageSize() const
    {
      return localData.size() * sizeof(Entry);
    }

  private:
    // This allows to hold a single array into which the local vectors are scattered, and
    // thus improves memory locality (for better cache hit rates and less false sharing).
    // Contains the local vectors one after the other.
    std::vector<Entry> localData;
    // Storage for local vectors.
    std::vector<value_type> localVectors;

    Vector* globalVector = nullptr;
    std::vector<std::mutex>* scatterMutexes = nullptr;
  };
  
}

#endif
