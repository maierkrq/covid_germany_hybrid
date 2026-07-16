/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2012-2019 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <new>
#include <numeric>  
#include <iostream>
#include <signal.h>
#include <sys/wait.h>
#include <unistd.h>

#ifdef HAVE_NUMA
#include "numa.h"
#endif

#include "utilities/detailed_exception.hh"
#include "utilities/threading.hh"

namespace Kaskade
{
  std::mutex DuneQuadratureRulesMutex;
  boost::mutex refElementMutex;
  
  //----------------------------------------------------------------------------
  
  void equalWeightRanges(std::vector<size_t>& w, size_t nRequest)
  {
    assert(nRequest > 0);
    size_t const totalWeight = std::accumulate(w.begin(),w.end(),0);
    size_t const N = w.size();
    size_t n = std::min(nRequest,N);
    
    // Initialize with a partitioning such that the first partition
    // covers the whole range.
    std::vector<size_t> x(nRequest+1,N);
    x[0] = 0;
    
    // Now step through the range and set (i.e. lower) the separation points
    // such that the weight of all partitions is roughly equal to N/nRequest.
    // This is achieved by setting x_k to x_{k-1}+1 and increasing further
    // until the covered weight reaches k*totalWeight/nRequest.
    size_t weight = 0; // count the weight covered so far
    int k = 1;
    for ( ; weight<totalWeight; ++k)
    {
      assert(k<=nRequest);

      x[k] = x[k-1]+1;           // at least one entry per partition
      weight += w[x[k-1]];       // and its weight is counted
      
      // Add more entries to the current partition until we're at the desired weight.
      // Zero weights are consumed immediately.
      // Note that if we run out of partitions (k==nRequest), the last one consumes
      // all the remaining weight (i.e. until weight=k*totalWeight/nRequest=totalWeight).
      // thus, the for-loop body is never executed with k>nRequest.
      while (weight < (k*totalWeight)/nRequest || (x[k] < N && w[x[k]]==0))
      {
        weight += w[x[k]];
        ++x[k];
        assert(x[k]<=N);
      }

      // Check that if we have consumed the whole range (i.e. x[k]==N),
      // then we also have consumed the whole weight (weight==totalWeight).
      assert(x[k]<N || weight==totalWeight);
    }

    // It may happen that the loop above stops early (or is never executed at all).
    // In that case, the last partition considered reaches x[k]=N, and the remaining
    // ones are empty (x[j]=N=x[k] for j>k) due to their initialization.

    // return the partitions in x
    std::swap(w,x);
  }
  
 
  //----------------------------------------------------------------------------

  namespace { 
   

    //----------------------------------------------------------------------------
    // Code to be executed by worker threads
    class Worker 
    {
    public:
      // Create thread pinned on given node (unless node is negative, then no pinning is done)
      Worker(int node_, ConcurrentQueue<Task>& tasks_)
      : node(node_), tasks(tasks_)
      {
        tasks.running(+1);
      }
       
      void operator()()
      {
	
#ifdef HAVE_NUMA
        if (numa_available()>=0 && node>=0)
        {
          numa_run_on_node(node);   // bind thread to the NUMA node, leaving choice of CPU to the operating system
          numa_set_preferred(node); // bind memory allocation from this thread to the memory of this node
        }
#endif
	
        Task t = tasks.pop_front();  //  get commission (blocks if there is no task to be done)
        
        while (true)                 // run forever, i.e. until explicit interruption
        {
          t();                       // do the work
          t = tasks.pop_front();     // get new commission (blocks if there is no task to be done)
        }
      }
      
    private:
      int node;                      // the node this worker thread runs on
      ConcurrentQueue<Task>& tasks;  // the queue from which tasks are extracted
    };
    
  } // end of anonymous namespace

  
  //----------------------------------------------------------------------------
  
  namespace ThreadingDetail 
  {

    NumaAllocatorBase::NumaAllocatorBase(int node_): nod(node_) 
    {
#ifdef HAVE_NUMA
      if (nod >= 0)
        allocator = &NumaThreadPool::instance().allocator(nod);
      else
        allocator = nullptr;
#endif
    }
    
    size_t NumaAllocatorBase::max_size() const
    {
#ifdef HAVE_NUMA      
      long free;
      if (nod>=0)
        return numa_node_size(nod,&free);
      else
      {
        std::allocator<char> a;
        return a.max_size();
      }
#else
      std::allocator<char> a;
      return a.max_size();
#endif 
    }
    
    void* NumaAllocatorBase::allocate(size_t n) 
    {
#ifdef HAVE_NUMA    
      void* mem = nod>=0? allocator->alloc(n): numa_alloc_interleaved(n);
      if (mem==nullptr) 
      {
        std::cerr << "cannot serve memory request (size=" << n << " B) on node " << nod << ".\n";
        throw std::bad_alloc();
      }
#else
      void* mem = std::malloc(n);
      if (mem==nullptr) 
      {
        std::cerr << "cannot serve memory request (size=" << n << " B) via malloc()\n";
        throw std::bad_alloc();
//         abort();
//         throw DetailedException("cannot serve memory request",__FILE__,__LINE__);
      }
#endif 

      return mem;
    }
    
    void NumaAllocatorBase::deallocate(void* p, size_t n)
    {

#ifdef HAVE_NUMA      
      if (nod>=0)
        allocator->free(p,n);
      else   
        numa_free(p,n);
#else
      std::free(p);
#endif 
    }
  
  }
  
  //----------------------------------------------------------------------------
    

  NumaThreadPool& NumaThreadPool::instance(int maxThreads)
  {
    static NumaThreadPool pool(maxThreads);
    return pool;
  }
  

  Ticket NumaThreadPool::run(Task&& task)
  {
    return runOnQueue(globalQueue,std::move(task));
  }
  
  Ticket NumaThreadPool::runOnNode(int node, Task&& task)
  {
    assert(0<=node && node<nodes());
    return runOnQueue(nodeQueue[node],std::move(task));
  }
  
  Ticket NumaThreadPool::runOnQueue(ConcurrentQueue<Task>& queue, Task&& task)
  {
    Ticket tick = task.get_future();
    if (sequential)
      task();                               // run on main thread
    else
      queue.push_back(std::move(task));     // submit it to the worker thread
    return tick;
  }
      

  NumaThreadPool::NumaThreadPool(int maxThreads)
  : sequential(maxThreads==1)
  {
    nCpu = std::max(1u,std::thread::hardware_concurrency());
    
#if defined(HAVE_NUMA) 
    // initialize libnuma
    if (numa_available()<0 || sequential) // no numa available
    {
      nNode = 1;
      
      // associate cpus and nodes trivially
      nodeByCpu.resize(nCpu,0);
      
      cpuByNode.resize(nNode);
      for (int c=0; c<nCpu; ++c)
        cpuByNode[0].push_back(c);
    }
    else // we have numa - look out for number of nodes and cpus association to nodes
    {
      nNode = numa_max_node()+1;
      nCpu = numa_num_configured_cpus(); 
      
      nodeByCpu.resize(nCpu);
      cpuByNode.resize(nNode);

      // for each node, obtain the cpus in this node
      for (int c=0; c<nCpu; ++c)
      {
        int n = numa_node_of_cpu(c);
        nodeByCpu[c] = n;
        cpuByNode[n].push_back(c);
      }
      
    }
#else // we don't know nothing about NUMA - probably we're on a single-socket-multi-core system or pureley sequential
      nNode = 1;
      
      // associate cpus and nodes trivially: only one node, containing all the cpus
      nodeByCpu.resize(nCpu,0);
      
      cpuByNode.resize(nNode);
      for (int c=0; c<nCpu; ++c)
        cpuByNode[0].push_back(c);
#endif
    
    nodeMemory.reserve(nNode);
    for (int i=0; i<nNode; ++i)
      nodeMemory.push_back(Kalloc(i,64,4*1024*1024,false));

    nodeQueue.resize(nNode);
    
    // On each node create one thread per cpu listening for the node task queue.
    for (int n=0; n<nNode; ++n)
      for (int i=0; i<cpuByNode[n].size(); ++i)
        threads.create_thread(Worker(n,nodeQueue[n]));
    
    // Create as many unpinned global threads as allowed. 
    for (int i=0; i<std::min(maxThreads,std::max(4,2*nCpu)); ++i)
      threads.create_thread(Worker(-1,globalQueue)); 
      
    // compute maximal number of CPUs on any node
    maxCpusPerNode = 0;
    for (auto const& cpus: cpuByNode)
      maxCpusPerNode = std::max(maxCpusPerNode,static_cast<int>(cpus.size()));
  }  
  
  Kalloc& NumaThreadPool::allocator(int node)
  {
    assert(0<=node && node<nNode);
    return nodeMemory[node];
  }
  
  void* NumaThreadPool::allocate(size_t len, int node)
  {
    assert(0<=node && node<nNode);
    void* mem = nodeMemory[node].alloc(len);
    
    if (mem==nullptr)
      throw std::bad_alloc();
    
    return mem;
  }
  
  void NumaThreadPool::deallocate(void* mem, size_t n, int node)
  {
    // look for memory block to get its length (required by numa)
    nodeMemory[node].free(mem,n);
  }
  
  size_t NumaThreadPool::alignment(int node) const
  {
    assert(0<=node && node<nodeMemory.size());
    return nodeMemory[node].alignment();
  }
  
  void NumaThreadPool::reserve(size_t n, size_t k, int node)
  {
    assert(0<=node && node<nodeMemory.size());
    nodeMemory[node].reserve(n,k);
  }

  NumaThreadPool::~NumaThreadPool()
  {
    // tell all worker threads to stop
    threads.interrupt_all();
    // wait for worker threads to terminate
    threads.join_all();
    
    // clean up the memory
    for (auto mem: memBlocks)
      nodeMemory[mem.second.second].free(mem.first,mem.second.first);
  }

  //----------------------------------------------------------------------------
  
  
  void runInBackground(std::function<void()>& f)
  {
    
    // use double fork to avoid zombie processes
    pid_t pid;
    if ((pid=fork()) == 0) // we're child
    {
      if ((pid=fork()) == 0)    // we're grandchild
        f();                    // do the work 

      kill(getpid(),SIGTERM);   // For some unknown reason, the child and grandchild do not terminate
                                // via exit() alone - hence we send an explicit SIGTERM. Funny is, this helps...
      exit(0);                  // terminate immediately if child, after work if grandchild
    }
    else if (pid>0)             // we're parent 
    {
      int status;               // wait for child to terminate (which it does immediately,
      waitpid(pid,&status,0);   // orphaning the grand child)
    }
    else
      throw DetailedException("fork() failed.",__FILE__,__LINE__);
  }

  
}
