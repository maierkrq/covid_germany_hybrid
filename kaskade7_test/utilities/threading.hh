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

#ifndef THREADING_HH
#define THREADING_HH

#include <fstream>
#include <functional>
#include <future>
#include <map>
#include <mutex>
#include <queue>
#include <utility>
#include <vector>

#include <boost/thread/condition_variable.hpp>
#include <boost/thread/locks.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/thread.hpp>

#include "utilities/kalloc.hh"
#include "utilities/timing.hh"

namespace Kaskade 
{
  /**
   * \ingroup threading
   * \brief A global lock for the Dune::QuadratureRules factory, which is not thread-safe as of 2015-01-01.
   */
  extern std::mutex DuneQuadratureRulesMutex;
  
  /**
   * \ingroup threading
   * \brief A global lock for the Dune::GenericReferenceElement singletons, which are not thread-safe.
   */
  extern boost::mutex refElementMutex;
  
  //---------------------------------------------------------------------------
  
  /**
   * \ingroup threading
   * \brief Computes partitioning points such that the sum of weights in each partition is roughly the same.
   * 
   * Let \f$ w_i, 0\le i < N \f$ denote the weights x[i]. On exit, x has size \f$ k=n+1 \f$ with values \f$ x_i \f$
   * such that \f$ x_0 = 0 \f$, \f$ x_{k-1} = N\f$, and 
   * \f[ \sum_{j=x_i}^{x_{i+1}-1} w_j \approx \frac{1}{k} \sum_{j=0}^{N-1} w_j. \f]
   * 
   * \param[in,out] x the array of (nonnegative) weights
   * \param[in]     n the desired number of partitions (positive). 
   */
  void equalWeightRanges(std::vector<size_t>& x, size_t n);
  
  /**
   * \ingroup threading
   * \brief Computes partitioning points of ranges for uniform weight distributions.
   * 
   * The functionality is similar to \ref equalWeightRanges, but for uniform weight distribution, the 
   * partitioning points can be computed directly.
   * 
   * \tparam Index an integral type
   * \param i the number of the range for which the starting point is to be computed. Has to be in [0,n]
   * \param n the number of ranges
   * \param m the total number of entries
   */
  template <class BlockIndex, class Index>
  Index uniformWeightRangeStart(BlockIndex i, BlockIndex n, Index m)
  {
    assert(i>=0 && i<=n && n>0);
    return (i*m)/n;
  }
  
  /**
   * \ingroup threading
   * \brief Computes the range in which an index is to be found when partitioned for uniform weights.
   *
   * \tparam Index an integral type
   * \param j the index for which the range number is to be computed. Has to be in [0,m[
   * \param n the number of ranges
   * \param m the total number of entries
   */
  template <class Index>
  Index uniformWeightRange(Index j, Index n, Index m)
  {
    assert(j>=0 && j<m && n>0 && m>0);
    
    // We're looking for i such that floor(i*m/n) <= j < floor((i+1)*m/n), i.e. index j has
    // to be in the half-open range given by the partitioning points of the range. Now this 
    // implies i*m/n-1 <= j < (i+1)*m/n and also j*n/m-1 < i <= (j+1)*n/m. As i is integer, 
    // this implies floor(j*n/m) <= i <= floor((j+1)*n/m). Typically, n/m is small, in which 
    // case this closed interval is likely to contain just one natural number - the result.
    Index low = (j*n)/m;       // floor is implied by integer arithmetics
    Index high = ((j+1)*n)/m;  // floor is implied by integer arithmetics
    
    if (low==high)
    {
      assert(uniformWeightRangeStart(low,n,m)<=j && j<uniformWeightRangeStart(low+1,n,m));
      return low;
    }
    
    // The implied inequality chains above are not sharp estimates, therefore the interval
    // [low,high] can contain several natural numbers. This is always the case if n>m, i.e., we have 
    // more ranges than entries and several ranges are empty. But it can also happen if j*n/m is
    // just below the next integral number. Now we have to walk through all natural numbers in 
    // the interval to find the correct range. TODO: is there a direct way?
    Index i = low;
    while (i<high && !(uniformWeightRangeStart(i,n,m)<=j && j<uniformWeightRangeStart(i+1,n,m)))
      ++i;
    
    assert(uniformWeightRangeStart(i,n,m)<=j && j<uniformWeightRangeStart(i+1,n,m));
    return i;
  }
  
  //---------------------------------------------------------------------------

  /**
   * \ingroup threading
   * \brief A concurrent fifo queue.
   * 
   * The queue is not copyable as it is intended to hold non-copyable std:packaged_task objects.
   */
  template <class T>
  class ConcurrentQueue {
  public:
    /**
     * \brief Constructs an empty queue.
     */
    ConcurrentQueue()
    : runningWorkers(0)
    {
    }
    
    /**
     * \brief Moves a queue.
     * 
     * This transferres all tasks waiting in q to the newly constructed queue, but leaves q intact.
     */
    ConcurrentQueue(ConcurrentQueue&& q)
    {
      boost::lock_guard<boost::mutex> lock(q.mutex);
      queue = std::move(q.queue);
    }
    
    /**
     * \brief Assignment.
     * 
     * This copies all tasks waiting in q to the newly constructed queue, but leaves q intact.
     */
    ConcurrentQueue& operator=(ConcurrentQueue const& q)
    {
      boost::lock_guard<boost::mutex> lock(q.mutex);
      queue = q.queue;
      return *this;
    }
    
    bool empty() const
    {
      // we need a lock here because someone may add / remove tasks right now
      boost::lock_guard<boost::mutex> lock(mutex); 
      return queue.empty();
    }
    
    /**
     * \brief Returns the number of tasks waiting.
     */
    size_t size() const
    {
      // we need a lock here because someone may add / remove tasks right now
      boost::lock_guard<boost::mutex> lock(mutex);
      return queue.size();
    }
    
    /**
     * \brief Stores an element at the end of the queue.
     */
    void push_back(T&& t)
    {
      {
      boost::lock_guard<boost::mutex> lock(mutex);
      queue.push(std::move(t));
      } // release lock before waking up consumers
      filled.notify_one();
    }
    
    /**
     * \brief Retrieves the foremost element.
     * 
     * This method blocks if the queue is empty and waits for data to become available.
     */
    T pop_front()
    {
      boost::unique_lock<boost::mutex> lock(mutex);
      
      // Wait for data to become available. 
      while (queue.empty())
        filled.wait(lock);
      
      // extract and remove data. When move assignment throws an exception, the queue
      // remains unmodified
      T t = std::move(queue.front());
      queue.pop();
      return t;
    }
    
    /**
     * \brief Change the number of running worker threads.
     * 
     * \param n if 1, increment the number of running workers, if -1, decrement it.
     * \return the number of running worker threads
     */
    int running(int n)
    {
      boost::lock_guard<boost::mutex> lock(mutex);
      runningWorkers += n;
      return runningWorkers;
    }
    
    /**
     * \brief Get the number of running worker threads.
     * \return the number of running worker threads
     */
    int running() const
    {
      return runningWorkers;
    }
    
  private:
    std::queue<T> queue;
    mutable boost::mutex mutex; // has to be mutable because of empty/size
    boost::condition_variable filled;
    int runningWorkers;         // number of worker threads running on this queue
  };
  
  //----------------------------------------------------------------------------

  /**
   * \ingroup threading
   * \brief Abstract interface for tasks to be scheduled for concurrent execution.
   */
  typedef std::packaged_task<void()> Task;
  
  /**
   * \ingroup threading
   * \brief Abstract waitable job ticket for submitted tasks.
   */
  typedef std::future<void> Ticket;
  
  

  //----------------------------------------------------------------------------

  /**
   * \ingroup threading
   * \brief Implementation of thread pools suitable for parallelization of (more or less) 
   *        memory-bound algorithms (not only) on NUMA machines.
   * 
   * This class maintains two thread pools satisfying different needs.
   * 
   * The threads in the first (global) pool can be moved by the operating system freely between 
   * nodes and CPUs, and should be used (by submitting tasks via \ref run)
   * whenever there is no particular need to execute the task on a particular NUMA node, 
   * i.e. if the task is not memory-bandwidth-bound or works on data that is not located 
   * on a particular NUMA node. The number of global threads defaults to twice the number 
   * of available CPUs (such that all CPUs can be busy even if some threads wait on a mutex), 
   * but is at least 4 unless limited to a smaller number on construction of the thread pool, 
   * see \ref instance.
   * 
   * The threads in the second (NUMA) pool are pinned to nodes (the OS is allowed to move them between CPUS on the same node), and should be used 
   * (by submitting tasks via \ref runOnNode) whenever the tasks are memory-bandwidth-bound and the data resides on a particular NUMA node. 
   * Note that as the threads are locked to the nodes, the operating system cannot move the threads. This is intended to keep
   * the thread close to its data in order to have local memory access. On the other hand, on multi-user machines it can lead to
   * several threads competing for the same CPU while other CPUs are idle. Use the node-locked threads only if locality of memory 
   * access is of top priority.
   * 
   * Relying on the first-touch policy for controlling memory locality is not guaranteed to keep threads and data close to each other.
   * First, the allocator may re-use memory blocks previously touched and released by a different thread, leading to remote data
   * access. Second, the operating system may decide to move a thread to a different node without knowledge of which thread
   * is memory-bound or compute-bound.
   * 
   * Caveat: When submitting tasks recursively to the task pool, it is easy to create a deadlock. While the top level tasks wait for the 
   *         completion of the lower level tasks, the lower level task is not processed as all threads in the pool are occupied by the 
   *         top level tasks. Take care not to submit more recursive tasks than there are worker threads available.
   */
  class NumaThreadPool 
  {
  public:
    /**
     * \brief Returns a globally unique thread pool instance.
     * 
     * On the very first call of this method, the singleton thread pool is created with a limitation of the number of global (unpinned)
     * threads limited by the given number of maxThreads. Later calls with different value of maxThreads do not change the number of threads.
     * In case the number of global threads shall be limited throughout, call this method at the very start of the program.
     * 
     * \param maxThreads an upper bound for the number of global threads to create.
     * 
     * As it makes little sense to have multiple thread pools fight for physical ressources, a single instance should be employed.
     */
    static NumaThreadPool& instance(int maxThreads = std::numeric_limits<int>::max());
    
    /**
     * \name System information
     * @{
     */
    
    /**
     * \brief Reports the number of NUMA nodes (i.e., memory interfaces/CPU sockets)
     */
    int nodes() const
    {
      return nNode;
    }
    
    /**
     * \brief Reports the total number of CPUs (usually a multiple of nodes).
     * 
     * This is based on what std::thread::hardware_concurrency() returns, and is not guaranteed to be what hardware is really available.
     * Precisely, it max(1,hardware_concurrency()).
     */
    int cpus() const
    {
      return nCpu;
    }
    
    /**
     * \brief Reports how many worker threads are running to work on the global task queue
     */
    int runningOnGlobalQueue() const
    {
      return globalQueue.running();
    }
        
    /**
     * \brief Reports the number of CPUs on the given node (usually the same for all nodes).
     */
    int cpus(int node) const
    {
      return cpuByNode[node].size();
    }
    
    /**
     * \brief Reports the maximal number of CPUs on one node.
     */
    int maxCpusOnNode() const
    {
      return maxCpusPerNode;
    }

    /**
     * \brief Returns true if tasks are executed sequentially.
     * sequential execution can be enforced by calling NumaThreadPool::instance(1) at program start. Note that nevertheless
     * tasks are executed in a std::packaged_task context, which means that exceptions are caught and rethrown in future:get().
     * Thus, stack context is lost when debugging even for sequential execution.
     */
    bool isSequential() const
    {
      return sequential;
    }
    
    /**
     * @}
     */
    
    /**
     * \name Task submission
     * @{
     */
        
    /**
     * \brief Schedules a task to be executed on an arbitrary CPU.
     * 
     * \param task the task object to be executed. 
     * 
     * \return a waitable ticket. Call wait() on the ticket in order to wait for the task to be completed.
     * 
     * Note that waiting for tasks blocks a thread. Waiting for other tasks within 
     * a task may use up threads from the thread pool and lead to deadlocks.
     */
    Ticket run(Task&& task);
    
    /**
     * \brief Schedules a task to be executed on a CPU belonging to the given NUMA node.
     * 
     * \param node the number of the node on which to execute the task. 0 <= node < nodes().
     * \param task the task object to be executed. The object has to live at least until the call to its call operator returns.
     * 
     * \return a waitable ticket. Call wait() on the ticket in order to wait for the task to be completed.
     * 
     * Note that waiting for tasks blocks a thread. Waiting for other tasks within 
     * a task may use up threads from the thread pool and lead to deadlocks.
     */
    Ticket runOnNode(int node, Task&& task);
    
    /**
     * @}
     */
    
    /**
     * \name NUMA-aware memory management
     * @{
     */
    
    /**
     * \brief Returns the allocator used for the given node.
     */
    Kalloc& allocator(int node);
    
    /**
     * \brief Allocates memory on a specific node.
     * 
     * Note: this is comparatively slow and should only be used for allocating large chunks of memory
     * to be managed locally. The memory has to be released by a subsequent call to deallocate.
     */
    void* allocate(size_t n, int node);
    
    /**
     * \brief frees a chunk of memory previously allocated
     */
    void deallocate(void* p, size_t n, int node);
    
    /**
     * \brief Reports the alignment size of allocator at given NUMA node.
     */
    size_t alignment(int node) const;
    
    /**
     * \brief Tells the allocator to prepare for subsequent allocation of several memory blocks of same size.
     * \param n the requested size of the memory blocks
     * \param k the number of memory blocks that will be requested
     * \param node on which NUMA node
     *
     * Use this as a hint for the allocator that allows to improve its performance.
     */
    void reserve(size_t n, size_t k, int node);
    
    /**
     * @}
     */
    
  private:
    
    // private constructor to be called by instance()
    NumaThreadPool(int maxThreads);
    ~NumaThreadPool();
    
    Ticket runOnQueue(ConcurrentQueue<Task>& queue, Task&& task);
    
    
    int nCpu, nNode, maxCpusPerNode;
    std::vector<ConcurrentQueue<Task>>    nodeQueue;    // task queues for passing to worker threads on nodes
    ConcurrentQueue<Task>                 globalQueue;  
    boost::thread_group                   threads;      // worker threads
    std::vector<int>                      nodeByCpu;    // NUMA topology info
    std::vector<std::vector<int>>         cpuByNode;    // NUMA topology info
    std::map<void*,std::pair<size_t,int>> memBlocks;    // NUMA memory allocations
    std::vector<Kalloc>                   nodeMemory;   // NUMA memory management
    bool sequential;                                    // if true, execute all tasks immediately
  };
  
  //-----------------------------------------------------------------------------------------------------
  
  /**
   * \ingroup threading
   * \brief A parallel for loop that executes the given functor in parallel on different CPUs
   *
   * \tparam Func a functor with an operator ()(int i, int n)
   * 
   * \param f the functor to call
   * \param maxTasks the maximal number of tasks to create
   *
   * The given function object is called in parallel for values (i,n) with i ranging from 0 to n-1, i.e.
   * it is the i-th call of total n calls. The function object is responsible for partitioning the work 
   * appropriately, e.g., using uniformWeightRange(), and perform locking if needed.
   * The number n is the same for all calls and guaranteed not to exceed maxTasks (but may be smaller). 
   * 
   * The function returns after the last task has been completed, therefore the computational 
   * effort for tasks 0,...,n-1 should be roughly equal, no matter which value n has. For an optimal performance,
   * maxTasks should be either a small multiple of the number of CPUs in the system (such that ideally all CPUs are 
   * busy the whole time) or much larger than that (such that an imbalance has only a small impact). 
   */
  template <class Func>
  void parallelFor(Func const& f, int maxTasks = std::numeric_limits<int>::max())
  {
    NumaThreadPool& pool = NumaThreadPool::instance();

    // If tasks are to be executed sequentially, we shortcut the task pool. Less for avoiding overhead,
    // but more for avoiding task packaging, which looses call stack info if exceptions are thrown.
    if (pool.isSequential())
    {
      f(0,1);
      return;
    }

    int nTasks = std::min(4*pool.cpus(),maxTasks);
    
    std::vector<Ticket> tickets(nTasks);
    for (int i=0; i<nTasks; ++i)
      tickets[i] = pool.run(Task([i,&f,nTasks] { f(i,nTasks); }));
    
    for (auto& ticket: tickets)   // wait for the tasks to be completed
      ticket.get();               // we use get() instead of wait() to rethrow any exceptions
  }
  
  /**
   * \ingroup threading
   * \brief A parallel for loop that executes the given functor in parallel on different CPUs
   *
   * \tparam Func a functor with an operator ()(int i)
   * 
   * \param f the functor to call
   * \param first the first iteration index, usually 0
   * \param last one after the last iteration index, i.e. the number of iterations.
   *             Precondition: last >= first.
   * 
   * The given function is called last-first times, exactly once for each iteration index,
   * with the only argument being the iteration index. The granularity can be quite fine.
   * The computational work for different iteration indices can be rather different without
   * sacrificing parallel efficiency.
   *
   * min(last-first,nCpus,nTasks) tasks are created. Each tasks grabs the next not yet
   * executed iteration index and calls the functor with that index, then repeats.
   */
  template <class Func>
  void parallelFor(size_t first, size_t last, Func const& f, size_t nTasks=std::numeric_limits<size_t>::max())
  {
    assert(last>=first);
    NumaThreadPool& pool = NumaThreadPool::instance();
    nTasks = std::min(last-first,std::min(nTasks,size_t(pool.cpus())));
    
    std::vector<Ticket> tickets(nTasks);
    std::mutex mutex;
    
    TaskTiming tt(nTasks+2);
    
    tt.start(nTasks);
    for (size_t i=0; i<nTasks; ++i)
    {
      tickets[i] = pool.run(Task([&f,&first,last,&mutex,i,&tt] 
                            { 
                              tt.start(i);
                              while (true)
                              {
                                std::unique_lock<std::mutex> lock(mutex);
                                size_t myPos = first;
                                ++first;
                                lock.unlock();
                                
                                if (myPos>=last)
                                {
                                  tt.stop(i);
                                  return;
                                }
                                
                                f(myPos); 
                              }
                            }));
      
      // Waking up threads appears to be quite time-consuming. If we have many threads in the pool and not so much work, 
      // the last threads are woken up when all the work has already been done. Then tasks are created without benefit.
      // In comparison, locking on a mutex appears to be blazingly fast, so we check whether there is still work to be
      // done, and if not, leave the loop. This should be particularly beneficial in hyperthreaded systems.
      std::lock_guard<std::mutex> lock(mutex);
      if (first>=last)
      {
        tickets.resize(i+1);
        break;
      }
    }
    tt.stop(nTasks);
    
    tt.start(nTasks+1);
    for (auto& ticket: tickets)   // wait for the tasks to be completed
      ticket.get();               // we use get() instead of wait() to rethrow any exceptions
    tt.stop(nTasks+1);
    
    std::ofstream out("timing.gnu");
    out << tt;
  }
  
  
  /**
   * \ingroup threading
   * \brief A parallel for loop that executes the given functor in parallel on different NUMA nodes
   *
   * \tparam Func a functor with an operator ()(int i, int n)
   * 
   * \param f the functor to call
   * \param maxTasks the maximal number of tasks to create
   *
   * The given function object is called in parallel for values (i,n) with i ranging from 0 to n-1. The number n is
   * min(nNodes,maxTasks) for all calls to f. Every task is executed on the node given by arguemnt i.
   * 
   * The function returns after the last task has been completed, therefore the computational 
   * effort for tasks 0,...,n-1 should be roughly equal, no matter which value n has.
   */
  template <class Func>
  void parallelForNodes(Func const& f, int maxTasks = std::numeric_limits<int>::max())
  {
    NumaThreadPool& pool = NumaThreadPool::instance();
    int nTasks = std::min(pool.nodes(),maxTasks);
    std::vector<Ticket> tickets(nTasks);
    for (int i=0; i<nTasks; ++i)
      tickets[i] = pool.runOnNode(i,Task([i,&f,&nTasks] { f(i,nTasks); }));
    for (auto& ticket: tickets) // wait for the tasks to be completed
      ticket.get();             // we use get() instead of wait() to rethrow any exceptions
  }
  
  //----------------------------------------------------------------------------
  
  /**
   * \cond internals
   */
  namespace ThreadingDetail 
  {
    class NumaAllocatorBase
    {
    public:
      NumaAllocatorBase(int node);
      
      size_t max_size() const;
      
      void* allocate(size_t n);
      
      void deallocate(void* p, size_t n);
      
      int node() const
      {
        return nod;
      }
      
    private:
      int nod;    
      Kalloc* allocator;
    };
  }
  /**
   * \endcond
   */
  
  /**
   * \ingroup threading
   * \brief An STL allocator that uses memory of a specific NUMA node only.
   */
  template <class T>
  class NumaAllocator
  {
  public:
    typedef T value_type;
    typedef T* pointer;
    typedef T const* const_pointer;
    typedef T& reference;
    typedef T const& const_reference;
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;
    
    // make sure that on copying, the target copy has its data in the same 
    // NUMA memory region as the source by default -- this improves data locality.
    typedef std::true_type propagate_on_container_copy_assignment;
    typedef std::true_type propagate_on_container_move_assignment;
    typedef std::true_type propagate_on_container_swap;
    
    template <class U>
    struct rebind
    {
      typedef NumaAllocator<U> other;
    };
    
    /**
     * \brief Construct an allocator for allocating on the given NUMA node.
     * 
     * \param node the NUMA node on which to allocate the memory. This has to be less than NumaThreadPool::instance().nodes(). 
     *             For negative values, the memory is allocated in interleaved mode.
     */
    NumaAllocator(int node): alloc(node)
    {
    }
    
    /**
     * \brief Reports the node on which we allocate.
     */
    int node() const
    {
      return alloc.node();
    }
    
    pointer address(reference x) const
    {
      return &x;
    }
    
    const_pointer address( const_reference x ) const
    {
      return &x;
    }
    
    size_type max_size() const 
    {
      return alloc.max_size() / sizeof(T);
    }
    
    /**
     * \brief Allocates the requested amount of memory.
     * 
     * \param n number of objects of type T
     * 
     * If \arg n == 0, a null pointer is returned. If the request cannot be fulfilled,
     * std::bad_alloc is thrown.
     */
    pointer allocate(size_type n, std::allocator<void>::const_pointer /* hint */ = 0)
    {
      if (n>0)
        return static_cast<pointer>(alloc.allocate(n*sizeof(T)));
      else
        return nullptr;
    }
    
    void deallocate(pointer p, size_type n)
    {
      if (p)
        alloc.deallocate(static_cast<void*>(p),n*sizeof(T));
    }
        
    template< class U, class... Args >
    void construct(U* p, Args&&... args)
    {
      ::new((void*)p) U(std::forward<Args>(args)...);
    }
    
    template <class U>
    void destroy(U* p)
    {
      p->~U();
    }
    
    /**
     * \brief comparison for equality
     * 
     * Allocators compare equal, if they allocate on the same NUMA node.
     */
    template <class U>
    bool operator==(NumaAllocator<U> const& other) const
    {
      return node()==other.node();
    }
    
    template <class U>
    bool operator!=(NumaAllocator<U> const& other) const
    {
      return !(node() == other.node());
    }
    
  private:
    ThreadingDetail::NumaAllocatorBase alloc;
  };
  
  //----------------------------------------------------------------------------
  
  /**
   * \ingroup threading
   * \brief A utility class implementing appropriate copy semantics for boost mutexes.
   * 
   * \todo: design appropriate semantics for move construction/assignment
   */
  class Mutex
  {
  public:
    /**
     * \brief Default constructor.
     * 
     * The constructed object's mutex is unlocked.
     */
    Mutex() = default;
    
    /**
     * \brief Copy constructor.
     * 
     * The constructed object's mutex is unlocked, independently of the state of m.
     */
    Mutex(Mutex const& m) {}
    
    /**
     * \brief Assignment.
     * 
     * Our mutex doesn't change its state.
     */
    Mutex& operator=(Mutex const& m) { return *this; }
    
    /**
     * \brief provides access to the mutex to perform the locking.
     */
    boost::mutex& get() { return mutex; }
    
  private:
    boost::mutex mutex; 
  };
  
  //-------------------------------------------------------------------------------------
  
  /**
   * \ingroup threading
   * \brief Executes a function in a child process.
   * 
   * The given function is executed independently of and in parallel to the caller in a child process.
   * Due to the address space separation, the caller may proceed to modify its data structures, while
   * the child is not affected. Uses this to execute expensive "fire and forget" tasks that rely on the 
   * data structures not to change and are difficult to parallelize. 
   * 
   * The typical use case is to write output files in the background, while the computation continues.
   */
  void runInBackground(std::function<void()>& f);
}

#endif
