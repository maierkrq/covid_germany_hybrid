/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2016-2022 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef TIMING_HH
#define TIMING_HH

#include <chrono>
#include <memory>
#include <mutex>
#include <ostream>
#include <stack>
#include <vector>


namespace Kaskade
{
  /**
   * \ingroup utilities
   * \brief Supports gathering and reporting execution times information for nested program parts.
   * 
   * Entering and leaving timed regions incurs a (small but measurable) overhead. Thus, the timed regions
   * should not be too quick. E.g., local assembly on a single cell is too fast for providing reasonably
   * accurate timings. 
   * 
   * While calling the methods is thread-safe, using it from within parallelized loops will lead to errors such
   * as stopping a section that is currently not running (because in a different thread some other section has
   * just been started). If needed, one timer per thread can be used.
   * \todo Implement merging of timer statistics
   * 
   * Example:
   * \code
   * auto& timer = Timings::instance();
   * timer.start("my expensive calculation");
   * f();   // if f throws an exception, the timer will not be stopped! (Use ScopedTimingSection)
   * timer.stop("my expensive calculation");
   * \endcode
   * 
   * \see ScopedTimingSection
   */
  class Timings
  {
    struct Times;

  public:
    /**
     * \brief Returns a reference to a single default instance.
     */
    static Timings& instance();
    
    /**
     * \brief Constructor
     */
    Timings();
    
    /**
     * \brief Prints a timing report to the given stream.
     * \param maxlevel report only timings up to a nesting level given here
     */
    std::ostream& report(std::ostream& out, int maxlevel=1000) const;
    
    /**
     * \brief Starts or continues the timing of given section.
     * 
     * If the section has not been started before, it is created.
     *
     * \return a ticket that can be used to check on stopping that the section to be stopped
     *         is the currently running one. Using tickets instead of names for section stopping
     *         prevents run time errors due to typos in section names (and may be slightly faster).
     *
     * Example:
     * \code
     * auto& timer = Timings::instance();
     * auto ticket = timer.start("my expensive calculation");
     * f();   // if f throws an exception, the timer will not be stopped! (Use ScopedTimingSection)
     * timer.stop(ticket);
     * \endcode
     */
    struct Times const* start(std::string const& name);
    
    /**
     * \brief Stops the timing of given section.
     * \param name if the running section's name does not coincide with the given one, a LookupException is thrown
     */
    void stop(std::string const& name);

    /**
     * \brief Stops the timing of given section.
     * \param ticket if the running section's ticket does not coincide with the given one,
     *               a LookupException is thrown. If no ticket is provided, no checking
     *               takes place (which can save a few CPU cycles, in case that's
     *               relevant).
     */
    void stop(struct Times const* ticket = nullptr);

    
    /**
     * \brief Resets the timer to an empty state.
     */
    void clear();
    
    // needed for use of unique_ptr with incomplete type
    ~Timings();
    
  private:
    std::unique_ptr<Times> all;
    std::stack<Times*> stack;
    mutable std::mutex mutex;
  };
  
  /**
   * \ingroup io
   * \brief write a timing report to a stream
   * \relates Timings 
   */
  std::ostream& operator<<(std::ostream& out, Timings const& timings);
  
  /**
   * \ingroup utilities
   * \brief A scope guard object that automatically closes a timing section on destruction.
   * 
   * This closes timing sections in case an exception is thrown.
   * 
   * Example:
   * \code
   * {
   *   ScopedTimingSection tsec("my expensive calculation");
   *   f();   // f may throw an exception, timing will be stopped
   * }
   * \endcode

   * \see Timings
   */
  class ScopedTimingSection
  {
  public:
    /**
     * \brief Constructor using the default timer provided by Timings::instance().
     * 
     * \param name the name of the section to be timed.
     * \param gather if false, do *not* gather any timings. This can be used to
     *               switch off timing collection in potentially multithreaded environments,
     *               while at the same time use ScopedTimingSection guards.
     * The timer is started immediately.
     */
    ScopedTimingSection(std::string const& name, bool gather=true);
    
    /**
     * \brief Constructor using provided timer.
     * 
     * \param timer a timer. This has to exist during the lifetime of the scoped timing section.
     * 
     * The timer is started immediately.
     */
    ScopedTimingSection(std::string const& name, Timings& timer, bool gather=true);
    
    ~ScopedTimingSection();
    
    /**
     * \brief Ensures the section is timed. 
     * 
     * If the timer is not running, it is started. Otherwise nothing is done.
     */
    void start();

    /**
     * \brief Ensures the section is not timed. 
     * 
     * If the timer is running, it is stopped. Otherwise nothing is done.
     */
    void stop();
    
    /**
     * \brief Reports whether the section is currently timed.
     */
    bool isRunning() const;

    /**
     * \brief Returns the used timer.
     */
    Timings& timer();
    
  private:
    std::string name;
    Timings& mytimer;
    bool running;
    bool gather;        // if false, do not call mytimer at all
  };
  
  // ----------------------------------------------------------------------------------------------
  
  /**
   * \ingroup utilities
   * \brief A class that gathers data on task timing and provides gnuplot visualization.
   * 
   * This is particularly useful as a performance investigation tool in multithreaded computations, where
   * one can investigate which task is performed during which time interval.
   * 
   * \tparam n the maximal number of segments for each task
   *
   * Example usage:
   * \code
   * TaskTiming tt(2);
   * start(0); doSomeWork(); stop(0);
   * start(1); doSomeOtherWork(); stop(1);
   * std::ofstream timings("timings.gnu"); timings << tt; 
   * \endcode
   * For visualization in form of horizontal task durations in second scale, use gnuplot as
   * \code
   * plot "timings.gnu" with linepoints
   * \endcode
   */
  class TaskTiming
  {
  public:
    /**
     * \brief Constructor.
     * \param nTasks number of tasks to measure. Has to be nonnegative.
     * \param m maximal number of segments per task to measure. Has to be nonnegative.
     */
    TaskTiming(int nTasks, int m=1);
    
    /**
     * \brief defines the start of given task.
     * \param task the task number in [0,nTasks[
     * If called multiple times for a task, only the last call counts.
     */
    void start(int task);

    /**
     * \brief defines the start of given task.
     * \param task the task number in [0,nTasks[
     * If called multiple times for a task, only the last call counts.
     */
    void stop(int task);
    
  private:
    using TimePoint = std::chrono::time_point<std::chrono::high_resolution_clock>;
    TimePoint zero;
    std::vector<std::array<TimePoint,2>> times;
    std::vector<int> segment;
    int nSegments;
    
    friend std::ostream& operator <<(std::ostream& out, TaskTiming const& tt);
  };
  
  /**
   * \ingroup io
   * \brief output of task timing to gnuplot file
   * \relates TaskTiming
   */
  std::ostream& operator <<(std::ostream& out, TaskTiming const& tt);
}



#endif
