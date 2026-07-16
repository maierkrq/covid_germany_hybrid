/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2016-2019 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <cassert>
#include <iomanip>

#include <boost/timer/timer.hpp>

#include "utilities/detailed_exception.hh"
#include "utilities/timing.hh"

namespace Kaskade
{

  struct Timings::Times
  {
    std::string             name;
    size_t                  calls;
    std::vector<Times>      children;
    boost::timer::cpu_timer timer;

    int maxNameLength(int indent) const
    {
      int len = name.size();
      for (auto const& c: children)
        len = std::max(len,indent+c.maxNameLength(2));
      return len;
    }

    // level: nesting depth, used for indentation
    // maxlevel: report only up to this detail
    // pos: column where to position the number of section calls
    // full: time spent in the 
    // returns: true if details were omitted, false if complete
    bool report(std::ostream& out, int level, int maxlevel, int pos, boost::timer::nanosecond_type full) const
    {
      out << std::setw(2*level) << ""
          << std::setw(3) << std::right << 100*timer.elapsed().wall/full << "% "
          << std::setw(pos-2*level) << std::left << name
          << std::setw(6) << std::right << calls << "  ";
      out << timer.format();

      full = timer.elapsed().wall;
      bool omittedDetails = false;
      for (auto const& c: children)
        omittedDetails = omittedDetails || c.report(out,level+1,maxlevel,pos,full);
      
      return omittedDetails;
    }
  };

  // ----------------------------------------------------------------------------------------------

  Timings& Timings::instance()
  {
    static Timings globalTimings;
    return globalTimings;
  }

  Timings::Timings()
  : all(new Times{"all",0,{},{}})
  {
    stack.push(all.get());
    all->timer.start();
  }

  std::ostream& Timings::report(std::ostream& out, int maxlevel) const
  {
    std::unique_lock<std::mutex> lock(mutex);
    int pos = all->maxNameLength(2)+4;
    out << std::setw(pos+4) << " " << "  calls   " << "\n";
    //~ out << std::setw(pos+30) << " " << "\n";
    bool omittedDetails = all->report(out,0,maxlevel,pos,all->timer.elapsed().wall);
    if (omittedDetails)
      out << "[deeper nested details omitted]\n";
    return out;
  }

  Timings::Times const* Timings::start(std::string const& name)
  {
    std::unique_lock<std::mutex> lock(mutex);
    
    // Find an child of the current section with the given name - then
    // this is a subsequent call.
    Times* section = nullptr;
    for (auto& c: stack.top()->children)
      if (name == c.name)
      {
        section = &c;
        section->timer.resume();
        break;
      }

    if (! section)
    {
      // Section not found - create a new one. Note that this can invalidate pointers into
      // the children container. Fortunately, there are no pointers/iterators around: The
      // stack contains only pointers into higher level containers, and any active tickets
      // handed out (i.e. pointers to Times structures) point to higher level containers.
      stack.top()->children.push_back(Times{name,0,{},{}});
      section = &(stack.top()->children.back());
      section->timer.start();
    }

    ++section->calls;
    stack.push(section);
    return section;
  }

  void Timings::stop(std::string const& name)
  {
    std::unique_lock<std::mutex> lock(mutex);
    if (stack.top()->name != name)
      throw LookupException("Timing section to be closed (" + name
                            + ") and running section (" + stack.top()->name
                            + ") have different names.",__FILE__,__LINE__);
    stack.top()->timer.stop();
    stack.pop();
  }

  void Timings::stop(Times const* ticket)
  {
    std::unique_lock<std::mutex> lock(mutex);
    if (stack.top() != ticket)
      throw LookupException("Timing section to be closed and running section have different tickets.",__FILE__,__LINE__);
    stack.top()->timer.stop();
    stack.pop();
  }

  void Timings::clear()
  {
    std::unique_lock<std::mutex> lock(mutex);
    while (stack.size()>1)  // remove everything (except the root)
      stack.pop();
    all->timer.start();      // reset timer to zero
    all->calls = 0;
    all->children.clear();
  }
  
  Timings::~Timings() = default;

  // ------------------------------------------------------------------------------------

  std::ostream& operator<<(std::ostream& out, Timings const& timings)
  {
    return timings.report(out);
  }
  
  // ------------------------------------------------------------------------------------
  
  ScopedTimingSection::ScopedTimingSection(std::string const& name_, bool gather_)
  : ScopedTimingSection(name_,Timings::instance(),gather_)
  {}
  
  ScopedTimingSection::ScopedTimingSection(std::string const& name_, Timings& timer_, bool gather_)
  : name(name_), mytimer(timer_), running(gather_), gather(gather_)
  {
    if (gather)
      mytimer.start(name);
  }

  ScopedTimingSection::~ScopedTimingSection()
  {
    stop();
  }

  void ScopedTimingSection::start()
  {
    if (gather)
    {
      if (!running)
        mytimer.start(name);
      running = true;
    }
  }
  
  void ScopedTimingSection::stop()
  {
    if (gather)
    {
      if (running)
        mytimer.stop(name);
      running = false;
    }
  }
  
  bool ScopedTimingSection::isRunning() const
  {
    return running;
  }

  Timings& ScopedTimingSection::timer()
  {
    return mytimer;
  }
 
  
  // ------------------------------------------------------------------------------------

  TaskTiming::TaskTiming(int nTasks, int m)
  : zero(std::chrono::high_resolution_clock::now())
  , times(nTasks*m), segment(nTasks,0), nSegments(m)
  { 
    assert(nTasks>=0);
    assert(m>=0);
  }    
  
  void TaskTiming::start(int task)
  {
    assert(0<=task && task<segment.size());
    times[segment[task]*segment.size()+task][0] = std::chrono::high_resolution_clock::now();
//     std::cerr << "starting task " << task << " at " << std::chrono::duration_cast<Seconds>(times[task][0]-zero).count() << "s\n";
  }
  
  void TaskTiming::stop(int task)
  {
    assert(0<=task && task<segment.size());
    times[segment[task]*segment.size()+task][1] = std::chrono::high_resolution_clock::now();
    ++segment[task];
//     std::cerr << "stopping task " << task << " at " << std::chrono::duration_cast<Seconds>(times[task][1]-zero).count() << "s\n";
  }
  
  std::ostream& operator <<(std::ostream& out, TaskTiming const& tt)
  {
    using Seconds = std::chrono::duration<double>;
    for (int i=0; i<tt.segment.size(); ++i)
      for (int j=0; j<tt.segment[i]; ++j)
        out << std::chrono::duration_cast<Seconds>(tt.times[j*tt.segment.size()+i][0]-tt.zero).count() << ' ' << i << '\n' 
            << std::chrono::duration_cast<Seconds>(tt.times[j*tt.segment.size()+i][1]-tt.zero).count() << ' ' << i << "\n\n"; 
    return out;
  }
  
}
