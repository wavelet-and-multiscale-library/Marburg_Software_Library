// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_FUNCTION_TIME_H
#define _MATHTL_FUNCTION_TIME_H

namespace MathTL
{
  /*!
    Base class for time-dependent functions

    (cf. deal.II library, version 5)
   */
  class FunctionTime
  {
  public:
    /*!
      default constructor, default initial time is zero
    */
    FunctionTime(const double initial_time = 0.0)
      : time(initial_time) {}

    /*!
      virtual destructor
    */
    virtual ~FunctionTime() {}
  
    /*!
      get current value of the time variable
    */
    inline double get_time() const { return time; }

    /*!
      set current time value
      (this could be overloaded by derived classes to perform
      additional calculations)
    */
    virtual void set_time(const double new_time) { time = new_time; }

    /*!
      advance current time value by some time step
    */
    virtual void advance_time(const double delta_t) { set_time(time + delta_t); }

  private:
    /*!
      current time value
    */
    double time;
  };
}

#endif
