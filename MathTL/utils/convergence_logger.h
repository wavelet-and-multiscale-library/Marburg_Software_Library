// -*- c++ -*-

// +--------------------------------------------------------------------+
// | convergence_logger.h, Copyright (c) 2018                           |
// | Henning Zickermann <zickermann@mathematik.uni-marburg.de>          |
// |                                                                    |
// | This file is part of MathTL - the Mathematical Template Library.   |
// |                                                                    |
// | Contact: AG Numerik, Philipps University Marburg                   |
// |          http://www.mathematik.uni-marburg.de/~numerik/            |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_CONVERGENCE_LOGGER_H
#define _MATHTL_CONVERGENCE_LOGGER_H

#include <string>
#include <map>
#include <vector>
#include <ctime>

#include "data_log.h"

namespace MathTL
{
  using std::map;

  /*!
    Abstract interface for logging the convergence of an iterative
    method, i.e. recording the (residual) approximation error in each
    iteration against the degrees of freedom and as a function of time.

    The iterative method is supposed to accept a non-const reference to
    an AbstractConvergenceLogger as input argument. Inside the method's
    body this reference should be used to call the virtual log
    function(s) once per iteration.
    
    This way it is left to the caller of the iterative method to decide
    in which way the data is actually "logged" (by supplying a suitable
    implementation of the AbstractConvergenceLogger interface).
    In console programs "logging" could mean to simply store the data
    inside a container object (as done by the MathTL::ConvergenceLogger
    class) so it can be exported as a Matlab plot file afterwards.
    In a GUI application "logging" could further imply live plotting of
    the convergence data.

    It is also possible to define additional (optional) data logs for
    logging other data related to the computation progress.
  */
  class AbstractConvergenceLogger
  {
  public:
    /*!
      virtual destructor
    */
    virtual ~AbstractConvergenceLogger() {}

    /*!
      Log a single convergence data point given by 'degrees_of_freedom'
      and 'approx_error'.
      Implementations are supposed to also record the elapsed time in
      seconds since the clock was started (see the function
      'getElapsedSeconds()' below).
      
      Inside an iterative method this should be called once per
      iteration.
    */
    virtual void logConvergenceData(double degrees_of_freedom,
                                    double approx_error) = 0;

    /*!
      Log a text message.
    */
    virtual void logMessage(const std::string& message) = 0;

    /*!
      Define an optional data log for logging data pairs of real-valued
      variables. The given 'logName' is used to distinguish between
      different optional data logs. The axis labels concern graphical
      representations of the data log. One can also specify if
      logarithmic scales are preferred in graphical representations.
    */
    virtual void defineOptionalDataLog(const std::string& logName,
                                       const std::string& x_axis_label,
                                       const std::string& y_axis_label,
                                       bool x_scale_logarithmic,
                                       bool y_scale_logarithmic) = 0;

    /*!
      Define an optional time log for logging a single, real-valued
      variable as a function of time. The given 'logName' is used to
      distinguish between different optional time logs. The axis label
      concerns graphical representations of the time log. One can also
      specify if logarithmic scales are preferred in graphical
      representations.
    */
    virtual void defineOptionalTimeLog(const std::string& logName,
                                       const std::string& y_axis_label,
                                       bool time_scale_logarithmic,
                                       bool y_scale_logarithmic) = 0;

    /*!
      Add a single data pair given by (x_value, y_value) to a previously
      defined optional data log with the given 'logName'.
    */
    virtual void addToDataLog(const std::string& logName,
                              double x_value, double y_value) = 0;

    /*!
      Add a single data pair given by (elapsedSeconds, y_value) to a
      previously defined optional time log with the given 'logName'.
      Here 'elapsedSeconds' is the elapsed time in seconds since the
      clock was started (see the function 'getElapsedSeconds()' below).
    */
    virtual void addToTimeLog(const std::string& logName,
                              double y_value) = 0;

    /*!
      Define an abort condition leading to termination of the
      computation in case the (residual) error gets greater than
      errorLimit (which could indicate that the computation might not
      converge at all).
    */
    virtual void setAbortCondition_ErrorExceeds(double errorLimit) = 0;

    /*!
      Define an abort condition leading to termination of the
      computation in case the number of iterations exceeds iterLimit.
    */
    virtual void setAbortCondition_IterationCountExceeds(int iterLimit) = 0;

    /*!
      Check if one of the previously defined abort conditions is met, in
      which case the computation is supposed to be terminated (f.e. by
      throwing an exception).
      In multi-threaded GUI applications it can also be used to
      terminate the computation gracefully when the user decides to
      cancel it.
      
      Because of the last-mentioned use case, this should be called very
      regularly (at least once per iteration).
    */
    virtual void checkAbortConditions() const = 0;

    /*!
      Start the clock which is used to measure time values for time
      logs.
      
      This should be called only once inside the iterative method's body
      (just before the first iteration).
    */
    virtual void startClock() = 0;

    /*!
      Pause the clock (f.e. during the computation of data for an
      optional log which is normally not computed as part of the
      iterative method).
    */
    virtual void pauseClock() = 0;

    /*!
      Continue the clock after it was paused.
    */
    virtual void continueClock() = 0;

    /*!
      Get the elapsed time (in seconds) since the clock was started by
      startClock(), not counting time intervals when the clock was
      paused (i.e. the time between a call of pauseClock() and
      continueClock()).
    */
    virtual double getElapsedSeconds() = 0;
  };




  /*!
    An implementation of the AbstractConvergenceLogger interface using
    the MathTL::DataLog class (as defined in MathTL/utils/data_log.h).

    For recording convergence data, the logger manages two different
    data logs: a "convergence log" and a "convergence time log" for
    recording the (residual) approximation error as a function of the
    degrees of freedom and as a function of time respectively.

    The convergence logs (as well as possible optional logs) can be
    written to a Matlab plot file when the logging process has finished.
  */
  class ConvergenceLogger : public AbstractConvergenceLogger
  {
  public:
    /*!
      Constructor. One can specify for both the convergence log and the
      convergence time log if logarithmic scales are preferred in
      graphical representations of these logs (in the form of Matlab
      plots).
    */
    ConvergenceLogger(bool convLog_x_scale_logarithmic = true,
                      bool convLog_y_scale_logarithmic = true,
                      bool time_scale_logarithmic = true,
                      bool timeLog_y_scale_logarithmic = true);

    /*!
      Log a single convergence data point. A data pair given by
      (degrees_of_freedom, approx_error) is added to the convergence log
      and a data pair given by (elapsed_seconds, approx_error) is added
      to the convergence time log, where 'elapsed_seconds' is the
      elapsed time in seconds since the clock was started (the value
      that 'getElapsedSeconds()' returns at the time of the logging).
    */
    void logConvergenceData(double degrees_of_freedom,
                            double approx_error);

    /*!
      Log a text message. Currently the message is just forwarded to the
      standard text output.
    */
    void logMessage(const std::string& message);

    /*!
      Define an optional data log for logging data pairs of real-valued
      variables. The given 'logName' is used to distinguish between
      different optional data logs. The axis labels concern graphical
      representations of the data log. One can also specify if
      logarithmic scales are preferred in graphical representations.
    */
    void defineOptionalDataLog(const std::string& logName,
                               const std::string& x_axis_label,
                               const std::string& y_axis_label,
                               bool x_scale_logarithmic,
                               bool y_scale_logarithmic);
    /*!
      Define an optional time log for logging a single, real-valued
      variable as a function of time. The given 'logName' is used to
      distinguish between different optional time logs. The axis label
      concerns graphical representations of the time log. One can also
      specify if logarithmic scales are preferred in graphical
      representations.
    */
    void defineOptionalTimeLog(const std::string& logName,
                               const std::string& y_axis_label,
                               bool time_scale_logarithmic,
                               bool y_scale_logarithmic);

    /*!
      Add a single data pair given by (x_value, y_value) to a previously
      defined optional data log with the given 'logName'.
    */
    void addToDataLog(const std::string& logName,
                      double x_value, double y_value);

    /*!
      Add a single data pair given by (elapsedSeconds, y_value) to a
      previously defined optional time log with the given 'logName'.
      Here 'elapsedSeconds' is the elapsed time in seconds since the
      clock was started (the value that 'getElapsedSeconds()' returns at
      the time of the call).
    */
    void addToTimeLog(const std::string& logName, double y_value);

    /*!
      Define an abort condition leading to termination of the
      computation in case the (residual) error gets greater than
      errorLimit (which could indicate that the computation might not
      converge at all).
    */
    void setAbortCondition_ErrorExceeds(double errorLimit);

    /*!
      Define an abort condition leading to termination of the
      computation in case the number of iterations exceeds iterLimit.
    */
    void setAbortCondition_IterationCountExceeds(int iterLimit);

    /*!
      Check if one of the previously defined abort conditions is met,
      in which case the computation is terminated by throwing an
      std::runtime_error.
    */
    void checkAbortConditions() const;

    /*!
      Start the clock which is used to measure time values for the
      convergence time log and optional time logs.
    */
    void startClock();

    /*!
      Pause the clock (f.e. during the computation of data for an
      optional log which is normally not computed as part of the
      iterative method).
    */
    void pauseClock();

    /*!
      Continue the clock after it was paused.
    */
    void continueClock();

    /*!
      Get the elapsed time (in seconds) since the clock was started by
      startClock(), not counting time intervals when the clock was
      paused (i.e. the time between a call of pauseClock() and
      continueClock()).
    */
    double getElapsedSeconds();

    /*!
      Read access to the convergence log which records the (residual)
      approximation error as a function of the degrees of freedom that
      are used for the approximation.
    */
    const DataLog& getConvergenceLog() const;
    
    /*!
      Read access to the convergence time log which records the
      (residual) approximation error as a function of time.
    */
    const TimeLog& getConvergenceTimeLog() const;
    
    /*!
      Read access to a previously defined optional data log
      named 'logName'.
    */
    const DataLog& getOptionalDataLog(const std::string& logName) const;
    
    /*!
      Read access to a previously defined optional time log
      named 'logName'.
    */
    const TimeLog& getOptionalTimeLog(const std::string& logName) const;
    
    /*!
      Get a list of the names of all previously defined optional
      data logs.
    */
    std::vector<std::string> getOptionalDataLogNames() const;

    /*!
      Get a list of the names of all previously defined optional
      time logs.
    */
    std::vector<std::string> getOptionalTimeLogNames() const;

    /*!
      Write the data of both the convergence log and the convergence
      time log to a given output stream (as Matlab plot data).
    */
    void writeConvergencePlots(std::ostream& os) const;
    
    /*!
      Write the data of a previously defined optional data log with the
      given logName to a given output stream (as Matlab plot data).
    */
    void writeOptionalDataPlot(const std::string& logName,
                               std::ostream& os) const;
    
    /*!
      Write the data of a previously defined optional time log with the
      given logName to a given output stream (as Matlab plot data).
    */
    void writeOptionalTimePlot(const std::string& logName,
                               std::ostream& os) const;

    /*!
      Reset the logger. This will delete all convergence data, optional
      logs and abort conditions resulting in the same state of the
      logger as it had after its construction.
    */
    void reset();

  protected:
    DataLog convergenceLog_;
    TimeLog convergenceTimeLog_;
    map<std::string, DataLog> optionalDataLogs_;
    map<std::string, TimeLog> optionalTimeLogs_;

    bool hasAbortCondition_Error_;
    double errorLimit_;
    double lastApproxError_;

    bool hasAbortCondition_Iterations_;
    int iterLimit_;
    int iterations_;

  private:
    std::clock_t clockTicksAtLastStart_;
    std::clock_t totalclockTicksUntilLastPause_;
    bool clockIsRunning_;
  };
  
  
  
  /*!
    A dummy implementation of the AbstractConvergenceLogger interface
    (doing nothing).
  */
  class DummyLogger : public AbstractConvergenceLogger
  {
  public:
    void logConvergenceData(double degrees_of_freedom,
                            double approx_error) {}

    void logMessage(const std::string& message) {}

    void defineOptionalDataLog(const std::string& logName,
                               const std::string& x_axis_label,
                               const std::string& y_axis_label,
                               bool x_scale_logarithmic,
                               bool y_scale_logarithmic) {}

    void defineOptionalTimeLog(const std::string& logName,
                               const std::string& y_axis_label,
                               bool time_scale_logarithmic,
                               bool y_scale_logarithmic) {}

    void addToDataLog(const std::string& logName,
                      double x_value, double y_value) {}

    void addToTimeLog(const std::string& logName, double y_value) {}

    void setAbortCondition_ErrorExceeds(double errorLimit) {}

    void setAbortCondition_IterationCountExceeds(int iterLimit) {}

    void checkAbortConditions() const {}

    void startClock() {}

    void pauseClock() {}

    void continueClock() {}

    double getElapsedSeconds()
    {
      return 0.0;
    }
  };

}

#include "convergence_logger.cpp"

#endif // _MATHTL_CONVERGENCE_LOGGER_H
