// -*- c++ -*-

// +--------------------------------------------------------------------+
// | convergence_logger.cpp, Copyright (c) 2018                         |
// | Henning Zickermann <zickermann@mathematik.uni-marburg.de>          |
// |                                                                    |
// | This file is part of MathTL - the Mathematical Template Library.   |
// |                                                                    |
// | Contact: AG Numerik, Philipps University Marburg                   |
// |          http://www.mathematik.uni-marburg.de/~numerik/            |
// +--------------------------------------------------------------------+


// implementation for convergence_logger.h

#include <iostream>
#include <stdexcept>

namespace MathTL
{

  inline
  ConvergenceLogger::ConvergenceLogger(bool convLog_x_scale_logarithmic,
                                       bool convLog_y_scale_logarithmic,
                                       bool time_scale_logarithmic,
                                       bool timeLog_y_scale_logarithmic)
    : convergenceLog_("Convergence", "# degrees of freedom", "(residual) error", convLog_x_scale_logarithmic, convLog_y_scale_logarithmic),
      convergenceTimeLog_("Temporal convergence", "(residual) error", time_scale_logarithmic, timeLog_y_scale_logarithmic),
      hasAbortCondition_Error_(false),
      errorLimit_(0.0),
      lastApproxError_(0.0),
      hasAbortCondition_Iterations_(false),
      iterLimit_(0),
      iterations_(0),
      totalclockTicksUntilLastPause_(0),
      clockIsRunning_(false)
  {

  }



  inline
  void ConvergenceLogger::logConvergenceData(double degrees_of_freedom, double approx_error)
  {
    convergenceLog_.addData(degrees_of_freedom, approx_error);
    convergenceTimeLog_.addData(getElapsedSeconds(), approx_error);

    lastApproxError_ = approx_error;
    iterations_++;
  }



  inline
  void ConvergenceLogger::logMessage(const std::string& message)
  {
    std::cout << message << std::endl;
  }



  inline
  void ConvergenceLogger::defineOptionalDataLog(const std::string& logName,
                                                const std::string& x_axis_label,
                                                const std::string& y_axis_label,
                                                bool x_scale_logarithmic,
                                                bool y_scale_logarithmic)
  {
    optionalDataLogs_.insert(std::make_pair(logName, DataLog(logName, x_axis_label, y_axis_label, x_scale_logarithmic, y_scale_logarithmic)));
  }



  inline
  void ConvergenceLogger::defineOptionalTimeLog(const std::string& logName,
                                                const std::string& y_axis_label,
                                                bool time_scale_logarithmic,
                                                bool y_scale_logarithmic)
  {
    optionalTimeLogs_.insert(std::make_pair(logName, TimeLog(logName, y_axis_label, time_scale_logarithmic, y_scale_logarithmic)));
  }



  inline
  void ConvergenceLogger::addToDataLog(const std::string& logName, double x_value, double y_value)
  {
    map<std::string, DataLog>::iterator search = optionalDataLogs_.find(logName);
    if (search != optionalDataLogs_.end())
    {
      search->second.addData(x_value, y_value);
    }
    else
    {
      std::cout << "Error: optional data log '" << logName << "' was not defined!" << std::endl;
    }
  }



  inline
  void ConvergenceLogger::addToTimeLog(const std::string& logName, double y_value)
  {
    map<std::string, TimeLog>::iterator search = optionalTimeLogs_.find(logName);
    if (search != optionalTimeLogs_.end())
    {
      search->second.addData(getElapsedSeconds(), y_value);
    }
    else
    {
      std::cout << "Error: optional time log '" << logName << "' was not defined!" << std::endl;
    }
  }



  inline
  void ConvergenceLogger::setAbortCondition_ErrorExceeds(double errorLimit)
  {
    errorLimit_ = errorLimit;
    hasAbortCondition_Error_ = true;
  }



  inline
  void ConvergenceLogger::setAbortCondition_IterationCountExceeds(int iterLimit)
  {
    iterLimit_ = iterLimit;
    hasAbortCondition_Iterations_ = true;
  }



  inline
  void ConvergenceLogger::checkAbortConditions() const
  {
    if (hasAbortCondition_Error_ && (lastApproxError_ > errorLimit_))
    {
      std::stringstream sstream;
      sstream << "Computation was aborted because the (residual) approximation error exceeded "
                 "the specified limit of " << errorLimit_ << std::endl;
      throw std::runtime_error(sstream.str());
    }
    if (hasAbortCondition_Iterations_ && (iterations_ > iterLimit_))
    {
      std::stringstream sstream;
      sstream << "Computation was aborted because the number of iterations exceeded the "
                 "specified limit of " << iterLimit_ << " iterations!" << std::endl;
      throw std::runtime_error(sstream.str());
    }
  }



  inline
  void ConvergenceLogger::startClock()
  {
    clockTicksAtLastStart_ = std::clock();
    totalclockTicksUntilLastPause_ = 0;
    clockIsRunning_ = true;
  }



  inline
  void ConvergenceLogger::pauseClock()
  {
    if (clockIsRunning_)
    {
      totalclockTicksUntilLastPause_ += std::clock() - clockTicksAtLastStart_;
      clockIsRunning_ = false;
    }
  }



  inline
  void ConvergenceLogger::continueClock()
  {
    if (!clockIsRunning_)
    {
      clockTicksAtLastStart_ = std::clock();
      clockIsRunning_ = true;
    }
  }



  inline
  double ConvergenceLogger::getElapsedSeconds()
  {
    if (clockIsRunning_)
    {
      return (totalclockTicksUntilLastPause_ + std::clock() - clockTicksAtLastStart_) / (double) CLOCKS_PER_SEC;
    }
    else
    {
      return totalclockTicksUntilLastPause_ / (double) CLOCKS_PER_SEC;
    }
  }



  inline
  const DataLog& ConvergenceLogger::getConvergenceLog() const
  {
    return convergenceLog_;
  }



  inline
  const TimeLog& ConvergenceLogger::getConvergenceTimeLog() const
  {
    return convergenceTimeLog_;
  }



  inline
  const DataLog& ConvergenceLogger::getOptionalDataLog(const std::string& logName) const
  {
    map<std::string, DataLog>::const_iterator search = optionalDataLogs_.find(logName);
    if (search != optionalDataLogs_.end())
    {
      return search->second;
    }
    else
    {
      std::stringstream sstream;
      sstream << "ConvergenceLogger: Error while trying to access data log '" << logName
              << "'. No optional data log defined with this name." << std::endl;
      throw std::out_of_range(sstream.str());
    }
  }
  
  
  
  inline
  const TimeLog& ConvergenceLogger::getOptionalTimeLog(const std::string& logName) const
  {
    map<std::string, TimeLog>::const_iterator search = optionalTimeLogs_.find(logName);
    if (search != optionalTimeLogs_.end())
    {
      return search->second;
    }
    else
    {
      std::stringstream sstream;
      sstream << "ConvergenceLogger: Error while trying to access time log '" << logName
              << "'. No optional time log defined with this name." << std::endl;
      throw std::out_of_range(sstream.str());
    }
  }



  inline
  std::vector<std::string> ConvergenceLogger::getOptionalDataLogNames() const
  {
    std::vector<std::string> names;

    for (map<std::string, DataLog>::const_iterator it = optionalDataLogs_.begin();
         it != optionalDataLogs_.end(); ++it)
    {
      names.push_back(it->first);
    }

    return names;
  }
  
  
  
  inline
  std::vector<std::string> ConvergenceLogger::getOptionalTimeLogNames() const
  {
    std::vector<std::string> names;

    for (map<std::string, TimeLog>::const_iterator it = optionalTimeLogs_.begin();
         it != optionalTimeLogs_.end(); ++it)
    {
      names.push_back(it->first);
    }

    return names;
  }



  inline
  void ConvergenceLogger::writeConvergencePlots(std::ostream& os) const
  {
    os << "subplot(2,1,1)" << std::endl;
    convergenceLog_.writePlotData(os);
    os << std::endl << "subplot(2,1,2)" << std::endl;
    convergenceTimeLog_.writePlotData(os);
  }



  inline
  void ConvergenceLogger::writeOptionalDataPlot(const std::string& logName, std::ostream& os) const
  {
    getOptionalDataLog(logName).writePlotData(os);
  }
  
  
  
  inline
  void ConvergenceLogger::writeOptionalTimePlot(const std::string& logName, std::ostream& os) const
  {
    getOptionalTimeLog(logName).writePlotData(os);
  }



  inline
  void ConvergenceLogger::reset()
  {
    totalclockTicksUntilLastPause_ = 0;
    clockIsRunning_ = false;

    convergenceLog_.clearData();
    convergenceTimeLog_.clearData();
    optionalDataLogs_.clear();
    optionalTimeLogs_.clear();

    hasAbortCondition_Error_ = false;
    lastApproxError_ = 0.0;
    hasAbortCondition_Iterations_ = false;
    iterations_ = 0;
  }

}
