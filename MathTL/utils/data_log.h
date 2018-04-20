// -*- c++ -*-

// +--------------------------------------------------------------------+
// | data_log.h, Copyright (c) 2018                                     |
// | Henning Zickermann <zickermann@mathematik.uni-marburg.de>          |
// |                                                                    |
// | This file is part of MathTL - the Mathematical Template Library.   |
// |                                                                    |
// | Contact: AG Numerik, Philipps University Marburg                   |
// |          http://www.mathematik.uni-marburg.de/~numerik/            |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_DATA_LOG_H
#define _MATHTL_DATA_LOG_H


#include <string>
#include <map>

namespace MathTL
{

  using std::map;

  /*!
    A class for recording ("logging") data pairs of real-valued
    variables that arise during a computation (e.g. convergence data).
    After the computation a graphical representation of the logged data
    can be obtained by writting the data to a Matlab plot file.
  */
  class DataLog
  {
  public:
    /*!
      Constructor. The logTitle and the axis labels affect the graphical
      representation of the data log as a Matlab plot.
      One can further specify if logarithmic scales are preferably used
      in that Matlab plot.
    */
    DataLog(const std::string& logTitle,
            const std::string& x_axis_label,
            const std::string& y_axis_label,
            bool x_scale_logarithmic,
            bool y_scale_logarithmic);

    /*!
      Add a single data pair given by (xValue, yValue) to the data log.
    */
    void addData(double xValue, double yValue);

    /*!
      Delete all logged data pairs.
    */
    void clearData();

    /*!
      Get a copy of the logged data pairs.
    */
    map<double,double> getData() const;

    /*!
      Write the logged data to a given output stream (as Matlab plot
      data). Logarithmic scales are used according to the preferences
      that were specified during the construction of the log.
    */
    void writePlotData(std::ostream& os) const;
    
    /*!
      Write the logged data to a given output stream (as Matlab plot
      data). Logarithmic scales are used according to the passed
      arguments (independently of the preferences that were specified
      during the construction of the log).
    */
    void writePlotData(std::ostream& os,
                       bool x_scale_logarithmic,
                       bool y_scale_logarithmic) const;

  protected:
    map<double,double> logData_;
    const std::string logTitle_;
    const std::string xAxisLabel_;
    const std::string yAxisLabel_;
    const bool x_scale_preferred_logarithmic_;
    const bool y_scale_preferred_logarithmic_;
  };



  /*!
    A data log for which the first variable of the logged data pairs
    represents computation time (in seconds).
  */
  class TimeLog : public DataLog
  {
  public:
    TimeLog(const std::string& logTitle,
            const std::string& y_axis_label,
            bool time_scale_logarithmic,
            bool y_scale_logarithmic);
  };

}


#include "data_log.cpp"

#endif // _MATHTL_DATA_LOG_H
