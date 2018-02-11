// -*- c++ -*-

// +--------------------------------------------------------------------+
// | data_log.cpp, Copyright (c) 2018                                   |
// | Henning Zickermann <zickermann@mathematik.uni-marburg.de>          |
// |                                                                    |
// | This file is part of MathTL - the Mathematical Template Library.   |
// |                                                                    |
// | Contact: AG Numerik, Philipps University Marburg                   |
// |          http://www.mathematik.uni-marburg.de/~numerik/            |
// +--------------------------------------------------------------------+


// implementation for data_log.h

#include <cmath>

#include "plot_tools.h"


namespace MathTL
{

  inline
  DataLog::DataLog(const std::string& logTitle,
                   const std::string& x_axis_label,
                   const std::string& y_axis_label,
                   bool x_scale_logarithmic,
                   bool y_scale_logarithmic)
    : logTitle_(logTitle), xAxisLabel_(x_axis_label), yAxisLabel_(y_axis_label),
      x_scale_preferred_logarithmic_(x_scale_logarithmic),
      y_scale_preferred_logarithmic_(y_scale_logarithmic)
  {

  }



  inline
  void DataLog::addData(double xValue, double yValue)
  {
    logData_[xValue] = yValue;
  }



  inline
  void DataLog::clearData()
  {
    logData_.clear();
  }



  inline
  map<double,double> DataLog::getData() const
  {
    return logData_;
  }



  inline
  void DataLog::writePlotData(std::ostream& os) const
  {
    writePlotData(os, x_scale_preferred_logarithmic_, y_scale_preferred_logarithmic_);
  }



  inline
  void DataLog::writePlotData(std::ostream& os,
                              bool x_scale_logarithmic,
                              bool y_scale_logarithmic) const
  {
    matlab_output(logData_, os);
    
    if (x_scale_logarithmic)
    {
      if (y_scale_logarithmic)
      {
        os << "loglog(x,y,'.-', 'MarkerSize', 10)" << std::endl;
      }
      else
      {
        os << "semilogx(x,y,'.-', 'MarkerSize', 10)" << std::endl;
      }
    }
    else
    {
      if (y_scale_logarithmic)
      {
        os << "semilogy(x,y,'.-', 'MarkerSize', 10)" << std::endl;
      }
      else
      {
        os << "plot(x,y,'.-', 'MarkerSize', 10)" << std::endl;
      }
    }

    os << "xlabel('" << xAxisLabel_ << "')" << std::endl;
    os << "ylabel('" << yAxisLabel_ << "')" << std::endl;
    os << "title('" << logTitle_ << "')" << std::endl;
  }



  inline
  TimeLog::TimeLog(const std::string& logTitle,
                   const std::string& y_axis_label,
                   bool time_scale_logarithmic,
                   bool y_scale_logarithmic)
    : DataLog(logTitle, "CPU time (in seconds)", y_axis_label, time_scale_logarithmic, y_scale_logarithmic)
  {

  }


}
