// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_PLOT_TOOLS_H
#define _MATHTL_PLOT_TOOLS_H

#include <iostream>
#include <map>

using std::map;

namespace MathTL {

  void matlab_output(const map<double,double>& log_10_residual_norms,
		     std::ostream& os);

  /*!
    enum type for the builtin Matlab color maps
  */
  enum MatlabColorMap {
//     autumn,
//     bone,
//     colorcube,
//     cool,
//     copper,
//     flag,
    gray,
//     hot,
//     hsv,
    jet,
//     lines,
//     pink,
//     prism,
//     spring,
//     summer,
//     white,
//     winter
  };
  
  /*!
    for a given x in [-1,1] and a given Matlab colormap,
    compute the corresponding RGB values
  */
  void get_color(const double x,
		 const MatlabColorMap colormap,
		 double& red, double& green, double& blue);
}

#include <utils/plot_tools.cpp>

#endif
