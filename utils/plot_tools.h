// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_PLOT_TOOLS_H
#define _MATHTL_PLOT_TOOLS_H

#include<map>
#include<geometry/sampled_mapping.h>

using std::map;

namespace MathTL {

  void matlab_output(const map<double,double>& log_10_residual_norms,
		     std::ostream& os)
  {
    Array1D<double> x(log_10_residual_norms.size());
    Array1D<double> values(log_10_residual_norms.size());
    
    map<double,double>::const_iterator it = log_10_residual_norms.begin();
    unsigned int i = 0;
    for (; it != log_10_residual_norms.end(); it++) {
      x[i] = it->first;
      values[i] = it->second;
      i++;
    }

    Grid<1> grid(x);
    SampledMapping<1> s_out(grid,values);
    s_out.matlab_output(os);
    
  }
}
#endif
