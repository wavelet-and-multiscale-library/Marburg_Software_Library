// implementation for plot_tools.h

#include <cassert>
#include <geometry/sampled_mapping.h>

namespace MathTL
{
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

  void get_color(const double x,
		 const MatlabColorMap colormap,
		 double& red, double& green, double& blue)
  {
    assert( x >= -1 && x <= 1);

    switch(colormap) {
    case gray:
      red = green = blue = (x+1.0)/2.0;
      break;
    default:
      break;
    }
  }
}
