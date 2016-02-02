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
//     case autumn:
//       break;
//     case bone:
//       break;
//     case colorcube:
//       break;
    case cool:
      red   = (x+1.0)/2.0;
      green = 1.0-red;
      blue  = 1.0;
      break;
//     case copper:
//       break;
//     case flag:
//       break;
    case gray:
      red = green = blue = (x+1.0)/2.0;
      break;
//     case hot:
//       break;
//     case hsv:
//       break;
    case jet:
      if (x <= -0.25)
	red = 0;
      else {
	if (x <= 0.25)
	  red = 2.0 * (x + 0.25);
	else {
	  if (x <= 0.75)
	    red = 1.0;
	  else
	    red = 2.5 - 2.0 * x;
	}
      }
      if (x <= -0.75)
	green = 0;
      else {
	if (x <= -0.25)
	  green = 2.0 * (x + 0.75);
	else {
	  if (x <= 0.25)
	    green = 1.0;
	  else {
	    if (x <= 0.75)
	      green = 1.5 - 2.0 * x;
	    else
	      green = 0;
	  }
	}
      }
      if (x <= -0.75)
	blue = 2.5 + 2.0 * x;
      else {
	if (x <= -0.25)
	  blue = 1.0;
	else {
	  if (x <= 0.25)
	    blue = 0.5 - 2.0 * x;
	  else
	    blue = 0;
	}
      }
      break;
//     case lines:
//       break;
//     case pink:
//       break;
//     case prism:
//       break;
//     case spring:
//       break;
//     case summer:
//       break;
//     case white:
//       break;
//     case winter:
//       break;
    default:
      break;
    }
  }
}
