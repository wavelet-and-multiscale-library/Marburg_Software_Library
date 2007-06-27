// implementation for biharmonic_1d_testcase.h

#include <cmath>

using namespace std;


namespace FrameTL
{
  //class Biharmonic1D_RHS_Integrand; // forward declaration for the compiler


  // members of class Biharmonic1D_Solution
  double
  Biharmonic1D_Solution::value(const Point<1>& p, const unsigned int component) const
  {
    const double x = p(0);
    double y;

    y = -cos(2*M_PI*x)+1.0; // smooth part
    if (x <= 0.5)
      y += 48.0*x*x*x*x - 63.0*x*x*x + 47.0/2.0*x*x;
    else
      y += 48.0*x*x*x*x - 127.0*x*x*x + 235.0/2.0*x*x - 46.0*x + 15.0/2.0;

    return y;
  }


  // members of Biharmonic1D_RHS_Integrand
  double
  Biharmonic1D_RHS_Integrand::value(const Point<1>& p, const unsigned int component) const
  {
    const double x = p(0); // the x coordinate

    return cos(2.0*M_PI*x)*M_PI*M_PI*M_PI*M_PI - 72.0;
  }


  // members of class Biharmonic1D_RHS
  template<class IBASIS>
  double
  Biharmonic1D_RHS<IBASIS>::evaluate(const typename AggregatedFrame<IBASIS, 1>::Index& lambda) const
  {
    double result;

    result = -16.0 * Functional<IBASIS,1,1>::evaluate(lambda); // use super class for computing the integral part

    result += 4.0 * Functional<IBASIS,1,1>::frame_->evaluate(true, 1, lambda, Point<1>(0.5)); // first derivative on 1/2
    result -= 384.0 * Functional<IBASIS,1,1>::frame_->evaluate(true, 0, lambda, Point<1>(0.5)); // frame evaluated on 1/2

    return result;
  }
}
