// implementation for ring_chart.h

#include <cmath>
#include <sstream>

namespace MathTL
{
  RingChart::RingChart(const double r0, const double r1)
    : r0_(r0), r1_(r1)
  {
  }
  
  RingChart::~RingChart()
  {
  }

  void
  RingChart::map_point(const Point<2>& x, Point<2>& y) const
  {
    const double arg = 2*M_PI*x[0];
    const double rofs = r0_+x[1]*(r1_-r0_);
    y[0] = rofs*cos(arg);
    y[1] = rofs*sin(arg);
  }

  void
  RingChart::map_point_inv(const Point<2>& x, Point<2>& y) const
  {
    y[0] = atan2(x[1],x[0])/(2*M_PI);
    if (y[0] < 0) y[0] += 1;
    y[1] = (sqrt(x[0]*x[0]+x[1]*x[1])-r0_)/(r1_-r0_);
  }
  
  const double
  inline
  RingChart::Gram_factor(const Point<2>& x) const
  {
    return sqrt(2*M_PI*(r1_-r0_)*(r0_+x[1]*(r1_-r0_)));
  }
  
  const double
  RingChart::Gram_D_factor(const unsigned int i,
			   const Point<2>& x) const
  {
    if (i==1) {
      return pow(r1_-r0_,1.5)/(M_2_SQRTPI*sqrt(r0_+x[1]*(r1_-r0_)));
    }
    return 0; // det Dkappa does not depend on phi
  }

  const double
  RingChart::Dkappa_inv(const unsigned int i,
			const unsigned int j,
			const Point<2>& x) const
  {
    // todo
    return 0;
  }

  const bool
  RingChart::in_patch(const Point<2>& x) const
  {
    // todo
    return true;
  }

  const string
  RingChart::to_string() const
  {
    std::stringstream strs;
 
    strs << "RingChart: "
	 << "r0=" << r0_
	 << ", r1=" << r1_
	 << endl;

    return strs.str();
  }
  
}
