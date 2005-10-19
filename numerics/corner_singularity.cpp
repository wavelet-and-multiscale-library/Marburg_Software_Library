// implementation for corner_singularity.h

#include <cmath>
#include <utils/tiny_tools.h>

namespace MathTL
{
  CornerSingularity::CornerSingularity(const Point<2>& x,
				       const double w0,
				       const double w,
				       const double t0,
				       const double t1)
    : Function<2>(), x0(x), theta0(w0), omega(w), r0(t0), r1(t1)
  {
  }
  
  double
  CornerSingularity::value(const Point<2>& p,
			   const unsigned int component) const
  {
    const Point<2> x(p-x0);
    const double r = hypot(x[0],x[1]);

    if (r >= r1) return 0.0;

    double theta = atan2(x[1],x[0]);

    // shift theta to [0,2*pi]
    if (theta < 0) theta += 2.0 * M_PI;
    theta -= theta0 * M_PI;
    if (theta < 0) theta += 2.0 * M_PI;

    if (theta >= omega * M_PI) return 0.0;

    return zeta(r) * pow(r, 1.0/omega) * sin(theta/omega);
  }

  inline
  void
  CornerSingularity::vector_value(const Point<2> &p,
				  Vector<double>& values) const
  {
    values[0] = value(p);
  }

  double
  CornerSingularity::zeta(const double r) const {
    if (r <= r0)
      return 1.0;
    else {
      if (r >= r1)
	return 0.0;
      else {
	const double help1(r1-r);
	const double help0(r-r0);
	return exp(-1.0/(help1*help1))/(exp(-1.0/(help0*help0))+exp(-1.0/(help1*help1)));
      }
    }
  }
}
