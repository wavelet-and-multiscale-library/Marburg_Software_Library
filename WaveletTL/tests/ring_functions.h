#ifndef _RING_FUNCTIONS
#define _RING_FUNCTIONS

#include <cmath>
#include <utils/function.h>

using MathTL::Function;
using std::max;

namespace WaveletTL
{
  // the constant function f(x)=23, divided by sqrt(2*Pi*(r1-r0)*r)
  class RingFunction1
    : public Function<2>
  {
  public:
    RingFunction1(const double r0, const double r1) : r0_(r0), r1_(r1) {}
    inline double value(const Point<2>& p, const unsigned int component = 0) const {
      const double r = sqrt(p[0]*p[0]+p[1]*p[1]);
      return 23/sqrt(2*M_PI*(r1_-r0_)*r);
    }
    void vector_value(const Point<2> &p, Vector<double>& values) const {
      values.resize(1, false);
      values[0] = value(p);
    }
  private:
    double r0_, r1_;
  };

  // a linear polynomial in r, divided by sqrt(2*Pi*(r1-r0)*r)
  class RingFunction2
    : public Function<2>
  {
  public:
    RingFunction2(const double r0, const double r1) : r0_(r0), r1_(r1) {}
    inline double value(const Point<2>& p, const unsigned int component = 0) const {
      const double r = sqrt(p[0]*p[0]+p[1]*p[1]);
      return (2*r-1)/sqrt(2*M_PI*(r1_-r0_)*r);
    }
    void vector_value(const Point<2> &p, Vector<double>& values) const {
      values.resize(1, false);
      values[0] = value(p);
    }
  private:
    double r0_, r1_;
  };

  // the hat function in r, divided by sqrt(2*Pi*(r1-r0)*r)
  class RingFunction3
    : public Function<2>
  {
  public:
    RingFunction3(const double r0, const double r1) : r0_(r0), r1_(r1) {}
    inline double value(const Point<2>& p, const unsigned int component = 0) const {
      const double r = sqrt(p[0]*p[0]+p[1]*p[1]);
      const double s = (r-r0_)/(r1_-r0_);
      return max(1-2*fabs(s-0.5),0.)/sqrt(2*M_PI*(r1_-r0_)*r);
    }
    void vector_value(const Point<2> &p, Vector<double>& values) const {
      values.resize(1, false);
      values[0] = value(p);
    }
  private:
    double r0_, r1_;
  };

  // a quadratic polynomial in r, divided by sqrt(2*Pi*(r1-r0)*r)
  class RingFunction4
    : public Function<2>
  {
  public:
    RingFunction4(const double r0, const double r1) : r0_(r0), r1_(r1) {}
    inline double value(const Point<2>& p, const unsigned int component = 0) const {
      const double r = sqrt(p[0]*p[0]+p[1]*p[1]);
      return (r1_-r)*(r-r0_)/sqrt(2*M_PI*(r1_-r0_)*r);
    }
    void vector_value(const Point<2> &p, Vector<double>& values) const {
      values.resize(1, false);
      values[0] = value(p);
    }
  private:
    double r0_, r1_;
  };

  // minus Laplace of RingFunction4
  class RingFunction4_RHS
    : public Function<2>
  {
  public:
    RingFunction4_RHS(const double r0, const double r1) : r0_(r0), r1_(r1) {}
    inline double value(const Point<2>& p, const unsigned int component = 0) const {
      const double r = sqrt(p[0]*p[0]+p[1]*p[1]);
      return -(-9*r*r+(r0_+r1_)*r-r0_*r1_)/(4*r*r*sqrt(2*M_PI*(r1_-r0_)*r));
    }
    void vector_value(const Point<2> &p, Vector<double>& values) const {
      values.resize(1, false);
      values[0] = value(p);
    }
  private:
    double r0_, r1_;
  };
  
  // minus Laplace plus Identity of RingFunction4
  class RingFunction4_RHS_Helmholtz
    : public Function<2>
  {
  public:
    RingFunction4_RHS_Helmholtz(const double r0, const double r1) : r0_(r0), r1_(r1) {}
    inline double value(const Point<2>& p, const unsigned int component = 0) const {
      const double r = sqrt(p[0]*p[0]+p[1]*p[1]);
      return -(-9*r*r+(r0_+r1_)*r-r0_*r1_)/(4*r*r*sqrt(2*M_PI*(r1_-r0_)*r))
	+ (r1_-r)*(r-r0_)/sqrt(2*M_PI*(r1_-r0_)*r);
    }
    void vector_value(const Point<2> &p, Vector<double>& values) const {
      values.resize(1, false);
      values[0] = value(p);
    }
  private:
    double r0_, r1_;
  };
}

#endif
