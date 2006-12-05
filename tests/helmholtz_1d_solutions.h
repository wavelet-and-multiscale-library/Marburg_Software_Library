#ifndef _HELMHOLTZ_1D_SOL
#define _HELMHOLTZ_1D_SOL

#include <utils/function.h>

using namespace MathTL;

namespace WaveletTL
{
  
  // a smooth solution of -u''+u=f, u(0)=u(1)=0
  class Solution1 : public Function<1> {
  public:
    inline double value(const Point<1>& p, const unsigned int component = 0) const {
      return M_SQRT2*sin(M_PI*p[0]);
    }
  
    void vector_value(const Point<1> &p, Vector<double>& values) const {
      values.resize(1, false);
      values[0] = value(p);
    }
  };

  // right-hand side for Solution1
  class RHS1 : public Function<1> {
  public:
    inline double value(const Point<1>& p, const unsigned int component = 0) const {
      return M_SQRT2*(1+M_PI*M_PI)*sin(M_PI*p[0]);
    }
  
    void vector_value(const Point<1> &p, Vector<double>& values) const {
      values.resize(1, false);
      values[0] = value(p);
    }
  };

  // polynomial solution of -u''=f, u(0)=u(1)=0
  class Solution2 : public Function<1> {
  public:
    inline double value(const Point<1>& p, const unsigned int component = 0) const {
      return p[0]*(1-p[0]);
    }
  
    void vector_value(const Point<1> &p, Vector<double>& values) const {
      values.resize(1, false);
      values[0] = value(p);
    }
  };

  // right-hand side for Solution2
  class RHS2 : public Function<1> {
    inline double value(const Point<1>& p, const unsigned int component = 0) const {
      return 2+p[0]*(1-p[0]);
    }
  
    void vector_value(const Point<1> &p, Vector<double>& values) const {
      values.resize(1, false);
      values[0] = value(p);
    }
  };

  // kink at 0<a<1
  class Solution3 : public Function<1> {
  public:
    Solution3(const double a = 0.5) : a_(a) {}

    inline double value(const Point<1>& p, const unsigned int component = 0) const {
      if (0. <= p[0] && p[0] < a_)
	return 1/(2*a_*a_)*p[0]*p[0];
 
      if (a_ <= p[0] && p[0] <= 1.0)
	return 0.5*(1-(p[0]-a_)/(1-a_))*(1-(p[0]-a_)/(1-a_));
 
      return 0.;
    }

    void vector_value(const Point<1> &p, Vector<double>& values) const {
      values.resize(1, false);
      values[0] = value(p);
    }

  protected:
    double a_;
  };

  // smooth part of corresponding right-hand side
  class RHS3_part : public Function<1> {
  public:
    RHS3_part(const double a = 0.5) : a_(a) {}

    inline double value(const Point<1>& p, const unsigned int component = 0) const {
      if (0. <= p[0] && p[0] < a_)
	return -1/(a_*a_) + 1/(2*a_*a_)*p[0]*p[0];
 
      if (a_ <= p[0] && p[0] <= 1.0)
	return -1/((1-a_)*(1-a_)) + 0.5*(1-(p[0]-a_)/(1-a_))*(1-(p[0]-a_)/(1-a_));
 
      return 0.;
    }

    void vector_value(const Point<1> &p, Vector<double>& values) const {
      values.resize(1, false);
      values[0] = value(p);
    }

  protected:
    double a_;
  };

}

#endif
