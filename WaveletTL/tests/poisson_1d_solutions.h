#ifndef _POISSON_1D_SOL
#define _POISSON_1D_SOL

#include <utils/function.h>

using namespace MathTL;

namespace WaveletTL
{
  // exact solution of -u''=1, u(0)=u(1)=0
  class Solution1 : public Function<1> {
  public:
    inline double value(const Point<1>& p, const unsigned int component = 0) const {
      return 0.5*p[0]*(1-p[0]);
    }
  
    void vector_value(const Point<1> &p, Vector<double>& values) const {
      values.resize(1, false);
      values[0] = value(p);
    }
  };

  // exact solution of -u''=f, u(0)=u(1)=0 with kink at x=0.5
  class Solution2 : public Function<1> {
  public:
    inline double value(const Point<1>& p, const unsigned int component = 0) const {
      if (0. <= p[0] && p[0] < 0.5)
	return -sin(3.*M_PI*p[0]) + 2.*p[0]*p[0];
    
      if (0.5 <= p[0] && p[0] <= 1.0)
	return -sin(3.*M_PI*p[0]) + 2.*(1-p[0])*(1-p[0]);
    
      return 0.;
    }
  
    void vector_value(const Point<1> &p, Vector<double>& values) const {
      values.resize(1, false);
      values[0] = value(p);
    }
  };

  // smooth part of corresponding right-hand side
  class RHS2_part : public Function<1> {
    inline double value(const Point<1>& p, const unsigned int component = 0) const {
      return -sin(3.*M_PI*p[0])*9.*M_PI*M_PI - 4.;
    }
  
    void vector_value(const Point<1> &p, Vector<double>& values) const {
      values.resize(1, false);
      values[0] = value(p);
    }
  };

  // smooth part of solution 2
  class Solution3 : public Function<1> {
  public:
    inline double value(const Point<1>& p, const unsigned int component = 0) const {
      return -sin(3.*M_PI*p[0]);
    }
  
    void vector_value(const Point<1> &p, Vector<double>& values) const {
      values.resize(1, false);
      values[0] = value(p);
    }
  };

  // smooth part of corresponding right-hand side
  class RHS3 : public Function<1> {
    inline double value(const Point<1>& p, const unsigned int component = 0) const {
      return -sin(3.*M_PI*p[0])*9.*M_PI*M_PI;
    }
  
    void vector_value(const Point<1> &p, Vector<double>& values) const {
      values.resize(1, false);
      values[0] = value(p);
    }
  };

  // kink at 0<a<1
  class Solution4 : public Function<1> {
  public:
    Solution4(const double a = 0.5) : a_(a) {}
  
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
  class RHS4_part : public Function<1> {
  public:
    RHS4_part(const double a = 0.5) : a_(a) {}
  
    inline double value(const Point<1>& p, const unsigned int component = 0) const {
      if (0. <= p[0] && p[0] < a_)
	return -1/(a_*a_);
    
      if (a_ <= p[0] && p[0] <= 1.0)
	return -1/((1-a_)*(1-a_));
    
      return 0.;
    }
  
    void vector_value(const Point<1> &p, Vector<double>& values) const {
      values.resize(1, false);
      values[0] = value(p);
    }
  
  protected:
    double a_;
  };

  // exact solution of -u''=1, u(0)=0, u'(1)=0
  class Solution5 : public Function<1> {
  public:
    inline double value(const Point<1>& p, const unsigned int component = 0) const {
      return 0.5*p[0]*(2-p[0]);
    }
  
    void vector_value(const Point<1> &p, Vector<double>& values) const {
      values.resize(1, false);
      values[0] = value(p);
    }
  };

  // exact solution of -u''=1, u'(0)=0, u(1)=0
  class Solution6 : public Function<1> {
  public:
    inline double value(const Point<1>& p, const unsigned int component = 0) const {
      return 0.5*(1-p[0]*p[0]);
    }
  
    void vector_value(const Point<1> &p, Vector<double>& values) const {
      values.resize(1, false);
      values[0] = value(p);
    }
  };

  // hat function
  class Hat : public Function<1>
  {
  public:
    inline double value(const Point<1>& p,
			const unsigned int component = 0) const
    {
      const double x = p[0];
      if (0. <= x && x < 0.5)
	return x;
      if (0.5 <= x && x <= 1.0)
	return 1-x;
      return 0;
    }
    
    void vector_value(const Point<1> &p,
		      Vector<double>& values) const
    {
      values.resize(1, false);
      values[0] = value(p);
    }
  };

  // "left" hat function
  class LeftHat : public Function<1>
  {
  public:
    inline double value(const Point<1>& p,
			const unsigned int component = 0) const
    {
      const double x = p[0];
      if (x >= 0 && x <= 1) {
	if (x >= 0.5)
	  return 0.5;
	else
	  return x;
      }
      return 0;
    }
    
    void vector_value(const Point<1> &p,
		      Vector<double>& values) const
    {
      values.resize(1, false);
      values[0] = value(p);
    }
  };

  // "right" hat function
  class RightHat : public Function<1>
  {
  public:
    inline double value(const Point<1>& p,
			const unsigned int component = 0) const
    {
      const double x = p[0];
      if (x >= 0 && x <= 1) {
	if (x >= 0.5)
	  return 1-x;
	else
	  return 0.5;
      }
      return 0;
    }
    
    void vector_value(const Point<1> &p,
		      Vector<double>& values) const
    {
      values.resize(1, false);
      values[0] = value(p);
    }
  };

  

}

#endif
