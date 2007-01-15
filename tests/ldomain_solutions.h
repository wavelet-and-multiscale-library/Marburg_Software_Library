#ifndef _LDOMAIN_SOL
#define _LDOMAIN_SOL

#include <utils/function.h>

using namespace MathTL;

namespace WaveletTL
{
  /*
    A test problem for the Poisson equation on the L--shaped domain with homogeneous Dirichlet b.c.'s:
      -Delta u(x,y) = 2*pi^2*sin(pi*x)*sin(pi*y)
    with exact solution
      u(x,y) = sin(pi*x)*sin(pi*y)
  */
  class EigenRHS
    : public Function<2,double>
  {
  public:
    virtual ~EigenRHS() {};
    double value(const Point<2>& p, const unsigned int component = 0) const {
      return 2*M_PI*M_PI*sin(M_PI*p[0])*sin(M_PI*p[1]);
    }
    void vector_value(const Point<2>& p, Vector<double>& values) const {
      values[0] = value(p);
    }
  };
  
  class EigenSolution
    : public Function<2,double>
  {
  public:
    virtual ~EigenSolution() {};
    double value(const Point<2>& p, const unsigned int component = 0) const {
      return sin(M_PI*p[0])*sin(M_PI*p[1]);
    }
    void vector_value(const Point<2>& p, Vector<double>& values) const {
      values[0] = value(p);
    }
  };

  class EigenSum
    : public Function<2,double>
  {
  public:
    virtual ~EigenSum() {};
    double value(const Point<2>& p, const unsigned int component = 0) const {
      return (1.0+2*M_PI*M_PI)*sin(M_PI*p[0])*sin(M_PI*p[1]);
    }
    void vector_value(const Point<2>& p, Vector<double>& values) const {
      values[0] = value(p);
    }
  };

  /*!
    u(x,y) = x*(1-x)*(1+x)*y*(1-y)*(1+y)
    f(x,y) = -Delta u(x,y) = 12*x*y-6*x*y^3-6*x^3*y
  */
  class PolyRHS
    : public Function<2,double>
  {
  public:
    virtual ~PolyRHS() {};
    double value(const Point<2>& p, const unsigned int component = 0) const {
      return 12*p[0]*p[1]-6*p[0]*p[1]*p[1]*p[1]-6*p[0]*p[0]*p[0]*p[1];
    }
    void vector_value(const Point<2>& p, Vector<double>& values) const {
      values[0] = value(p);
    }
  };

  class PolySolution
    : public Function<2,double>
  {
  public:
    virtual ~PolySolution() {};
    double value(const Point<2>& p, const unsigned int component = 0) const {
      return p[0]*(1-p[0])*(1+p[0])*p[1]*(1-p[1])*(1+p[1]);
    }
    void vector_value(const Point<2>& p, Vector<double>& values) const {
      values[0] = value(p);
    }
  };

  class PolySum
    : public Function<2,double>
  {
  public:
    virtual ~PolySum() {};
    double value(const Point<2>& p, const unsigned int component = 0) const {
      return p[0]*(1-p[0])*(1+p[0])*p[1]*(1-p[1])*(1+p[1])
	+12*p[0]*p[1]-6*p[0]*p[1]*p[1]*p[1]-6*p[0]*p[0]*p[0]*p[1];
    }
    void vector_value(const Point<2>& p, Vector<double>& values) const {
      values[0] = value(p);
    }
  };
}

#endif
