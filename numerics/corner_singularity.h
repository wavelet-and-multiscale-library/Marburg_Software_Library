// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2008                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_CORNER_SINGULARITY_H
#define _MATHTL_CORNER_SINGULARITY_H

#include <geometry/point.h>
#include <utils/function.h>

namespace MathTL
{
  /*!
    Corner singularity function w.r.t. a reentrant corner at x0 in a
    polygonal domain in R^2. In polar coordinates w.r.t. the corner,
    the function looks like

      s(r,theta) = zeta(r)*r^{pi/omega}*sin((pi/omega)*theta),

    where zeta:[0,infty)->R is the smooth cutoff function

      zeta(r) = w(r1-r)/(w(r-r0)+w(r1-r))
      w(t) = exp(-1/t^2) for t>0 and 0 elsewhere

    This yields exactly the function S_{l,1} from [D].

    References:
    [D] Dahlke, Besov regularity for elliptic boundary value problems on
        polygonal domains, Appl. Math. Lett. 12(1999), 31-36
  */
  class CornerSingularity
    : public Function<2>
  {
  public:
    /*!
      constructor from a corner x0,
      a starting angle theta0 (times pi, against positive x-axis)
      and an inner angle omega (times pi);
    */
    CornerSingularity(const Point<2>& x0,
		      const double theta0,
		      const double omega,
		      const double r0 = 0.01,
		      const double r1 = 0.99);

    //! destructor
    virtual ~CornerSingularity() {}

    //! point value at x
    double value(const Point<2>& x, const unsigned int component = 0) const;

    //! vector-valued value at x (for compatibility with Function)
    void vector_value(const Point<2>& p, Vector<double>& values) const;
    
  protected:
    //! corner
    Point<2> x0;

    //! starting angle
    double theta0;

    //! inner angle
    double omega;

    //! cutoff parameters
    double r0, r1;

    //! the cutoff function
    double zeta(const double r) const;
  };

  /*!
    This class models
      -\Delta s(x),
    where s(x) is the corner singularity function from above.
  */
  class CornerSingularityRHS
    : public Function<2>
  {
  public:
    /*!
      constructor from a corner x0,
      a starting angle theta0 (times pi, against positive x-axis)
      and an inner angle omega (times pi);
    */
    CornerSingularityRHS(const Point<2>& x0,
			 const double theta0,
			 const double omega,
			 const double r0 = 0.01,
			 const double r1 = 0.99);
    
    //! destructor
    virtual ~CornerSingularityRHS() {}

    //! point value at x
    double value(const Point<2>& x, const unsigned int component = 0) const;

    //! vector-valued value at x (for compatibility with Function)
    void vector_value(const Point<2>& p, Vector<double>& values) const;
    
  protected:
    //! corner
    Point<2> x0;

    //! starting angle
    double theta0;

    //! inner angle
    double omega;

    //! cutoff parameters
    double r0, r1;

    //! the cutoff function
    double zeta(const double r) const;

    //! first derivative of the cutoff function
    double zeta_prime(const double r) const;

    //! second derivative of the cutoff function
    double zeta_primeprime(const double r) const;
  };

  /*!
    A time--dependent corner singularity function w.r.t. to the reentrant corner at x0 in a
    polygonal domain in R^2. In polar coordinates w.r.t. the corner,
    the function looks like

      s(r,theta) = t^{3/4}*r^{pi/omega}*sin((pi/omega)*theta)*(1-x^2)*(1-y^2)

    with x=r*cos(theta0+theta), y=r*sin(theta0+theta).

    This function was used in [GSS+] as a test function for the heat equation
    on the L--shaped domain in R^2.

    References:
    [GSS+]
  */
  class CornerTimeSingularity
    : public Function<2>
  {
  public:
    /*!
      constructor from a corner x0,
      a starting angle theta0 (times pi, against positive x-axis)
      and an inner angle omega (times pi);
    */
    CornerTimeSingularity(const Point<2>& x0,
			  const double theta0,
			  const double omega);

    //! destructor
    virtual ~CornerTimeSingularity() {}

    //! point value at x
    double value(const Point<2>& x, const unsigned int component = 0) const;

    //! vector-valued value at x (for compatibility with Function)
    void vector_value(const Point<2>& p, Vector<double>& values) const;
    
  protected:
    //! corner
    Point<2> x0;

    //! starting angle
    double theta0;

    //! inner angle
    double omega;
  };

  /*!
    u'(t,x,y)-Delta u(t,x,y) of CornerTimeSingularity
  */
  class CornerTimeSingularityRHS
    : public Function<2>
  {
  public:
    /*!
      constructor from a corner x0,
      a starting angle theta0 (times pi, against positive x-axis)
      and an inner angle omega (times pi);
    */
    CornerTimeSingularityRHS(const Point<2>& x0,
			     const double theta0,
			     const double omega);

    //! destructor
    virtual ~CornerTimeSingularityRHS() {}

    //! point value at x
    double value(const Point<2>& x, const unsigned int component = 0) const;

    //! vector-valued value at x (for compatibility with Function)
    void vector_value(const Point<2>& p, Vector<double>& values) const;
    
  protected:
    //! corner
    Point<2> x0;

    //! starting angle
    double theta0;

    //! inner angle
    double omega;
  };

  /*!
    temporal derivative of CornerTimeSingularityRHS
  */
  class CornerTimeSingularityRHSt
    : public Function<2>
  {
  public:
    /*!
      constructor from a corner x0,
      a starting angle theta0 (times pi, against positive x-axis)
      and an inner angle omega (times pi);
    */
    CornerTimeSingularityRHSt(const Point<2>& x0,
			      const double theta0,
			      const double omega);

    //! destructor
    virtual ~CornerTimeSingularityRHSt() {}
    
    //! point value at x
    double value(const Point<2>& x, const unsigned int component = 0) const;
    
    //! vector-valued value at x (for compatibility with Function)
    void vector_value(const Point<2>& p, Vector<double>& values) const;
    
  protected:
    //! corner
    Point<2> x0;

    //! starting angle
    double theta0;

    //! inner angle
    double omega;
  };

}

// implementations of inline functions
#include "numerics/corner_singularity.cpp"

#endif
