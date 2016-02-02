// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_CORNER_SINGULARITY_BIHARMONIC_H
#define _MATHTL_CORNER_SINGULARITY_BIHARMONIC_H

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
  class CornerSingularityBiharmonic
    : public Function<2>
  {
  public:
    /*!
      constructor from a corner x0,
      a starting angle theta0 (times pi, against positive x-axis)
      and an inner angle omega (times pi);
    */
    CornerSingularityBiharmonic(const Point<2>& x0,
		      const double theta0,
		      const double omega,
		      const double r0 = 0.01,
		      const double r1 = 0.99);

    //! destructor
    virtual ~CornerSingularityBiharmonic() {}

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
  class CornerSingularityBiharmonicRHS
    : public Function<2>
  {
  public:
    /*!
      constructor from a corner x0,
      a starting angle theta0 (times pi, against positive x-axis)
      and an inner angle omega (times pi);
    */
    CornerSingularityBiharmonicRHS(const Point<2>& x0,
			 const double theta0,
			 const double omega,
			 const double r0 = 0.01,
			 const double r1 = 0.99);
    
    //! destructor
    virtual ~CornerSingularityBiharmonicRHS() {}

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
 
    //! third derivative of the cutoff function
    double zeta_third_der(const double r) const;

    //! fourth derivative of the cutoff function
    double zeta_fourth_der(const double r) const;
 };

}

// implementations of inline functions
#include "numerics/corner_singularity_biharmonic.cpp"

#endif
