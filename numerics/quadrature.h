// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_QUADRATURE_H
#define _MATHTL_QUADRATURE_H

#include <utils/array1d.h>
#include <utils/function.h>
#include <geometry/point.h>

namespace MathTL
{
  /*!
    A base class for N-point quadrature rules in DIM space dimensions
    on a box spanned by a,b\in\mathbb R^{DIM}
      Q(f) = \sum_{k=1}^N w_k f(x_k)
    with nodes x_k and weights w_k.
    Internally, we store the nodes and weights for quadrature on [0,1]^{DIM}
    and scale the nodes appropriately when quadrature has to be performed

    (this class is _not_ intended to be used in practice,
    please use derived special versions below instead!)
   */
  template <unsigned int DIM>
  class QuadratureRule
  {
  public:
    /*!
      typedef for lower-dimensional quadrature rule
    */
    typedef QuadratureRule<DIM-1> SubQuadratureRule;

    /*!
      default constructor: yields empty quadrature rule
    */
    QuadratureRule();

    /*!
      copy constructor: copy points and weights
    */
    QuadratureRule(const QuadratureRule<DIM>& Q);

    /*!
      construct quadrature rule as a tensor product from
      a lower dimensional one and a one-dimensional one
    */
    QuadratureRule(const SubQuadratureRule& Q,
 		   const QuadratureRule<1>& Q1);
    
    /*!
      virtual destructor
    */
    virtual ~QuadratureRule();

    /*!
      return number of quadrature points
    */
    unsigned int get_N() const;

    /*!
      return all quadrature points
    */
    void get_points(Array1D<Point<DIM> >& points) const;

    /*!
      return quadrature weights
    */
    void get_weights(Array1D<double>& weights) const;

    /*!
      evaluate quadrature rule on [0,1]^{DIM}
      (we assume that the function is real-valued)
    */
    double integrate(const Function<DIM, double>& f) const;

    /*!
      Evaluate quadrature rule on [a,b] where a,b are points in \mathbb R^d.
      We do so by rescaling the integral to one over [0,1]^{DIM}
      (we assume that the function is real-valued)
    */
    double integrate(const Function<DIM, double>& f,
		     const Point<DIM>& a, const Point<DIM>& b) const;

  protected:
    /*!
      nodes of the quadrature rule
    */
    Array1D<Point<DIM> > points_;

    /*!
      weights of the quadrature rule
    */
    Array1D<double> weights_;
  };

  /*!
    1D midpoint rule Q(f) = (b-a) * f((a+b)/2)
  */
  class MidpointRule : public QuadratureRule<1>
  {
  public:
    /*!
      construct midpoint rule
    */
    MidpointRule();
  };

  /*!
    closed 1D Newton-Cotes rule with N subintervals, i.e., N+1 points
  */
  class ClosedNewtonCotesRule : public QuadratureRule<1>
  {
  public:
    /*!
      construct closed Newton Cotes rule with N subintervals
    */
    ClosedNewtonCotesRule(const unsigned int N);
  };

  /*!
    trapezoidal rule Q(f) = (b-a)/2 * (f(a)+f(b))
    (equivalent to ClosedNewtonCotesRule<1>)
   */
  class TrapezoidalRule: public QuadratureRule<1>
  {
  public:
    /*!
      construct trapezoidal rule
    */
    TrapezoidalRule();
  };

  /*!
    Simpson rule Q(f) = (b-a)/6 * (f(a)+4*f(a+b)/2+f(b))
    (equivalent to ClosedNewtonCotesRule<2>)
   */
  class SimpsonRule: public QuadratureRule<1>
  {
  public:
    /*!
      construct Simpson rule
    */
    SimpsonRule();
  };

  /*!
    construct composite quadrature rule from a simple one
   */
  template <unsigned int DIM>
  class CompositeRule
    : public QuadratureRule<DIM>
  {
  public:
    /*!
      construct composite quadrature rule on [0,1]^{DIM} with
      N subintervals from a simple quadrature rule
    */
    CompositeRule(const QuadratureRule<1>& Q,
		  const unsigned int N = 1);
    
  protected:
    /*!
      basic 1D quadrature rule
    */
    const QuadratureRule<1> Q_;

    /*!
      number of subintervals
    */
    unsigned int N_;

  private:
    /*!
      check whether a given 1D quadrature rule uses both endpoints
      (needed to glue N quadrature rules neatly together,
      this is only called for DIM==1)
    */
    static bool uses_both_endpoints(const QuadratureRule<1>& Q);
  };
}


// include implementation of inline functions and specializations
#include <numerics/quadrature.cpp>

#endif
