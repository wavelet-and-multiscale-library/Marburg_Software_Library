// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_SCHOENBERG_SPLINES_H
#define _MATHTL_SCHOENBERG_SPLINES_H

#include <cmath>
#include <utils/function.h>
#include <numerics/splines.h>
#include <numerics/cardinal_splines.h>

namespace MathTL
{
  /*!
    The following routines provide the pointwise evaluation of
    d-th order Schoenberg B-splines N_{k,d}, k >= -d+1,
    which use the special infinite knot sequence

      t^0_{-d+1} = ... = t^0_0 = 0          (knot with multiplicity d at x=0)
      t^0_k = k, k >= 1

    i.e., t^0_k = max(0,k) for k >= -d+1.
  */
  template <int d>
  class SchoenbergKnotSequence
    : public KnotSequence
  {
  public:
    /*
      index of the first knot
    */
    int k0() const { return -d+1; }
    
    /*!
      compute the k-th knot
    */
    double knot(const int k) const { return std::max(0,k); }
  };
  
  /*!
    evaluate an arbitrary Schoenberg B-spline N_{k,d}(x)
  */
  template <int d>
  double EvaluateSchoenbergBSpline(const int k, const double x)
  {
    double r(0);

    // take care of the multiple knot t_{-d+1} = ... = t_0
#if 0
    const double diff1 = std::max(0,k+d-1) - std::max(0,k);
    if (diff1 > 0) r += (x - std::max(0,k)) * EvaluateSchoenbergBSpline<d-1>(k, x) / diff1;

    const double diff2 = std::max(0,k+d) - std::max(0,k+1);
    if (diff2 > 0) r += (std::max(0,k+d) - x) * EvaluateSchoenbergBSpline<d-1>(k+1, x) / diff2;
#else
    if (k <= 0)
      {
	const double diff1 = std::max(0,k+d-1);
	if (diff1 > 0) r += x * EvaluateSchoenbergBSpline<d-1>(k, x) / diff1;
	
	const double diff2 = std::max(0,k+d) - std::max(0,k+1);
	if (diff2 > 0) r += (std::max(0,k+d) - x) * EvaluateSchoenbergBSpline<d-1>(k+1, x) / diff2;
      }
    else
      return ((x-k) * EvaluateSchoenbergBSpline<d-1>(k, x)
	      +(k+d-x) * EvaluateSchoenbergBSpline<d-1>(k+1, x)) / (d - 1);
#endif
    
    return r;
  }
  
  /*!
    evaluate an arbitrary Schoenberg B-spline N_{k,1}(x) = \chi_{[t^0_k,t^0_{k+1})}
  */
  template <>
  inline
  double EvaluateSchoenbergBSpline<1>(const int k, const double x)
  {
#if 0
    return (x >= std::max(0,k) && x < std::max(0,k+1) ? 1.0 : 0.0);
#else
    if (k <= 0)
      return (x >= 0 && x < std::max(0,k+1) ? 1.0 : 0.0);
    else
      return (x >= k && x < k+1 ? 1.0 : 0.0);
#endif
  }

  /*!
    evaluate an arbitrary Schoenberg B-spline N_{k,2}(x) = \chi_{[0,infty)}(x) * N_2(x)
  */
  template <>
  inline
  double EvaluateSchoenbergBSpline<2>(const int k, const double x)
  {
    return (x < 0 ? 0 : EvaluateCardinalBSpline<2>(k, x));
  }

  /*!
    evaluate a primal [P] function
      phi_{j,k}(x) = 2^{j/2}N_{k-d/2,d}(2^jx)
  */
  template <int d>
  inline
  double EvaluateSchoenbergBSpline_td(const int j, const int k, const double x)
  {
#if 0
    const double factor(1 << j);
    return sqrt(factor) * EvaluateSchoenbergBSpline<d>(k-(d/2), factor * x);
#else
    return twotothejhalf(j) * EvaluateSchoenbergBSpline<d>(k-(d/2), (1<<j) * x);
#endif
  }
  /*!
    evaluate the first derivative N_{k,d}'(x) of an arbitrary Schoenberg B-spline
  */
  template <int d>
  inline
  double EvaluateSchoenbergBSpline_x(const int k, const double x)
  {
    double r(0);

    if (k == 1-d) {
      r = (1-d) * EvaluateSchoenbergBSpline<d-1>(2-d, x);
    } else {
      if (k >= 0)
	r = EvaluateSchoenbergBSpline<d-1>(k, x) - EvaluateSchoenbergBSpline<d-1>(k+1, x);
      else
	r = (d-1) * (EvaluateSchoenbergBSpline<d-1>(k, x) / (k+d-1)
		     - EvaluateSchoenbergBSpline<d-1>(k+1, x) / (k+d));
    }

    return r;
  }
  
  /*!
    evaluate the first derivative N_{k,1}'(x) of an arbitrary Schoenberg B-spline
  */
  template <>
  inline
  double EvaluateSchoenbergBSpline_x<1>(const int k, const double x)
  {
    return 0.;
  }

  /*!
    evaluate the second derivative N_{k,d}''(x) of an arbitrary Schoenberg B-spline
  */
  template <int d>
  inline
  double EvaluateSchoenbergBSpline_xx(const int k, const double x)
  {
    double r(0);
    
    if (k == 1-d) {
      r = (1-d) * EvaluateSchoenbergBSpline_x<d-1>(2-d, x);
    } else {
      if (k >= 0)
	r = EvaluateSchoenbergBSpline_x<d-1>(k, x) - EvaluateSchoenbergBSpline_x<d-1>(k+1, x);
      else
	r = (d-1) * (EvaluateSchoenbergBSpline_x<d-1>(k, x) / (k+d-1)
		     - EvaluateSchoenbergBSpline_x<d-1>(k+1, x) / (k+d));
    }

    return r;
  }
  
  /*!
    evaluate the first derivative N_{k,1}'(x) of an arbitrary Schoenberg B-spline
  */
  template <>
  inline
  double EvaluateSchoenbergBSpline_xx<1>(const int k, const double x)
  {
    return 0.;
  }

  /*!
    evaluate the first derivative of a primal [P] function
      phi_{j,k}'(x) = 2^{3*j/2}N_{k-d/2,d}'(2^jx)
  */
  template <int d>
  inline
  double EvaluateSchoenbergBSpline_td_x(const int j, const int k, const double x)
  {
#if 0
    const double factor(1 << j);
    return factor * sqrt(factor) * EvaluateSchoenbergBSpline_x<d>(k-(d/2), factor * x);
#else
    return twotothejhalf(3*j) * EvaluateSchoenbergBSpline_x<d>(k-(d/2), (1<<j) * x);
#endif
  }

  /*!
    evaluate the second derivative of a primal [P] function
      phi_{j,k}''(x) = 2^{5*j/2}N_{k-d/2,d}''(2^jx)
  */
  template <int d>
  inline
  double EvaluateSchoenbergBSpline_td_xx(const int j, const int k, const double x)
  {
#if 0
    const double factor(1 << j);
    return  factor * factor * sqrt(factor) * EvaluateSchoenbergBSpline_xx<d>(k-(d/2), factor * x);
#else
    return twotothejhalf(5*j) * EvaluateSchoenbergBSpline_xx<d>(k-(d/2), (1<<j) * x);
#endif
  }


/*
   expand a Schoenbergspline as a Picewise Polynomial funktion 
   f = 1 falls es sich um einen Generator an der rechten Intervall Grenze handelt
   f = 0 sonst
*/
template <int d>
  Piecewise<double> ExpandSchoenbergBspline(const int j, const int k, const int f)
{
  Piecewise<double> r(j);
  Polynomial<double> p, q;
  if (f==0){
    p.set_coefficient(0, floor((double) d/2.0) -k);  //floor((double) d/2.0)
    p.set_coefficient(1, ldexp(1.0, j)); // p(x)=2^jx-k+floor(d/2)
  }
  else{
    p.set_coefficient(0, floor((double) d/2.0) + ldexp(1.0, j) -k);
    p.set_coefficient(1, -ldexp(1.0, j));
  }
  double factor = sqrt(ldexp(1.0, j));
  int m, n;
  SparseMatrix<double> coeffs(d, d);
  switch(d) // precalculated entries for 1<=d<=4
    {
    case 1:
      {
	if (k>=0) {coeffs.set_entry(0, 0, 1);} // k-floor((double) d/2.0) >=0
	else {coeffs.set_entry(0, 0, 0);}
      }
      break;
    case 2:
      {
        if (k>=1) {   // k-floor((double) d/2.0) >=0
	coeffs.set_entry(1, 0, 1);

	coeffs.set_entry(0, 1, 2);
	coeffs.set_entry(1, 1, -1);}
        else if (k>= 0) {  // k-floor((double) d/2.0) >=-1
	coeffs.set_entry(0, 1, 2);
	coeffs.set_entry(1, 1, -1);}
      }
      break;
    case 3:
      {
        if (k>=1) { // k-floor((double) d/2.0) >=0
	coeffs.set_entry(2, 0, 0.5);

	coeffs.set_entry(0, 1, -1.5);
	coeffs.set_entry(1, 1, 3);
	coeffs.set_entry(2, 1, -1);

	coeffs.set_entry(0, 2, 4.5);
	coeffs.set_entry(1, 2, -3);
	coeffs.set_entry(2, 2, 0.5);
	}
	else if (k>=0){ // k-floor((double) d/2.0) >=-1
	coeffs.set_entry(0, 1, -3.5);
	coeffs.set_entry(1, 1, 5.0);
	coeffs.set_entry(2, 1, -1.5);

	coeffs.set_entry(0, 2, 4.5);
	coeffs.set_entry(1, 2, -3);
	coeffs.set_entry(2, 2, 0.5);
	}
	else if (k>=-1){ // k-floor((double) d/2.0) >=-2
	coeffs.set_entry(0, 2, 9);
	coeffs.set_entry(1, 2, -6);
	coeffs.set_entry(2, 2, 1);}
      }
      break;
    case 4:
      {
	if (k>=2) { // k-floor((double) d/2.0) >=0
	coeffs.set_entry(3, 0, 1.0/6.0);

	coeffs.set_entry(0, 1, 2.0/3.0);
	coeffs.set_entry(1, 1, -2);
	coeffs.set_entry(2, 1, 2);
	coeffs.set_entry(3, 1, -0.5);

	coeffs.set_entry(0, 2, -22.0/3.0);
	coeffs.set_entry(1, 2, 10);
	coeffs.set_entry(2, 2, -4);
	coeffs.set_entry(3, 2, 0.5);

	coeffs.set_entry(0, 3, 32.0/3.0);
	coeffs.set_entry(1, 3, -8);
	coeffs.set_entry(2, 3, 2);
	coeffs.set_entry(3, 3, -1.0/6.0);
	}
	else if (k>=1){ // k-floor((double) d/2.0) >=-1
	coeffs.set_entry(0, 1, 29.0/12.0);
	coeffs.set_entry(1, 1, -5.75);
	coeffs.set_entry(2, 1, 4.25);
	coeffs.set_entry(3, 1, -5.5/6.0);

	coeffs.set_entry(0, 2, -115.0/12.0); 
	coeffs.set_entry(1, 2, 12.25); 
	coeffs.set_entry(2, 2, -4.75);
	coeffs.set_entry(3, 2, 7.0/12.0);

	coeffs.set_entry(0, 3, 32.0/3.0);
	coeffs.set_entry(1, 3, -8);
	coeffs.set_entry(2, 3, 2);
	coeffs.set_entry(3, 3, -1.0/6.0);
	}
	else if (k>=0){ // k-floor((double) d/2.0) >=-2
	coeffs.set_entry(0, 2, -38.0); 
	coeffs.set_entry(1, 2, 42.0); 
	coeffs.set_entry(2, 2, -15.0);
	coeffs.set_entry(3, 2, 21.0/12.0);

	coeffs.set_entry(0, 3, 16.0);
	coeffs.set_entry(1, 3, -12.0);
	coeffs.set_entry(2, 3, 3.0);
	coeffs.set_entry(3, 3, -0.25);
	}
	else if (k>=-1){ // k-floor((double) d/2.0) >=-3
	coeffs.set_entry(0, 3, 64.0);
	coeffs.set_entry(1, 3, -48.0);
	coeffs.set_entry(2, 3, 12.0);
	coeffs.set_entry(3, 3, -1.0);
	}
      }
      break;
    default:
    {
      return r;
      //coeffs = localBSplineCoeffs(d);
    }
}
  for (n=0; n<=d-1; n++)
    {
     // coeffs.findc(n, ind);
     // CoeffsType col;
      for (m=0; m<=d-1; m++)
	q.set_coefficient(m, coeffs.get_entry(m, n));
      if (f==0){
        r.set_local_expansion(n + k - floor((double) d/2.0), q.substitute_into(p));  //(int)floor((double) d/2.0)
      } 
      else {
        r.set_local_expansion((int)ldexp(1.0, j) - n - k + floor((double) d/2.0) -1, q.substitute_into(p));
      }
    }
  r *= factor;

  return r;
   
}





  /*!
    a translated and dilated Schoenberg B-spline as a function object
  */
  template <int d>
  class SchoenbergBSpline_td : public Function<1>
  {
  public:
    //! constructor from j, k
    SchoenbergBSpline_td(const int j, const int k)
      : j_(j), k_(k) {}
    
    //! virtual destructor
    virtual ~SchoenbergBSpline_td() {}

    //! point value
    inline double value(const Point<1>& p,
			const unsigned int component = 0) const
    {
      return EvaluateSchoenbergBSpline_td<d>(j_, k_, p[0]);
    }
  
    //! point value
    void vector_value(const Point<1> &p,
		      Vector<double>& values) const
    {
      values.resize(1, false);
      values[0] = value(p);
    }
   
    //! set j
    void set_j(const int j) { j_ = j; }
    
    //! set k
    void set_k(const int k) { k_ = k; }
    
  protected:
    int j_, k_;
  };

}

#endif
