// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_CARDINAL_SPLINES_H
#define _MATHTL_CARDINAL_SPLINES_H

#include <cmath>
#include <utils/function.h>
#include <numerics/splines.h>
#include <utils/tiny_tools.h>
#include <algebra/piecewise.h>
#include <algebra/polynomial.h>
#include <algebra/sparse_matrix.h>

namespace MathTL
{
  /*!
    (semi)infinite knot sequence for cardinal splines supported in [0,\infty)
  */
  template <int d>
  class CardinalKnotSequence
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
    double knot(const int k) const { return k; }
  };

  /*!
    evaluate a cardinal B-spline N_d(x) via recursion
    (fastest possible variant)
  */
  template <int d>
  inline
  double EvaluateCardinalBSpline(const double x)
  {
    return (x < 0 || x >= d ? 0
	    : (x * EvaluateCardinalBSpline<d-1>(x)
	       + (d-x) * EvaluateCardinalBSpline<d-1>(x-1)) / (d-1));

  }

  /*!
    evaluate a cardinal B-spline N_3(x)
    (fastest possible variant)
  */
  template <>
  double EvaluateCardinalBSpline<3>(const double x)
  {
    if (x < 0 || x >= 3)
      return 0;
    if (x < 1)
      return x*x/2;
    else
      return (x < 2 ? x*(3-x)-1.5 : (x*(x-6)+9)/2);
  }
  
  /*!
    evaluate a cardinal B-spline N_2(x)
    (fastest possible variant)

  */
  template <>
  inline
  double EvaluateCardinalBSpline<2>(const double x)
  {
    return (x < 0 || x >= 2 ? 0 : 1-abs(x-1));
  }
  
  /*!
    evaluate a cardinal B-spline N_1(x)
    (fastest possible variant)
  */
  template <>
  inline
  double EvaluateCardinalBSpline<1>(const double x)
  {
    return (x < 0 || x >= 1 ? 0 : 1);
  }

  /*!
    evaluate a shifted cardinal B-spline N_d(x-k) via recursion
  */
  template <int d>
  inline
  double EvaluateCardinalBSpline(const int k, const double x)
  {
    if (x < k || x >= k+d)
      return 0.;
    return ((x-k) * EvaluateCardinalBSpline<d-1>(k, x)
 	    + (k+d-x) * EvaluateCardinalBSpline<d-1>(k+1, x)) / (d-1);
  }
  
  /*!
    evaluate a shifted cardinal B-spline N_1(x-k)
  */
  template <>
  inline
  double EvaluateCardinalBSpline<1>(const int k, const double x)
  {
    if (x < k || x >= k+1)
      return 0.0;
    return 1.0;
  }

  /*!
    evaluate a primal CDF function
      phi_{j,k}(x) = 2^{j/2}N_d(2^jx-k+d/2)
  */
  template <int d>
  inline
  double EvaluateCardinalBSpline_td(const int j, const int k, const double x)
  {
#if 0
    return sqrt(factor) * EvaluateCardinalBSpline<d>(pow(2,j) * x - k + d/2); Funktioniert auch für j<0
#else
    return twotothejhalf(j) * EvaluateCardinalBSpline<d>((1<<j) * x - k + d/2);
#endif
  }
  
  /*!
    evaluate the first derivative N_d'(x) of a cardinal B-spline
    (fastest possible variant)
  */
  template <int d>
  inline
  double EvaluateCardinalBSpline_x(const double x)
  {
    return (x < 0 || x >= d ? 0
	    : EvaluateCardinalBSpline<d-1>(x) - EvaluateCardinalBSpline<d-1>(x-1));
  }
  
  /*!
    evaluate the first derivative N_2'(x) of a cardinal B-spline
    (fastest possible variant)
  */
  template <>
  inline
  double EvaluateCardinalBSpline_x<2>(const double x)
  {
    if (x < 0 || x >= 2)
      return 0;
    return (x < 1 ? 1: -1);
  }

  /*!
    evaluate the first derivative N_1'(x) of a cardinal B-spline
  */
  template <>
  inline
  double EvaluateCardinalBSpline_x<1>(const double x)
  {
    return 0;
  }

  /*!
    evaluate the first derivative N_d'(x-k) of a shifted cardinal B-spline
  */
  template <int d>
  inline
  double EvaluateCardinalBSpline_x(const int k, const double x)
  {
    if (x < k || x >= k+d)
      return 0;
    return EvaluateCardinalBSpline<d-1>(k, x) - EvaluateCardinalBSpline<d-1>(k+1, x);
  }

  /*!
    evaluate the second derivative N_d''(x-k) of a shifted cardinal B-spline
  */
  template <int d>
  inline
  double EvaluateCardinalBSpline_x_sq(const int k, const double x)
  {
    return EvaluateCardinalBSpline_x<d-1>(k, x) - EvaluateCardinalBSpline_x<d-1>(k+1, x);
  }

  /*!
    evaluate the first derivative N_1'(x-k) of a shifted cardinal B-spline
  */
  template <>
  inline
  double EvaluateCardinalBSpline_x<1>(const int k, const double x)
  {
    return 0;
  }
  
  /*!
    evaluate the first derivative of a primal CDF function
      phi_{j,k}'(x) = 2^{3*j/2}N_d'(2^jx-k+d/2)
  */
  template <int d>
  inline
  double EvaluateCardinalBSpline_td_x(const int j, const int k, const double x)
  {
#if 0
    const double factor(1<<j);
    return factor * sqrt(factor) * EvaluateCardinalBSpline_x<d>(k, factor * x + d/2);
#else
    return twotothejhalf(3*j) * EvaluateCardinalBSpline_x<d>((1<<j) * x - k + d/2);
#endif
  }

  /*!
    evaluate the first derivative of a primal CDF function
      phi_{j,k}''(x) = 2^{5*j/2}N_d''(2^jx-k+d/2)
  */
  template <int d>
  inline
  double EvaluateCardinalBSpline_td_x_sq(const int j, const int k, const double x)
  {
#if 0
    const double factor(1<<j);
    return factor * factor * sqrt(factor) * EvaluateCardinalBSpline_x_sq<d>(k, factor * x + d/2);
#else
    return twotothejhalf(5*j) * EvaluateCardinalBSpline_x_sq<d>(k, (1<<j) * x + d/2);
#endif
  }


  /*!
    evaluate a shifted cardinal B-spline N_d(x-k) via recursion
    (remark: only use this version if the spline order d is unknown at compile time)
  */

  double EvaluateCardinalBSpline(const int d, const int k, double x)
  {
    if (x < k)
      return 0.;
    else
      {
	if (x >= k+d)
	  return 0.;
	else
	  {
	    // we know that x\in\supp N_d(.-k) 
	    if (d == 1)
	      return 1.;
	    else
	      {
		// hard-encode case d=2 
		if (d == 2)
		  {
		    if (x < k+1)
		      return x-k;
		    else
		      return 2.-(x-k);
		  }
		else
		  return ((x-k) * EvaluateCardinalBSpline(d-1, k, x)
			  + (k+d-x) * EvaluateCardinalBSpline(d-1, k+1, x)) / (d-1);
	      }
	  }
      }
  }

  
  /*!
    evaluate a primal CDF function
      phi_{j,k}(x) = 2^{j/2}N_d(2^jx-k+d/2)
  */
  inline double EvaluateCardinalBSpline_td(const int d, const int j, const int k, const double x)
  {
    const double factor(1<<j);
    return sqrt(factor) * EvaluateCardinalBSpline(d, k, factor * x + d/2);
  }
  
  /*!#include <algebra/polynomial.h>
    evaluate the first derivative N_d'(x-k) of a shifted cardinal B-spline
  */
  inline double EvaluateCardinalBSpline_x(const int d, const int k, const double x)
  {
    if (d == 1)
      return 0.;
    else
      return EvaluateCardinalBSpline(d-1, k, x) - EvaluateCardinalBSpline(d-1, k+1, x);
  }
  
  /*!
    evaluate the first derivative of a primal CDF function
      phi_{j,k}'(x) = 2^{j/2}N_d(2^jx-k+d/2)
  */
  inline double EvaluateCardinalBSpline_td_x(const int d, const int j, const int k, const double x)
  {
    const double factor(1<<j);
    return factor * sqrt(factor) * EvaluateCardinalBSpline_x(d, k, factor * x + d/2);
  }


/*------------------------------------------------------------------*
 *                        n       n   n - k + 1     n       n       *
 *  Binomialkoeffizient (   ) = (   ) --------- , (   ) = (   ) = 1 *
 *                        k      k-1      k         0       n       *
 *                                                                  *
 *  Rekursive Berechnung, damit keine zu grossen Zahlen entstehen   *
 *------------------------------------------------------------------*/

double binomd (long int n, long int k)
   {
   long int i;
   long int binomialkoeffizient = 1;

   if (n <  0 || k < 0 || k > n) return (-1);

   if (k == 0 || k == n) return (binomialkoeffizient);

   for (i = 1; i <= k; i++)
      binomialkoeffizient = binomialkoeffizient * (n - i + 1) / i;
   return (double)(binomialkoeffizient);
   } 


// helper function for cBspline
SparseMatrix<double> localBSplineCoeffs(const int d)
{
  int k, mu, nu;

  SparseMatrix<double> A(d,d);
  if(d==1)
    {
      A.set_entry(0, 0, 1);
      return A;
    }
  else 
    {
      SparseMatrix<double> B = localBSplineCoeffs(d-1);
      // first row (nu=0)
//       A.set_entry(0, 0, 0);
      for(k=1; k<d-1; k++)
	{
// 	  A.set_entry(0, k, 0);
	  for(mu=0; mu<=d-2; mu++)
	    A.set_entry(0, k, A.get_entry(0, k) + (pow(k, mu+1)*(B.get_entry(mu, k-1)-B.get_entry(mu, k)) + pow(-1, mu)*B.get_entry(mu, k-1))/(mu+1));
	}
//       A.set_entry(0, d-1, 0);
      for(mu=0; mu<=d-2; mu++)
	A.set_entry(0, d-1, (A.get_entry(0, d-1) + (pow(k, mu+1)*B.get_entry(mu, d-2) + pow(-1, mu)*B.get_entry(mu, d-2))/(mu+1)));	  
      // second to last row (nu>0)
      for(nu=1; nu<d; nu++)
	{
	  A.set_entry(nu, 0, B.get_entry(nu-1, 0)/nu);
	  for(k=1; k<d-1; k++)
	    {
	      A.set_entry(nu, k, B.get_entry(nu-1, k)/nu);
	      for(mu=nu; mu<=d-1; mu++)
		A.set_entry(nu, k, A.get_entry(nu, k) - (pow(-1, mu-nu)*binomd(mu, nu)*B.get_entry(mu-1, k-1))/mu);
	    }
// 	  A.put(nu, d-1) = 0;
	  for(mu=nu; mu<=d-1; mu++)
	    A.set_entry(nu, d-1, A.get_entry(nu, d-1) - (pow(-1, mu-nu)*binomd(mu, nu)*B.get_entry(mu-1, d-2))/mu);
	}

      return A;
    }
}



  /*!
     Expand a Bspline as a Picewise only primal CDF function
      phi_{j,k}(x) = 2^{j/2}N_d(2^jx-k+d/2)
  */
template <int d>
  Piecewise<double> ExpandBspline(const int j, const int k)
{
  Piecewise<double> r(j);
  Polynomial<double> p, q;
  p.set_coefficient(0, floor((double) d/2.0) -k);  //floor((double) d/2.0)
  p.set_coefficient(1, ldexp(1.0, j)); // p(x)=2^jx-k+floor(d/2)
  double factor = sqrt(ldexp(1.0, j));
  int m, n;
  SparseMatrix<double> coeffs(d, d);
  switch(d) // precalculated entries for 1<=d<=4
    {
    case 1:
      {
	coeffs.set_entry(0, 0, 1);
      }
      break;
    case 2:
      {
	coeffs.set_entry(1, 0, 1);

	coeffs.set_entry(0, 1, 2);
	coeffs.set_entry(1, 1, -1);
      }
      break;
    case 3:
      {
	coeffs.set_entry(2, 0, 0.5);

	coeffs.set_entry(0, 1, -1.5);
	coeffs.set_entry(1, 1, 3);
	coeffs.set_entry(2, 1, -1);

	coeffs.set_entry(0, 2, 4.5);
	coeffs.set_entry(1, 2, -3);
	coeffs.set_entry(2, 2, 0.5);
      }
      break;
    case 4:
      {
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
      break;
    default:
    {
      //return r;
      coeffs = localBSplineCoeffs(d);
    }
}
  for (n=0; n<=d-1; n++)
    {
     // coeffs.findc(n, ind);
     // CoeffsType col;
      for (m=0; m<=d-1; m++)
	q.set_coefficient(m, coeffs.get_entry(m, n));

      r.set_local_expansion(n + k - floor((double) d/2.0), q.substitute_into(p));  //(int)floor((double) d/2.0)
    }
  r *= factor;

  return r;
}




  /*!
    cardinal B-spline N_d(x) as Function object
  */
template <int d>
  class CardinalBSpline : public Function<1>
  {
  public:
    /*!
      default constructor: B-splines are real-valued
    */
    CardinalBSpline() : Function<1>(1) {}

    /*!
      virtual destructor
    */
    virtual ~CardinalBSpline() {}

    /*!
      value of a B-spline
    */
    inline double value(const Point<1>& p,
			const unsigned int component = 0) const
    {
      return EvaluateCardinalBSpline<d>(0, p(0));
    }
  
    /*!
      value of a B-spline
    */
    void vector_value(const Point<1> &p,
		      Vector<double>& values) const
    {
      values.resize(1, false);
      values[0] = value(p);
    }
  };


template <int d, int j, int k>
  class CardinalBSpline_td : public Function<1>
  {
  public:
    /*!
      default constructor: B-splines are real-valued
    */
    CardinalBSpline_td() : Function<1>(1) {}

    /*!
      virtual destructor
    */
    virtual ~CardinalBSpline_td() {}

    /*!
      value of a B-spline
    */
    inline double value(const Point<1>& p,
			const unsigned int component = 0) const
    {
      return EvaluateCardinalBSpline_td<d>(j, k, p(0));
    }
  
    /*!
      value of a B-spline
    */
    void vector_value(const Point<1> &p,
		      Vector<double>& values) const
    {
      values.resize(1, false);
      values[0] = value(p);
    }
  };
}

#endif
