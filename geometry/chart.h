// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_CHART_H
#define _MATHTL_CHART_H

#include <iostream>
#include <string>
#include <algebra/matrix.h>
#include <algebra/vector.h>
#include <geometry/point.h>

using std::string;

namespace MathTL
{
  /*!
    Abstract base class for smooth parametrizations
      kappa: (0,1)^d -> R^m
    of single "patches" in R^m.
   */
  template <unsigned int DIM_d, unsigned int DIM_m = DIM_d>
  class Chart
  {
  public:
    //! virtual destructor
    virtual ~Chart() {}

    /*!
      map a point (forward) y = kappa(x)
     */
    virtual void map_point(const Point<DIM_d>& x, Point<DIM_m>& y) const = 0;

    /*!
      inverse mapping y = kappa^{-1}(x)
     */
    virtual void map_point_inv(const Point<DIM_d>& x, Point<DIM_m>& y) const = 0;

    /*!
      square root of the Gram determinant sqrt(det(Dkappa(x)^T * Dkappa(x)))
      (additional factor for integration over "plain" functions)
    */
    virtual const double Gram_factor(const Point<DIM_d>& x) const = 0;

    /*!
      i-th partial derivative of Gram factor
      (additional factor for integration over first derivatives)
    */
    virtual const double Gram_D_factor(const unsigned int i,
				       const Point<DIM_d>& x) const = 0;

    /*!
      (i,j)-th element of (D (kappa^{-1}))(x)
    */
    virtual const double Dkappa_inv(const unsigned int i,
				    const unsigned int j,
				    const Point<DIM_d>& x) const = 0;

    /*!
      checks whether a special point x lies in the patch represented by this
      parametrization
    */
    virtual const bool in_patch(const Point<DIM_m>& x) const = 0;

    /*!
      returns a string representation of this object
     */
    virtual const string to_string() const = 0;
  };
  
  /*!
    stream output for arbitrary charts
  */
  template <unsigned int DIM_d, unsigned int DIM_m>
  std::ostream& operator << (std::ostream& s, const Chart<DIM_d,DIM_m>&);

  //
  // Some examples:

  /*!
    affine linear mapping y = A*x+b
   */
  template <unsigned int DIM>
  class AffineLinearMapping
    : public Chart<DIM,DIM>
  {
  public:
    //! default constructor, yields the identity mapping
    AffineLinearMapping();

    //! constructor from A and b (dimensions should fit)
    AffineLinearMapping(const Matrix<double>& A, const Point<DIM>& b);

    void map_point(const Point<DIM>&, Point<DIM>&) const;
    void map_point_inv(const Point<DIM>&, Point<DIM>&) const;
    const double Gram_factor(const Point<DIM>&) const;
    const double Gram_D_factor(const unsigned int i, const Point<DIM>& x) const;
    const double Dkappa_inv(const unsigned int i, const unsigned int j,
			    const Point<DIM>& x) const;
    const bool in_patch(const Point<DIM>& x) const;

    const string to_string() const;

    /*!
      static field to store the name of the class
     */
    static const string className;

    //! read access to A
    const Matrix<double>& A() const { return A_; }
    
    //! read access to b
    const Point<DIM>& b() const { return b_; }

  protected:
    Matrix<double> A_, A_inv;
    double det_A;
    Point<DIM> b_;
  };

  /*!
    This class models parametrizations for arbitrary quadrangles in \mathbb R^2.
    It's crucial functionality is to map a single point, lying in (0,1)^2 to a
    point in the quadrangle and vice versa. The involved mapping explicitely looks like:

        k(s,t) := (1-s)*(1-t)b_00 + s*(1-t)b_10 + (1-s)*t*b_01 + s*t*b_11,

    where the b_ij are the vertices of the qudrangle at hand.
    
   */
  class LinearBezierMapping : public Chart<2,2>
  {

  public:

    /*!
      pureyl virtual destructor
     */
    ~LinearBezierMapping () {};
    
    /*!
      default constructor:
    */
    LinearBezierMapping ();
    /*!
      copy constructor
     */
    LinearBezierMapping (const LinearBezierMapping&);

    /*!
      constructor for initialization of the four vertices
     */
    LinearBezierMapping (const Point<2> &, const Point<2> &,
			 const Point<2> &, const Point<2> &);

    /*!
      assignment operator
    */
    LinearBezierMapping& operator = (const LinearBezierMapping& x);

    /*!
      setup routine, called by preceding constructor,
      sets up generic qudrangle, needed for beeingable to invert the mapping.
      idea: by appropriate shifting rotation and shearing, the qudrangle can be
      transformed into another qudrangle for which it is clear how its parametrization
      be inverted.
     */
    void setup ();

    /*!
      access to vertices
    */
    const Point<2>& get_b_00() const;
    const Point<2>& get_b_01() const;
    const Point<2>& get_b_10() const;
    const Point<2>& get_b_11() const;
 
    void map_point(const Point<2>&, Point<2>&) const;

    void map_point_inv(const Point<2>&, Point<2>&) const;

    const double Gram_factor(const Point<2>& x) const;
    const double Gram_D_factor(const unsigned int i, const Point<2>& x) const;
    const double Dkappa_inv(const unsigned int i, const unsigned int j,
			    const Point<2>& x) const;
    const bool in_patch(const Point<2>& x) const;
    
    const string to_string() const;
   
    /*!
      static field to store the name of the class
    */
    static const string className;

    //vertices of the quadrangle to parametrize
    Point<2> b_00;    
    Point<2> b_01;
    Point<2> b_10;
    Point<2> b_11;

  private:
    //the following two routines are helpers that are only of interest
    //in this special 2D BezierMapping case!
    //only for the special case of LinearBezierMapping
    const double det_DKappa(const Point<2>& p) const;

    //partial derivative with respect to i-th component
    //of j-th component of Kappa
    const double partial_i_Kappa_j(const unsigned int i,
				   const unsigned int j,
				   const Point<2>& x) const;

    //vertices of generic qudrangle, needed for inverting this mapping,
    //initialized in constructor
    Point<2> b_gen_00;
    Point<2> b_gen_10;
    Point<2> b_gen_01;
    Point<2> b_gen_11;

    Vector<double> d_ds_d_dt_kappa_r;
    Vector<double> d_dt_d_ds_kappa_r;

    Vector<double> min_b00_plus_b10;
    Vector<double> min_b00_plus_b01;

    //some quantities, worth storing to prevent dispensable recomutation
    double cos_rot_angle;
    double sin_rot_angle;
    double rot_angle;
    double shearing_param;
    double scaleX;
    double scaleY;
    
    //signum of det( D (mapPoint(x,y) ) )
    //identical for all (x,y)!
    bool sgn_det_D;


  };
}

#include "geometry/chart.cpp"

#endif
