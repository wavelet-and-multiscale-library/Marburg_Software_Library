// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of FrameTL - the Wavelet Template Library        |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Manuel Werner                                                      |
// +--------------------------------------------------------------------+

#ifndef _FRAMETL_PARAMETRIZATION_H
#define _FRAMETL_PARAMETRIZATION_H

#include <iostream>
#include <geometry/point.h>
#include <utils/array1d.h>
#include <utils/function.h>
#include <algebra/vector.h>
#include <algebra/matrix.h>

using std::cout;
using std::endl;
using MathTL::Point;
using MathTL::Array1D;
using MathTL::Function;
using MathTL::Vector;
using MathTL::Matrix;

namespace FrameTL
{
  /*!
    Abstract base class for parametrization of single patches forming
    a bounded domain in \matbb R.
   */
  template <unsigned int DIM_m, unsigned int DIM_d>
  class Parametrization
  {
  public:
    /*!
      virtual destructor
     */
    virtual ~Parametrization () {};

    /*!
      signum of det( D (mapPoint(x,y) ) )
      identical for all (x,y)!
     */
    virtual const unsigned short int get_sgn_det_D() const = 0;

    /*!
      maps a point
     */
    virtual void mapPoint(Point<DIM_d>&, const Point<DIM_m>&) const = 0;

    /*!
      inverse mapping
     */
    virtual void mapPointInv(Point<DIM_m>&, const Point<DIM_d>&) const = 0;

    /*!
      det( D (mapPoint(x,y)) )
     */
    virtual const double det_D(const Point<DIM_m>&) const = 0;
    
    /*!
      |det( D (mapPoint(x,y)) )|
     */
    virtual const double abs_Det_D(const Point<DIM_m>&) const = 0;

    /*!
      \partial/\partial x (det(D kappa))(s,t)
     */
    virtual const double d_x_det_D(const Point<DIM_m>&) const = 0;

    /*!
     \partial/\partial y (det(D kappa))(s,t)
     */
    virtual const double d_y_det_D(const Point<DIM_m>&) const = 0;
   
    /*!
      \partial / \partial_dim \kappa^(direc)
    */    
    virtual const double d_dim_kappa_direc(const unsigned short int& dim,
					   const unsigned short int& direc,
					   const Point<DIM_m>&) const = 0;

  };

  /*!
    This class models parametrizations for arbitrary quadrangles in \mathbb R^2.
    It's crucial functionality is to map a single point, lying in (0,1)^2 to a
    point in the quadrangle and vice versa. The involved mapping explicitely looks like:

        k(s,t) := (1-s)*(1-t)b_00 + s*(1-t)b_10 + (1-s)*t*b_01 + s*t*b_11,

    where the b_ij are the vertices of the qudrangle at hand.
    
   */
  class LinearBezierMapping : public Parametrization<2,2>
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
    
    /*!
      signum of det( D (mapPoint(x,y) ) )
      identical for all (x,y)!
     */
    const unsigned short int get_sgn_det_D() const;

    /*!
      maps a point
     */
    void mapPoint(Point<2>&, const Point<2>&) const;

    /*!
      inverse mapping
     */
    void mapPointInv(Point<2>&, const Point<2>&) const;

    /*!
      det( D (mapPoint(x,y)) )
     */
    const double det_D(const Point<2>&) const;
    
    /*!
      |det( D (mapPoint(x,y)) )|
     */
    const double abs_Det_D(const Point<2>&) const;

    /*!
      \partial/\partial x (det(D kappa))(s,t)
     */
    const double d_x_det_D(const Point<2>&) const;

    /*!
     \partial/\partial y (det(D kappa))(s,t)
     */
    const double d_y_det_D(const Point<2>&) const;
   
    /*!
      \partial / \partial_dim \kappa^(direc)
    */    
    const double d_dim_kappa_direc(const unsigned short int& dim,
				   const unsigned short int& direc,
				   const Point<2>&) const;

  protected:
    //vertices of the quadrangle to parametrize
    Point<2> b_00;    
    Point<2> b_01;
    Point<2> b_10;    
    Point<2> b_11;

  private:
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

  /*!
    stream output for LinearBezierMapping
   */
  std::ostream& operator << (std::ostream&, const LinearBezierMapping&);






  /*!
    This class models parametrizations for patches forming
    cuboids in \mathbb R^d.
    This mapping is of the simple form

    k( point ) = A * point + b, where A \in \mathbb K^{d x d} and b \in \mathbb K^d.


   */
  template <unsigned int DIM>
  class AffinLinearMapping : public Parametrization<DIM,DIM>
  {
  public:

    /*!
      size type (cf. STL containers)
    */
    typedef size_t size_type;

    /*!
      pureyl virtual destructor
     */
    ~AffinLinearMapping () {};

    /*!
      default constructor:
    */
    AffinLinearMapping();
   
    /*!
      copy constructor
     */
    AffinLinearMapping (const AffinLinearMapping&);

    /*!
      constructor for initialization of the four vertices
     */
    AffinLinearMapping (const Matrix<double> &, const Vector<double> &);

    /*!
      assignment operator
    */
    AffinLinearMapping& operator = (const AffinLinearMapping&);

    const Matrix<double>& get_A() const;
    const Vector<double>& get_b() const;

   /*!
      signum of det( D (mapPoint(x,y) ) )
      identical for all (x,y)!
     */
    const unsigned short int get_sgn_det_D() const;

    /*!
      maps a point
     */
    void mapPoint(Point<DIM>&, const Point<DIM>&) const;

    /*!
      inverse mapping
     */
    void mapPointInv(Point<DIM>&, const Point<DIM>&) const;
    /*!
      det( D (mapPoint(x,y)) )
     */
    const double det_D(const Point<DIM>&) const;
    
    /*!
      |det( D (mapPoint(x,y)) )|
     */
    const double abs_Det_D(const Point<DIM>&) const;

    /*!
      \partial/\partial x (det(D kappa))(s,t)
     */
    const double d_x_det_D(const Point<DIM>&) const;

    /*!
     \partial/\partial y (det(D kappa))(s,t)
     */
    const double d_y_det_D(const Point<DIM>&) const;
   
    /*!
      \partial / \partial_dim \kappa^(direc)
    */    
    const double d_dim_kappa_direc(const unsigned short int& dim,
				   const unsigned short int& direc,
				   const Point<DIM>&) const;
   
  protected:
    Matrix<double> A;
    Vector<double> b;

  private:
    //signum of det( D (mapPoint(x,y) ) )
    //identical for all (x,y)!
    bool sgn_det_D;

  };

  /*!
    stream output for AffinLinearMapping
  */
  template <unsigned int DIM>
  std::ostream& operator << (std::ostream&, const AffinLinearMapping<DIM>&);

}

// include implementation
#include "parametrization.cpp"

#endif
