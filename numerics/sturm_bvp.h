// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_STURM_BVP_H
#define _MATHTL_STURM_BVP_H

#include <numerics/bvp.h>

namespace MathTL
{
  /*!
    Abstract base class for a (real valued) Sturm boundary value problem on [0,1]
    
      -(py')'(t) + q(t)y(t) = g(t), 0 <= t <= 1 
      
    with third order homogeneous boundary conditions
   
      R_0y := alpha_0*y(0) + alpha_1*p(0)*y'(0) = 0
      R_1y := beta_0 *y(1) + beta_1 *p(1)*y'(1) = 0

    You have to take care that for a fundamental system {y_1,y_2} of the diff. eq.,
    one has a nontrivial determinant R_0y_1*R_1y_2-R_1y_1*R_0y_2 != 0.
  */
  class SturmBVP
    : public BVP<2>
  {
  public:
    /*!
      virtual destructor
     */
    virtual ~SturmBVP ();
    
    /*!
      diffusion coefficient
    */
    virtual double p(const double t) const = 0;

    /*!
      first derivative of the diffusion coefficient
    */
    virtual double p_prime(const double t) const = 0;

    /*!
      reaction coefficient
    */
    virtual double q(const double t) const = 0;

    /*!
      right-hand side
    */
    virtual double g(const double t) const = 0;
  
    /*!
      left boundary condition coefficient 0
     */
    virtual double alpha0() const = 0;

    /*!
      left boundary condition coefficient 1
     */
    virtual double alpha1() const = 0;

    /*!
      right boundary condition coefficient 0
     */
    virtual double beta0() const = 0;

    /*!
      left boundary condition coefficient 1
    */
    virtual double beta1() const = 0;
    
    //! right-hand side as a general BVP (inherited)
    void apply_f(const double t, const Point<2>& v, Point<2>& result) const;

    //! left boundary condition matrix (inherited)
    void apply_A(const Point<2>& v, Point<2>& result) const;
    
    //! right boundary condition matrix (inherited)
    void apply_B(const Point<2>& v, Point<2>& result) const;
  };
}

#include <numerics/sturm_bvp.cpp>

#endif
