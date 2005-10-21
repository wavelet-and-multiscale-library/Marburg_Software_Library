// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_RECURSION_H
#define _MATHTL_RECURSION_H

#include <cmath>
#include <algebra/polynomial.h>
#include <algebra/vector.h>

namespace MathTL
{
  // some tools for three-term recursions
  // reference: Deuflhard/Bornemann, Numerische Mathematik 1

  /*!
    base class for inhomogeneous three-term recursions of the form
      p_k = a_k * p_{k-1} + b_k * p_{k-2} + c_k,   k=2,3,...
   */
  class Recursion
  {
  public:
    virtual ~Recursion() {}

    virtual double a(const unsigned int k) const = 0;
    virtual double b(const unsigned int k) const = 0;
    virtual double c(const unsigned int k) const = 0;

    //! evaluation of p_n
    double operator () (const unsigned int n) const;

    //! (trivial) forward summation
    double forwardSummation(const Vector<double>& coeffs) const;

    //! adjoint summation
    /*!
      adjoint summation of 
        \sum_{k=0}^n \alpha_k * p_k
      remarks:
      - stable for dominant solutions of the three-term recursion
      - potentially unstable for minimal solutions
    */
    double adjointSummation(const Vector<double>& coeffs) const;
    
  protected:
    //! first two entries of the recursion
    double p0_, p1_;
  };

  //! base class for homogeneous three-term recursions (c_k=0)
  class HomogeneousRecursion : public Recursion
  {
  public:
    double c(const unsigned int k) const { return 0.0; }
  };

  //! monomial recursion
  class MonomialRecursion: public HomogeneousRecursion
  {
  public:
    MonomialRecursion(const double x) : x_(x) { p0_ = 1.0; p1_ = x_; }

    double a(const unsigned int k) const { return x_; }
    double b(const unsigned int k) const { return 0.0; }

  protected:
    double x_;
  };

  //! recursion for the chebyshev polynomials P_n
  class ChebyshevRecursion : public HomogeneousRecursion
  {
  public:
    ChebyshevRecursion(const double x) : x_(x) { p0_ = 1.0; p1_ = x_; }

    double a(const unsigned int k) const { return 2*x_; }
    double b(const unsigned int k) const { return -1.0; }

  protected:
    double x_;
  };

  //! recursion for P_{n}=\cos(n.)
  class CosineRecursion : public HomogeneousRecursion
  {
  public:
    CosineRecursion(const double x) : x_(x) { p0_ = 1.0; p1_ = cos(x_); }

    double a(const unsigned int k) const { return 2*cos(x_); }
    double b(const unsigned int k) const { return -1.0; }

  protected:
    double x_;
  };

  //! recursion for P_{n}=\sin(n.)
  class SineRecursion : public HomogeneousRecursion
  {
  public:
    SineRecursion(const double x) : x_(x) { p0_ = 0.0; p1_ = sin(x_); }

    double a(const unsigned int k) const { return 2*cos(x_); }
    double b(const unsigned int k) const { return -1.0; }

  protected:
    double x_;
  };
  
  //! recursion for Jacobi polynomials
  class JacobiRecursion : public HomogeneousRecursion
  {
  public:
    JacobiRecursion(const double a, const double b, const double x);
  
    double a(const unsigned int k) const;
    double b(const unsigned int k) const;
    
  protected:
    double a_, b_, x_;
  };

  //! recursion for Legendre polynomials
  class LegendreRecursion : public JacobiRecursion
  {
  public:
    LegendreRecursion(const double x) : JacobiRecursion(0.0,0.0,x) {}
  };
}

// include implementation of inline functions
#include <numerics/recursion.cpp>

#endif
