// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_BVP_H
#define _MATHTL_BVP_H

#include <geometry/point.h>
#include <geometry/atlas.h>
#include <algebra/matrix.h>
#include <utils/array1d.h>
#include <utils/fixed_array1d.h>
#include <utils/function.h>

namespace MathTL
{
  /*!
    Abstract base class for a general (vector valued)
    two-point boundary value problem on [0,1]
    
      y'(t) = f(t,y(t))

    with linear, homogeneous boundary conditions

      r(y(0),y(1)) = Ay(0) + By(1) = 0
  */
  template <unsigned int DIM>
  class TwoPointBVP
  {
  public:
    /*!
      virtual destructor
     */
    virtual ~TwoPointBVP();

    /*!
      apply the left boundary condition matrix to v
    */
    virtual void apply_A(const Point<DIM>& v, Point<DIM>& result) const = 0;
    
    /*!
      apply the right boundary condition matrix to v
    */
    virtual void apply_B(const Point<DIM>& v, Point<DIM>& result) const = 0;
    
    /*!
      apply the right--hand side f to (t,v)
    */
    virtual void apply_f(const double t, const Point<DIM>& v, Point<DIM>& result) const = 0;
  };

  /*!
    Base class for a symmetric, second-order elliptic
    boundary value problem in divergence form over some domain
    Omega in R^d with boundary Gamma=dOmega

      -div(a(x)grad u(x)) + q(x)u(x) = f(x) in Omega

    with some boundary conditions (Dirichlet, Neumann, Robin).
    However, we only specify the parameters a, q and f in this class
    and postpone the boundary condition treatment to the
    discretization process. It will be implicitly assumed that
    the wavelet bases or frames used in a wavelet-Galerkin scheme
    fulfill the appropriate boundary conditions.
  */
  template <unsigned int DIM>
  class EllipticBVP
  {
  public:
    /*!
      constructor with given (scalar) coefficients
    */
    EllipticBVP(const Function<DIM>* a,
		const Function<DIM>* q,
		const Function<DIM>* f);

    /*!
      virtual destructor
    */
    virtual ~EllipticBVP() {}

    /*!
      diffusion coefficient a
     */
    virtual const double a(const Point<DIM>& x) const
    {
      return a_->value(x);
    }

    /*!
      reaction coefficient q
    */
    virtual const double q(const Point<DIM>& x) const
    {
      return q_->value(x);
    }

    /*!
      flag to indicate whether all coefficients are constants
      (speeds up quadrature a bit)
    */
    virtual const bool constant_coefficients() const
    {
        return false;
    }
    
    /*!
      right-hand side f
    */
    virtual const double f(const Point<DIM>& x) const
    {
      return f_->value(x);
    }

    /*!
      set the right-hand side to another function
    */
    void set_f(const Function<DIM>* f);
    
  protected:
    //! diffusion coefficient
    const Function<DIM>* a_;

    //! reaction coefficient
    const Function<DIM>* q_;
    
    //! right-hand side
    const Function<DIM>* f_;
  };

  /*!
    The Poisson equation on a d-dimensional domain
      -Delta u(x) = f(x)
  */
  template <unsigned int DIM>
  class PoissonBVP
    : public EllipticBVP<DIM>
  {
  public:
    /*!
      constructor with given right-hand side
    */
    PoissonBVP(const Function<DIM>* f);

    /*!
      diffusion coefficient a
     */
    const double a(const Point<DIM>& x) const { return 1.0; }

    /*!
      reaction coefficient q
    */
    const double q(const Point<DIM>& x) const { return 0.0; }

    /*!
      flag for constant coefficients
    */
    const bool constant_coefficients() const { return true; }
  };
  
  
  /*!
    The Poisson equation with non-constant diffusion coefficient a(x)
    on a d-dimensional domain
      -div(a(x)grad u(x)) = f(x) in Omega
  */
  template <unsigned int DIM>
  class PoissonBVP_Coeff
    : public EllipticBVP<DIM>
  {
  public:
    /*!
      constructor with given diffusion coefficient and right-hand side
    */
    PoissonBVP_Coeff(const Function<DIM>* a, const Function<DIM>* f);

    /*!
      constructor with given right-hand side, where also an alternative representation f_ar: \Omega -> \R^2 of the rhs is given in the sense of:
      rhs_{\lambda} = \int_{\Omega} a(x) nabla u(x) * nabla \Psi_{\lambda}(x) dx = \int_{\Omega} f(x) \Psi_{\lambda}(x) dx for all \lambda \in \Lambda
      alternative representation of rhs:
      rhs_{\lambda} = \int_{\Omega} a(x) nabla u(x) * nabla \Psi_{\lambda}(x) dx = \int_{\Omega} f_ar(x) * \nabla \Psi_{\lambda}(x) dx for all \lambda \in \Lambda
    */
    PoissonBVP_Coeff(const Function<DIM>* a, const Function<DIM>* f, const Function<DIM>* f_ar);

    /*!
      reaction coefficient q
    */
    const double q(const Point<DIM>& x) const { return 0.0; }

    /*!
      flag for constant coefficients
    */
    const bool constant_coefficients() const { return false; }

    /*!
      set the coefficient a to another function
    */
    void set_a(const Function<DIM>* a);

    /*!
      evaluate right-hand side f_ar
    */
    const void f_ar(const Point<DIM>& x, Vector<double>& values) const { f_ar_->vector_value(x, values); }

    protected:

    //! alternative representation of the right-hand side
    const Function<DIM>* f_ar_;

  };
  
  

  /*!
    The Identity equation on a d-dimensional domain
      u(x) = f(x).
    Actually the identity operator is noe second order operator.
    We neglect this point for the moment.
  */
  template <unsigned int DIM>
  class IdentityBVP
    : public EllipticBVP<DIM>
  {
  public:
    /*!
      constructor with given right-hand side
    */
    IdentityBVP(const Function<DIM>* f);

    /*!
      diffusion coefficient a
     */
    const double a(const Point<DIM>& x) const { return 0.0; }

    /*!
      reaction coefficient q
    */
    const double q(const Point<DIM>& x) const { return 1.0; }

    /*!
      flag for constant coefficients
    */
    const bool constant_coefficients() const { return true; }
  };



  /*!
    Base class for a symmetric, second-order elliptic
    boundary value problem in divergence form over some domain
    Omega in R^d with boundary Gamma=dOmega

      -delta u(x)= f(x) in Omega

    with some boundary conditions (Dirichlet, Neumann, Robin).
    However, we only specify the parameter f in this class
    and postpone the boundary condition treatment to the
    discretization process. It will be implicitly assumed that
    the wavelet bases or frames used in a wavelet-Galerkin scheme
    fulfill the appropriate boundary conditions.
  */
   template <unsigned int DIM>
   class BiharmonicBVP
   {
  public:
    /*!
      constructor with given (scalar) coefficients
    */
    BiharmonicBVP(const Function<DIM>* f);

    /*!
      virtual destructor
    */
    virtual ~BiharmonicBVP() {}

    /*!
      diffusion coefficient a
     */
    

    /*!
      flag to indicate whether all coefficients are constants
      (speeds up quadrature a bit)
    */
    const bool constant_coefficients() const;
    
    /*!
      right-hand side f
    */
    const double f(const Point<DIM>& x) const
    {
      return f_->value(x);
    }

    /*!
      set the right-hand side to another function
    */
    void set_f(const Function<DIM>* f);
    
  protected:
    
    //! right-hand side
    const Function<DIM>* f_;
  };


}

#include <numerics/bvp.cpp>

#endif
