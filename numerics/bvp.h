// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
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
      right-hand side f
    */
    virtual const double f(const Point<DIM>& x) const
    {
      return f_->value(x);
    }
    
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
  };

}

#include <numerics/bvp.cpp>

#endif
