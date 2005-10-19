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
    Omega in R^d with boundary Gamma=dOmega,
    with homogeneous Dirichlet/Neumann/Robin boundary conditions

      -div(a(x)grad u(x)) + q(x)u(x) = f(x) in Omega
                                u(x) = 0 on Gamma_D
                            du/dn(x) = 0 on Gamma\Gamma_D

    The entire problem of course depends on an atlas of Omega,
    which is given as a template parameter ATLAS.
    For each patch kappa_i(\Box) of the atlas, you have to specify
    2*d Dirichlet boundary condition orders
    (0 <-> no b.c., 1 <-> Dirichlet b.c.).
    However, the instance of ATLAS will be accessed only when it
    comes to a discretization, e.g., in a wavelet-Galerkin scheme.
    For examples concerning atlas, cf. geometry/atlas.h
  */
  template <unsigned int DIM>
  class EllipticBVP
  {
  public:
    /*!
      constructor with a given atlas, boundary conditions
      and scalar coefficients
    */
    EllipticBVP(const Atlas<DIM>* atlas,
		const Array1D<FixedArray1D<int,2*DIM> >& bc,
		const Function<DIM>* a,
		const Function<DIM>* q,
		const Function<DIM>* f);

    //! virtual destructor
    virtual ~EllipticBVP();

    /*!
      diffusion coefficient a
     */
    const double a(const Point<DIM>& x) const
    {
      return a_->value(x);
    }

    /*!
      reaction coefficient q
    */
    const double q(const Point<DIM>& x) const
    {
      return q_->value(x);
    }

    /*!
      right-hand side f
    */
    const double f(const Point<DIM>& x) const
    {
      return f_->value(x);
    }

  protected:
    //! pointer to the underlying atlas
    const Atlas<DIM>* atlas_;

    //! flag for deletion of the atlas
    bool delete_atlas;

    //! boundary conditions
    Array1D<FixedArray1D<int,2*DIM> > bc_;

    //! flag for deletion of the functions
    bool delete_functions;
    
    //! diffusion coefficient
    const Function<DIM>* a_;

    //! reaction coefficient
    const Function<DIM>* q_;
    
    //! right-hand side
    const Function<DIM>* f_;
  };
}

#include <numerics/bvp.cpp>

#endif
