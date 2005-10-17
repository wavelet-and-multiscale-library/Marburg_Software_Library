// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_IVP_H
#define _MATHTL_IVP_H

#include <geometry/point.h>
#include <algebra/matrix.h>

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
    Abstract base class for a two--dimensional second-order elliptic
    boundary value problem in divergence form over some domain
    Omega in R^d with boundary Gamma=dOmega (reaction-diffusion),
    with homogeneous Dirichlet/Neumann/Robin boundary conditions

      div(a(x)grad u(x)) + q(x)u(x) = f(x) in Omega
      u(x)     = 0 on Gamma_D
      du/dn(x) = 0 on Gamma_N

    The entire problem of course depends on an atlas of Omega,
    which is given as a template parameter ATLAS. The atlas is
    responsible for the specification of the various boundary conditions
    on Gamma. However, the instance of ATLAS will be accessed only when it
    comes to a discretization, e.g., in a wavelet-Galerkin scheme.
    For examples concerning atlas, cf. geometry/atlas.h
  */
  template <class ATLAS>
  class EllipticBVP
  {
  public:
    /*!
      default constructor
    */
    EllipticBVP(const ATLAS& atlas);

    //! virtual destructor
    virtual ~EllipticBVP();
  protected:
    //! reference to the underlying atlas
    const ATLAS& atlas_;
  };
}

#include <numerics/bvp.cpp>

#endif
