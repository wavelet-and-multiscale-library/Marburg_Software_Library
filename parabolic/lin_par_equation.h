// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_LIN_PAR_EQUATION_H
#define _WAVELETTL_LIN_PAR_EQUATION_H

namespace WaveletTL
{
  /*!
    This class models a linear parabolic equation of the form

      u'(t) = Au(t) + f(t) =: F(t,u(t)),  0 < t <= T
       u(0) = u_0

    where -A: \ell_{2,D} -> \ell_{2,D^{-1}} is a positive isomorphism
    and f: (0,T] -> \ell_{2,D^{-1}}, D being a diagonal preconditioner (see below).
    An equation of this form arises when equivalently reformulating
    a problem of the analogous form

      v'(t) = Bv(t) + g(t) =: G(t,v(t)), 0 < t <= T
       v(0) = v_0

    with isomorphism -B: H -> H' and right-hand side g: (0,T] -> H'
    by means of a biorthogonal wavelet basis Psi=\{\psi_\lambda\},
    setting F(t,v):=<G(t,v^\top\Psi),\tilde\Psi> for all v in \ell_{2,D}.
    
    The (unpreconditioned) stiffness matrix

      -A = (a(\psi_\lambda',\psi_\lambda))_{\lambda,\lambda'}

    is modeled in the template parameter class ELLIPTIC_EQ.
    
  */
  template <class ELLIPTIC_EQ>
  class LinearParabolicEquation
  {
  public:
  };
}

#include <parabolic/lin_par_equation.cpp>

#endif
