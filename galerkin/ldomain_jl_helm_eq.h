// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_LDOMAIN_JL_HH_EQ_H
#define _WAVELETTL_LDOMAIN_JL_HH_EQ_H

#include <Ldomain/ldomain_jl_basis.h>
#include <galerkin/infinite_preconditioner.h>

namespace WaveletTL
{

  /*!
    This class models the (preconditioned) infinite-dimensional matrix problem
    
      Au = D^{-1}LD^{-1}u = D^{-1}F

    when reformulating a Helmholtz equation on the L--shaped domain
    
      -Delta u(x)) + alpha*u(x) = f(x)

    with first order b.c.'s as an equivalent operator equation
    within \ell_2 by means of an LDomainJLBasis wavelet basis.

    The corresponding bilinear form in

      L = (a(\psi_\nu,\psi_\lambda))_{\lambda,\nu}

    is

      a(u,v) = \int_Omega [a(x)grad u(x)grad v(x)+q(x)u(x)v(x)] dx
      
    and the right-hand side is
 
      F(v) = \int_0^1 f(x)v(x) dx
   
    The evaluation of a(.,.) and f is possible for arguments \psi_\lambda
    which stem from a wavelet basis \Psi=\{\psi_\lambda\} of the corresponding
    function space over Omega.
  */
  class LDomainJLHelmholtzEquation
    : public FullyDiagonalEnergyNormPreconditioner<LDomainJLBasis::Index>
  {
  public:
    /*!
      the wavelet basis class
    */
    typedef LDomainJLBasis WaveletBasis;

    /*!
      wavelet index class
    */
    typedef WaveletBasis::Index Index;

    /*!
      constructor from a right-hand side, a basis, alpha and specified b.c.'s
    */
    LDomainJLHelmholtzEquation(const WaveletBasis& basis,
			       const double alpha,
			       const InfiniteVector<double,Index>& y);
  };
}

#include <galerkin/ldomain_jl_helm_eq.cpp>

#endif
