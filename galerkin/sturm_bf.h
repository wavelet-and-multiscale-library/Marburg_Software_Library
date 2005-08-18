// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_STURM_BILINEAR_FORM_H
#define _WAVELETTL_STURM_BILINEAR_FORM_H

namespace WaveletTL
{
  /*!
    This class models the bilinear form

      a(u,v) = \int_0^1 [p(t)u'(t)v'(t)+q(t)u(t)v(t)] dt

    corresponding to the Sturm boundary value problem on [0,1]
    
      -(py')'(t) + q(t)y(t) = g(t), 0 <= t <= 1 
      
    with first (Dirichlet) or second (Neumann) order b.c.'s as modeled in
    the class simpleSturmBVP.
   
    The evaluation of a(.,.) is possible for arguments \psi_\lambda
    which stem from a wavelet basis \Psi=\{\psi_\lambda\} of the corresponding
    function space over [0,1].
    To achieve independence from the concrete choice of \Psi, the wavelet basis
    class is given as a template parameter WBASIS and should provide a constructor of
    the form
       WBASIS::WBASIS(const bool bc_left, const bool bc_right)
    where the parameters bc_* indicate where to enforce homogeneous Dirichlet
    boundary conditions.
    Of course a natural concrete value for WBASIS is the template class DKUBasis<d,dT>.
   */
  template <class WBASIS>
  class SturmBilinearForm
  {
  public:
    SturmBilinearForm(const simpleSturmBVP& bvp);
    
    /*!
      evaluate the bilinear form
    */
    double operator () (const typename WBASIS::Index& lambda,
			const typename WBASIS::Index& nu) const;

  protected:
    const simpleSturmBVP& bvp_;
    WBASIS wbasis_;
  };
}

#include <galerkin/sturm_bf.cpp>

#endif
