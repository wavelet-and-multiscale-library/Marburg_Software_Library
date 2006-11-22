// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_FULL_LAPLACIAN_H
#define _WAVELETTL_FULL_LAPLACIAN_H

#include <algebra/qs_matrix.h>

namespace WaveletTL
{
  /*!
    The following class models the discretization of the Laplacian operator
    in a primal spline wavelet basis on a level j.
    We use the decomposition

      <A Psi_j,Psi_j>^T = T_j * <A Phi_j,Phi_j>^T * T_j^T

    where T_j models the multiscale transformation.
    The system is diagonally preconditioned.
  */
  template <int d, int dT>
  class FullLaplacian
  {
  public:
    /*!
      constructor taking an information object on some spline wavelet basis
    */
    FullLaplacian(const SplineBasisData<d,dT>& sd);

  protected:
    const SplineBasisData<d,dT>& sd_;
  };
}

#include <galerkin/full_laplacian.cpp>

#endif
