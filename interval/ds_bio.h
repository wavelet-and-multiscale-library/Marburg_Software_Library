// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_DS_BIO_H
#define _WAVELETTL_DS_BIO_H

namespace WaveletTL
{
  /*!
    biorthogonalization methods for the DS basis, see, e.g., [B]

    The methods 'partialSVD' and 'BernsteinSVD' enable the [DKU]/[DS] boundary treatment:
    at the boundary, exactly one generator and one wavelet does not vanish,
    which can then be modified to satisfy homogeneous boundary conditions.
  */
  enum DSBiorthogonalizationMethod
    {
      none,         // method #1 in IGPMlib, C_L = I
      SVD,          // method #2 in IGPMlib, Gamma_L = U*S*V, C_L = S^{-1/2}U^T, C_L_T = S^{-1/2}V
      Bernstein,    // method #3 in IGPMlib, transformation to Bernstein basis on [0,b]
      partialSVD,   // method #4 in IGPMlib, partial SVD
      BernsteinSVD  // method #5 in IGPMlib, transformation to Bernstein basis plus SVD
    };
}

#endif
