// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_SPLINE_BASIS_FLAVOR_H
#define _WAVELETTL_SPLINE_BASIS_FLAVOR_H

namespace WaveletTL
{
  /*!
    different flavors for SplineBasis
  */
  enum SplineBasisFlavor {
    P_construction,        // [Pr06]
    DS_construction_bio1,  // [DKU99,DS98], IGPMlib method #1, C_L = I					       
    DS_construction_bio2,  // ", IGPMlib method #2, Gamma_L = U*S*V, C_L = S^{-1/2}U^T, C_L_T = S^{-1/2}V
    DS_construction_bio3,  // ", IGPMlib method #3, transformation to Bernstein basis on [0,b]	       
    DS_construction_bio4,  // ", IGPMlib method #4, partial SVD				       
    DS_construction_bio5,  // ", IGPMlib method #5, transformation to Bernstein basis plus partial SVD
    DS_construction_bio5e  // ", IGPMlib method #5 + energy norm orth. of boundary wavelets [Ba01]
  };
}

#endif
