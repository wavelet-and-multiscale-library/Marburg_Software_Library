// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_DKU_H
#define _WAVELETTL_DKU_H

#include <iostream>
#include <cmath>
#include <Rd/cdf_utils.h>
#include <Rd/cdf_basis.h>
#include <algebra/matrix.h>
#include <utils/array1d.h>

using MathTL::Matrix;

namespace WaveletTL
{
  /*!
    biorthogonalization methods for the DKU basis, see, e.g., [B]
  */
  enum DKUBiorthogonalizationMethod
    {
      none,         // C_L = I
      SVD,          // Gamma_L = U*S*V, C_L = S^{-1/2}U^T, C_L_T = S^{-1/2}V
      Bernstein,    // transformation to Bernstein basis on [0,b]
      partialSVD,   // partial SVD
      BernsteinSVD  // transformation to Bernstein basis plus partial SVD
    };

  /*!
    Template class for the wavelet bases on the interval introduced in [DKU], [DS].
    All formulas refer to [DKU], except those explicitly denoted by [DS].

    References:
    [B]   Barsch:
          Adaptive Multiskalenverfahren fuer elliptische partielle Dgln. - Realisierung,
	  Umsetzung und numerische Ergebnisse
    [DKU] Dahmen, Kunoth, Urban:
          Biorthogonal spline-wavelets on the interval - Stability and moment conditions
    [DS]  Dahmen, Schneider:
          Wavelets with complementary boundary conditions - Function spaces on the cube
  */
  template <int d, int dT>
  class DKUBasis
  {
  public:
    /*!
      constructor
      (ellT = ell2T)
    */
    explicit DKUBasis(DKUBiorthogonalizationMethod bio = none);

    //! DKU parameter ell-tilde (3.2.10)
    inline const int ellT() const { return ellT_; }

    //! DKU abbreviation (3.2.16)
    inline const int ell() const { return ell_; }

    //! coarsest possible level
    inline const int j0() const { return (int) ceil(log(ellT_+ell2T_-1.)/log(2.0)+1); }

  protected:
    int ell1_, ell2_, ell1T_, ell2T_, ell_, ellT_;

    /*!
      an instance of the CDF basis on R
    */
    CDFBasis<d, dT> cdf_;

    Matrix<double> Alpha_, AlphaT_;
    Matrix<double> BetaL_, BetaLT_;
    Matrix<double> GammaL_;
    Matrix<double> CL_, CLT_;

  private:
    // private helper routines
    double alpha_(const int m, const int r) const;
    double alphaT_(const int m, const int r) const;
    double betaL_(const int m, const int r) const;
    double betaLT_(const int m, const int r) const;
  };
}

#include <interval/dku.cpp>

#endif
