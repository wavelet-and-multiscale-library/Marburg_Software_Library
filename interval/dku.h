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
#include <Rd/cdf_basis.h>
#include <algebra/matrix.h>

using MathTL::Matrix;

namespace WaveletTL
{
  /*!
    Template class for the wavelet bases on the interval introduced in [DKU].

    References:
    [DKU]
  */
  template <int d, int dt>
  class DKUBasis
  {
  public:
    //! default constructor
    DKUBasis();

    //! left support bound for a primal CDF function
    inline const int ell1() const { return ell1_; }

    //! right support bound for a primal CDF function
    inline const int ell2() const { return ell2_; }

    //! left support bound for a dual CDF function
    inline const int ell1T() const { return ell1T_; }

    //! right support bound for a dual CDF function
    inline const int ell2T() const { return ell2T_; }

    //! DKU parameter ell-tilde (3.2.10)
    inline const int ellT() const { return ellT_; }

    //! DKU abbreviation (3.2.16)
    inline const int ell() const { return ell_; }

    //! coarsest possible level
    inline const int j0() const { return (int) ceil(log(ellT()+ell2T()-1.)/log(2.0)+1); }

  protected:
    int ell1_, ell2_, ell1T_, ell2T_, ellT_, ell_;

    /*!
      an instance of the CDF basis on R
    */
    CDFBasis<d, dt> cdf_;

    Matrix<double> Alpha_, AlphaT_;
    Matrix<double> BetaL_, BetaLT_;

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
