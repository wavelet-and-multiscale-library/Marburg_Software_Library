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
  template <unsigned int d, unsigned int dt>
  class DKUBasis
  {
  public:
    //! default constructor
    DKUBasis();

    //! left support bound for a primal CDF function
    inline const int ell1() const { return -d/2; }

    //! right support bound for a primal CDF function
    inline const int ell2() const { return d - d/2; }

    //! left support bound for a dual CDF function
    inline const int ell1T() const { return ell1()-dt+1; }

    //! right support bound for a dual CDF function
    inline const int ell2T() const { return ell2()+dt-1; }

    //! DKU parameter ell-tilde (3.2.10)
    inline const int ellT() const { return ellT_; }

    //! DKU abbreviation (3.2.16)
    inline const int ell() const { return ellT()-((int)dt-(int)d); }

    //! coarsest possible level
    inline const int j0() const { return (int) ceil(log(ellT()+ell2T()-1.)/log(2.0)+1); }

    /*!
      \alpha_{\tilde\theta,r}(y) = \int x^r\tilde\theta(x-y)\,dx
     */
    double alpha(unsigned int r, int y) const;

  protected:
    /*!
      an instance of the CDF basis on R
    */
    CDFBasis<d, dt> cdf_;

    int ellT_;
    Matrix<double> Alpha_, AlphaT_;

  private:
    double alpha_(const int m, const unsigned int r) const;
    double alphaT_(const int m, const unsigned int r) const;
  };
}

#include <interval/dku.cpp>

#endif
