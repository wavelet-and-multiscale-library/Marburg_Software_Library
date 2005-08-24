// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_DS_BASIS_H
#define _WAVELETTL_DS_BASIS_H

#include <algebra/matrix.h>
#include <algebra/sparse_matrix.h>
#include <algebra/infinite_vector.h>
#include <utils/array1d.h>

#include <Rd/cdf_utils.h>
#include <Rd/cdf_basis.h>
#include <interval/i_index.h>

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
      BernsteinSVD  // method #5 in IGPMlib, transformation to Bernstein basis plus partial SVD
    };

  /*!
    Template class for the wavelet bases on the interval as introduced in [DS].
    All formulas refer to the preprint versions of [DS] (and [DKU], where indicated).

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
  class DSBasis
  {
  public:
    /*!
      constructor
      
      You can specify the order of either the primal (s) or the dual (sT) boundary conditions at
      the left and right end of the interval [0,1]. The corresponding dual basis will then be
      constructed to fulfill the corresponding complementary boundary conditions.
    */
    DSBasis(const int s0 = 1, const int s1 = 1, const int sT0 = 0, const int sT1 = 0,
	    DSBiorthogonalizationMethod bio = SVD);

    //! freezing parameters, (4.11)
    inline const int ellT_l() const { return ell2T<d,dT>() + s0 + sT0; }
    inline const int ellT_r() const { return ell2T<d,dT>() + s1 + sT1; }
    inline const int ell_l()  const { return ellT_l() + d - dT; }
    inline const int ell_r()  const { return ellT_r() + d - dT; }
    
    //! coarsest possible level, (4.20)
    inline const int j0() const { return (int) ceil(log(std::max(ellT_l(),ellT_r())+ell2T<d,dT>()-1.)/log(2.0)+1); }
    
    /*!
      wavelet index class
    */
    typedef IIndex<DSBasis<d,dT> > Index;

    /*!
      boundary indices in \Delta_j^X and \tilde\Delta_j^X (4.10),(4.14),(4.26)
     */
    inline const int DeltaLmin() const { return ell_l()-d; }
    inline const int DeltaLmax() const { return ell_l()-1-s0; }
    inline const int Delta0min() const { return DeltaLmax()+1; }
    inline const int Delta0max(const int j) const { return DeltaRmin(j)-1; }
    inline const int DeltaRmin(const int j) const { return (1<<j)-(d%2)-(ell_r()-1-s1); }
    inline const int DeltaRmax(const int j) const { return (1<<j)-(d%2)-(ell_r()-d); }
    
    inline const int DeltaLTmin() const { return ellT_l()-dT; } // == DeltaLmin()
    inline const int DeltaLTmax() const { return ellT_l()-1-sT0; }
    inline const int Delta0Tmin() const { return DeltaLTmax()+1; }
    inline const int Delta0Tmax(const int j) const { return DeltaRTmin()-1; }
    inline const int DeltaRTmin(const int j) const { return (1<<j)-(d%2)-(ellT_r()-1-sT1); }
    inline const int DeltaRTmax(const int j) const { return (1<<j)-(d%2)-(ellT_r()-dT); } // == DeltaRmax()

    //! size of Delta_j
    inline static const int Deltasize(const int j) { return DeltaRmax(j)-DeltaLmin()+1; }
    
    /*!
      boundary indices in \nabla_j
    */
    inline static const int Nablamin() { return 0; }
    inline static const int Nablamax(const int j) { return (1<<j)-1; }

  protected:
    //! boundary condition orders at 0 and 1
    int s0, s1, sT0, sT1;

    //! the generator biorthogonalization method
    DSBiorthogonalizationMethod bio;

    //! one instance of a CDF basis (for faster access to the primal and dual masks)
    CDFBasis<d,dT> cdf;
    
    //! single moments \alpha_{m,r} := \int_{\mathbb R} x^r\phi(x-m)\,dx
    const double alpha(const int m, const unsigned int r) const;

    //! single moments \alphaT_{m,r} := \int_{\mathbb R} x^r\phiT(x-m)\,dx
    const double alphaT(const int m, const unsigned int r) const;

    //! refinement coeffients of left dual boundary generators
    const double betaL(const int m, const unsigned int r) const;

    //! refinement coeffients of left dual boundary generators
    const double betaLT(const int m, const unsigned int r) const;

    //! refinement coeffients of left dual boundary generators (m reversed)
    const double betaR(const int m, const unsigned int r) const;

    //! refinement coeffients of left dual boundary generators (m reversed)
    const double betaRT(const int m, const unsigned int r) const;


    //! compute Gramian of left and right unbiorthogonalized primal boundary functions
    void setup_GammaLR();

    //! storage for these Gramians
    Matrix<double> GammaL, GammaR;

    //! setup the boundary blocks for the generator biorthogonalization
    void setup_CX_CXT();

    //! storage for these blocks
    Matrix<double> CL, CLT, inv_CL, inv_CLT;
    Matrix<double> CR, CRT, inv_CR, inv_CRT;

    //! setup expansion coefficients w.r.t. the (restricted) CDF basis
    void setup_CXA_CXAT();

    //! storage for these coefficients
    Matrix<double> CLA, CRA, CLAT, CRAT;
  };
}

#include <interval/ds_basis.cpp>

#endif
