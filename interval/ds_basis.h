// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_DS_BASIS_H
#define _WAVELETTL_DS_BASIS_H

#include <Rd/cdf_utils.h>
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
    the different possibilities for the boundary condition set Z
    (dual b.c.: complementary)
  */
  enum DSZ
    {
      empty,  // no b.c. for the primal basis
      Zero,   // primal b.c. at 0, no b.c. at 1
      One,    // no b.c. at 1, primal b.c. at 1
      ZeroOne // b.c.'s at 0 and 1
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
  template <int d, int dT, int ellTl = d-(d/2)+dT-1, int ellTr = d-(d/2)+dT-1>
  class DSBasis
  {
  public:
    /*!
      constructor
      
      You can specify the order of either the primal (s) or the dual (sT) boundary conditions at
      the left and right end of the interval [0,1]. The corresponding dual basis will then be
      constructed to fulfill the corresponding complementary boundary conditions.
    */
    DSBasis(DSZ bc = ZeroOne, const int s = 1, const int sT = 1,
	    DSBiorthogonalizationMethod bio = none);

    //! coarsest possible level, (4.20)
    inline static const int j0() { return (int) ceil(log(std::max(ellTl,ellTr)+ell2T<d,dT>()-1.)/log(2.0)+1); }
    
    //! the other freezing parameter ell, (4.11)
    inline static const int elll() { return ellTl + d - dT; }
    inline static const int ellr() { return ellTr + d - dT; }
    
    /*!
      wavelet index class
    */
    typedef IIndex<DSBasis<d,dT,ellTl,ellTr> > Index;

    /*!
      boundary indices in \Delta_j^X and \tilde\Delta_j^X (4.10),(4.14),(4.26)
     */
    inline static const int DeltaLmin() { return elll() - d; }
//     inline static const int DeltaLmax() { return elll-1-s0; }
//     inline static const int Delta0min() { return DeltaLmax()+1; }
//     inline static const int Delta0max(const int j) { return DeltaRmin(j)-1; }
//     inline static const int DeltaRmin(const int j) { return (1<<j)-ell1_-ell2_-(ellr_-1-Z[1]); }
    inline static const int DeltaRmax(const int j) { return (1<<j)-(d%2)-(ellr()-d); }

    inline static const int DeltaLTmin() { return ellTl - dT; } // == DeltaLmin()
//     inline const int DeltaLTmax() const { return ellTl_-1-ZT[0]; }
//     inline const int Delta0Tmin() const { return DeltaLTmax()+1; }
//     inline const int Delta0Tmax(const int j) const { return DeltaRTmin()-1; }
//     inline const int DeltaRTmin(const int j) const { return (1<<j)-ell1_-ell2_-(ellTr_-1-ZT[1]); }
    inline static const int DeltaRTmax(const int j) { return (1<<j)-(d%2)-(ellTr-dT); } // == DeltaRmax()

    //! size of Delta_j
    inline static const int Deltasize(const int j) { return DeltaRmax(j)-DeltaLmin()+1; }

    /*!
      boundary indices in \nabla_j
    */
    inline static const int Nablamin() { return 0; }
    inline static const int Nablamax(const int j) { return (1<<j)-1; }

  protected:
    // order of b.c. at 0 and 1
    int s0, s1, sT0, sT1;
  };
}

#include <interval/ds_basis.cpp>

#endif
