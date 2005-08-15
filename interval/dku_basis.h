// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_DKU_BASIS_H
#define _WAVELETTL_DKU_BASIS_H

#include <iostream>
#include <cmath>

#include <algebra/matrix.h>
#include <algebra/sparse_matrix.h>
#include <algebra/infinite_vector.h>
#include <utils/array1d.h>

#include <Rd/cdf_utils.h>
#include <Rd/cdf_basis.h>
#include <interval/dku_index.h>

using MathTL::Matrix;

namespace WaveletTL
{
  /*!
    biorthogonalization methods for the DKU basis, see, e.g., [B]

    The methods 'partialSVD' and 'BernsteinSVD' enable the [DKU]/[DS] boundary treatment:
    at the boundary, exactly one generator and one wavelet does not vanish,
    which can be modified to satisfy homogeneous boundary conditions for the primal
    basis.
  */
  enum DKUBiorthogonalizationMethod
    {
      none,         // method #1 in IGPMlib, C_L = I
      SVD,          // method #2 in IGPMlib, Gamma_L = U*S*V, C_L = S^{-1/2}U^T, C_L_T = S^{-1/2}V
      Bernstein,    // method #3 in IGPMlib, transformation to Bernstein basis on [0,b]
      partialSVD,   // method #4 in IGPMlib, partial SVD
      BernsteinSVD  // method #5 in IGPMlib, transformation to Bernstein basis plus partial SVD
    };

  /*!
    Template class for the wavelet bases on the interval as introduced in [DKU], [DS].
    All formulas refer to the preprint version of [DKU],
    except those explicitly denoted by [DS].

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

      You can toggle Dirichlet boundary conditions for the primal basis with the parameters
      bc_left/bc_right. The dual basis will then be constructed to fulfill complementary boundary
      conditions, see [DS].
    */
    explicit DKUBasis(bool bc_left = false,
		      bool bc_right = false,
		      DKUBiorthogonalizationMethod bio = BernsteinSVD);

    //! coarsest possible level
    inline const int j0() const { return (int) ceil(log(ellT_l+ell2T_-1.)/log(2.0)+1); }

    /*!
      boundary indices in \Delta_j^X and \tilde\Delta_j^X (3.2.17)
     */
    inline const int DeltaLmin() const { return ell_l-d; }
    inline const int DeltaLmax() const { return ell_l-1-Z[0]; }
    inline const int Delta0min() const { return DeltaLmax()+1; }
    inline const int Delta0max(const int j) const { DeltaRmin(j)-1; }
    inline const int DeltaRmin(const int j) const { return (1<<j)-ell1_-ell2_-(ell_r-1-Z[1]); }
    inline const int DeltaRmax(const int j) const { return (1<<j)-ell1_-ell2_-(ell_r-d); }

    inline const int DeltaLTmin() const { return ellT_l-dT; } // == DeltaLmin()
    inline const int DeltaLTmax() const { return ellT_l-1-ZT[0]; }
    inline const int Delta0Tmin() const { return DeltaLTmax()+1; }
    inline const int Delta0Tmax(const int j) const { return DeltaRTmin()-1; }
    inline const int DeltaRTmin(const int j) const { return (1<<j)-ell1_-ell2_-(ellT_r-1-ZT[1]); }
    inline const int DeltaRTmax(const int j) const { return (1<<j)-ell1_-ell2_-(ellT_r-dT); } // == DeltaRmax()

    //! size of Delta_j
    inline const int Deltasize(const int j) const { return DeltaRmax(j)-DeltaLmin()+1; }

    /*!
      boundary indices in \nabla_j
    */
    inline const int Nablamin() const { return 0; }
    inline const int Nablamax(const int j) const { return (1<<j)-1; }

    /*!
      wavelet index class
    */
    typedef DKUIndex<d, dT> Index;

    /*!
      first (leftmost) generator on scale j >= j0
    */
    Index firstGenerator(const int j) const;

    /*!
      last (rightmost) generator on scale j >= j0
    */
    Index lastGenerator(const int j) const;

    /*!
      first (leftmost) wavelet on scale j >= j0
    */
    Index firstWavelet(const int j) const;

    /*!
      last (rightmost) wavelet on scale j >= j0
    */
    Index lastWavelet(const int j) const;

    //! DECOMPOSE routine, simple version
    /*!
      Constructs for a given single wavelet index lambda a coefficient set c,
      such that
        \psi_lambda = \sum_{\lambda'}c_{\lambda'}\psi_{\lambda'}
      where the multiscale decomposition starts with the coarsest
      generator level jmin.
     */
    void decompose_1(const Index& lambda, const int jmin,
		     InfiniteVector<double, Index>& c) const;

    //! dual DECOMPOSE routine, simple version
    /*!
      Constructs for a given single wavelet index lambda a coefficient set c,
      such that
        \tilde\psi_lambda = \sum_{\lambda'}c_{\lambda'}\tilde\psi_{\lambda'}
      where the multiscale decomposition starts with the coarsest
      generator level jmin.
     */
    void decompose_t_1(const Index& lambda, const int jmin,
		       InfiniteVector<double, Index>& c) const;

    //! DECOMPOSE routine, full version
    /*!
      constructs for a given coefficient set c another one v with level >= jmin,
      such that
        \sum_{\lambda}c_\lambda\psi_lambda = \sum_{\lambda'}v_{\lambda'}\psi_{\lambda'}
    */
    void decompose(const InfiniteVector<double, Index>& c, const int jmin,
		   InfiniteVector<double, Index>& v) const;

    //! dual DECOMPOSE routine, full version
    /*!
      constructs for a given coefficient set c another one v with level >= jmin,
      such that
        \sum_{\lambda}c_\lambda\tilde\psi_lambda = \sum_{\lambda'}d_{\lambda'}\tilde\psi_{\lambda'}
    */
    void decompose_t(const InfiniteVector<double, Index>& c, const int jmin,
		     InfiniteVector<double, Index>& v) const;

    //! RECONSTRUCT routine, simple version
    /*!
      Constructs for a given single wavelet index lambda a coefficient set c,
      such that
        \psi_lambda = \sum_{\lambda'}c_{\lambda'}\psi_{\lambda'}
      where always |\lambda'|>=j
     */
    void reconstruct_1(const Index& lambda, const int j,
		       InfiniteVector<double, Index>& c) const;

    //! RECONSTRUCT routine, full version
    /*!
      Constructs for a given coefficient set c another one v,
      such that
        \sum_{\lambda}c_\lambda\psi_lambda = \sum_{\lambda'}v_{\lambda'}\psi_{\lambda'}
      where always |\lambda'|>=j
    */
    void reconstruct(const InfiniteVector<double, Index>& c, const int j,
		     InfiniteVector<double, Index>& v) const;

    //! dual RECONSTRUCT routine, simple version
    /*!
      Constructs for a given single wavelet index lambda a coefficient set c,
      such that
        \tilde\psi_lambda = \sum_{\lambda'}c_{\lambda'}\tilde\psi_{\lambda'}
      where always |\lambda'|>=j
     */
    void reconstruct_t_1(const Index& lambda, const int j,
			 InfiniteVector<double, Index>& c) const;

    //! dual RECONSTRUCT routine, full version
    /*!
      Constructs for a given coefficient set c another one v,
      such that
        \sum_{\lambda}c_\lambda\tilde\psi_\lambda = \sum_{\lambda'}v_{\lambda'}\tilde\psi_{\lambda'}
      where always |\lambda'|>=j
    */
    void reconstruct_t(const InfiniteVector<double, Index>& c, const int j,
		       InfiniteVector<double, Index>& v) const;

    //! setup the refinement matrix M_{j,0} for a given level j
    void assemble_Mj0(const int j, SparseMatrix<double>& mj0) const;

    //! setup the refinement matrix \tilde M_{j,0} for a given level j
    void assemble_Mj0T(const int j, SparseMatrix<double>& mj0T) const;

    //! setup the refinement matrix M_{j,1} for a given level j
    void assemble_Mj1(const int j, SparseMatrix<double>& mj1) const;

    //! setup the refinement matrix \tilde M_{j,1} for a given level j
    void assemble_Mj1T(const int j, SparseMatrix<double>& mj1T) const;

    //! setup the transposed refinement matrix M_{j,0} for a given level j
    void assemble_Mj0_t(const int j, SparseMatrix<double>& mj0_t) const;

    //! setup the transposed refinement matrix \tilde M_{j,0} for a given level j
    void assemble_Mj0T_t(const int j, SparseMatrix<double>& mj0T_t) const;

    //! setup the transposed refinement matrix M_{j,1} for a given level j
    void assemble_Mj1_t(const int j, SparseMatrix<double>& mj1_t) const;

    //! setup the transposed refinement matrix \tilde M_{j,1} for a given level j
    void assemble_Mj1T_t(const int j, SparseMatrix<double>& mj1T_t) const;

    //! compute single rows of these matrices on higher levels than j0
    void Mj0_get_row   (const int j, const Vector<double>::size_type row,
			InfiniteVector<double, Vector<double>::size_type>& v) const;
    void Mj0T_get_row  (const int j, const Vector<double>::size_type row,
			InfiniteVector<double, Vector<double>::size_type>& v) const;
    void Mj1_get_row   (const int j, const Vector<double>::size_type row,
			InfiniteVector<double, Vector<double>::size_type>& v) const;
    void Mj1T_get_row  (const int j, const Vector<double>::size_type row,
			InfiniteVector<double, Vector<double>::size_type>& v) const;
    void Mj0_t_get_row (const int j, const Vector<double>::size_type row,
			InfiniteVector<double, Vector<double>::size_type>& v) const;
    void Mj0T_t_get_row(const int j, const Vector<double>::size_type row,
			InfiniteVector<double, Vector<double>::size_type>& v) const;
    void Mj1_t_get_row (const int j, const Vector<double>::size_type row,
			InfiniteVector<double, Vector<double>::size_type>& v) const;
    void Mj1T_t_get_row(const int j, const Vector<double>::size_type row,
			InfiniteVector<double, Vector<double>::size_type>& v) const;
    
    /*!
      Evaluate a single primal/dual generator or wavelet \psi_\lambda
      on a dyadic subgrid of [0,1].
     */
    SampledMapping<1> evaluate(const Index& lambda,
			       const bool primal,
			       const int resolution) const;

    /*!
      Evaluate an arbitrary linear combination of primal or dual
      wavelets on a dyadic subgrid of [0,1].
    */
    SampledMapping<1> evaluate(const InfiniteVector<double, Index>& coeffs,
			       const bool primal,
			       const int resolution) const;

  protected:
    int ell1_, ell2_, ell1T_, ell2T_;
    int ell_l, ell_r, ellT_l, ellT_r;
    DKUBiorthogonalizationMethod bio_;
    Array1D<int> Z, ZT;

    CDFBasis<d, dT> cdf_;

    Matrix<double> Alpha_, AlphaT_;
    Matrix<double> BetaL_, BetaLT_, BetaR_, BetaRT_;
    Matrix<double> GammaL_, GammaR_;
    Matrix<double> CL_, CLT_, inv_CL_, inv_CLT_;
    Matrix<double> CR_, CRT_, inv_CR_, inv_CRT_;
    
    /*!
      coefficients combining (3.2.25), (3.2.26) and the biorthogonalization
      (represent biorthogonalized boundary generators as linear combinations
      of restricted ones from the line)
    */
    Matrix<double> CLA_;  // left primal boundary generators
    Matrix<double> CRA_;  // right primal boundary generators
    Matrix<double> CLAT_; // left dual boundary generators
    Matrix<double> CRAT_; // right dual boundary generators

    //! storage for transformation matrices on level j0 and j0+1
    SparseMatrix<double> Cj_, CjT_, Cjp_, CjpT_;
    SparseMatrix<double> inv_Cj_, inv_CjT_, inv_Cjp_, inv_CjpT_;

    SparseMatrix<double> Mj0_, Mj0T_, Mj1_, Mj1T_;     // refinement matrices on the coarsest level j0()
    SparseMatrix<double> Mj0_t, Mj0T_t, Mj1_t, Mj1T_t; // transposed versions, for performance reasons

    void setup_Alpha();
    void setup_AlphaT();
    void setup_BetaL();
    void setup_BetaLT();
    void setup_BetaR();
    void setup_BetaRT();
    void setup_GammaLR();
    void setup_CX_CXT();
    void setup_CXA_CXAT();

    // Cj, CjT, Cjp, CjpT (5.2.5)
    void setup_Cj();

    // ML, MR (3.5.2)
    Matrix<double> ML() const;
    Matrix<double> MR() const;

    // MLTs, MRTs (3.5.6)
    Matrix<double> MLTp() const;
    Matrix<double> MRTp() const;

    // Mj0, Mj0Tp (3.5.1), (3.5.5)
    void setup_Mj0  (const Matrix<double>& ML,   const Matrix<double>& MR,   SparseMatrix<double>& Mj0  );
    void setup_Mj0Tp(const Matrix<double>& MLTp, const Matrix<double>& MRTp, SparseMatrix<double>& Mj0Tp);

    // routines for the stable completion, [DKU section 4.1]
    void F(SparseMatrix<double>& FF); // (4.1.11), (4.1.14)
    void P(const Matrix<double>& ML, const Matrix<double>& MR, SparseMatrix<double>& PP); // (4.1.22)
    void GSetup(SparseMatrix<double>& A, SparseMatrix<double>& H, SparseMatrix<double>& Hinv); // (4.1.1), (4.1.13)
    void GElim (SparseMatrix<double>& A, SparseMatrix<double>& H, SparseMatrix<double>& Hinv); // elimination/factorization
    void InvertP(const SparseMatrix<double>& PP, SparseMatrix<double>& PPinv);
    void BT(const SparseMatrix<double>& A, SparseMatrix<double>& BB); // (4.1.9), (4.1.13)

    // boundary generator/wavelet modifications [CTU],[DS]
    void boundary_modifications(SparseMatrix<double>& Mj1, SparseMatrix<double>& Mj1T);
  };
}

#include <interval/dku_basis.cpp>

#endif
