// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner                                     |
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
#include <interval/ds_bio.h>

// for convenience, include also some functionality
#include <interval/ds_support.h>
#include <interval/ds_evaluate.h>

namespace WaveletTL
{
  /*!
    Template class for the wavelet bases on the interval as introduced in [DKU] and [DS].
    All formulas refer to the preprint versions of [DS] (and [DKU], where indicated).

    The biorthogonalization method is used as template parameter, to control some basis properties
    also in those cases where an instance of DSBasis is constructed later (e.g., in LDomainBasis).

    References:
    [B]   Barsch:
    Adaptive Multiskalenverfahren fuer elliptische partielle Dgln. - Realisierung,
    Umsetzung und numerische Ergebnisse
    [DKU] Dahmen, Kunoth, Urban:
    Biorthogonal spline-wavelets on the interval - Stability and moment conditions
    [DS]  Dahmen, Schneider:
    Wavelets with complementary boundary conditions - Function spaces on the cube
  */
  template <int d, int dT, DSBiorthogonalizationMethod BIO = Bernstein>
  class DSBasis
  {
  public:
    /*!
      constructor
      
      You can specify the order of either the primal (s) or the dual (sT) boundary conditions at
      the left and right end of the interval [0,1]. Several combinations are possible:
      
      si=sTi=0  : no b.c., original [DKU] construction
      si=s>0=sTi: primal basis with b.c., dual basis without b.c.
      si=0<s=sTi: primal basis without b.c., dual basis with b.c.
      si=sTi>0  : b.c. for primal and dual basis, like [CTU] (not recommended, loss of approximation order)

      The default basis is the original [DKU] construction without any boundary conditions.
    */
    DSBasis(const int s0 = 0, const int s1 = 0, const int sT0 = 0, const int sT1 = 0);

    /*!
      alternative constructor, you can specify whether first order homogeneous Dirichlet b.c.'s
      for the primal functions are set or not. The dual functions have free b.c.'s.
    */
    DSBasis(const bool bc_left, const bool bc_right);

    //! freezing parameters, (4.11)
    inline const int ellT_l() const { return ell2T<d,dT>() + s0 + sT0; }
    inline const int ellT_r() const { return ell2T<d,dT>() + s1 + sT1; }
    inline const int ell_l()  const { return ellT_l() + d - dT; }
    inline const int ell_r()  const { return ellT_r() + d - dT; }
    
    //! coarsest possible level, (4.20)
    inline const int j0() const { return j0_; }
    
    /*!
      wavelet index class
    */
    typedef IntervalIndex<DSBasis<d,dT,BIO> > Index;

    /*!
      size_type, for convenience
    */
    typedef Vector<double>::size_type size_type;

    /*!
      geometric type of the support sets
    */
    typedef struct {
      int j;
      int k1;
      int k2;
    } Support;

    /*!
      Compute an interval 2^{-j}[k1,k2] which contains the support of a
      single primal DS generator or wavelet \psi_\lambda.
      (j == lambda.j()+lambda.e() is neglected for performance reasons)
    */
    void support(const Index& lambda, int& k1, int& k2) const;
    
    /*!
      space dimension of the underlying domain
    */
    static const int space_dimension = 1;

    /*!
      critical Sobolev regularity for the primal generators/wavelets
    */
    static double primal_regularity() { return d - 0.5; }

    /*!
      degree of polynomial reproduction for the primal generators/wavelets
    */
    static unsigned int primal_polynomial_degree() { return d; }

    /*!
      number of vanishing moments for the primal wavelets
    */
    static unsigned int primal_vanishing_moments() { return dT; }

    //! read access to the primal b.c. order at x=0
    const int get_s0() const { return s0; }

    //! read access to the primal b.c. order at x=1
    const int get_s1() const { return s1; }

    //! read access to the primal b.c. order at x=0
    const int get_sT0() const { return sT0; }

    //! read access to the primal b.c. order at x=1
    const int get_sT1() const { return sT1; }

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
    inline const int Deltasize(const int j) const { return DeltaRmax(j)-DeltaLmin()+1; }
    
    /*!
      boundary indices in \nabla_j
    */
    inline const int Nablamin() const { return 0; }
    inline const int Nablamax(const int j) const { return (1<<j)-1; }

    //! size of Nabla_j
    inline const int Nablasize(const int j) const { return 1<<j; }

    //! index of first (leftmost) generator on level j >= j0
    Index first_generator(const int j) const;

    //! index of last (rightmost) generator on level j >= j0
    Index last_generator(const int j) const;

    //! index of first (leftmost) wavelet on level j >= j0
    Index first_wavelet(const int j) const;

    //! index of last (rightmost) wavelet on level j >= j0
    Index last_wavelet(const int j) const;

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

    /*!
      read access to the internal instance of the CDF basis
    */
    const CDFBasis<d,dT>& get_CDF_basis() const { return cdf; }

    /*!
      read access to the boundary generator expansion coefficients
    */
    const Matrix<double>& get_CLA() const { return CLA; }
    const Matrix<double>& get_CRA() const { return CRA; }
    const Matrix<double>& get_CLAT() const { return CLAT; }
    const Matrix<double>& get_CRAT() const { return CRAT; }

    /*!
      read access to the diverse refinement matrices on level j0
    */
    const SparseMatrix<double>& get_Mj0()  const { return Mj0; }
    const SparseMatrix<double>& get_Mj0T() const { return Mj0T; }
    const SparseMatrix<double>& get_Mj1()  const { return Mj1; }
    const SparseMatrix<double>& get_Mj1T() const { return Mj1T; }
    const SparseMatrix<double>& get_Mj0_t()  const { return Mj0_t; }
    const SparseMatrix<double>& get_Mj0T_t() const { return Mj0T_t; }
    const SparseMatrix<double>& get_Mj1_t()  const { return Mj1_t; }
    const SparseMatrix<double>& get_Mj1T_t() const { return Mj1T_t; }

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

  protected:
    //! coarsest possible level
    int j0_;

    //! boundary condition orders at 0 and 1
    int s0, s1, sT0, sT1;

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

    //! general setup routine which is shared by the different constructors
    void setup();

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

    //! generator biorthogonalization matrices on level j0 and j0+1 Cj, CjT, Cjp, CjpT (5.2.5)
    void setup_Cj();

    //! those matrices
    SparseMatrix<double> Cj, CjT, Cjp, CjpT;
    SparseMatrix<double> inv_Cj, inv_CjT, inv_Cjp, inv_CjpT;

    //! setup refinement matrix blocks ML, MR (3.5.2)
    Matrix<double> ML() const;
    Matrix<double> MR() const;

    //! setup refinement matrix blocks MLTs, MRTs (3.5.6)
    Matrix<double> MLTp() const;
    Matrix<double> MRTp() const;

    //! setup initial refinement matrices Mj0, Mj0Tp (3.5.1), (3.5.5)
    void setup_Mj0  (const Matrix<double>& ML,   const Matrix<double>& MR,   SparseMatrix<double>& Mj0  );
    void setup_Mj0Tp(const Matrix<double>& MLTp, const Matrix<double>& MRTp, SparseMatrix<double>& Mj0Tp);

    //! refinement matrices on the coarsest level j0 and their transposed versions
    SparseMatrix<double> Mj0, Mj0T, Mj1, Mj1T;     
    SparseMatrix<double> Mj0_t, Mj0T_t, Mj1_t, Mj1T_t;

    // routines for the stable completion, [DKU section 4.1]
    void F(SparseMatrix<double>& FF); // (4.1.11), (4.1.14)
    void P(const Matrix<double>& ML, const Matrix<double>& MR, SparseMatrix<double>& PP); // (4.1.22)
    void GSetup(SparseMatrix<double>& A, SparseMatrix<double>& H, SparseMatrix<double>& Hinv); // (4.1.1), (4.1.13)
    void GElim (SparseMatrix<double>& A, SparseMatrix<double>& H, SparseMatrix<double>& Hinv); // elimination/factorization
    void InvertP(const SparseMatrix<double>& PP, SparseMatrix<double>& PPinv);
    void BT(const SparseMatrix<double>& A, SparseMatrix<double>& BB); // (4.1.9), (4.1.13)

    // wavelet symmetrization from [DS]
    void DS_symmetrization(SparseMatrix<double>& Mj1, SparseMatrix<double>& Mj1T);
  };
}

#include <interval/ds_basis.cpp>

#endif
