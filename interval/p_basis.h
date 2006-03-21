// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_P_BASIS_H
#define _WAVELETTL_P_BASIS_H

#include <algebra/vector.h>
#include <algebra/matrix.h>
#include <algebra/sparse_matrix.h>
#include <Rd/cdf_utils.h>
#include <Rd/cdf_basis.h>
#include <interval/i_index.h>

using MathTL::Vector;
using MathTL::Matrix;

namespace WaveletTL
{
  /*!
    Template class for the wavelet bases on the interval as introduced in [P].

    The primal generators are exactly those B-splines associated with the
    Schoenberg knot sequence

      t^j_{-d+1} = ... = t_0 = 0          (knot with multiplicity d at x=0)
      t^j_k = k * 2^{-j}, 1 <= k <= 2^j-1
      t^j_{2^j} = ... = t_{2^j+d-1} = 1   (knot with multiplicity d at x=1)

    i.e.

      B_{j,k}(x) = (t^j_{k+d}-t^j_k)[t^j_k,...,t^j_{k+d}](t-x)^{d-1}_+

    with supp(B_{j,k}(x) = [t^j_k, t^j_{k+d}].
    In other words, we have exactly

     d-1     left boundary splines  (k=-d+1,...,-1),
     2^j-d+1 inner splines          (k=0,...,2^j-d),
     d-1     right boundary splines (k=2^j-d+1,...,2^j-1)

    Since the primal CDF generators are centered around floor(d/2)=-ell_1,
    we perform an index shift by ell_1, i.e., we use the generators
 
      phi_{j,k}(x) = 2^{j/2} B_{j,k-ell_1}

    So, if no boundary conditions are imposed, the index of the leftmost generator
    will be 1-d+floor(d/2). See default constructor for possible b.c.'s
    
    References:
    [P] Primbs:
        Stabile biorthogonale Wavelet-Basen auf dem Intervall
	Dissertation, Univ. Duisburg-Essen, 2006
  */
  template <int d, int dT>
  class PBasis
  {
  public:
    /*!
      constructor
      
      At the moment, you may (only) specify the order of the primal (s) boundary conditions
      at the left and right end of the interval [0,1].
      For technical reasons (easing the setup of the dual generators),
      only the values si >= d-2 are currently implemented.
      In other words, linear spline wavelets may or may not fulfill Dirichlet b.c.'s.
      Higher order wavelet bases will always fulfill b.c.'s of order at least d-2.
      The dual wavelet basis will have no b.c.'s in either case and can reproduce the
      full range of polynomials of order dT.
    */
    PBasis(const int s0 = 0, const int s1 = 0);

    //! coarsest possible level
    inline const int j0() const { return j0_; }

    //! freezing parameters
    inline const int ellT_l() const { return -ell1T<d,dT>()+s0+2-d; }
    inline const int ellT_r() const { return -ell1T<d,dT>()+s1+2-d; }
    inline const int ell_l()  const { return ellT_l() + d - dT; }
    inline const int ell_r()  const { return ellT_r() + d - dT; }
    
    //! wavelet index class
    typedef IntervalIndex<PBasis<d,dT> > Index;
    
    //! size_type, for convenience
    typedef Vector<double>::size_type size_type;

    //! extremal generator indices
    inline const int DeltaLmin() const { return 1-d-ell1<d>()+s0; }
    inline const int DeltaLmax() const { return -ell1<d>(); }
    inline const int Delta0min() const { return DeltaLmax()+1; }
    inline const int Delta0max(const int j) const { return DeltaRmin(j)-1; }
    inline const int DeltaRmin(const int j) const { return (1<<j)-(d%2)-(ell_r()-1-s1); }
    inline const int DeltaRmax(const int j) const { return (1<<j)-1-ell1<d>()-s1; }

    //! size of Delta_j
    inline const int Deltasize(const int j) const { return DeltaRmax(j)-DeltaLmin()+1; }

    //! boundary indices in \nabla_j
    inline const int Nablamin() const { return 0; }
    inline const int Nablamax(const int j) const { return (1<<j)-1; }
    
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
    int s0, s1;

    //! general setup routine which is shared by the different constructors
    void setup();

    //! one instance of a CDF basis (for faster access to the primal and dual masks)
    CDFBasis<d,dT> cdf;

    //! single CDF moments \alpha_{m,r} := \int_{\mathbb R} x^r\phi(x-m)\,dx
    const double alpha(const int m, const unsigned int r) const;

     //! refinement coeffients of left dual boundary generators
    const double betaL(const int m, const unsigned int r) const;

    //! refinement coeffients of left dual boundary generators (m reversed)
    const double betaR(const int m, const unsigned int r) const;

    //! boundary blocks in Mj0
    Matrix<double> ML_, MR_;

    //! boundary blocks in Mj0T
    Matrix<double> MLT_, MRT_;

    //! Gramian matrices for the left and right generators (primal against unbiorth. dual)
    Matrix<double> GammaL, GammaR;

    //! refinement matrices on the coarsest level j0 and their transposed versions
    SparseMatrix<double> Mj0, Mj0T, Mj1, Mj1T;
    SparseMatrix<double> Mj0_t, Mj0T_t, Mj1_t, Mj1T_t;

    //! setup initial refinement matrices Mj0, Mj0Tp [DKU, (3.5.1), (3.5.5)]
    void setup_Mj0  (const Matrix<double>& ML,   const Matrix<double>& MR,   SparseMatrix<double>& Mj0  );
    void setup_Mj0Tp(const Matrix<double>& MLTp, const Matrix<double>& MRTp, SparseMatrix<double>& Mj0Tp);

    //! generator biorthogonalization matrices on level j0 and j0+1 CjT, CjpT (5.2.5)
    void setup_Cj();

    //! those matrices
    SparseMatrix<double> CjT, CjpT;
    SparseMatrix<double> inv_CjT, inv_CjpT;

    // routines for the stable completion, cf. [DKU section 4.1]
    void F(SparseMatrix<double>& FF); // (4.1.11), (4.1.14)
    void P(const Matrix<double>& ML, const Matrix<double>& MR, SparseMatrix<double>& PP); // (4.1.22)
    void GSetup(SparseMatrix<double>& A, SparseMatrix<double>& H, SparseMatrix<double>& Hinv); // (4.1.1), (4.1.13)
    void GElim (SparseMatrix<double>& A, SparseMatrix<double>& H, SparseMatrix<double>& Hinv); // elimination/factorization
    void InvertP(const SparseMatrix<double>& PP, SparseMatrix<double>& PPinv);
    double BT(const SparseMatrix<double>& A, SparseMatrix<double>& BB); // (4.1.9), (4.1.13)

    // wavelet symmetrization from [DS]
    void DS_symmetrization(SparseMatrix<double>& Mj1, SparseMatrix<double>& Mj1T);
  };
}

#include <interval/p_basis.cpp>

#endif
