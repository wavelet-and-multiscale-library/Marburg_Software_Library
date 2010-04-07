

// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_P_BASIS_H
#define _WAVELETTL_P_BASIS_H

#define _PP_AUSWERTUNG_DER_WAVELETS   //Auskommentieren von dieser Zeile führt zum alten code

#include <iostream>
#include <algebra/vector.h>
#include <algebra/matrix.h>
#include <algebra/sparse_matrix.h>
#include <Rd/cdf_utils.h>
#include <Rd/cdf_basis.h>
#include <interval/i_index.h>

// for convenience, include also some functionality
#include <interval/p_support.h>
#include <interval/p_evaluate.h>

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
      The dual wavelet basis will have no b.c.'s in either case and can reproduce the
      full range of polynomials of order dT.

    */
    PBasis(const int s0 = 0, const int s1 = 0);

    /*!
      alternative constructor, you can specify whether first order homogeneous Dirichlet b.c.'s
      for the primal functions are set or not. The dual functions have free b.c.'s
    */
    PBasis(const bool bc_left, const bool bc_right);

    //! coarsest possible level
    inline const int j0() const { return j0_; }

    //! freezing parameters
    inline const int ellT_l() const { return (s0 >= (d-2)) ? -ell1T<d,dT>()+s0+2-d : -ell1T<d,dT>(); }
    inline const int ellT_r() const { return (s1 >= (d-2)) ? -ell1T<d,dT>()+s1+2-d : -ell1T<d,dT>(); }
    inline const int ell_l()  const { return ellT_l() + d - dT; }
    inline const int ell_r()  const { return ellT_r() + d - dT; }
    
    //! wavelet index class
    typedef IntervalIndex<PBasis<d,dT> > Index;
    
    //! size_type, for convenience
    typedef Vector<double>::size_type size_type;

    //! geometric type of the support sets
    typedef struct {
      int j;
      int k1;
      int k2;
    } Support;

    /*!
      Compute an interval 2^{-j}[k1,k2] which contains the support of a
      single primal [P] generator or wavelet \psi_\lambda.
      (j == lambda.j()+lambda.e() is neglected for performance reasons)
    */
    void support(const Index& lambda, int& k1, int& k2) const;

    //! space dimension of the underlying domain
    static const int space_dimension = 1;

    //! critical Sobolev regularity for the primal generators/wavelets
    static double primal_regularity() { return d - 0.5; }

    //! degree of polynomial reproduction for the primal generators/wavelets
    static unsigned int primal_polynomial_degree() { return d; }

    //! number of vanishing moments for the primal wavelets
    static unsigned int primal_vanishing_moments() { return dT; }

    //! read access to the primal b.c. order at x=0
    const int get_s0() const { return s0; }

    //! read access to the primal b.c. order at x=1
    const int get_s1() const { return s1; }

    //! read access to the primal b.c. order at x=0
    const int get_sT0() const { return 0; }

    //! read access to the primal b.c. order at x=1
    const int get_sT1() const { return 0; }

    void set_jmax(const int jmax) {
      jmax_ = jmax;
      setup_full_collection();
    }


    //! extremal generator indices
    inline const int DeltaLmin() const { return 1-d-ell1<d>()+s0; }
    inline const int DeltaLmax() const { return -ell1<d>(); }
    inline const int Delta0min() const { return DeltaLmax()+1; }
    inline const int Delta0max(const int j) const { return DeltaRmin(j)-1; }
    inline const int DeltaRmin(const int j) const { return (1<<j)-(d%2)-(ell_r()-1-s1); }
    inline const int DeltaRmax(const int j) const { return (1<<j)-1-ell1<d>()-s1; }

    inline const int DeltaLTmin() const
    { 
      return (s0 >= (d-2)) ? ellT_l()-dT : ellT_l()-(dT+d-2-s0);
    } // == DeltaLmin()

    inline const int DeltaLTmax() const
    { 
      return (s0 >= (d-2)) ? DeltaLTmin()+(dT-1) : DeltaLTmin()+(dT+d-3-s0);
    }

    inline const int DeltaRTmin(const int j) const
    {
      return (s1 >= (d-2)) ? DeltaRTmax(j)-(dT-1) : DeltaRTmax(j)-(dT+d-3-s1);
    }

    inline const int DeltaRTmax(const int j) const
    {
      return (s1 >= (d-2)) ? (1<<j)-(d%2)-(ellT_r()-dT) : (1<<j)-(d%2) - (ellT_r()-(dT+d-2-s1));
    } // == DeltaRmax()


    //! size of Delta_j
    inline const int Deltasize(const int j) const { return DeltaRmax(j)-DeltaLmin()+1; }

    //! boundary indices in \nabla_j
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


#ifdef _PP_AUSWERTUNG_DER_WAVELETS
    /*
      Compute the Picewiese expansion of all wavelets and Generatoren
      for given j, d
    */
    void waveletPP(const int j,Array1D<Piecewise<double> >& wavelets) const;
#endif

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
      point evaluation of (derivatives) of a single primal or dual
      generator or wavelet \psi_\lambda or \tilde\psi_\lambda
    */
    double evaluate(const unsigned int derivative, const Index& lambda, const double x) const;

    /*!
      point evaluation of (derivatives) of a single primal generator
      or wavelet \psi_\lambda at several points simultaneously
    */
    void
    evaluate
    (const unsigned int derivative,
     const Index& lambda,
     const Array1D<double>& points, Array1D<double>& values) const;


    //! read access to the internal instance of the CDF basis
    const CDFBasis<d,dT>& get_CDF_basis() const { return cdf; }

    /*!
      read access to the boundary generator expansion coefficients
      (CLA == CRA == I)
    */
    const Matrix<double>& get_CLAT() const { return CLAT; }
    const Matrix<double>& get_CRAT() const { return CRAT; }

    //! read access to the diverse refinement matrices on level j0
    const Matrix<double>& get_MLT_BB()  const { return MLT_BB; }
    const Matrix<double>& get_MRT_BB()  const { return MRT_BB; }

    const Matrix<double>& get_CLAT_BB()  const { return CLAT_BB; }
    const Matrix<double>& get_CRAT_BB()  const { return CRAT_BB; }

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

    //! get the wavelet index corresponding to a specified number
    const inline Index* get_wavelet (const int number) const {
      return &full_collection[number];
    }

    //! number of wavelets between coarsest and finest level
    const int degrees_of_freedom() const { return full_collection.size(); };

    //! Wavelets eines bestimmten levels

#ifdef _PP_AUSWERTUNG_DER_WAVELETS

    Array1D<Array1D<Piecewise<double> > > wavelets;

    void setWavelets(){
    int jmax = 12;  // wählen des max levels   12: 2,41 sek
                    //                         16: 37,0 sek
                    //  Pro Level erhähung verdopelt sich die Zeit
    wavelets.resize(jmax+1);
    for(int i = 3; i<=jmax; i++){
       waveletPP(i,wavelets[i]);
    }
    }

#endif



  protected:
    //! coarsest possible level
    int j0_;

    //! finest possible level
    int jmax_;

    //! boundary condition orders at 0 and 1
    int s0, s1;

    //! general setup routine which is shared by the different constructors
    void setup();

    //! setup full collectin of wavelets between j0_ and jmax_ as long as a jmax_ has been specified
    void setup_full_collection();

    //! collection of all wavelets between coarsest and finest level
    Array1D<Index> full_collection;

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

    //! refinement matrices for dual boundary generators BEFORE BIORTHOGONALIZATION
    Matrix<double> MLT_BB, MRT_BB;

    //! boundary blocks in Mj0T
    Matrix<double> MLT_, MRT_;

    //! Gramian matrices for the left and right generators (primal against unbiorth. dual)
    Matrix<double> GammaL, GammaR;

    //! setup the boundary blocks for the dual biorthogonalized generators
    void setup_CXT();

    //! storage for these blocks
    Matrix<double> CLT, inv_CLT, CRT, inv_CRT;

    //! refinement matrices on the coarsest level j0 and their transposed versions
    SparseMatrix<double> Mj0, Mj0T, Mj1, Mj1T;
    SparseMatrix<double> Mj0_t, Mj0T_t, Mj1_t, Mj1T_t;

    //! setup initial refinement matrices Mj0, Mj0Tp [DKU, (3.5.1), (3.5.5)]
    void setup_Mj0  (const Matrix<double>& ML,   const Matrix<double>& MR,   SparseMatrix<double>& Mj0  );
    void setup_Mj0Tp(const Matrix<double>& MLTp, const Matrix<double>& MRTp, SparseMatrix<double>& Mj0Tp);

    //! setup expansion coefficients w.r.t. the (restricted) dual CDF basis
    void setup_CXAT();

    //! setup refinement coefficients for additional dual generators,
    //! see [P], Section 4.2
    void setup_additional_duals_L(Matrix<double>& MLTp);
    void setup_additional_duals_R(Matrix<double>& MRTp);

    //! storage for these coefficients
    Matrix<double> CLAT, CRAT;
    //! representation of the (polynomial reproducing) dual boundary generators BEFORE BIORTHOGONALIZATION
    Matrix<double> CLAT_BB, CRAT_BB;

    
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
