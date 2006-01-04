// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_CDD1_H
#define _WAVELETTL_CDD1_H

#include <set>
#include <algebra/infinite_vector.h>
#include <adaptive/compression.h>

using std::set;

namespace WaveletTL
{
  /*
    An adaptive, residual-based solver for the infinite-dimensional problem

      Au = F,

    as developed in [CDD1] and [BB+], where A is assumed to be s.p.d.
    Given the problem and a target accuracy epsilon,
    the algorithm constructs a coefficient vector u_epsilon, such that
    
      ||u-u_epsilon|| <= epsilon.

    You can specify a maximal level jmax for the internal APPLY calls.

    References:
    [BB+]  Barinka/Barsch/Charton/Cohen/Dahlke/Dahmen/Urban:
           Adaptive Wavelet Schemes For Elliptic Problems: Implementation and Numerical Experiments
    [CDD1] Cohen/Dahmen/DeVore:
           Adaptive Wavelet Methods II - Beyond the Elliptic Case
  */

  /*!
    the routine ALGORITHMc from [BB+]
   */
  template <class PROBLEM>
  void CDD1_SOLVE(const PROBLEM& P, const double epsilon,
		  InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& u_epsilon,
		  const int jmax = 99,
		  const CompressionStrategy strategy = St04a);
  /*!
    the routine ALGORITHMc from [BB+],
    with a given initial guess for u_epsilon
   */
  template <class PROBLEM>
  void CDD1_SOLVE(const PROBLEM& P, const double epsilon,
		  const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& guess,
		  InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& u_epsilon,
		  const int jmax = 99,
		  const CompressionStrategy strategy = St04a);

  /*!
    the parameters chosen or computed in the INIT phase of ALGORITHMc
   */
  typedef struct {
    double c1, c2;
    double kappa;
    double gamma;
    double F;
    double q0, q1, q2, q3, q4;
    unsigned int K;
    double theta, theta_bar;
  } CDD1Parameters;
  
  /*!
    NPROG:
    Given an approximation v (the support of which is contained in Lambda)
    to the exact Galerkin solution u of Au = F with
    ||u-v||_2 <= delta, compute a new approximation v_hat supported in Lambda_hat,
    such that ||u-v_hat||_2 <= delta/2.
    An approximate residual r_hat as well as the last iterand ubar before the final thresholding
    are also returned.
   */
  template <class PROBLEM>
  void NPROG(const PROBLEM& P, const CDD1Parameters& params,
	     const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& F,
	     const set<typename PROBLEM::WaveletBasis::Index>& Lambda,
	     const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& v,
	     const double delta,
	     InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& v_hat,
	     set<typename PROBLEM::WaveletBasis::Index>& Lambda_hat,
	     InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& r_hat,
	     InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& u_Lambda_k,
	     const int jmax = 99,
	     const CompressionStrategy strategy = St04a);

  /*!
    GALERKIN:
    Given an approximation v to the exact Galerkin solution u_Lambda of Au = F w.r.t. the
    index set Lambda, such that ||u_Lambda-v||_2 <= delta, and a target accuracy eta,
    compute an approximation u_bar to u_Lambda which is supported on Lambda and satisfies
    ||u_bar-u_Lambda||_2 <= eta.
  */
  template <class PROBLEM>
  void GALERKIN(const PROBLEM& P, const CDD1Parameters& params,
		const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& F,
		const set<typename PROBLEM::WaveletBasis::Index>& Lambda,
 		const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& v,
 		const double delta,
		const double eta,
 		InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& ubar,
		const int jmax = 99,
		const CompressionStrategy strategy = St04a);

  /*!
    NGROW:
    Given a set Lambda, an initial approximation ubar (supported in Lambda) to the
    Galerkin solutin u_Lambda of Au = F, calculate an approximate residual r with
      ||r-r_Lambda||_2 <= xi_1 + xi_2 + c_2 * ||ubar-u_Lambda||_2
    and a new index set Lambda_tilde\supset Lambda as small as possible such that
      ||P_{Lambda_tilde\setminus Lambda}r||_2 >= gamma * ||r||_2
  */
  template <class PROBLEM>
  void NGROW(const PROBLEM& P, const CDD1Parameters& params,
	     const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& F,
	     const set<typename PROBLEM::WaveletBasis::Index>& Lambda,
 	     const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& ubar,
 	     const double xi1,
 	     const double xi2,
 	     set<typename PROBLEM::WaveletBasis::Index>& Lambda_tilde,
 	     InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& r,
	     const int jmax = 99,
	     const CompressionStrategy strategy = St04a);

  /*!
    INRESIDUAL:
    Given an index set Lambda, an approximation v to the exact Galerkin solution
    u_Lambda of Au = F, calculate an approximate INternal residual r, such that
      ||r - (A_Lambda v - P_Lambda f)||_2 <= eta_1 + eta_2
  */
  template <class PROBLEM>
  void INRESIDUAL(const PROBLEM& P, const CDD1Parameters& params,
		  const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& F,
		  const set<typename PROBLEM::WaveletBasis::Index>& Lambda,
		  const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& v,
		  const double eta1,
		  const double eta2,
		  InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& r,
		  const int jmax = 99,
		  const CompressionStrategy strategy = St04a);

  /*!
    NRESIDUAL:
    Given an index set Lambda, an approximation v to the exact Galerkin solution
    u_Lambda of Au = F, calculate an approximate residual r (not necessarily supported
    in J\Lambda), such that
      ||r - r_Lambda||_2 <= eta_1 + eta_2 + c_2 * ||v-u_Lambda||_2
    The routine also returns the support set Lambda_tilde of the approximate residual r.
  */
  template <class PROBLEM>
  void NRESIDUAL(const PROBLEM& P, const CDD1Parameters& params,
		 const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& F,
		 const set<typename PROBLEM::WaveletBasis::Index>& Lambda,
		 const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& v,
		 const double eta1,
		 const double eta2,
		 InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& r,
		 set<typename PROBLEM::WaveletBasis::Index>& Lambda_tilde,
		 const int jmax = 99,
		 const CompressionStrategy strategy = St04a);
}

#include <adaptive/cdd1.cpp>

#endif
