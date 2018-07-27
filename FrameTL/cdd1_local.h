// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2007                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _FRAMETL_CDD1_LOCALH
#define _FRAMETL_CDD1_LOCALH

#include <set>
#include <algebra/infinite_vector.h>
#include <adaptive/compression.h>

using std::set;
using WaveletTL::CompressionStrategy;
using WaveletTL::St04a;
using WaveletTL::CDD1;

/*!\namespace FrameTL
  The namespace FrameTL.
 */
namespace FrameTL
 {

   /*! \file cdd1_local.h
     \brief (The routines are taken from WaveletTL/adaptive/cdd1.{cpp,h}.
     They are here adapted in order to be able to use them as a local solver for a frame domain decomposition method.)

     An adaptive, residual-based solver for the infinite-dimensional problem

     \f$Au = F\f$,

     as developed in [CDD1] and [BB+], where A is assumed to be s.p.d.
     Given the problem and a target accuracy epsilon,
     the algorithm constructs a coefficient vector u_epsilon, such that

     \f$||u-u_\epsilon|| \leq \epsilon\f$.

     You can specify a maximal level jmax for the internal APPLY calls.

     References:
     [BB+]  Barinka/Barsch/Charton/Cohen/Dahlke/Dahmen/Urban:
     Adaptive Wavelet Schemes For Elliptic Problems: Implementation and Numerical Experiments
     [CDD1] Cohen/Dahmen/DeVore:
     Adaptive Wavelet Methods II - Beyond the Elliptic Case.
   */

  /*! \fn template <class PROBLEM> void CDD1_LOCAL_SOLVE(const PROBLEM& P, const int patch, const double epsilon,
    const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& guess,
    InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& u_epsilon,
    const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& v_k,
    const int jmax = 99,
    const CompressionStrategy strategy = St04a)
    \brief The routine ALGORITHMc from [BB+], with a given initial guess for u_epsilon.

    \param P The global discrete (cached) problem.
    \param patch The patch on which the local auxiliary problem has to be solved.
    \param epsilon The target accuracy.
    \param gues The initial gues from which the adaptive Galerkin scheme is started.
    \param u_epsilon The final approximation.
    \param v_k For the case of the multiplicative Schwarz method a la Stevenson/Werner 2009,
    v_k have to be the coefficients u_k^{(j)}, j\neq i f from eq. (6.1.23) in
    Manuel's PhD thesis. In the additive case, v_k are the frame coeffients of the function in the first
    argument of the bilinear form a(.,.) in the right-hand side of the local problem in algorithm AddSchw
    on page 165 of Manuel's thesis.
    param jmax The maximal alevel of resolution.
    param strategy The compression strategy (the way to create the matrix A_j from Definition 4.2 in Manuels
    PhD thesis).}
   */

  /*! \fn template <class PROBLEM> void CDD1_LOCAL_SOLVE(const PROBLEM& P, const int patch,
    const double epsilon,
    const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& guess,
    InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& u_epsilon,
    const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& u_k_very_sparse,
    const double c1,
    const double c2,
    const int jmax = 99,
    const CompressionStrategy strategy = St04a)
    \brief The routine ALGORITHMc from [BB+], with a given initial guess for u_epsilon and for the parameters c1,c2.

  */

  /*! \struct  typedef struct CDD1Parameters
    \brief The parameters chosen or computed in the INIT phase of ALGORITHMc.
  */

  /*! \fn  FrameTL::template <class PROBLEM> void NPROG(const PROBLEM& P,
    const int patch,
    const CDD1Parameters& params,
    const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& F,
    const set<typename PROBLEM::WaveletBasis::Index>& Lambda,
    const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& v,
    const double delta,
    InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& v_hat,
    set<typename PROBLEM::WaveletBasis::Index>& Lambda_hat,
    InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& r_hat,
    InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& u_Lambda_k,
    const int jmax = 99,
    const CompressionStrategy strategy = St04a)
    \brief NPROG.

    Given an approximation v (the support of which is contained in Lambda)
    to the exact Galerkin solution u of Au = F with
    ||u-v||_2 <= delta, compute a new approximation v_hat supported in Lambda_hat,
    such that ||u-v_hat||_2 <= delta/2.
    An approximate residual r_hat as well as the last iterand ubar before the final thresholding
    are also returned.
  */

  /*! \fn  template <class PROBLEM> void GALERKIN(const PROBLEM& P,
    const int patch,
    const CDD1Parameters& params,
    const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& F,
    const set<typename PROBLEM::WaveletBasis::Index>& Lambda,
    const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& v,
    const double delta,
    const double eta,
    InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& ubar,
    const int jmax = 99,
    const CompressionStrategy strategy = St04a)
    \brief GALERKIN.

    Given an approximation v to the exact Galerkin solution u_Lambda of Au = F w.r.t. the
    index set Lambda, such that ||u_Lambda-v||_2 <= delta, and a target accuracy eta,
    compute an approximation u_bar to u_Lambda which is supported on Lambda and satisfies
    ||u_bar-u_Lambda||_2 <= eta.
  */

  /*! \fn  template <class PROBLEM> void NGROW(const PROBLEM& P,
    const int patch,
    const CDD1Parameters& params,
    const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& F,
    const set<typename PROBLEM::WaveletBasis::Index>& Lambda,
    const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& ubar,
    const double xi1,
    const double xi2,
    set<typename PROBLEM::WaveletBasis::Index>& Lambda_tilde,
    InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& r,
    const int jmax = 99,
    const CompressionStrategy strategy = St04a);
    \brief NGROW.

    Given a set Lambda, an initial approximation ubar (supported in Lambda) to the
    Galerkin solutin u_Lambda of Au = F, calculate an approximate residual r with
    ||r-r_Lambda||_2 <= xi_1 + xi_2 + c_2 * ||ubar-u_Lambda||_2
    and a new index set Lambda_tilde\supset Lambda as small as possible such that
    ||P_{Lambda_tilde\setminus Lambda}r||_2 >= gamma * ||r||_2
  */

  /*! \fn  template <class PROBLEM> void INRESIDUAL(const PROBLEM& P,
    const int patch,
    const CDD1Parameters& params,
    const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& F,
    const set<typename PROBLEM::WaveletBasis::Index>& Lambda,
    const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& v,
    const double eta1,
    const double eta2,
    InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& r,
    const int jmax = 99,
    const CompressionStrategy strategy = St04a)
    \brief INRESIDUAL.

    Given an index set Lambda, an approximation v to the exact Galerkin solution
    u_Lambda of Au = F, calculate an approximate INternal residual r, such that
    ||r - (A_Lambda v - P_Lambda f)||_2 <= eta_1 + eta_2
  */

  /*! \fn  template <class PROBLEM> void NRESIDUAL(const PROBLEM& P,
    const int patch,
    const CDD1Parameters& params,
    const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& F,
    const set<typename PROBLEM::WaveletBasis::Index>& Lambda,
    const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& v,
    const double eta1,
    const double eta2,
    InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& r,
    set<typename PROBLEM::WaveletBasis::Index>& Lambda_tilde,
    const int jmax = 99,
    const CompressionStrategy strategy = St04a)
    \brief NRESIDUAL.

    Given an index set Lambda, an approximation v to the exact Galerkin solution
    u_Lambda of Au = F, calculate an approximate residual r (not necessarily supported
    in J\Lambda), such that
    ||r - r_Lambda||_2 <= eta_1 + eta_2 + c_2 * ||v-u_Lambda||_2
    The routine also returns the support set Lambda_tilde of the approximate residual r.
  */

  template <class PROBLEM>
  void CDD1_LOCAL_SOLVE(const PROBLEM& P,
			const int patch,
			const double epsilon,
			const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& guess,
			InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& u_epsilon,
			const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& v_k,
			const int jmax = 99,
			const CompressionStrategy strategy = St04a);
  template <class PROBLEM>
  void CDD1_LOCAL_SOLVE(const PROBLEM& P,
            const int patch,
			const double epsilon,
			const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& guess,
			InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& u_epsilon,
			const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& u_k_very_sparse,
			const double c1,
			const double c2,
			const int jmax = 99,
			const CompressionStrategy strategy = St04a);
  typedef struct {
    double c1, c2;
    double kappa;
    double gamma;
    double F;
    double q0, q1, q2, q3, q4;
    unsigned int K;
    double theta, theta_bar;
  } CDD1Parameters;

  template <class PROBLEM>
  void NPROG(const PROBLEM& P,
	     const int patch,
	     const CDD1Parameters& params,
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

  template <class PROBLEM>
  void GALERKIN(const PROBLEM& P,
		const int patch,
		const CDD1Parameters& params,
		const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& F,
		const set<typename PROBLEM::WaveletBasis::Index>& Lambda,
 		const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& v,
 		const double delta,
		const double eta,
 		InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& ubar,
		const int jmax = 99,
		const CompressionStrategy strategy = St04a);

  template <class PROBLEM>
  void NGROW(const PROBLEM& P,
	     const int patch,
	     const CDD1Parameters& params,
	     const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& F,
	     const set<typename PROBLEM::WaveletBasis::Index>& Lambda,
 	     const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& ubar,
 	     const double xi1,
 	     const double xi2,
 	     set<typename PROBLEM::WaveletBasis::Index>& Lambda_tilde,
 	     InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& r,
	     const int jmax = 99,
	     const CompressionStrategy strategy = St04a);

  template <class PROBLEM>
  void INRESIDUAL(const PROBLEM& P,
		  const int patch,
		  const CDD1Parameters& params,
		  const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& F,
		  const set<typename PROBLEM::WaveletBasis::Index>& Lambda,
		  const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& v,
		  const double eta1,
		  const double eta2,
		  InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& r,
		  const int jmax = 99,
		  const CompressionStrategy strategy = St04a);

  template <class PROBLEM>
  void NRESIDUAL(const PROBLEM& P,
		 const int patch,
		 const CDD1Parameters& params,
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

 #include <cdd1_local.cpp>

 #endif
