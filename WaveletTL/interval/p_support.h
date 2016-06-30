// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_P_SUPPORT_H
#define _WAVELETTL_P_SUPPORT_H

#include <list>
#include <set>


namespace WaveletTL
{
  template <int d, int dT> class PBasis;

  /*!
    Compute an interval 2^{-j}[k1,k2] which contains the support of a
    single primal [P] generator or wavelet \psi_\lambda.
    (j == lambda.j()+lambda.e() is neglected for performance reasons)
  */
  template <int d, int dT>
  void support(const PBasis<d,dT>& basis,
	       const typename PBasis<d,dT>::Index& lambda,
	       int& k1, int& k2);
  template <int d, int dT>
  inline
  void support(const PBasis<d,dT>& basis,
               const int j_, const int e_, const int k_,
               int& k1, int& k2);
  
  /*!
    Decide whether the support of a generator/wavelet \psi_\lambda
    has a nontrivial intersection with a given dyadic interval 2^{-m}[a,b]
    and compute it in the form 2^{-j}[k1,k2]. If the return value is false,
    k1 and k2 will have no meaningful value for performance reasons.
  */
  template <int d, int dT>
  bool intersect_supports(const PBasis<d,dT>& basis,
			  const typename PBasis<d,dT>::Index& lambda,
			  const int m, const int a, const int b,
			  int& j, int& k1, int& k2);
  
  /*!
    Decide whether the supports of two generators/wavelets \psi_\lambda and
    \psi_\nu have an intersection of positive measure and compute it
    in the form 2^{-j}[k1,k2]. If the return value is false, the computed
    support intersection will have no meaningful values, for performance reasons.
  */
  template <int d, int dT>
  bool intersect_supports(const PBasis<d,dT>& basis,
			  const typename PBasis<d,dT>::Index& lambda,
			  const typename PBasis<d,dT>::Index& nu,
			  typename PBasis<d,dT>::Support& supp);



  /*!
    For a given wavelet \psi_\lambda, compute all generators/wavelets
    \psi_\nu with level |\nu|=j, such that the respective supports
    have a nontrivial intersection.
  */
  template <int d, int dT>
  void intersecting_wavelets(const PBasis<d,dT>& basis,
			     const typename PBasis<d,dT>::Index& lambda,
			     const int j, const bool generators,
			     std::list<typename PBasis<d,dT>::Index>& intersecting);
  
  /*!
    For a given wavelet \psi_\lambda, compute all generators/wavelets
    \psi_\nu with level |\nu|=j, such that the respective supports
    have a nontrivial intersection;
    return those intersections

PERFORMANCE:
   * computation of the intersections is inefficient! 
   * It needs to be replaced by a direct approach!
   * For known k the support of a generator/ wavelet is immideately known. 
   * As a result for a given range [firstk, lastk] of intersecting wavelets the respective support overlapps with supp(lambda) are directly accessible.
  */
    template <int d, int dT>
    void intersecting_wavelets(const PBasis<d,dT>& basis,
    			       const typename PBasis<d,dT>::Index& lambda,
			       const int j, const bool generators,
			       std::list<std::pair<typename PBasis<d,dT>::Index, typename PBasis<d,dT>::Support> >& intersecting);

  /*!
    Decide whether the support of a given (primal) generator/wavelet \psi_\lambda
    intersects the singular support of a union 2^{-m}[a,b] of dyadic unit intervals.
    If this is the case, return true and the intersection of \supp\psi_\lambda
    and \supp\psi_\nu in the form 2^{-j}[k1,k2].
    Otherwise, return false (in this case, j, k1 and k2 will have
    no meaningful values, for performance reasons).
    The routine should only be called in the case lambda.j()+lambda.e() >= m.
  */
  template <int d, int dT>
  bool intersect_singular_support(const PBasis<d,dT>& basis,
				  const typename PBasis<d,dT>::Index& lambda,
				  const int m, const int a, const int b,
				  int& j, int& k1, int& k2);

  /*!
    Decide whether the support of a given (primal) generator/wavelet \psi_\lambda
    intersects the singular support of another (primal) generator/wavelet \psi_\nu.
    If this is the case, return true and the intersection of \supp\psi_\lambda
    and \supp\psi_\nu in the form 2^{-j}[k1,k2].
    Otherwise, return false (in this case, j, k1 and k2 will have
    no meaningful values, for performance reasons).
  */
  template <int d, int dT>
  bool intersect_singular_support(const PBasis<d,dT>& basis,
				  const typename PBasis<d,dT>::Index& lambda,
				  const typename PBasis<d,dT>::Index& nu,
				  int& j, int& k1, int& k2);
  /*
   * This version seems useless. ToDo: check usage.
   */
  template <int d, int dT>
  bool intersect_singular_support(const PBasis<d,dT>& basis,
				  const typename PBasis<d,dT>::Index& lambda,
				  const typename PBasis<d,dT>::Index& nu);
  /*
   * new version: 
   * speedup 1.3 compared to the old code for the first method
   * speedup 1.1 compared to the new code for the first method (!!)
   * observe the different meaning of (muj,mue,muk) and (m,a,b) (in the above version of this method
   */
  template <int d, int dT>
  bool intersect_singular_support(const PBasis<d,dT>& basis,
				  const int lamj, const int lame, const int lamk,
                                  const int muj, const int mue, const int muk,
				  int& j, int& k1, int& k2);

  /*
   * This method gives the minimal and maximal shiftparameter k for such that the wavelets or generators on level j 
   * intersect with the wavelet or generator lambda.
   * 
   * Moved from interval/support.h
  */
  template <int d, int dT>
  void get_intersecting_wavelets_on_level(const PBasis<d,dT>& basis,
                                          const typename PBasis<d,dT>::Index& lambda,
                                          const int j, 
                                          const bool generators,
                                          int& mink, 
                                          int& maxk);
  template <int d, int dT>
  void get_intersecting_wavelets_on_level(const PBasis<d,dT>& basis,
                                          const int lamj, const int lame, const int lamk,
                                          const int j, 
                                          const bool generators,
                                          int& mink, 
                                          int& maxk);
  /*
   * In this variant k1_lambda, k2_lambda correspond to the output of
   * support(lambda_basis, lamj, lame, lamk, k1_lambda, k2_lambda)
   * 
   * this is useful if lambda_basis fulfills other boundary conditions than basis
   */
// PERFORMANCE :: this method only uses (lam_j+lam_e)! However, if you replace the 2 arguments lam_j and lam_e with only one, you end up with the same signature as the lamj,lame,lamk version above :(
  template <int d, int dT>
  void get_intersecting_wavelets_on_level(const PBasis<d,dT>& basis,
                                            const int lamj, const int lame, const int k1_lambda, const int k2_lambda,
                                            const int j, 
                                            const bool generators,
                                            int& mink, 
                                            int& maxk);
  
}

#include <interval/p_support.cpp>

#endif
