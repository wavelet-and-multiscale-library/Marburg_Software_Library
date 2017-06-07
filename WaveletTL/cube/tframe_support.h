// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner, Ulrich Friedrich                   |
// +--------------------------------------------------------------------+
#ifndef _WAVELETTL_TFRAME_SUPPORT_H
#define	_WAVELETTL_TFRAME_SUPPORT_H

#include <list>
#include <set>

namespace WaveletTL
{

    template <class IFRAME, unsigned int DIM> class TensorFrame;


    /*
     * For a given interval frame IFRAME, compute a cube
     * 2^{-j_}<a_,b_> = 2^{-j_1}[a_1,b_1]x...x2^{-j_n}[a_n,b_n]
     * which contains the support of a single primal cube quark
     * or quarklet psi_lambda.
     */

    template <class IFRAME, unsigned int DIM>
    void support(const TensorFrame<IFRAME,DIM>& frame,
                 const typename TensorFrame<IFRAME,DIM>::Index& lambda,
                 typename TensorFrame<IFRAME,DIM>::Support& supp);

    /*
     * For a given interval frame, compute a cube
     * 2^{-j_}<a_,b_> = 2^{-j_1}[a_1,b_1]x...x2^{-j_n}[a_n,b_n]
     * representing the intersection of the support cubes
     * corresponding to the indices lambda and mu.
     * Function returns true if a nontrivial intersection
     * exists, false otherwise. In the latter case 'supp'
     * has no meaningful value.
     */
    template <class IFRAME, unsigned int DIM>
    bool intersect_supports(const TensorFrame<IFRAME,DIM>& frame,
  			  const typename TensorFrame<IFRAME,DIM>::Index& lambda,
			  const typename TensorFrame<IFRAME,DIM>::Index& mu,
			  typename TensorFrame<IFRAME,DIM>::Support& supp);

  
    /*
     * For a given quarklet \psi_\lambda, compute all quarks OR quarklets
     * \psi_\nu with level |\nu|=j (multiindex), such that the respective supports
     * have a nontrivial intersection.
     * OUTPUT IS SORTED
     * 
     * output std::list<int>& intersecting is easy to implement and may be faster
     */
    template <class IFRAME, unsigned int DIM>
    void intersecting_quarklets(const TensorFrame<IFRAME,DIM>& frame,
                               const typename TensorFrame<IFRAME,DIM>::Index& lambda,
                               const MultiIndex<int,DIM> j, const bool generators,
                               std::list<typename TensorFrame<IFRAME,DIM>::Index>& intersecting,
                               const MultiIndex<int,DIM> p);

    /*
     * For a given quarklet \psi_\lambda, compute all quarks AND quarklets
     * \psi_\nu with level |\nu|=j (multiindex), such that the respective supports
     * have a nontrivial intersection.
     * 
     * Routine is slow. Use is disadvised. However, it propably produces correct sets, so it might be useful for debugging
     */
    /*
    template <class IFRAME, unsigned int DIM>
    void intersecting_elements(const TensorFrame<IFRAME,DIM>& frame,
                               const typename TensorFrame<IFRAME,DIM>::Index& lambda,
                               const MultiIndex<int,DIM> j,
                               std::list<typename TensorFrame<IFRAME,DIM>::Index>& intersecting);
*/
  /*!
   * Decide whether the support of a given (primal) quark/quarklet \psi_\lambda
   * intersects the singular support of another (primal) quark/quarklet \psi_\nu..
   *
   * Note: We have intersection of the singular supports if and only if:
   * (cube_support:)   one of the components has this property in one dimension
   * (tframe_support:) all of the components have this property
  */
  template <class IFRAME, unsigned int DIM>
  bool intersect_singular_support(const TensorFrame<IFRAME,DIM>& frame,
				  const typename TensorFrame<IFRAME,DIM>::Index& lambda,
				  const typename TensorFrame<IFRAME,DIM>::Index& nu);
}

#include <cube/tframe_support.cpp>


#endif	/* _WAVELETTL_TFRAME_SUPPORT_H */

