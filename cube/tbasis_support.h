// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner, Ulrich Friedrich                   |
// +--------------------------------------------------------------------+
#ifndef _TBASIS_SUPPORT_H
#define	_TBASIS_SUPPORT_H

#include <list>
#include <set>

namespace WaveletTL
{

    template <class IBASIS, unsigned int DIM> class TensorBasis;


    /*
     * For a given interval basis IBASIS, compute a cube
     * 2^{-j_}<a_,b_> = 2^{-j_1}[a_1,b_1]x...x2^{-j_n}[a_n,b_n]
     * which contains the support of a single primal cube generator
     * or wavelet psi_lambda.
     */

    template <class IBASIS, unsigned int DIM>
    void support(const TensorBasis<IBASIS,DIM>& basis,
                 const typename TensorBasis<IBASIS,DIM>::Index& lambda,
                 typename TensorBasis<IBASIS,DIM>::Support& supp);

    /*
     * For a given interval basis, compute a cube
     * 2^{-j_}<a_,b_> = 2^{-j_1}[a_1,b_1]x...x2^{-j_n}[a_n,b_n]
     * representing the intersection of the support cubes
     * corresponding to the indices lambda and mu.
     * Function returns true if a nontrivial intersection
     * exists, false otherwise. In the latter case 'supp'
     * has no meaningful value.
     */
    template <class IBASIS, unsigned int DIM>
    bool intersect_supports(const TensorBasis<IBASIS,DIM>& basis,
  			  const typename TensorBasis<IBASIS,DIM>::Index& lambda,
			  const typename TensorBasis<IBASIS,DIM>::Index& mu,
			  typename TensorBasis<IBASIS,DIM>::Support& supp);

  /*!
    For a given wavelet \psi_\lambda, compute all generators/wavelets
    \psi_\nu with level |\nu|=j (multiindex), such that the respective supports
    have a nontrivial intersection
  */
  template <class IBASIS, unsigned int DIM>
  void intersecting_wavelets(const TensorBasis<IBASIS,DIM>& basis,
			     const typename TensorBasis<IBASIS,DIM>::Index& lambda,
			     const MultiIndex<int,DIM> j, const bool generators,
			     std::list<typename TensorBasis<IBASIS,DIM>::Index>& intersecting);

  /*!
    Decide whether the support of a given (primal) generator/wavelet \psi_\lambda
    intersects the singular support of another (primal) generator/wavelet \psi_\nu.
  */
  template <class IBASIS, unsigned int DIM>
  bool intersect_singular_support(const TensorBasis<IBASIS,DIM>& basis,
				  const typename TensorBasis<IBASIS,DIM>::Index& lambda,
				  const typename TensorBasis<IBASIS,DIM>::Index& nu);
}

#include "tbasis_support.cpp"


#endif	/* _TBASIS_SUPPORT_H */

