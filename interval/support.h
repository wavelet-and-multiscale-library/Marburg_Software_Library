// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner, Gregor Kriwet                      |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_SUPPORT_H
#define _WAVELETTL_SUPPORT_H

#include <list>
#include <set>
#include <Rd/cdf_utils.h>

namespace WaveletTL
/*
  This class is needed for the new computation of the intersecting wavelets in cube_support, frame_support and tbasis_support.
  The method is testet for the P Basis and the DS Basis.
*/

{

  /*
     This method gives the minimal and maximal shiftparameter k for such that the wavelets or generators on level j 
     intersects with the wavelet or generator lambda.
  */
  template <class IBASIS>
  void get_intersecting_wavelets_on_level(const IBASIS& basis,
			     const typename IBASIS::Index& lambda,
			     const int j, const bool generators,
			     int& mink, int& maxk)
  { typedef typename IBASIS::Index Index;
    typedef typename IBASIS::Support Support;

    // compute support of \psi_\lambda
    const int j_lambda = lambda.j() + lambda.e();
    int k1_lambda, k2_lambda;
    support(basis, lambda, k1_lambda, k2_lambda);
    int d = basis.get_primalorder();
    int dT = basis.get_dualorder();
    

    // new code
    if (generators) {
      // the leftmost generator on the level j, s.th. its support intersects [a,b], fulfills
      //   2^{-j}(k+ell2) > a  but  2^{-j}(k-1+ell2) <= a,
      // so that ...
      mink = std::max(basis.DeltaLmin(), (int) floor(ldexp(1.0,j-j_lambda)*k1_lambda-(d-d/2))+1);
      
      // the rightmost generator on the level j, s.th. its support intersects [a,b], fulfills
      //   2^{-j}(k+ell1) < b  but  2^{-j}(k+1+ell1) >= b,
      // so that ...
      maxk = std::min(basis.DeltaRmax(j), (int) ceil(ldexp(1.0,j-j_lambda)*k2_lambda-(-d/2))-1);

    } else {
      // if the left boundary wavelets are not involved, then
      // the rightmost wavelet on the level j, s.th. its support intersects [a,b], fulfills
      //   2^{-j}(k+(d+dT)/2) > a  but  2^{-j}(k-1+(d+dT)/2) <= a,
      // so that ...
      mink = (ldexp(1.0,j-j_lambda)*k1_lambda < d+dT-1 // overestimate, TODO!
			  ? 0
			  : (int) floor(ldexp(1.0,j-j_lambda)*k1_lambda-(d+dT)/2)+1);

      // if the right boundary wavelets are not involved, then
      // the leftmost wavelet on the level j, s.th. its support intersects [a,b], fulfills
      //   2^{-j}(k-(d+dT)/2+1) < b  but  2^{-j}(k-(d+dT)/2+2) >= a,
      // so that ...
      maxk = (ldexp(1.0,-j_lambda)*k2_lambda > 1 - ldexp(1.0,-j)*(d+dT-1) // overestimate, TODO!
			 ? (1<<j)-1
			 : (int) ceil(ldexp(1.0,j-j_lambda)*k2_lambda+(d+dT)/2)-2);
    }
   }
}

#endif
