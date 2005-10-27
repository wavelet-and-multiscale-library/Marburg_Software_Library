// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_CACHED_PROBLEM_H
#define _WAVELETTL_CACHED_PROBLEM_H

#include <map>

namespace WaveletTL
{
  /*!
    This class provides a cache layer for generic (preconditioned)
    infinite-dimensional matrix problems of the form
    
      Au = D^{-1}LD^{-1}u = D^{-1}F.

    The operator equation is not assumed to be induced by a local operator,
    i.e., the cache class should also work in the case of integral operators.
    All evaluations of the bilinear form a(.,.) are cached.
    Internally, the cache is managed as follows. The nonzero values of the bilinear
    form a(.,.) are stored in a map.

    The template class CachedProblem implements the minimal signature to be
    used within the APPLY routine.
  */
  template <class PROBLEM>
    class CachedProblem: public PROBLEM
    {
    public:
      /*!
	default constructor
      */
      CachedProblem(const PROBLEM& P);
      
      /*!
	make wavelet basis type accessible
      */
      typedef typename PROBLEM::WaveletBasis WaveletBasis;

      /*!
	evaluate the (unpreconditioned) bilinear form a
	(cached)
      */
      double a(const typename WaveletBasis::Index& lambda,
	       const typename WaveletBasis::Index& nu) const;

      /*!
	estimate the spectral norm ||A||
      */
      double norm_A() const;
      
      /*!
	estimate the spectral norm ||A^{-1}||
      */
      double norm_Ainv() const;
      
      /*!
	estimate the compression constants alpha_k in
        ||A-A_k|| <= alpha_k * 2^{-s*k}
      */
      double alphak(const unsigned int k) const {
	return 2*norm_A(); // suboptimal
      }

    protected:
      // type of one column in the entry cache of A
      typedef std::map<typename WaveletBasis::Index, double> Column;
      
      // type of the entry cache of A
      typedef std::map<typename WaveletBasis::Index, Column> ColumnCache;
      
      // entries cache for A (mutable to overcome the constness of add_column())
      mutable ColumnCache entries_cache;

      // estimates for ||A|| and ||A^{-1}||
      mutable double normA, normAinv;
    };
}

#include <galerkin/cached_problem.cpp>

#endif
