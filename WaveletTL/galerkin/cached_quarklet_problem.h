// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner, Philipp Keding                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_CACHED_QUARKLET_PROBLEM_H
#define _WAVELETTL_CACHED_QUARKLET_PROBLEM_H

#include <map>
#include <algebra/infinite_vector.h>
#include <adaptive/compression.h>
#include <galerkin/infinite_preconditioner.h>

using MathTL::InfiniteVector;

namespace WaveletTL
{
  /*!
    This class provides a cache layer for generic (preconditioned, cf. precond.h)
    infinite-dimensional matrix problems of the form
    
      Au = P^{-1}LQ^{-1}u = P^{-1}F.

    The operator equation is not assumed to be induced by a local operator,
    i.e., the cache class should also work in the case of integral operators.
    All evaluations of the bilinear form a(.,.) are cached.
    Internally, the cache is managed as follows. The nonzero values of the bilinear
    form a(.,.) are stored in a map.

    The template class CachedQuarkletProblem implements the minimal signature to be
    used within the APPLY routine.
    The cache organization and the addlevel function is adapted to the DKOR strategy is different to CachedProblem class.
    
    
  */
  template <class PROBLEM>
  class CachedQuarkletProblem
   :public FullyDiagonalQuarkletPreconditioner<typename PROBLEM::Index>
  {
  public:
    /*!
      constructor from an uncached problem,
      you can specify the estimates for ||A|| and ||A^{-1}||
      (if zero, CachedProblem will compute the estimates)
    */
    CachedQuarkletProblem(const PROBLEM* P,
		  const double normA = 0.0,
		  const double normAinv = 0.0);
    
    /*!
      make wavelet basis type accessible
    */
    typedef typename PROBLEM::WaveletBasis WaveletBasis;
    
    /*!
      wavelet index class
    */
    typedef typename PROBLEM::Index Index;
    
    /*!
      read access to the basis
    */
    const WaveletBasis& basis() const { return problem->basis(); }
    
    /*!
      space dimension of the problem
    */
    static const int space_dimension = PROBLEM::space_dimension;
    
    /*!
      locality of the operator
    */
    static bool local_operator() { return PROBLEM::local_operator(); }
    
    /*!
      (half) order t of the operator
    */
    double operator_order() const { return problem->operator_order(); }
    
    /*!
      evaluate the diagonal preconditioner D
    */
    double D(const Index& lambda) const {
      return problem->D(lambda);
    }

    /*
     * access to the underlying problem
     */
    const PROBLEM* get_problem() const {return problem;}

    /*!
      evaluate the (unpreconditioned) bilinear form a
      (cached)
    */
    double a(const Index& lambda,
	     const Index& nu) const;
    
    /*!
      estimate the spectral norm ||A||
    */
    double norm_A() const;
    
    /*!
      estimate the spectral norm ||A^{-1}||
    */
    double norm_Ainv() const;
    
    /*!
      estimate compressibility exponent s^*
      (we assume that the coefficients a(x),q(x) are smooth)
    */
    double s_star() const {
      return problem->s_star();
    }
    
    /*!
      estimate the compression constants alpha_k in
      ||A-A_k|| <= alpha_k * 2^{-s*k}
    */
    double alphak(const unsigned int k) const {
      return 2*norm_A(); // pessimistic
    }
    
    /*!
      evaluate the (unpreconditioned) right-hand side f
    */
    double f(const Index& lambda) const {
      return problem->f(lambda);
    }
    
    /*!
      approximate the wavelet coefficient set of the preconditioned right-hand side F
      within a prescribed \ell_2 error tolerance
    */
    void RHS(const double eta,
	     InfiniteVector<double, Index>& coeffs) const {
      problem->RHS(eta, coeffs);
    }

    
    /*!
      compute (or estimate) ||F||_2
    */
    double F_norm() const { return problem->F_norm(); }
    
    /*!
      w += factor * (stiffness matrix entries in column lambda on level j, p)
    */
    void add_level (const Index& lambda,
		    //InfiniteVector<double, Index>& w,
		    Vector<double>& w,
                    const int p, 
		    const int j,
		    const double factor,
		    const int J,
		    const CompressionStrategy strategy = DKR,
                    const int jmax = 99,
                    const int pmax = 0,
                    const double a = 0,
                    const double b = 0) const;
    
    int number (const Index& lambda, const int jmax) const;
    
  protected:
    //! the underlying (uncached) problem
    const PROBLEM* problem;
    
    // type of one subblock in one column of stiffness matrix  A
    typedef std::map<int, double> Subblock;
   
    // type of one block in one column of stiffness matrix  A
    // the key codes the level, that data are the entries
    typedef std::map<int, Subblock> Block;
    
    // type of one column in the entry cache of A
    // the key codes the polynomial, that data are the Subblocks to different levels
    //typedef std::map<int[2], Block> Column;
    typedef std::map<int, Block> Column;
    
    // type of the entry cache of A
    //typedef std::map<Index, Column> ColumnCache;
    typedef std::map<int, Column> ColumnCache;

    // entries cache for A (mutable to overcome the constness of add_column())
    mutable ColumnCache entries_cache;
    
    // estimates for ||A|| and ||A^{-1}||
    mutable double normA, normAinv;
  };  
}

#include <galerkin/cached_quarklet_problem.cpp>

#endif
