// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_STURM_EQUATION_H
#define _WAVELETTL_STURM_EQUATION_H

#include <set>
#include <utils/array1d.h>

namespace WaveletTL
{
  /*!
    This class models the (preconditioned) infinite-dimensional matrix problem
    
      Au = D^{-1}LD^{-1}u = D^{-1}F

    when reformulating a Sturm boundary value problem on [0,1]
    
      -(py')'(t) + q(t)y(t) = g(t), 0 <= t <= 1 
      
    with first (Dirichlet) or second (Neumann) order b.c.'s as modeled in
    the class simpleSturmBVP as an equivalent operator equation
    within \ell_2 by means of a wavelet basis.

    The corresponding bilinear form in

      L = (a(\psi_\nu,\psi_\lambda))_{\lambda,\nu}

    is

      a(u,v) = \int_0^1 [p(t)u'(t)v'(t)+q(t)u(t)v(t)] dt
      
    and the right-hand side is
 
      f(v) = \int_0^1 g(t)v(t) dt
   
    The evaluation of a(.,.) and f is possible for arguments \psi_\lambda
    which stem from a wavelet basis \Psi=\{\psi_\lambda\} of the corresponding
    function space over [0,1].
    To achieve independence from the concrete choice of \Psi, the wavelet basis
    class is given as a template parameter WBASIS and should provide a constructor of
    the form

       WBASIS::WBASIS(const bool bc_left, const bool bc_right)

    where the parameters bc_* indicate where to enforce homogeneous Dirichlet
    boundary conditions.
    Of course a natural concrete value for WBASIS is the template class DSBasis<d,dT>.

    References:
    [St04a] Stevenson:
            On the compressibility of operators in wavelet coordinates
  */
  template <class WBASIS>
  class SturmEquation
  {
  public:
    SturmEquation(const simpleSturmBVP& bvp);

    /*!
      make template argument accessible
    */
    typedef WBASIS WaveletBasis;

    /*!
      read access to the basis
    */
    const WBASIS& basis() const { return basis_; }
    
    /*!
      evaluate the diagonal preconditioner D
    */
    double D(const typename WBASIS::Index& lambda) const;

    /*!
      rescale a coefficient vector by an integer power of D, c |-> D^{n}c
    */
    void rescale(InfiniteVector<double, typename WBASIS::Index>& coeffs, const int n) const;

    /*!
      evaluate the (unpreconditioned) bilinear form a;
      you can specify the order p of the quadrature rule, i.e.,
      (piecewise) polynomials of maximal degree p will be integrated exactly.
      Internally, we use an m-point composite Gauss quadrature rule adapted
      to the singular supports of the spline wavelets involved,
      so that m = (p+1)/2;
    */
    double a(const typename WBASIS::Index& lambda,
	     const typename WBASIS::Index& nu,
	     const unsigned int p = 4) const;

    /*!
      evaluate the (unpreconditioned) bilinear form a;
      you can give a hint about the support intersection
      (see the add_column())
    */
    double a(const typename WBASIS::Index& lambda,
	     const typename WBASIS::Index& nu,
	     const typename WBASIS::Support& supp,
	     const unsigned int p = 4) const;

    /*!
      given an index set Lambda, setup the corresponding  preconditioned stiffness matrix
     */
    void setup_stiffness_matrix(const std::set<typename WBASIS::Index>& Lambda,
				SparseMatrix<double>& A_Lambda) const;

    /*!
      estimate the spectral norm ||A||
    */
    double norm_A() const { return normA; }

    /*!
      estimate the spectral norm ||A^{-1}||
    */
    double norm_Ainv() const { return normAinv; }

    /*!
      estimate compressibility exponent s^*
    */
    double s_star() const {
      return 1.0 + WBASIS::primal_vanishing_moments(); // [St04a], Th. 2.3 for n=1
    }

    /*!
      estimate the compression constants alpha_k in
        ||A-A_k|| <= alpha_k * 2^{-s*k}
    */
    double alphak(const unsigned int k) const {
      return 0.1; // first quick hack, estimate the constants!
    }

    /*!
      add (a constant multiple of) the lambda-th column of A to a working coefficient set w,
      applying the J-th truncation rule from [CDD1, Prop. 3.4]
      (for the APPLY routine)
    */
    void add_column(const double factor,
		    const typename WBASIS::Index& lambda,
		    const int J,
		    InfiniteVector<double, typename WBASIS::Index>& w) const;
    
    /*!
      evaluate the (unpreconditioned) right-hand side f
    */
    double f(const typename WBASIS::Index& lambda) const;

    /*!
      given an index set Lambda, setup the corresponding preconditioned right-hand side
    */
    void setup_righthand_side(const std::set<typename WBASIS::Index>& Lambda,
			      Vector<double>& F_Lambda) const;
    
    /*!
      approximate the wavelet coefficient set of the preconditioned right-hand side F
      within a prescribed \ell_2 error tolerance
    */
    void RHS(const double eta, InfiniteVector<double, typename WBASIS::Index>& coeffs) const;

  protected:
    const simpleSturmBVP& bvp_;
    WBASIS basis_;
    double normA, normAinv;

    // type of one block in a given column of A
    class MatrixBlock {
    public:
      Array1D<typename WBASIS::Index> indices;
      Array1D<double> entries;
    };

#define _WAVELETTL_STURM_EQUATION_CACHE
#ifdef _WAVELETTL_STURM_EQUATION_CACHE
    // type of one column of A ((j0-1)-th entry corresponds to the generators on level j0)
    typedef std::map<int, MatrixBlock> MatrixBlockCache;

    // type of the cached columns of A
    typedef std::map<typename WBASIS::Index,
		     MatrixBlockCache,
		     std::less<typename WBASIS::Index> > MatrixColumnCache;

    // entry cache for A (mutable to overcome constness of add_column())
    mutable MatrixColumnCache cache_;
#endif

    // helper function to compute subblocks
    void compute_matrix_block(const typename WBASIS::Index& lambda,
			      const int level,
			      MatrixBlock& block) const;
  };
}

#include <galerkin/sturm_equation.cpp>

#endif
