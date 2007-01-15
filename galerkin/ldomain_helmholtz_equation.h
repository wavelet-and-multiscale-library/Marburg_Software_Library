// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2007                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_LDOMAIN_HELMHOLTZ_EQUATION_H
#define _WAVELETTL_LDOMAIN_HELMHOLTZ_EQUATION_H

#include <utils/function.h>
#include <utils/array1d.h>
#include <algebra/infinite_vector.h>
#include <interval/spline_basis.h>
#include <galerkin/galerkin_utils.h>
#include <galerkin/infinite_preconditioner.h>
#include <galerkin/ldomain_gramian.h>
#include <galerkin/ldomain_laplacian.h>
#include <galerkin/cached_problem.h>
#include <adaptive/compression.h>
#include <Ldomain/ldomain_basis.h>

using namespace MathTL;

namespace WaveletTL
{
  /*!
    This class models the (preconditioned) infinite-dimensional matrix problem
    
      Au = D^{-1}LD^{-1}u = D^{-1}F

    when reformulating a Helmholtz equation on the L--shaped domain
    
      -Delta u(x)) + alpha*u(x) = f(x)

    with first order b.c.'s as an equivalent operator equation
    within \ell_2 by means of a wavelet basis.

    The corresponding bilinear form in

      L = (a(\psi_\nu,\psi_\lambda))_{\lambda,\nu}

    is

      a(u,v) = \int_Omega [a(x)grad u(x)grad v(x)+q(x)u(x)v(x)] dx
      
    and the right-hand side is
 
      F(v) = \int_0^1 f(x)v(x) dx
   
    The evaluation of a(.,.) and f is possible for arguments \psi_\lambda
    which stem from a wavelet basis \Psi=\{\psi_\lambda\} of the corresponding
    function space over Omega.

    The wavelet bases useable for the discretization are those
    modeled by SplineBasis<d,dT,DS_construction>.

    Internally, HelmholtzEquation1D holds two distinct caches for the
    Gramian (identity) part of the stiffness matrix and the Laplacian part.
  */
  template <int d, int dT>
  class LDomainHelmholtzEquation
    : public FullyDiagonalEnergyNormPreconditioner<typename LDomainBasis<SplineBasis<d,dT,DS_construction> >::Index>
  {
  public:
    /*!
      make template argument accessible
    */
    typedef LDomainBasis<SplineBasis<d,dT,DS_construction> > WaveletBasis;

    /*!
      wavelet index class
    */
    typedef typename WaveletBasis::Index Index;

    /*!
      index type of vectors and matrices
    */
    typedef typename Vector<double>::size_type size_type;

    /*!
      constructor from a given wavelet basis,
      two filenames for the precomputed Gramian and Laplacian (plus corresponding jmax)
      and a right-hand side y
    */
    LDomainHelmholtzEquation(const WaveletBasis& basis,
			     const char* G_file,
			     const char* A_file,
			     const int jmax,
			     const double alpha,
			     const InfiniteVector<double, typename WaveletBasis::Index>& y);

    /*!
      set right-hand side y
    */
    void set_rhs(const InfiniteVector<double,Index>& y) const;

    /*!
      set reaction coefficient alpha
    */
    void set_alpha(const double alpha) const;

    /*!
      read access to the basis
    */
    const WaveletBasis& basis() const { return basis_; }
    
    /*!
      space dimension of the problem
    */
    static const int space_dimension = 2;

    /*!
      differential operators are local
    */
    static bool local_operator() { return true; }

    /*!
      (half) order t of the operator
    */
    double operator_order() const { return 1.; }
    
    /*!
      evaluate the diagonal preconditioner D
    */
    double D(const Index& lambda) const;

    /*!
      evaluate the (unpreconditioned) bilinear form a
      (inherited from EnergyNormPreconditioner)
    */
    double a(const Index& lambda, const Index& nu) const;

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
    */
    double s_star() const;
    
    /*!
      estimate the compression constants alpha_k in
        ||A-A_k|| <= alpha_k * 2^{-s*k}
    */
    double alphak(const unsigned int k) const {
      return 2*norm_A(); // suboptimal
    }

    /*!
      evaluate the (unpreconditioned) right-hand side f
    */
    double f(const Index& lambda) const {
      return y_.get_coefficient(lambda);
    }

    /*!
      approximate the wavelet coefficient set of the preconditioned right-hand side F
      within a prescribed \ell_2 error tolerance
    */
    void RHS(const double eta, InfiniteVector<double,Index>& coeffs) const
    {
      coeffs = y_precond; // dirty
    }

    /*!
      compute (or estimate) ||F||_2
    */
    double F_norm() const { return l2_norm(y_); }

    /*!
      w += factor * (stiffness matrix entries in column lambda on level j)
    */
    void add_level (const Index& lambda,
		    InfiniteVector<double, Index>& w, const int j,
		    const double factor,
		    const int J,
		    const CompressionStrategy strategy = St04a) const;

  protected:
    WaveletBasis basis_;
    mutable double alpha_;
    mutable InfiniteVector<double, typename WaveletBasis::Index> y_;
    mutable InfiniteVector<double, typename WaveletBasis::Index> y_precond;

  public:
    LDomainGramian<SplineBasis<d,dT,DS_construction> > G_;
    CachedProblemFromFile<LDomainGramian<SplineBasis<d,dT,DS_construction> > > GC_;
    LDomainLaplacian<SplineBasis<d,dT,DS_construction> > A_;
    CachedProblemFromFile<LDomainLaplacian<SplineBasis<d,dT,DS_construction> > > AC_;

    // estimates for ||A|| and ||A^{-1}||
    mutable double normA, normAinv;
  };
}

#include <galerkin/ldomain_helmholtz_equation.cpp>

#endif
