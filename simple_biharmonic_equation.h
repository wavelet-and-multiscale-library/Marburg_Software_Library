// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of FrameTL - the Frame Template Library          |
// |                                                                    |
// | Copyright (c) 2002-2007                                            |
// | Manuel Werner, Andreas Schneider                                   |
// +--------------------------------------------------------------------+

#ifndef _FRAMETL_SIMPLE_BIHARMONIC_EQUATION_H
#define _FRAMETL_SIMPLE_BIHARMONIC_EQUATION_H

#include <aggregated_frame.h>
#include <functional.h>
#include <adaptive/compression.h>
#include <galerkin/infinite_preconditioner.h>
#include <index1D.h>
#include <frame_support.h>

using FrameTL::AggregatedFrame;
using WaveletTL::CompressionStrategy;
using WaveletTL::FullyDiagonalEnergyNormPreconditioner;

namespace FrameTL
{

 
  /*!
    quadrature strategies for computation of
    stiffness matrix entries, use only TrivialAffine!
  */
  enum QuadratureStrategy
    {
      Composite,
      TrivialAffine,
      SplineInterpolation
    };

  /*!
    This class models the (preconditioned) infinite-dimensional matrix problem
    
    Au = D^{-1}LD^{-1}u = D^{-1}F

    when reformulating a symmetric, forth-order elliptic
    boundary value problem in divergence form over some domain
    Omega in R^d with boundary Gamma=dOmega,
    with homogeneous Dirichlet/Neumann/Robin boundary conditions

    -Delta^2 u(x) = f(x) in Omega
                             u(x) = 0 on Gamma_D
                         du/dn(x) = 0 on Gamma\Gamma_D.
    
    The corresponding bilinear form in

    L = (a(\psi_\nu,\psi_\lambda))_{\lambda,\nu}

    is

    a(u,v) = \int_Omega <Delta u(x), Delta v(x)>  dx 
     
    and the right-hand side is a functional
     
    v -> f(v).

    The evaluation of a(.,.) and f is possible for arguments \psi_\lambda
    which stem from an aggregated wavelet frame \Psi=\{\psi_\lambda\} of the corresponding
    function space over Omega.     
  */
  template <class IBASIS, unsigned int DIM>
  class SimpleBiharmonicEquation
  //  : public FullyDiagonalDyadicPreconditioner<typename AggregatedFrame<IBASIS,DIM>::Index>
    : public FullyDiagonalEnergyNormPreconditioner<typename AggregatedFrame<IBASIS,DIM>::Index>
  {
  public:	

    /*!
      constructor
     */
    SimpleBiharmonicEquation(const Functional<IBASIS,DIM>* rhs,
		       const AggregatedFrame<IBASIS,DIM>* frame,
		       const int jmax,
		       const QuadratureStrategy qsrtat = TrivialAffine,
		       const bool precompute_rhs = true);

    /*!
      make template argument accessible
    */
    typedef AggregatedFrame<IBASIS,DIM> Frame;

    /*!
      dummy typedef to be compatible with WaveletTL
      routines
     */
    typedef AggregatedFrame<IBASIS,DIM> WaveletBasis;
    
    /*!
      make template argument accessible
    */
    typedef typename Frame::Index Index;

    /*!
      read access to the frame
    */
    const AggregatedFrame<IBASIS,DIM>& frame() const { return *frame_; }

    /*!
      get the righthand side
    */
    const Functional<IBASIS,DIM>&  get_rhs() const { return *rhs_; }


    /*!
      read access to the frame but with a somewhat weird
      function name.
      this is just a first hack to be able to use
      routines in WaveletTL's compression.h
    */
    const AggregatedFrame<IBASIS,DIM>& basis() const { return *frame_; }  

    /*!
      space dimension of the problem
    */
    static const int space_dimension = DIM;

    /*!
      differential operators are local
    */
    static bool local_operator() { return true; }

    /*!
      order of the operator
    */
    static double operator_order() { return 2; }
    
    /*!
      evaluate the diagonal preconditioner D
    */
    double D(const typename AggregatedFrame<IBASIS,DIM>::Index& lambda) const;

    /*!
      evaluate the (unpreconditioned) bilinear form a;
    */
    double a(const typename AggregatedFrame<IBASIS,DIM>::Index& lambda,
	     const typename AggregatedFrame<IBASIS,DIM>::Index& nu) const;

    /*!
      estimate the spectral norm ||A||
    */
    double norm_A() const;

    /*!
      sets estimate for ||A||
    */
    void set_norm_A(const double _normA) { normA = _normA; }

    /*!
      sets estimate for ||A^{-1}||
    */
    void set_Ainv(const double nAinv) { normAinv = nAinv; };

    /*!
      estimate compressibility exponent s^*
    */   
     double s_star() const;

    /*!
      estimate the compression constants alpha_k in
      ||A-A_k|| <= alpha_k * 2^{-s*k}
    */
    double alphak(const unsigned int k) const {
      cout << "works" << endl;
      return 2*norm_A(); // suboptimal
    }
   
    /*!
      evaluate the (unpreconditioned) right-hand side f
    */
    double f(const typename AggregatedFrame<IBASIS,DIM>::Index& lambda) const;

    /*!
      approximate the wavelet coefficient set of the preconditioned right-hand side F
      within a prescribed \ell_2 error tolerance
    */
    void RHS(const double eta, InfiniteVector<double, 
	     typename AggregatedFrame<IBASIS,DIM>::Index>& coeffs) const;

   /*!
      compute (or estimate) ||F||_2
    */
    double F_norm() const { return sqrt(fnorm_sqr); }


   protected:
    
    //! corresponding righthand side
    const Functional<IBASIS,DIM>* rhs_;
    //! corresponding frame
    const AggregatedFrame<IBASIS,DIM>* frame_;
    //! maximal level to be used
    const int jmax_;



    //####################
    typedef std::map<Index1D<IBASIS>,double > Column1D;
    typedef std::map<Index1D<IBASIS>,Column1D> One_D_IntegralCache;
    
    /*!
      cache for one dimensional integrals
      ONLY USED TOGETHER WITH TrivialAffine QUADRATURE RULE OPTION
     */
    mutable One_D_IntegralCache one_d_integrals;
    
    

    //####################

  private:

     /*!
     */
    double integrate(const Index1D<IBASIS>& lambda,
		     const Index1D<IBASIS>& mu,
		     const FixedArray1D<Array1D<double>,DIM >& irregular_grid,
		     const int N_Gauss, const int dir) const;

   
    // precompute the right-hand side
    void compute_rhs(bool b);

    // precompute diagonal of stiffness matrix
    void compute_diagonal();


    // right-hand side coefficients on a fine level, sorted by modulus
    Array1D<std::pair<typename AggregatedFrame<IBASIS,DIM>::Index,double> > fcoeffs;

    // right-hand side coefficients on a fine level, sorted by modulus
    // square root of coefficients on diagonal of stiffness matrix
    Array1D<double> stiff_diagonal;

    // (squared) \ell_2 norm of the precomputed right-hand side
    double fnorm_sqr;

    mutable double normA, normAinv;

    QuadratureStrategy qstrat_; 

  };
}

#include <simple_biharmonic_equation.cpp>

#endif // _FRAMETL_SIMPLE_BIHARMONIC_EQUATION_H
