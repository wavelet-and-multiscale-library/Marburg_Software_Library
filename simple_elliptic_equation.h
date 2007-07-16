// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of FrameTL - the Frame Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Manuel Werner                                                      |
// +--------------------------------------------------------------------+

#ifndef _FRAMETL_SIMPLE_ELLIPTIC_EQUATION_H
#define _FRAMETL_SIMPLE_ELLIPTIC_EQUATION_H

#include <aggregated_frame.h>
#include <numerics/bvp.h>
#include <adaptive/compression.h>
#include <interval/i_index.h>
#include <galerkin/infinite_preconditioner.h>
#include <frame_support.h>
#include <index1D.h>

using FrameTL::AggregatedFrame;
using MathTL::EllipticBVP;
using WaveletTL::CompressionStrategy;
using WaveletTL::IntervalIndex;
using WaveletTL::FullyDiagonalEnergyNormPreconditioner;

namespace FrameTL
{

  /*!
    This class models the (preconditioned) infinite-dimensional matrix problem
    
    Au = D^{-1}LD^{-1}u = D^{-1}F

    when reformulating a symmetric, second-order elliptic
    boundary value problem in divergence form over some domain
    Omega in R^d with boundary Gamma=dOmega,
    with homogeneous Dirichlet/Neumann/Robin boundary conditions

    -div(a(x)grad u(x)) + q(x)u(x) = f(x) in Omega
                             u(x) = 0 on Gamma_D
                         du/dn(x) = 0 on Gamma\Gamma_D.
    
    The corresponding bilinear form in

    L = (a(\psi_\nu,\psi_\lambda))_{\lambda,\nu}

    is

    a(u,v) = \int_Omega <a(x)*grad u(x), grad v(x)>  dx +
              \int_Omega q(x) * u(x) * v(x) dx
     
    and the right-hand side is
     
    f(v) = \int_Omega f(x)*v(x)  dx.

    The evaluation of a(.,.) and f is possible for arguments \psi_\lambda
    which stem from an aggregated wavelet frame \Psi=\{\psi_\lambda\} of the corresponding
    function space over Omega.

    WE ASSUME THAT THE COEFFICIENTS OF THE ELLIPTIC PDE ARE SEPERABLE AND SMOOTH AND THAT THE PATCHES
    OF THE UNDERLYING DOMAIN DECOMPOSITION ARE RECATANGULAR AND ALIGNED WITH THE CARTESIAN
    COORDINATES.
    FOR THIS SPECIAL CASE, a(.,.) CAN BE EXACTLY COMPUTED AT UNIT COST AND TENSOR PRODUCT
    STRUCTURE CAN BE EXPLOITED.
  */
  template <class IBASIS, unsigned int DIM>
  class SimpleEllipticEquation
  //: public FullyDiagonalDyadicPreconditioner<typename AggregatedFrame<IBASIS,DIM>::Index>
      : public FullyDiagonalEnergyNormPreconditioner<typename AggregatedFrame<IBASIS,DIM>::Index>
  {
  public:

    /*!
      constructor
     */
    SimpleEllipticEquation(const EllipticBVP<DIM>* ell_bvp,
			   const AggregatedFrame<IBASIS,DIM>* frame,
			   const int jmax);

    
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
      get the boundary value problem
    */
    const EllipticBVP<DIM>&  get_bvp() const { return *ell_bvp_; }


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
    double operator_order() const { return 1; }
    
    /*!
      evaluate the diagonal preconditioner D
    */
    double D(const typename AggregatedFrame<IBASIS,DIM>::Index& lambda) const;

    /*!
      rescale a coefficient vector by an integer power of D, c |-> D^{n}c
    */
    void rescale(InfiniteVector<double, typename AggregatedFrame<IBASIS,DIM>::Index>& coeffs,
		 const int n) const;

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
      returns spectral norm ||A^{-1}||
      estimate for ||A^{-1}|| has to be
      externally computed and to be set
      during initialization of the program.
      We assume this because ||A^{-1}|| quantity is hardly
      available in the frame case and
      quite complicated eigenvalue/eigenvector
      methods have to applied that are not implemented so
      far.
    */
    double norm_Ainv() const { return normAinv; };

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
      //return pow(2,(-k))*norm_A(); // suboptimal
      return 2.*norm_A(); // suboptimal
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
      approximate the wavelet coefficient set of the preconditioned right-hand side F restricted
      to patch p within a prescribed \ell_2 error tolerance
    */
    void RHS(const double eta, const int p,
	     InfiniteVector<double, 
	     typename AggregatedFrame<IBASIS,DIM>::Index>& coeffs) const;


    /*!
      compute (or estimate) ||F||_2
    */
    double F_norm() const { return sqrt(fnorm_sqr); }

    /*!
      set the boundary value problem
    */
    void set_bvp(const EllipticBVP<DIM>*);


    /*!
      w += factor * (stiffness matrix entries of in column lambda on level j)
    */
    void add_level (const Index& lambda,
		    InfiniteVector<double, Index>& w, const int j,
		    const double factor,
		    const int J,
		    const CompressionStrategy strategy) const;

   protected:
    
    /*!
      corresponding elliptic boundary value problem
     */
    const EllipticBVP<DIM>* ell_bvp_;

    /*!
      underlying frame
     */
    const AggregatedFrame<IBASIS,DIM>* frame_;

    //#################### Caching ##################
    typedef std::map<Index1D<IBASIS>,double > Column1D;
    typedef std::map<Index1D<IBASIS>,Column1D> One_D_IntegralCache;
    
    /*!
      cache for one dimensional integrals
      ONLY USED TOGETHER WITH TrivialAffine QUADRATURE RULE OPTION
     */
    mutable One_D_IntegralCache one_d_integrals;
    //###############################################

  private:

    /*!
     */
    double integrate(const Index1D<IBASIS>& lambda,
		     const Index1D<IBASIS>& mu,
		     const int N_Gauss,
		     const int dir,
		     const typename CubeBasis<IBASIS,DIM>::Support* supp_lambda,
		     const typename CubeBasis<IBASIS,DIM>::Support* supp_mu) const;

    // precompute the right-hand side
    void compute_rhs();

    // precompute diagonal of stiffness matrix
    void compute_diagonal();


    const int jmax_;

    // right-hand side coefficients on a fine level, sorted by modulus
    Array1D<std::pair<typename AggregatedFrame<IBASIS,DIM>::Index,double> > fcoeffs;

    // patchwise right-hand side coefficients on a fine level, sorted by modulus
    Array1D<Array1D<std::pair<typename AggregatedFrame<IBASIS,DIM>::Index,double> > > fcoeffs_patch;


    //     // right-hand side coefficients on a fine level, sorted by modulus
    //     InfiniteVector<double,typename AggregatedFrame<IBASIS,DIM>::Index> stiff_diagonal;
    
    
    // square root of coefficients on diagonal of stiffness matrix
    Array1D<double> stiff_diagonal;


    // (squared) \ell_2 norm of the precomputed right-hand side
    double fnorm_sqr;
    Array1D<double> fnorms_sqr_patch;

    // reminder: The keyword mutable can only be applied to non-static
    // and non-const data members of a class. If a data member is declared mutable,
    // then it is legal to assign a value to this data member from
    // a const member function.
    mutable double normA, normAinv;
  };
}

#include <simple_elliptic_equation.cpp>

#endif
