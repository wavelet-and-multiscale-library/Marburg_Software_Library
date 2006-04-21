// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Manuel Werner                                                      |
// +--------------------------------------------------------------------+

#ifndef _FRAMETL_ELLIPTIC_EQUATION_H
#define _FRAMETL_ELLIPTIC_EQUATION_H

#include <aggregated_frame.h>
#include <numerics/bvp.h>
#include <adaptive/compression.h>
#include <interval/i_index.h>

#include <galerkin/infinite_preconditioner.h>

using FrameTL::AggregatedFrame;
using MathTL::EllipticBVP;
using WaveletTL::CompressionStrategy;
using WaveletTL::IntervalIndex;
using WaveletTL::FullyDiagonalEnergyNormPreconditioner;

namespace FrameTL
{

  template <class IBASIS>
  class Index1D
  {

  public:
    Index1D (const IntervalIndex<IBASIS>& ind,
	     const unsigned int p, const unsigned int dir,
	     const unsigned int der);

    bool operator < (const Index1D<IBASIS>& lambda) const;
    bool operator == (const Index1D<IBASIS>& lambda) const;
    bool operator != (const Index1D<IBASIS>& lambda) const;
    bool operator <= (const Index1D<IBASIS>& lambda) const;

    IntervalIndex<IBASIS> index() const { return ind_; };
    unsigned int derivative() const { return der_; };
    unsigned int p() const { return p_; };
    unsigned int direction() const { return dir_; };

  protected:
    IntervalIndex<IBASIS> ind_;
    unsigned int p_;
    unsigned int dir_;
    unsigned int der_;

  };

  /*!
    quadrature strategies for computation of
    stiffness matrix entries
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
  */
  template <class IBASIS, unsigned int DIM>
  class EllipticEquation
  //  : public FullyDiagonalDyadicPreconditioner<typename AggregatedFrame<IBASIS,DIM>::Index>
    : public FullyDiagonalEnergyNormPreconditioner<typename AggregatedFrame<IBASIS,DIM>::Index>
  {
  public:

//     /*!
//       constructor from a boundary value problem and specified b.c.'s
//     */
//     EllipticEquation(const EllipticBVP<DIM>* bvp,
// 		     const FixedArray1D<bool,2*DIM>& bc);

    /*!
      constructor
     */
    EllipticEquation(const EllipticBVP<DIM>* ell_bvp,
		     const AggregatedFrame<IBASIS,DIM>* frame,
		     const QuadratureStrategy qsrtat = Composite);

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
    static int operator_order() { return 1; }
    
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
    */
    double a_quad(const typename AggregatedFrame<IBASIS,DIM>::Index& lambda,
		  const typename AggregatedFrame<IBASIS,DIM>::Index& nu,
		  const unsigned int p, const unsigned int N) const;

    /*!
      evaluate the (unpreconditioned) bilinear form a;
      you can specify the order p of the quadrature rule, i.e.,
      (piecewise) polynomials of maximal degree p will be integrated exactly.
      Internally, we use an m-point composite Gauss quadrature rule adapted
      to the singular supports of the spline wavelets involved,
      so that m = (p+1)/2;
    */
//     double a(const typename AggregatedFrame<IBASIS,DIM>::Index& lambda,
// 	     const typename AggregatedFrame<IBASIS,DIM>::Index& nu,
// 	     const unsigned int p = 2) const;

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
      return pow(2,(-k))*norm_A(); // suboptimal
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
		    const CompressionStrategy strategy = St04a) const;

   protected:
    
    /*!
      corresponding elliptic boundary value problem
     */
    const EllipticBVP<DIM>* ell_bvp_;

    /*!
      underlying frame
     */
    const AggregatedFrame<IBASIS,DIM>* frame_;


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
      helper routines for a (.. , ..). Entries in diagonal and non-diagonal
      blocks of the stiffness matrix have to be treated differently.
     */
    double a_same_patches(const typename AggregatedFrame<IBASIS,DIM>::Index& lambda,
			  const typename AggregatedFrame<IBASIS,DIM>::Index& nu,
			  const unsigned int q_order = 2) const;

    /*!
     */
    double integrate(const Index1D<IBASIS>& lambda,
		     const Index1D<IBASIS>& mu,
		     const FixedArray1D<Array1D<double>,DIM >& irregular_grid,
		     const int N_Gauss) const;

    /*!
      helper routines for a (.. , ..). Entries in diagonal and non-diagonal
      blocks of the stiffness matrix have to be treated differently.
     */
    double a_different_patches(const typename AggregatedFrame<IBASIS,DIM>::Index& lambda,
			       const typename AggregatedFrame<IBASIS,DIM>::Index& nu,
			       const unsigned int q_order = 2, const unsigned int rank = 1) const;


    double a_different_patches_adaptive(const typename AggregatedFrame<IBASIS,DIM>::Index& lambda,
					const typename AggregatedFrame<IBASIS,DIM>::Index& nu) const;


    // precompute the right-hand side
    void compute_rhs();

    // precompute diagonal of stiffness matrix
    void compute_diagonal();


    // right-hand side coefficients on a fine level, sorted by modulus
    Array1D<std::pair<typename AggregatedFrame<IBASIS,DIM>::Index,double> > fcoeffs;

    // right-hand side coefficients on a fine level, sorted by modulus
    InfiniteVector<double,typename AggregatedFrame<IBASIS,DIM>::Index> stiff_diagonal;


    // (squared) \ell_2 norm of the precomputed right-hand side
    double fnorm_sqr;

    // reminder: This keyword can only be applied to non-static
    // and non-const data members of a class. If a data member is declared mutable,
    // then it is legal to assign a value to this data member from
    // a const member function.
    // estimates for ||A|| and ||A^{-1}||
    mutable double normA, normAinv;

    QuadratureStrategy qstrat_; 

  };
}

#include <elliptic_equation.cpp>

#endif
