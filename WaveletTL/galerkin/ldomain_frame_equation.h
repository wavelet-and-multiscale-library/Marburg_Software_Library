// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_LDOMAIN_FRAME_EQUATION_H
#define _WAVELETTL_LDOMAIN_FRAME_EQUATION_H

#include <set>
#include <utils/fixed_array1d.h>
#include <utils/array1d.h>
#include <numerics/bvp.h>
#include <Ldomain/ldomain_frame_support.h>

#include <galerkin/galerkin_utils.h>
#include <galerkin/infinite_preconditioner.h>
#include <interval/indexq1D.h>
#include <interval/i_q_index.h>

using MathTL::FixedArray1D;
using MathTL::EllipticBVP;

namespace WaveletTL
{
  template <class IFRAME> class LDomainFrame;

  /*!
    This class models the (preconditioned) infinite-dimensional matrix problem
    
      Au = D^{-1}LD^{-1}u = D^{-1}F

    when reformulating an poisson boundary value problem on the L--shaped domain
    
      -div(grad u(x)) = f(x)

    with first order b.c.'s as an equivalent operator equation
    within \ell_2 by means of a quarklet frame.

    The corresponding bilinear form in

      L = (a(\psi_\nu,\psi_\lambda))_{\lambda,\nu}

    is

      a(u,v) = \int_Omega [grad u(x)grad v(x)] dx
      
    and the right-hand side is
 
      F(v) = \int_0^1 f(x)v(x) dx
   
    The evaluation of a(.,.) and f is possible for arguments \psi_\lambda
    which stem from a quarklet frame \Psi=\{\psi_\lambda\} of the corresponding
    function space over Omega.

    (Currently, ) this class is designed to work with the quarklet type frame.
    However, in the future, this may be changed to allow also other wavelet frames
    over the L--shaped domain.
  */
  //template <class IFRAME>
  template <class IFRAME,  class LDOMAINFRAME = LDomainFrame<IFRAME> >
  class LDomainFrameEquation
//     : public FullyDiagonalDyadicPreconditioner<typename LDomainBasis<IBASIS>::Index>
#ifdef DYADIC  
  : public FullyDiagonalQuarkletPreconditioner<typename LDOMAINFRAME::Index>

    
#else
#ifdef TRIVIAL
    : public TrivialPreconditioner<typename LDOMAINFRAME::Index>
#else
    : public FullyDiagonalEnergyNormPreconditioner<typename LDOMAINFRAME::Index>
#endif
#endif
  {
  public:
      
    /*!
    make template argument accessible
    */
    typedef LDOMAINFRAME Frame;
//    typedef LDomainFrame<IFRAME> Frame;
    /*!
      constructor from a boundary value problem 
      not used anymore
    */
//    LDomainFrameEquation(const EllipticBVP<2>* bvp,
//		    const bool precompute_rhs = true);
    
    
    LDomainFrameEquation(const EllipticBVP<2>* bvp, const Frame* frame,
		    const bool precompute_rhs = true);
    /*!
      constructor from a boundary value problem and specified b.c.'s
    */
    //LDomainFrameEquation(const EllipticBVP<2>* bvp,
	//	    const FixedArray1D<bool,8>& bc,
	//	    const bool precompute_rhs = true);

    /*!
      copy constructor
    */
    LDomainFrameEquation(const LDomainFrameEquation&);
    
    
    /*!
      quarklet index class
    */
    typedef typename Frame::Index Index;

    /*!
      read access to the frame
    */
    const Frame& frame() const { return *frame_; }
    
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
      evaluate the (unpreconditioned) bilinear form a;
      you can specify the order p of the quadrature rule, i.e.,
      (piecewise) polynomials of maximal degree p will be integrated exactly.
      Internally, we use an m-point composite tensor product Gauss rule adapted
      to the singular supports of the spline wavelets involved,
      so that m = (p+1)/2;
    */
    double a(const Index& lambda, const Index& nu,
	     const unsigned int p) const; // p=4

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
    double f(const Index& lambda) const;

    /*!
      approximate the wavelet coefficient set of the preconditioned right-hand side F
      within a prescribed \ell_2 error tolerance
    */
    void RHS(const double eta,
	     InfiniteVector<double,Index>& coeffs) const;
    
    void RHS(const double eta,
	     InfiniteVector<double,int>& coeffs) const;

    /*!
      compute (or estimate) ||F||_2
    */
    double F_norm() const { return sqrt(fnorm_sqr); }

    /*!
      set the boundary value problem
    */
    void set_bvp(const EllipticBVP<2>*);
    
    /*
     * set the maximal wavelet level jmax.
     * modifies
     *   frame_->full_collection, fcoeffs, fnorm_sqr
     */
//    inline void set_jpmax(const unsigned int jmax, const unsigned int pmax, const bool computerhs = true)
//    {
//        cout << "hier4" << endl;
//        frame_->set_jpmax(jmax, pmax);
//        cout << "hier5" << endl;
//        if (computerhs) {
//#ifndef DYADIC
//            compute_diagonal();
//#endif
//            compute_rhs();
//        }
//    }

  protected:
    const EllipticBVP<2>* bvp_;
    const Frame* frame_;

    // right-hand side coefficients on a fine level, sorted by modulus
    Array1D<std::pair<Index,double> > fcoeffs;
    Array1D<std::pair<int,double> > fcoeffs_int;

    // precompute the right-hand side
    void compute_rhs();
    
    
    
    // #####################################################################################
    // Caching of appearing 1D integrals when making use of the tensor product structure of the wavelets
    // during the evaluation of the bilinear form.
    typedef std::map<IndexQ1D<IFRAME>,double > Column1D;
    typedef std::map<IndexQ1D<IFRAME>,Column1D> One_D_IntegralCache;
    
    mutable One_D_IntegralCache one_d_integrals;
    // #####################################################################################

    /*
    @param lambda Index of first wavelet or generator.
    @param mu Index of second wavelet or generator.
    @param N_Gauss Number of Gauss quadrature knots to be used. This should be cosen equal
    to the spline order in the constant coefficient case to be sure to integrate exactly.
    @param dir The spatial direction under considerattion.
    @param supp_lambda Support of the reference wavelet on the cube having the function given by lambda
    as component in the direction dir.
    @param supp_mu Support of the reference wavelet on the cube having the function given by mu
    as component in the direction dir.
     */
    double integrate(const IndexQ1D<IFRAME>& lambda,
		     const IndexQ1D<IFRAME>& mu,
		     const int N_Gauss,
		     const int dir,
                     typename Frame::Support supp) const;
    
    /*!
      precomputation of the right-hand side
      (constness is not nice but necessary to have RHS a const function)
    */
    void compute_diagonal();
    
    //! Square root of coefficients on diagonal of stiffness matrix.
    Array1D<double> stiff_diagonal;

    // (squared) \ell_2 norm of the precomputed right-hand side
    double fnorm_sqr;

    // estimates for ||A|| and ||A^{-1}||
    mutable double normA, normAinv;
  };
}

#include <galerkin/ldomain_frame_equation.cpp>

#endif
