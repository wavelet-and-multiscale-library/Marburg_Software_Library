// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner, Ulrich Friedrich                   |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_TFRAME_EQUATION_H
#define	_WAVELETTL_TFRAME_EQUATION_H

#include <set>
#include <utils/fixed_array1d.h>
#include <utils/array1d.h>
#include <numerics/bvp.h>
#include <cube/tframe_support.h>
#include <utils/function.h>
#include <geometry/point.h>

#include <galerkin/galerkin_utils.h>
#include <galerkin/infinite_preconditioner.h>
#include <interval/pq_expansion.h>
#include <interval/indexq1D.h>
#include <interval/i_q_index.h>

using MathTL::FixedArray1D;
using MathTL::EllipticBVP;

namespace WaveletTL
{
    template <class IFRAME, unsigned int DIM> class TensorFrame;


    /*!
    This class models the (preconditioned) infinite-dimensional matrix problem

      Au = D^{-1}LD^{-1}u = D^{-1}F

    when reformulating an elliptic boundary value problem on the cube [0,1]^d

      -div(a(x)grad u(x)) + q(x)u(x) = f(x)

    with first (Dirichlet) or second (Neumann) order b.c.'s (modeled in
    the class EllipticBVP) as an equivalent operator equation
    within \ell_2 by means of a wavelet frame.

    The corresponding bilinear form in

      L = (a(\psi_\nu,\psi_\lambda))_{\lambda,\nu}

    is

      a(u,v) = \int_Omega [a(x)grad u(x)grad v(x)+q(x)u(x)v(x)] dx

    and the right-hand side is

      F(v) = \int_0^1 f(x)v(x) dx

    The evaluation of a(.,.) and f is possible for arguments \psi_\lambda
    which stem from a wavelet frame \Psi=\{\psi_\lambda\} of the corresponding
    function space over Omega.
    \Psi is assumed to be of tensor structure. Unlike in cube_equation it is
    assumed that the minimal level of \Psi is of MultiIndex type, as modeled
    by tframe. To achieve independence from the concrete choice of \Psi,
    the wavelet frame class is given as a template parameter TENSORFRAME.
    It should have a constructor of the form

       TENSORFRAME::TENSORFRAME(const FixedArray1D<bool,2*DIM>& bc);

    where bc indicates the enforcement of homogeneous Dirichlet boundary conditions (true).
    A natural concrete value for TENSORFRAME is the TensorFrame<DSFrame<d,dT> >.
    */
 
    template <class IFRAME, unsigned int DIM, class TENSORFRAME = TensorFrame<IFRAME,DIM> >
    class TensorFrameEquation
    //     : public FullyDiagonalDyadicPreconditioner<typename TENSORFRAME::Index>
#ifdef DYADIC  
    : public FullyDiagonalQuarkletPreconditioner<typename TENSORFRAME::Index, DIM>   
#endif
#ifdef TRIVIAL
    : public TrivialPreconditioner<typename TENSORFRAME::Index>
#endif
#ifdef ENERGY
    : public FullyDiagonalEnergyNormPreconditioner<typename TENSORFRAME::Index>
#endif
#ifdef DYPLUSEN
    : public FullyDiagonalDyPlusEnNormPreconditioner<typename TENSORFRAME::Index, DIM>
#endif
    {
    public:
        /*!
    make template argument accessible
    */
         typedef TENSORFRAME Frame;
        /*
         * constructor from a boundary value problem and specified b.c.'s
         */

//        TensorFrameEquation(EllipticBVP<DIM>* bvp,
//                       const FixedArray1D<bool,2*DIM>& bc,
//                       const bool precompute_rhs = true);
        
//     /*!
//     */
//     TensorFrameEquation(const TENSORFRAME* frame,
// 		 const EllipticBVP<DIM>* bvp,
// 		 const FixedArray1D<bool,2*DIM>& bc);


//     /*!
//     */
//     TensorFrameEquation(const TENSORFRAME* frame,
// 		 const EllipticBVP<DIM>* bvp,
// 		 const FixedArray1D<int,2*DIM>& bc);

        /*
         * constructor from a boundary value problem and specified b.c.'s
         */
//        TensorFrameEquation(const EllipticBVP<DIM>* bvp,
//                       const FixedArray1D<int,2*DIM>& bc,
//                       const bool precompute_rhs = true);

		/*
		 * Constructor from a boundary value problem and a given tensor frame.
		 * The tensor frame is supposed to incorporate the desired boundary
		 * conditions and maximal levels 'jmax' and 'pmax' should have been set
		 * to it beforehand (generally by a call of set_jpmax(int,int)).
		 */
        TensorFrameEquation(const EllipticBVP<DIM>* bvp, const Frame* frame, const bool precompute_rhs = true);
        
        //constructor with given neumann-bc function
        TensorFrameEquation(const EllipticBVP<DIM>* bvp, const Frame* frame, const Function<DIM>* phi, const bool precompute_rhs = true);

        /*
         * copy constructor
         */
        TensorFrameEquation(const TensorFrameEquation&);     

        /*
         * quarklet index class
         */
        typedef typename Frame::Index Index;
        
        /*
         * read access to the frame
         */
        const Frame& frame() const { return *frame_; }

        /*
         * space dimension of the problem
         */
        static const int space_dimension = DIM;

        /*
         * differential operators are local
         */
        static bool local_operator() { return true; }
        

        /*
         * (half) order t of the operator
         * (inherited from FullyDiagonalEnergyNormPreconditioner)
         */
        double operator_order() const { return 1.; }

        /*
         * evaluate the diagonal preconditioner D
         */
        double D(const Index& lambda) const;

        /*
         * evaluate the (unpreconditioned) bilinear form a
         * (inherited from FullyDiagonalEnergyNormPreconditioner)
         */
        double a(const Index& lambda,
                 const Index& nu) const;

        /*
         * evaluate the (unpreconditioned) bilinear form a;
         * you can specify the order p of the quadrature rule, i.e.,
         * (piecewise) polynomials of maximal degree p will be integrated exactly.
         * Internally, we use an m-point composite tensor product Gauss rule adapted
         * to the singular supports of the spline wavelets involved,
         * so that m = (p+1)/2;
         */
        double a(const Index& lambda,
                 const Index& nu,
                 const unsigned int p) const;
        
                    
        /*
         * estimate the spectral norm ||A||
         * PERFORMANCE :: use setup_full_collection
         */
        double norm_A() const;

        /*
         * estimate the spectral norm ||A^{-1}||
         */
        double norm_Ainv() const;

        /*
         * estimate compressibility exponent s^*
         * (we assume that the coefficients a(x),q(x) are smooth)
         */
        double s_star() const {return 1.0;}

        /*
         * estimate the compression constants alpha_k in
         * ||A-A_k|| <= alpha_k * 2^{-s*k}
         */
        inline double alphak(const unsigned int k) const {
            return 2*norm_A(); // suboptimal
        }

        /*
         * evaluate the (unpreconditioned) right-hand side f
         */
        double f(const Index& lambda) const;

        /*
         * approximate the wavelet coefficient set of the preconditioned right-hand side F
         * within a prescribed \ell_2 error tolerance.
         * uses only entries from fcoeffs.
         */
        void RHS(const double eta,
                 InfiniteVector<double,Index>& coeffs) const;
        void RHS(const double eta,
                 InfiniteVector<double, int>& coeffs) const;

        /*
         * compute (or estimate) ||F||_2
         */
        double F_norm() const { return sqrt(fnorm_sqr); }

        /*
         * set the boundary value problem
         */
        void set_bvp(const EllipticBVP<DIM>*);

        /*
         * set or change the righthandside
         */
        void set_f(const Function<DIM>* fnew);
        
        /*!
         function for neumann-bc
        */
        const double phi(const Point<DIM>& x) const
        {
            return phi_->value(x);
        }
        
        /*
         * set the maximal wavelet level jmax.
         * modifies
         *   frame_.full_collection, fcoeffs, fnorm_sqr
         */
//        inline void set_jpmax(const unsigned int jmax, const unsigned int pmax, const bool computerhs = true)
//        {
////            frame_.set_jpmax(jmax, pmax);
////            if (computerhs) compute_rhs();
//        }
        
//        double integrate(const int lp, const int lj, const int le, const int lk,
//		   const int mp, const int mj, const int me, const int mk, const int derivative=0) const;

    protected:
        const bool nonzeroneumann_;
        const EllipticBVP<DIM>* bvp_;
        const Frame* frame_;
        
        //! function for neumann-bc
        const Function<DIM>* phi_;
        
        // right-hand side coefficients on a fine level, sorted by modulus
        Array1D<std::pair<Index,double> > fcoeffs;
        Array1D<std::pair<int,double> > fcoeffs_int;

        // precompute the right-hand side
        void compute_rhs();
        
        // #####################################################################################
        // Caching of appearing 1D integrals when making use of the tensor product structure of the wavelets
        // during the evaluation of the bilinear form.
//        typedef std::map<Index1D,double > Column1D;
//        typedef std::map<Index1D,Column1D> One_D_IntegralCache;
        typedef std::map<IndexQ1D<IFRAME>,double > Column1D;
        typedef std::map<IndexQ1D<IFRAME>,Column1D> One_D_IntegralCache;

        mutable One_D_IntegralCache one_d_integrals;
        // #####################################################################################
        
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

#include <galerkin/tframe_equation.cpp>

#endif	/* _WAVELETTL_TFRAME_EQUATION_H */

