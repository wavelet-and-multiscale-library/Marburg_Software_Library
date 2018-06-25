// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner, Ulrich Friedrich , Philipp Keding, |
// | Alexander Sieber                                                   |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_CACHED_QUARKLET_RECRING_PROBLEM_H
#define	_WAVELETTL_CACHED_QUARKLET_RECRING_PROBLEM_H

#include <map>
#include <cmath>
#include <adaptive/compression.h>
#include <algebra/infinite_vector.h>
#include <algebra/vector.h>
#include <galerkin/galerkin_utils.h>
#include <galerkin/infinite_preconditioner.h>
#include <numerics/eigenvalues.h>


using MathTL::InfiniteVector;
using std::max;
using std::min;
using std::set;

namespace WaveletTL
{
    /*
     * This class provides a cache layer for generic (preconditioned, cf. precond.h)
     * infinite-dimensional matrix problems of the form
     *
     *    Au = P^{-1}LQ^{-1}u = P^{-1}F.
     *
     * The operator equation is not assumed to be induced by a local operator,
     * i.e., the cache class should also work in the case of integral operators.
     * All evaluations of the bilinear form a(.,.) are cached. The underlying
     * frame is assumed to be of recring tensor type as modeled by recring_frame.h.
     * Internally, the cache is managed as follows. The nonzero values of the bilinear
     * form a(.,.) are stored in a map.
     * For a given sublevel j\in\N^d the sublevels {k\in\N^d with ||j-k||_1 = dist}
     * are computed with add_levelsphere. The 1d integrals needed for the computation
     * of sphere(dist) are stored, since they are needed to compute sphere(dist+1).
     *
     * The template class CachedQuarkletRecRingProblem implements the minimal signature to be
     * used within the APPLY routine.
     */
    template <class PROBLEM>
    class CachedQuarkletRecRingProblem
#ifdef DYADIC
      : public FullyDiagonalQuarkletPreconditioner<typename PROBLEM::Index, 2>
    
#endif
#ifdef TRIVIAL
    : public TrivialPreconditioner<typename PROBLEM::Index>
#endif
#ifdef ENERGY
    : public FullyDiagonalEnergyNormPreconditioner<typename PROBLEM::Index>
#endif
#ifdef DYPLUSEN
    : public FullyDiagonalDyPlusEnNormPreconditioner<typename PROBLEM::Index, 2>
#endif

    {
    public:
        /*
         * make quarklet frame type accessible
         */
        typedef typename PROBLEM::Frame QuarkletFrame;

        /*
         * quarklet index class
         */
        typedef typename PROBLEM::Index Index;

        /*
         * type corresponding to the levels of the quarklet frame
         */
        typedef typename Index::level_type level_type;
        typedef typename Index::polynomial_type polynomial_type;
        /*
         * constructor from an uncached problem,
         * you can specify the estimates for ||A|| and ||A^{-1}||
         * (if zero, CachedQuarkletRecRingProblem will compute the estimates)
         */
        CachedQuarkletRecRingProblem(PROBLEM* P,
                       const double normA = 0.0,
                       const double normAinv = 0.0);

        /*
         * read access to the frame
         */
        const QuarkletFrame& frame() const { return problem->frame(); }

        /*
         * space dimension of the problem
         */
        static const int space_dimension = 2;

        /*
         * locality of the operator
         */
        inline static bool local_operator() { return PROBLEM::local_operator(); }

        /*
         * (half) order t of the operator
         */
        inline double operator_order() const { return problem->operator_order(); }

        /*
         * evaluate the diagonal preconditioner D
         * Values are taken from the cache if computed previously, otherwise stored. 
         * Maybe this is too much effort if only a(mu,mu) is of interest and no
         * other element of the muth column. In such cases a mid caching strategy might be effective
         */
        inline double D(const Index& lambda) const {
//            cout << "D" << endl;
//            return sqrt(a(lambda,lambda));//this needs to be changed to dyadic preconditioner

            return problem->D(lambda);

        }

        /*
         * evaluate the bilinear form a of the underlying problem
         * caching: compute and store all entries a(lambda2,nu) with
         * lambda2 having the same sublevel as lambda. add_ball relies on this caching strategy
         */
        double a(const Index& lambda,
                 const Index& nu) const;

        /*
         * estimate the spectral norm ||A||
         */
        double norm_A() const;


        /*
         * estimate the spectral norm ||A|| and ||A^{-1}|| for given offset
         */
        void normtest(unsigned int offsetj, unsigned int offsetp = 0) const;

        /*
         * estimate the spectral norm ||A^{-1}||
         */
        double norm_Ainv() const;

        /*
         * set the spectral norm ||A||
         */
        void set_normA(double value) { normA = value; }

        /*
         * set the spectral norm ||A^{-1}||
         */
        void set_normAinv(double value) { normAinv = value; }

        /*
         * estimate compressibility exponent s^*
         * (we assume that the coefficients a(x),q(x) are smooth)
         *
         * useless since only s^*=\infty is considered.
         */
         
        inline double s_star() const {
            return problem->s_star();
        }
         

        /*
         * estimate the compression constants alpha_k in
         * ||A-A_k|| <= alpha_k * 2^{-rho*k}
         *
         * we assume they are independent of k and hardcode some results from the theory.
         *
         * rho tends to 1/2 for tframe
         *
         * influences jp_tilde in apply_tensor.cpp
         */
        inline double alphak(const unsigned int k) const {
            switch (space_dimension){
                case 1:
                    return 2.4;
                    break;
                case 2:
                    return 8.2;
                    break;
                default:
                    return 48.0;
            }
            //return 2*norm_A(); // pessimistic
        }

        /*
         * evaluate the (unpreconditioned) right-hand side f
         */
        inline double f(const Index& lambda) const {
            return problem->f(lambda);

        }

        /*
         * change the righthandside f
         */
        void set_f(const Function<PROBLEM::space_dimension>* fnew);


        /*
         * read access to the stored preconditioned coefficients of the righthand side
         */
        //inline double get_fcoeffs(unsigned int i) const
        //{
        //    return problem->fcoeffs[i].second;
        //}

        /*
         * approximate the quarklet coefficient set of the preconditioned right-hand side F
         * within a prescribed \ell_2 error tolerance
         */
        inline void RHS(const double eta,
                 InfiniteVector<double, Index>& coeffs) const
        {
            problem->RHS(eta, coeffs);
        }
        
        inline void RHS(const double eta,
                 InfiniteVector<double, int>& coeffs) const
        {
            problem->RHS(eta, coeffs);
        }

        

        /*
         * compute (or estimate) ||F||_2
         */
        inline double F_norm() const { return problem->F_norm(); }

        /*
         * Called by APPLY // add_compressed_column
         * w += factor * (stiffness matrix entries in column lambda with ||nu-lambda|| <= range && ||nu|| <= maxlevel)
         * Hack? : works with non dynamic vector w of size = degrees_of_freedom
         *
         * "strategy" argument ist needed for compatibility with compression.h .
         * There is only 1 strategy for TFrame at the moment
         */
        void add_ball(const Index& lambda,
                     //InfiniteVector<double, Index>& w,
                     Vector<double>& w,
                     const int radius,
                     const double factor,
                     const int maxlevel,
                     const CompressionStrategy strategy = tensor_simple,
                     const bool precond = true,
                     const int maxpolynomial = 0,
                     const double A = 0,
                     const double B = 0) const;
          protected:
        //! the underlying (uncached) problem
        PROBLEM* problem;

        // type of one subblock in one block of stiffness matrix  A
        // entries are indexed by the number of the quarklet.
        typedef std::map<int, double> Subblock;
        
        // type of one block in one column of stiffness matrix  A
        // the key codes the sublevel, that data are the entries
        // a sublevel is a MultiIndex. Entries are indexed by the number of the
        // sublevel as given by the ordering in MultiIndex, beginning at jmin
        typedef std::map<int, Subblock> Block;

        // type of one column in the entry cache of A
        // the key codes the subpolynomial, that data are the entries
        // a subpolynomial is a MultiIndex. Entries are indexed by the number of the
        // sublevel as given by the ordering in MultiIndex, beginning at (0,0)
        typedef std::map<int, Block> Column;

        // type of the entry cache of A
        //typedef std::map<Index, Column> ColumnCache;
        // entries are indexed by the number of the (Quarklet-) Index
        typedef std::map<int, Column> ColumnCache;

        // entries cache for A (mutable to overcome the constness of add_column())
        mutable ColumnCache entries_cache;

        // estimates for ||A|| and ||A^{-1}||
        mutable double normA, normAinv;
    };

}

#include "cached_quarklet_recring_problem.cpp"

#endif	/* _WAVELETTL__CACHED_QUARKLET_RECRINGPROBLEM_H */


