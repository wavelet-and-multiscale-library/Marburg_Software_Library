// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2013                                            |
// | Thorsten Raasch, Manuel Werner, Ulrich Friedrich                   |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_CACHED_QTPROBLEM_H
#define	_WAVELETTL_CACHED_QTPROBLEM_H


#include <adaptive/compression.h>
#include <map>
#include <utils/fixed_array1d.h>
#include <algebra/fixed_matrix.h>
#include <galerkin/infinite_preconditioner.h>
#include <interval/p_evaluate.h>
#include <iostream>
#include <fstream>
        
namespace WaveletTL
{
    /*
     * This class provides a cache layer for generic (preconditioned, cf. precond.h)
     * infinite-dimensional matrix problems of the form
     *
     *    Au = P^{-1}LQ^{-1}u = P^{-1}F,
     * 
     * where Lu = div a grad u + q u. 
     * The domain is assumed to be given by a quadrangulation, i.e., \Omega = \sum_k e_k (0,1)^d, e_k\in\Z^d.
     * The basis used for the discretization is given QTBasis. 
     * a and q are assumed to be given by Haar wavelet coefficients on each cube.
     * It is possible to change the coefficient vectors of a,q.
     * 
     * All evaluations of the bilinear form a(.,.) are cached.
     * Internally, the cache is managed as a map over the columns of A and then as a map over the level of the columns.
     * If an entry a(psi_lambda,psi_mu) is requested, the whole column a(psi_lambda,psi_nu) with |nu| == |mu| is computed and stored.
     * The needed integrals break down into one-dimensional integrals of a Haar wavelet and 2 wavelets from PBasis.
     * All these integrals are stored (gramian and for the first order derivatives).
     * Depending on the basis, there are only a finite number of different types of wavelets, depending on the boundary conditions, that need to be integrated against each other.
     * 
     * The template class CachedQTProblem implements the minimal signature to be
     * used within the APPLY routine. 
     * Here, it is assumed that A is s^*-compressible with s^* = \infty and consequently a certain compression strategy is possible.
     */

    template <class QTBASIS, unsigned int ONEDIMHAARCOUNT, unsigned int DER_ONEDIMHAARCOUNT = ONEDIMHAARCOUNT>
    class CachedQTProblem
    //: public FullyDiagonalEnergyNormPreconditioner<typename PROBLEM::Index>
    //: public FullyDiagonalEnergyNormPreconditioner<int>
    : public InfiniteDiagonalMatrix<double>
    {
        public:
        /*
         * for compatibility with cdd1: make QTBasis visible as 'WaveletBasis'
         */
        typedef QTBASIS WaveletBasis;
        /*
         * wavelet index class
         */
        typedef typename QTBASIS::Index Index;

        /*
         * type corresponding to the levels of the wavelet basis
         */
        typedef typename Index::level_type Index_lt;

        static const unsigned int DIM = QTBASIS::space_dimension;

// cleanup
        //const unsigned int haar_maxlevel;
        //static const int space_dimension = PROBLEM::space_dimension;

        /*
         * Constructor
         * Specify the basis and  a,q with respect to the same numbering of the cubes.
         * Specify the right-hand side f. 
         * Specify estimates for ||A|| and ||A^{-1}||
// TODO
         * (if zero, CachedQTProblem will compute estimates)
         */
        CachedQTProblem(QTBASIS* basis,
                const Array1D<FixedMatrix<double, DER_ONEDIMHAARCOUNT> > & agencoeffs,
                const Array1D<FixedMatrix<double, ONEDIMHAARCOUNT> > & qgencoeffs,
                const Function<DIM>* f,
                const char* rhs_filename,
                const double normA,
                const double normAinv);
        
        /* 
         * variant without initializing a,q,D or the rhs
         * is used in aff_lin_par_eq. There only the cache qualities of this class are needed
         */
        CachedQTProblem(QTBASIS* basis,
                const double normA,
                const double normAinv);

        /*
         * Evaluate the diagonal preconditioner D, i.e., the square root of the diagonal of A.
         * Values are taken from the cache if possible. 
         * If not, the whole column a_{mu,nu} with |mu|=|nu| is computed and stored.
         * 
         * Maybe this is too much effort!
         * Alternative: a separate cache for the diagonal of A.
         */
        inline double D(const int lambdanum) const
        {
            return diagonal_cache_[lambdanum];
        }
        

        /*
         * Evaluate the bilinear form a of the underlying problem
         * Caching: compute and store all onedim integrals needed for the entries a(mu2,nu) with
         * mu2 having the same sublevel as mu. 
         * 
         * We assume that the operator is local
         * 
         * add_ball (below) relies on this caching strategy.
         */
        double a(const int munum,
                 const int nunum);
        
        /*
         * only for compatibility with 
         * FullyDiagonalEnergyNormPreconditioner<int>
         * 
         * enables the use of InfiniteVector . scale (this,-1) <-- useful
         */
        inline 
        double 
        diag(const int& lambdanum) const
        {
            return diagonal_cache_[lambdanum];
        }
        
        inline double a(const int& munum,
                 const int& nunum) const
        {
            abort();
            return 0;
        }
        
        /* For compatibility with CDD1 */
        double F_norm() const { return sqrt(f_precond_norm_sqr_); }
        
        /*
         * read access to the basis
         */
        inline const QTBASIS* basis() const { return basis_; };
        
        /*
         * read access to the rhs function f
         */
        //inline const Function<DIM> * f() const { return f_; };
        
        /*
         * estimate the spectral norm ||A||
         * Set up stiffness matrix consisting of all wavelets up to level jmax.
         * also computes norm of A^{-1}
         */
        double norm_A(const unsigned int jmax = 3*DIM) const;


        /*
         * estimate the spectral norm ||A^{-1}||
         * Set up stiffness matrix consisting of all wavelets up to level jmax.
         * also computes norm of A^{-1}
         */
        double norm_Ainv(const unsigned int jmax = 2*DIM) const;
        
        /*
         * estimate the spectral norm ||A|| and ||A^{-1}|| for given offset
         */
        void normtest(unsigned int offset);

        /*
         * approximate the wavelet coefficient set of the preconditioned right-hand side F
         * within a prescribed \ell_2 error tolerance.
         * uses only entries from fcoeffs.
         */
        void RHS(const double eta,
                 InfiniteVector<double, Index>& coeffs) const;
        void RHS(const double eta,
                 InfiniteVector<double, int>& coeffs) const;
        
            // \| A -A^{jp} \| \leq D 2^{-\rho jp} 
        
        /*
         * estimate the compression constants alpha_k in
         * ||A-A_k|| <= alpha_k * 2^{-rho*k} (A is compressible with s* = \infty)
         * 
         * A_k has C*k^dim nontrivial entries per row/column which need C*k^dim many operations to be computed
         * 
         * we assume the alpha_k are independent of k. Values are hardcoded (and were computed for tbasis)
// TODO: compute alpha_k for qtbasis
         *
         * rho tends to 1/2 for tbasis
         *
         * method is used for computation of jp_tilde in apply_tensor.cpp
         */
        inline double alphak(const unsigned int k) const {
            switch (DIM){
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
         * Called by APPLY_TENSOR
         * Hack? : works with non dynamic vector w of size = degrees_of_freedom, this is probably faster than using an InfiniteVector
         *
         * "strategy" argument is needed for compatibility with compression.h .
         * There is only 1 strategy for TBasis at the moment: 
         * computes
         * w += factor * (stiffness matrix entries in column lambda with ||nu-lambda|| <= radius && ||nu|| <= jmax)
         * 
         * argumnent maxlevel is unused at the moment (build in for compatibility with tbasis)
         * Version with Index lambda instead of int lambdanum is partly implemented, but commented out. You can find it at the end of the cpp file
         */
#if 0 
        // index version is discontinued code!
        void add_ball(const Index& lambda,
                Vector<double>& w,
                const int radius,
                const double factor,
                const int maxlevel = 99,
                const CompressionStrategy strategy = tensor_simple,
                const bool precond = true);
#endif        
        void add_ball(const unsigned int& lambdanum,
                Vector<double>& w,
                const int radius,
                const double factor,
                const int maxlevel = 99,
                const CompressionStrategy strategy = tensor_simple,
                const bool precond = true);
        
        /*
         * called by add_ball for DIM >= 3
         * calls itself recursively
         * computes:
         * compute all levels with j[0] == center_j[0],..., j[fixed_dimensions] == center_j[fixed_dimensions]
         * and j[fixed_dimensions],...,j[DIM-1] inside the ball with radius radius around center_j[fixed_dimensions],...,center_j[DIM-1]
         * add factor*a(mu,lambda) to w for all mu inside such a level (divided by D(mu) if precond == true)
         * 
         * It may happen that the minimal levels for the patches are not the same. 
         * Therefore it may happen, that not all patches allow for a level with center_j[0],...,center_j[fixed_dimensions]
         * If this is the case it holds that set_is_restriction == true and j0_restriction_of_patches is a list of all patches that allow the current restrictions by center_j
         */
        void add_leveldisc_recurse(
                const unsigned int& lambdanum,
                Vector<double>& w,
                const MultiIndex<int,DIM> center_j,
                const int fixed_dimensions,
                const unsigned int radius,
                const double factor,
                const list<int> & j0_restriction_of_patches,
                const bool set_is_restriction,
                const bool precond);
        
        /*
         * access to the RHS : f(lambda)/D(Lambda)
         */
        inline double f(const unsigned int lambdanum) const
        {
            return fcoeffs_precond_unsorted_[lambdanum];
        }
        
#if 0        
        /*
         * applys the Galerkin system matrix corresponding to the given index set
         * 'window' to vector x. Missing entries will be computed&cached.
         * The k-th entry of x is associated to the k-th entry of window.
         * Longer Vectors x may be inserted into apply. res will have the same
         * length than x but all entries beyond length(window) will be 0
         */
        void apply(const std::set<int>& window,
                   const Vector<double>& x,
                   Vector<double>& res) const;
        
        template <class PROBLEM>
  void APPLY_TENSOR(const PROBLEM& P,
          const InfiniteVector<double, int>& v,
          const double eta,
          InfiniteVector<double, int>& w,
          const int jmax = 99,
          const CompressionStrategy strategy = tensor_simple,
          const bool preconditioning = true);
#endif
        
#if 0        
        inline
        set<int> get_patch_with_minimal_j0_(unsigned int i) const
        {
            return patch_with_minimal_j0_[i];
        }
        
//        inline
//        set<int> get_minimal_level_y() const
//        {
//            return minimal_level_y_;
//        }
        
        /* 
         * Same, but returns (D^-1 A D^-1)_{mu,nu}.
         * This is probably more efficient than calling a and D separately in setup_stiffness_matrix
         */
        /*
        double aprecond(const int munum,
                        const int nunum) const;
         * */

        /*
         * locality of the operator
         */
        inline static bool local_operator() { return PROBLEM::local_operator(); }

        /*
         * (half) order t of the operator
         */
        inline double operator_order() const { return problem->operator_order(); }
        
        

        
        
            
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
         * approximate the wavelet coefficient set of the preconditioned right-hand side F
         * within a prescribed \ell_2 error tolerance
         */
        inline void RHS(const double eta,
                 InfiniteVector<double, Index>& coeffs) const
        {
            problem->RHS(eta, coeffs);
        }
        
        

        /*
         * applys the galerkin system matrix corresponding to the given index set
         * 'window' to vector x. Missing entries will be computed
         */
        bool CG(const std::set<int>& window,
                const Vector<double> &b,
                Vector<double> &xk,
                const double tol,
                const unsigned int maxiter,
                unsigned int& iterations);

        /*
         * compute (or estimate) ||F||_2
         */
        inline double F_norm() const { return problem->F_norm(); }

        
        
        /*
         * Called by add_ball. Recursivly all levels with |J-|lambda||<=range are added.
         * All dimensions are visited with the current dimension denoted by current_dim.
         * Indices with the highest recursion depth (level = DIM) are added to w
         */
        void add_level_recurse(const Index& lambda,
                               Vector<double>& w,
                               const int radius,
                               const double factor,
                               const index_lt & current_level,
                               const int current_dim,
                               const int maxlevel,
                               const bool legal,
                               const double d1,
                               const bool precond) const;
#endif

       // protected:
        
        /*
         * w += factor * A_levelnum e_lambda
         * A_levelnum = matrix where all entries a_(mu,nu) with nu != lambda and |mu| != levelnum are 0
         * e_lambda = lambda_th unit vector
         */
        /*
        void add_level(const Index& lambda,
                Vector<double>& w,
                const int levelnum,
                const double factor) const;
        */
        void add_level(const unsigned int& lambdanum,
                Vector<double>& w,
                const unsigned int levelnum,
                const double factor,
                const bool precond);  
        
        
        
        // The cache:
        
        /* 
         * Idea 1: we never integrate a single wavelet lambda against a single wavelet mu. 
         * we always integrate a single wavelet lambda against ALL wavelets mu on some (dim-dimensional level) that intersect with lambda,
         * Idea 2: we do not store the values of dim-dimensional integrals. 
         * the DIM-dimensional integrals are decomposed into one-dimensional ones. 
         * For ALL one-dimensional pairs of wavelets lambda_i, mu_i that come out of the dim-dimensional indices:
         * we store the integrals of lambda_i and mu_i against ALL Haar-Generators (i.e. constant functions).
         * Idea 3: we only consider integrals on the unit interval.
         * 
         * Observation: lambda_i and mu_i might stem from different bases. Therefore different caches are needed -> more cases
         * 
         * Depending on the whole extension business, there are 3 cases of how lambda_i and mu_i are related to each other:
         * type I: LL, MM, RR: here all possible combinations of bases for lambda_i, mu_i may come together. 
         * If both wavelets are extended, the integral corresponds to a non extended one (i.e., the MM case)
         * type II: LM, RM: lambda_i is a boundary wavelet (BC on this side is known). 
         * we store lambda_i reflected integrated against the mu_i 
         * (supp_lambda_i is transformed into information w.r.t. mu_i's side)
         * type III: ML, MR: lambda_i is arbitrary. 
         * we only consider mu_is that are boundary wavelets. 
         * we store the integral of lambda_i reflected against mu_i.
         */ 
        
        /*
         * one-dim cache
         * type I:
         * mu_i_level != min_level -> wavelets:
         * [basistypenumber] -> lambda_i_num -> mu_i_level -> mu_i_k** -> entries
         * basistypenumber = 4*lambda_i_basistype+mu_i_basistype 
         * mu_i_k** is not the number that Index1D produces (unneccessary cost!) it is simply the k value
         * entries is a pair of FixedVectors of length #HaarGenerators containing the values of the integrals of
         * psi_mui psi_nui and psi_mui' psi_nui' against all Haar generators
         * mu_i_level == min_level -> additionally store generators:
         * [basistypenumber] -> lambda_i_num               -> mu_i_k** -> entries
         * type II:
         * same structure as type I. lambda_i_basisnum can only be 1,2,3 == BC 01 10 11. 
         * basistypenumber = 4*(lambda_i_basistype-1)+mu_i_basistype
         * type III:
         * same structure as type I. mu_i_bases can only be BC 01 10 11. 
         * Observe: if mui_basisnum corresponds to 11 and extinfo is ML (lambda is not extended, but mu is), only left boundary wavelets mu are considered. 
         * Similarly, if intinfo is MR, only mui that are right boundary wavelets are considered.
         * To store these Situations we need 2 caches for the case that mu[i] stems from a 11 basis
         * => mui_basiscode: BC 01 10 11 (ML) 11 (MR)  // (basiscode != basisnum)
         * basistypenumber = 4*lambda_i_basistype+mu_i_basiscode 
         * the number of different mu_i_k will be the same (i.e. since all boundary mu_is have the same support)
         * we nee to store integrals \int  mu.nu and \int mu'.nu' (with Haar generators)
         */
        
        //typedef FixedArray1D<double,ONEDIMHAARCOUNT> entries;
        typedef std::pair<FixedArray1D<double,ONEDIMHAARCOUNT>, FixedArray1D<double,DER_ONEDIMHAARCOUNT> > entries;
        typedef std::map<int, entries> Block;
        typedef std::map<int, Block> Column;
        typedef std::map<int, Column> ColumnCache;
        typedef FixedArray1D<ColumnCache,16> typeIcache;
        typedef FixedArray1D<ColumnCache,12> typeIIcache;
        typedef FixedArray1D<ColumnCache,12> typeIIIcache;
        typedef FixedArray1D<Column,16> typeIcachedGens; // there is only one level with generators, so we can omit one "map"
        typedef FixedArray1D<Column,12> typeIIcachedGens;
        typedef FixedArray1D<Column,16> typeIIIcachedGens;
        typeIcache typeIcache_;
        typeIIcache typeIIcache_;
        typeIIIcache typeIIIcache_;
        typeIcachedGens typeIcachedGens_;
        typeIIcachedGens typeIIcachedGens_;
        typeIIIcachedGens typeIIIcachedGens_;
        
        Vector<double> diagonal_cache_;
        /*
         * Input: 
         * patch (lambda's patch)
         * specialized_intinfo (needs to include info from both wavelets! not only one wavelet and (j,p) of the other! get_LMR_info is NOT sufficient)
         * integralshares (one dimensional integrals; lambda_i is reflected to mu's patch if necessary)
         * 
         * output: a(mu, lambda)
         * 
         * specialized_intinfo :
         * this now contains information about both wavelets!
         * eg. "LR" means "nu is extended left and mu is extended right. They are both boundary wavelets/generators"
         * The output of get_LMR_info is different: "LR" means:
         * "nu is extended left and the patch of mu is extended right, e.g., there exists some mu that are extended right, but there also exist some that are NOT extended right and that intersect nu.
         * nu is a boundary wavelets/generators. No information about the type of the other wavelet/generator is included.
         */  
        double compute_sum_over_patches (const int nu_p, 
                const MultiIndex<unsigned int, DIM> intinfo, 
                const FixedArray1D<entries,DIM> integralshares) const;
        /*
         * compose the entries a(mu,nu).
         * For i=0,...,DIM-2: use the onedim-integral values given by block_it.
         * For i=DIM-1: take all entries in the Block starting at block_it[DIM-1] (the remaining dimension)
         * blocksize = size of the block with first element block_it[DIM-1]
         * mu_adapted_info_it: analog to block_it
         * add factor times those entries to w[number],...,w[number+block_it.size()-1]
         * The patches intersecting the support of mu are computed using nu_p and the intinfo provided by mu_adapted_intinfo1 and _2
         * if (precond == true) devide w[number] by D(number)
         * 
         * ATTENTION: header is for arbitrary dimensions!
         * Internally we assume that agencoefs and qgencoeffs are matrices, iff DIM = 2
         * Otherwise we assume ONEDIMHAARCOUNT == 1 and compute the Poisson equation
         */
        void compose_wavelets (Vector<double>& w,
                const unsigned int start,
                const unsigned int blocksize,
                const double factor,
                const unsigned int nu_p,
                const FixedArray1D< typename Array1D<unsigned int>::const_iterator ,DIM> mu_adapted_intinfo_it,
                const FixedArray1D< typename Block::iterator, DIM> block_it,
                const bool precond) const;
        
        /*
         * compute one dimensional integrals 
         * \int g_\eta \psi_{\mu_i} \psi_{\nu_i}, \eta =0,...,ONEDIMHAARCOUNT-1
         * \int g_\eta \psi_{\mu_i}' \psi_{\nu_i}', \eta =0,...,ONEDIMHAARCOUNT-1
         * \mu_i and \nu_i are onedimensional wavelets with respect to the bases
         * nui_basisnum and mui_basisnum, respectively.
         */
        void compute_onedim_haar_integrals(const bool reflected,
                const int nui_j,
                const int nui_e,
                const int nui_k,
                const unsigned int nui_basisnum,
                const int mui_j,
                const int mui_e,
                const int mui_k,
                const unsigned int mui_basisnum,
                FixedArray1D<double,ONEDIMHAARCOUNT>& gram,
                FixedArray1D<double,DER_ONEDIMHAARCOUNT>& der) const;
        
        /*
         * compute square root of the diagonal entries of A.
         */
        void compute_D();
        
        
        /*
         * set up the fcoeffs_, f_precond_norm_sqr_
         * store f_coeffs_precond_unsorted_ to disc
         */
        void compute_rhs(const char* rhs_filename);

        /*
         * set up the fcoeffs_, f_precond_norm_sqr_
         * to this end: load f_coeffs_precond_unsorted_ from disc
         */
        void load_rhs(const char* rhs_filename);
        
        /*
         * return \int f psi_lambda, i.e.,
         * the unpreconditioned right-hand side f, integrated against a primal wavelet
         */
//        inline double compute_f(const unsigned int lambdanum) const
//        {
//            return basis_->integrate(f_,lambdanum);
//        }
        
        //! the underlying (uncached) problem
        QTBASIS* basis_;

        //! the coefficients of the partial differential equation
        // agencoefs_[patch][i][j] == generator coeff on patch 'patch', ith in x direction, jth in y direction
        // dont get confused with the directions: i: left-right; j: bottom-up
        //const InfiniteVector<double, int> acoeffs_;
        //const InfiniteVector<double, int> qcoeffs_;
        const Array1D<FixedMatrix<double, DER_ONEDIMHAARCOUNT> > agencoeffs_;
        const Array1D<FixedMatrix<double, ONEDIMHAARCOUNT> > qgencoeffs_;
        
        // estimates for ||A|| and ||A^{-1}||
        //mutable 
        double normA_, normAinv_;
        
        // right-hand-site f
        const Function<DIM> * f_;
        
        // \int f \psi_lambda coefficient vector
        // right-hand side coefficients on a fine level, sorted by modulus
        Array1D<std::pair<int,double> > fcoeffs_precond_sorted_;
        Vector<double> fcoeffs_precond_unsorted_;
        //Array1D<std::pair<typename WaveletBasis::Index,double> > fcoeffs;
        
        // (squared) \ell_2 norm of the precomputed and preconditioned(!) right-hand side
        double f_precond_norm_sqr_;
        
        
        
        // store the patches where j0[p][0] (j0[p][1]) is minimal
        // designed for the 2d case
        //list<int> minimal_level_x_, minimal_level_y_;
        FixedArray1D<list<int>,DIM> patch_with_minimal_j0_;
        
#if 0        
         /*
         * type IV: LR RL: here both lambda_i and mu_i are extended. 
         * However, with Idea2, we integrate against Haar-Generators.
         * Those only live on one patch. 
         * So is is sufficient so store: lambda_i reflected against mu_i. 
         * The integral values are contained in those computed for the type II and type III case.
         * Depending on the domain decomposition the data cached in type IV cache might be redundant or not.
        // // type IV:
        // // LR: possible BC for lambda_i and mu_i bases: {10,11} x {01, 11} 0,1,2,3 = 1001 1011 1101 1111
        // // RL: possible BC for lambda_i and mu_i bases: {01,11} x {10, 11} 0,1,2,3 = 0110 0111 1110 1111
        // // the number of different mu_i_k will be the same (i.e. since all boundary mu_is have the same support)
         */
        
        //typedef FixedArray1D<ColumnCache,8> typeIVcache;
        //typedef FixedArray1D<Column,8> typeIVcachedGens;
        //typeIVcache typeIVcache_;
        //typeIVcachedGens typeIVcachedGens_;
        
        // The cache:
        
        // type of one block in one column of stiffness matrix  A
        // entries are indexed by the number of the wavelet.
        typedef std::map<int, double> Block;

        // type of one column in the entry cache of A
        // the key codes the sublevel, that data are the entries
        // a sublevel is a MultiIndex. Entries are indexed by the number of the
        // sublevel as given by the ordering in MultiIndex, beginning at j_min
        typedef std::map<int, Block> Column;

        // type of the entry cache of A
        //typedef std::map<Index, Column> ColumnCache;
        // entries are indexed by the number of the (Wavelet-) Index
        typedef std::map<int, Column> ColumnCache;

        // entries cache for A (i.e. a,q) 
        //(mutable to overcome the constness of add_column())
        // mutable 
        ColumnCache entries_cache_;
        
        // Matrix cache for the Haar matrices
        //typedef std::map<int, double> Block; // entries on one level for one fixed Haar wavelet
        //typedef std::map<int, Block> BlockHaar; // entries on one level for all Haar wavelets
        //typedef std::map<int, BlockHaar> ColumnHaar; // a column containing all entries for all Haar wavelets
        //typedef std::map<int, ColumnHaar> ColumnCacheHaar; // all column containers together
        // Haar cache: sort by : Haar index -> column index -> level number -> entries!
        typedef std::map<int, ColumnCache> HaarMatrixCache;
        HaarMatrixCache haar_cache_a_, haar_cache_q_;
        
        // 1d cache: sort by : Haar index -> column index -> level number -> type of 1d basis involved -> entries!
        // one additional layer
        typedef std::map<int, HaarMatrixCache> Haar1dCache;
        Haar1dCache haar_1d_diff_cache_, haar_1d_gram_cache;
#endif        
#if 0
        /*
         * read access to the basis
         */
        inline const QTBASIS& basis() const { return basis_; }

        /*
         * Estimate the spectral norm ||A||
         * Recomputes estimate if stored value is 0.
         */
        double norm_A() const;

// CLEANUP
        /*
         * Estimate the spectral norm ||A|| and ||A^{-1}|| for given offset
         */
        void normtest(unsigned int offset) const;

        /*
         * Estimate the spectral norm ||A^{-1}||
         * Recomputes estimate if stored value is 0.
         */
        double norm_Ainv() const;

        /*
         * estimate the compression constants alpha_k in
         * ||A-A_k|| <= alpha_k * 2^{-rho*k}
         *
         * We assume they are independent of k and hardcode (HACK!) some numbers.
         *
         * rho is close to 1/2 for qtbasis
         *
         * influences jp_tilde in apply_tensor.cpp
         */
        inline double alphak(const unsigned int k) const {
            switch (DIM){
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
         * change/set the righthandside f.
         * fcoeffs will be recomputed.
         */
        void set_f(const Function<DIM>* fnew);
        /*
         * approximate the wavelet coefficient set of the preconditioned right-hand side Q^-1F
         * within a prescribed \ell_2 error tolerance
         */
        void RHS(const double eta,
                 InfiniteVector<double, Index>& coeffs) const;
#endif
// CLEANUP
        // nice routines, but not needed for apply - only for the non adaptive case
#if 0
        /*
         * applys the galerkin system matrix corresponding to the given index set
         * 'window' to vector x. Missing entries will be computed
         * The k-th entry of x is associated to the k-th entry of window.
         * Longer Vectors x may be inserted into apply. res will have the same
         * length than x but all entries beyond length(window) will be 0
         */
        void apply(const std::set<int>& window,
                   const Vector<double>& x,
                   Vector<double>& res) const;

        /*
         * applys the galerkin system matrix corresponding to the given index set
         * 'window' to vector x. Missing entries will be computed
         */
        bool CG(const std::set<int>& window,
                const Vector<double> &b,
                Vector<double> &xk,
                const double tol,
                const unsigned int maxiter,
                unsigned int& iterations);
        /*
         * compute (or estimate) ||F||_2
         */
        inline double F_norm() const { return sqrt(f_precond_norm_sqr_); }

        /*
         * Called by APPLY // add_compressed_column
         * w += factor * (stiffness matrix entries in column lambda with ||nu-lambda|| <= radius && ||nu|| <= maxlevel)
         * Hack? : works with non dynamic vector w of size = degrees_of_freedom
         *
         * "strategy" argument ist needed for compatibility with compression.h .
         * There is only 1 strategy for TBasis at the moment
         */
        void add_ball(const Index& lambda,
                     //InfiniteVector<double, Index>& w,
                     Vector<double>& w,
                     const int radius,
                     const double factor,
                     const int maxlevel,
                     const CompressionStrategy strategy = tensor_simple,
                     const bool precond = true) const;
        /*
         * Called by add_ball. Add all levels with |J-|lambda||<=range by recursion.
         * All dimensions are visited with the current dimension denoted by current_dim.
         * Indices with the highest recursion depth (level = DIM) are added to w
         * 
         * Specialization to 1,2,3 space dimensions may be faster
         */
        void add_level_recurse(const Index& lambda,
                               Vector<double>& w,
                               const int radius,
                               const double factor,
                               const Index_lt & current_level,
                               const int current_dim,
                               const int maxlevel,
                               const bool legal,
                               const double d1,
                               const bool precond) const;
#endif
    };
}
#include "cached_qtproblem.cpp"
#endif	/* _WAVELETTL_CACHED_QTPROBLEM_H */

