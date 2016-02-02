// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2015                                            |
// | Thorsten Raasch, Manuel Werner, Ulrich Friedrich                   |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_QTBasis_H
#define _WAVELETTL_QTBasis_H

#include <list>

#include <algebra/infinite_vector.h>
#include <utils/fixed_array1d.h>
#include <utils/multiindex.h>
#include <utils/function.h>
#include <geometry/point.h>
#include <utils/array1d.h>
#include <algebra/matrix.h>

#include <general_domain/qtbasis_index.h>

// for convenience, include also some functionality
// #include "qtbasis_support.h"

using std::list;
using std::pair;
using MathTL::Point;
using MathTL::Function;
using MathTL::FixedArray1D;
using MathTL::MultiIndex;
using MathTL::InfiniteVector;
using MathTL::Array1D;

namespace WaveletTL
{

    /*
     * Implementation of the generalized tensor wavelet basis.
     * 
     * [CDFS13] Piecewise Tensor Product Wavelet Bases by Extensions and Approximation Rates
     *          N. Chegini, S. Dahlke, U. Friedrich, R. Stevenson
     *          Mathematics of Computation 82 (2013), 2157-2190. 
     * 
     * Whenever wavelet intersections are computed the result is stored in a cache.
     * 
     * Code may contain bugs when bases with different minimal levels are combines
     * Code contains several assert commands
     * Code contains several lines beginning with "CLEANUP". What follows are often comments for understanding/debugging. Sometimes assertions.
     * 
     */
    template <class IBASIS, unsigned int DIM = 2>
    class QTBasis {
    public:
        //! default constructor. 1 patch. Homogeneous Dirichlet BCs
        QTBasis();

 	/*
    	 * Constructor. For each patch specify its offset (spatial translation
         * from the origin), its neighbors and is boundary conditions. These infos
         * are given per face of the patches. Dirichlet boundary conditions for
    	 * the primal functions are specified, the dual functions will be constructed
         * to fulfill free b.c.'s.
         * An entry of -1 in neighbors means no neighbor in this direction
    	 * i-th direction at x=0 <-> index 2*i
    	 * i-th direction at x=1 <-> index 2*i+1
         * Which basis is extended is determined by the boundary conditions imposed at the interfaces.
    	 */
        QTBasis(const Array1D<Point<DIM, int> >& corners,
                const Array1D<FixedArray1D<int,2*DIM> >& neighbours,
                const Array1D<FixedArray1D<bool,2*DIM> >& bc) ;

        /*
    	 * Constructor with specified boundary condition orders
         * Warning: at the moment only boundary conditions of first order can be handled (size of bases_infact is static!)
         * So only use values of 0 and 1 as entries for bc.
    	 */
    	QTBasis(const Array1D<Point<DIM, int> >& corners,
                const Array1D<FixedArray1D<int,2*DIM> >& neighbours,
                const Array1D<FixedArray1D<int,2*DIM> >& bc) ;

        virtual ~QTBasis();


        //! Coarsest level j0 for each patch
    	inline const Array1D < MultiIndex<int,DIM> > j0() const { return j0_; }

        inline void set_jmax(const int jmax) {
            jmax_ = jmax;
            precompute_firstlast_wavelets();
            setup_full_collection();
    	}
        
        inline const unsigned int get_jmax() const {
            return jmax_;
        }

        //! Return number of patches
        inline const unsigned int get_nop() const {
            return j0_.size();
        }


        /*
         * Compute the indices of the first and last wavelet for all wavelet levels up to
         * \|\lambda\| == jmax. Indices are stored. 
         * first/last_wavelet routines do not compute indices any more
         */
        void precompute_firstlast_wavelets();
        
    	/*!
         * Setup full collection of wavelets between j0_ and jmax_.
         * This codes the mapping \N -> wavelet_indices
         * jmax_ has to be specified previously.
         */
    	void setup_full_collection();


        //! Read access to the bases
    	inline const Array1D<FixedArray1D<IBASIS*,DIM> > & bases() const { return bases_; }
        
        inline const IBASIS* get_bases(const int p, const unsigned int i) const { return bases_[p][i]; }

        /*
    	 * Geometric type of the support sets
    	 */
    	typedef struct {
      		int j[DIM];
      		int a[DIM];
      		int b[DIM];
    	} Support;

    	//! Wavelet index class
    	typedef QTIndex<IBASIS,DIM,QTBasis<IBASIS,DIM> > Index;
        
        //! make template argument accessible
        typedef IBASIS IntervalBasis;
        
        /*
         * space dimension of the problem
         */
        static const int space_dimension = DIM;

        /*
         * There are many ways the information whether a certain basis function is extended can be stored.
         * Let [j,a,b] be the support on the main patch, [j,a*,b*] the support after extending and [j,a**,b**] be the support where the offset of the patch is added to a,b.
         * Then you could store for instance
         * A: {j,a,b} offset and extension information can be extracted from somewhere else
         * B: {j,a,b,offset,exinfo}
         * C: {j,a*,b*,offset}
         * D: {j,a**,b**}
         * D gives the real support (human readable), as does C (where we have a little more effort).
         * A, B are tailored to save the program from evaluating too much, where in case A the minimal amount of redundant information is stored
         * At the moment variant A is implemented. evaluate, integrate ... are influenced. this means that the method does not differ from the tbasis case.
         */
        void support_local(const Index& lambda, Support& supp) const; // variant A
        void support_full(const Index& lambda, Support& supp) const; // variant D
        
    	//! Point evaluation of (derivatives of) primal generators or wavelets \psi_\lambda
        // first version checks whether x is contained in the support of \psi_\lambda
    	double evaluate_check(const unsigned int derivative,
                        const int lambdanum,
                        const Point<DIM> x) const;
        
        //! assumes x lies in the support of lambda
        double evaluate_simple(const unsigned int derivative,
                               const int lambdanum,
                               const Point<DIM> x) const;

    	/*
    	 * Helper function, integrate a smooth function f against a
    	 * primal generator or wavelet
    	 */
    	double integrate(const Function<DIM>* f,
                         const unsigned int lambdanum) const;
        
        /*
    	 * The inner product of two wavelets is not implemented here.
         * Instead, it is realized in cachedQTProblem
         * Reason: a cache for the 1D integrals makes sense!
         * However, it doesn't seem that it fits into this class.
    	double integrate(const unsigned int lambdanum,
                         const unsigned int munum) const;
    	 */
        
        /*
    	 * For a given function, compute all integrals w.r.t. the primal
    	 * or dual generators/wavelets \psi_\lambda with |\lambda|\le jmax.
    	 * - When integrating against the primal functions, the integrand has to be smooth
    	 *   to be accurately reproduced by the dual basis.
    	 * - When integration against dual functions is specified,
    	 *   we integrate against the primal ones instead and multiply the resulting
    	 *   coefficients with the inverse of the primal gramian.
    	 *
    	 * Maybe a thresholding of the returned coefficients is helpful (e.g. for
    	 * expansions of spline functions).
    	 */
    	void expand(const Function<DIM>* f,
                    const bool primal,
                    //const MultiIndex<int,DIM> jmax,
                    InfiniteVector<double,Index>& coeffs) const;
        void expand(const Function<DIM>* f,
                    const bool primal,
                    InfiniteVector<double,int>& coeffs) const;
        
        //! Index of first generator
    	Index first_generator() const;

        //! Index of last generator
        // doesnt make much sense since gen and wav are thrown together in one set
    	//Index last_generator() const;

    	//! Index of first wavelet on level j >= j0
    	//Index first_wavelet(const MultiIndex<int,DIM> j) const;
        
        //! Index of first wavelet on the level (j,p) coded by levelnum (as in get_level)
        // Caution: first_wavelet(0) == first_generator !! (Generators are treated as wavelets)
        inline Index first_wavelet(const unsigned int levelnum) const
        {
            return first_wavelets_[levelnum];
        }

        //! Index of first wavelet on level j >= ||j0||_1
        //Index first_wavelet(const int level) const;

    	//! Index of last wavelet on sublevel j >= j0
        // TODO handling of invalid argument
    	// Index last_wavelet(const MultiIndex<int,DIM> j) const;

        //! Index of last wavelet on the level (j,p) coded by levelnum (as in get_level)
        inline Index last_wavelet(const unsigned int levelnum) const
        {
            return last_wavelets_[levelnum];
        }
        
        //! Index of last wavelet on level j >= ||j0||_1
        // TODO handling of invalid argument
    	// Index last_wavelet(const unsigned int level) const;
        
        //! Number of the last wavelet with levelnorm level
        //unsigned int last_wavelet_num(const unsigned int level) const;
        
        //! Number of wavelets between coarsest and finest level
    	const inline int degrees_of_freedom() const { return full_collection_.size(); };

    	//! Get the wavelet index corresponding to a specified number
    	const inline Index* get_wavelet (const unsigned int number) const {
      		return &full_collection_[number];
    	}

        inline unsigned int get_levelnum(const MultiIndex<int,DIM> & j, const unsigned int p) const
        {
            MultiIndex<int,DIM> temp_mi(j);
            for (unsigned int i=0; i<DIM; ++i)
                temp_mi[i] -= j0_[p][i];
            return level_to_num_[p][temp_mi.number()];
        }
        
        inline void get_level(const unsigned int levelnum, MultiIndex<int,DIM>& j, unsigned int& p) const
        {
            j = num_to_level_[levelnum].first;
            p = num_to_level_[levelnum].second;
        }
        
        inline unsigned int get_numberoflevels() const
        {
            return num_to_level_.size();
        }
        
    	/*
    	 * Critical Sobolev regularity for the primal generators/wavelets.
    	 * We assume the same regularity in each dimension.
    	 */
    	inline static double primal_regularity() { return IBASIS::primal_regularity(); }

    	/*
    	 * Degree of polynomial reproduction for the primal generators/wavelets
    	 * We assume the same polynomial degree in each dimension.
    	 * */
    	inline static unsigned int primal_polynomial_degree() { return IBASIS::primal_polynomial_degree(); }

    	/*
    	 * Number of vanishing moments for the primal wavelets.
    	 * We assume the same number of vanishing moments in each dimension.
    	 */
    	inline static unsigned int primal_vanishing_moments() { return IBASIS::primal_vanishing_moments(); }

    
        /*
         * compute all intersecting mus for a given lambda.
         * lambda = (lambda_j,lambda_e,lambda_k) stems from the basis lambda_basisnum.
         * we are only interested in mus on level mu_level from the basis mu_basisnum
         * the boundary conditions of the bases imply their number: 00,01,10,11 == 0,1,2,3
         * 
         * reflected == true : reflect support of lambda, i.e.,
         * 
         *              scale_lambda = lambda.j()[i] + lambda.e()[i];
         *              temp_i1 = (1 << scale_lambda)-lam_k2;
         *              lam_k2 = (1 << scale_lambda)-lam_k1;
         *              lam_k1 = temp_i1;
         * (useful to reduce all needed integral cases in the extension business to an integral over [0,1])
         * 
         * generator == true : consider generators mu instead of wavelets (makes only sense for mu_level == minimal level j0 !)
         * 
         * returns false if there is no intersection
         * stores in kmin,kmax the k-values of the intersecting mus
         * kmin,kmax contain meaningless data if false is returned!
         * 
         * computed values are cached, e.g.,
         * nonreflected_wavs[lambda_basisnum*4+mu_basisnum][lambda_num][mu_level-j0] = (kmin,kmax)
         * reflected/gens similar
         * 
         * Caution:
         * If lambda is by consstruction not extended, 
         * then this method does NOT yield the wavelets on the ,e.g., 
         * left patch that are exteded to lambdas patch and that 
         * intersect its support.
         * It yields ALL wavelets/generators that intersect the reflected 
         * support of lambda. This includes NONboundary wavelets that by 
         * construction are not extended to lambdas patch!
         * Therefore in get_onedim_intersections the extension types 3 and 5 
         * are considered seperately and only in other cases 
         * get_cached_onedim_intersections is called.
         */
        void get_cached_onedim_intersections(const unsigned int lambda_basisnum,
                const unsigned int lambda_j,
                const unsigned int lambda_e,
                const unsigned int lambda_k,
                const bool reflected,
                const unsigned int mu_basisnum,
                const unsigned int mu_j,
                const bool generators,
                int& kmin,
                int& kmax);
        
        /*
         * Compute the type of intersection that a wavelet lambda and a wavelet from patch mu_p can have.
         * In each dimension their intersection is of one of the types {L,M,R}^2 described
         * in "intersecting_wavelets".
         * This method takes into account
         * - whether lambda is extended in a direction or not
         * - the extension info of the basis on patch mu_p
         * NOT taken into account:
         * - the actual position of lambda (if lambda is an inner wavelet).
         * 
         * Observe:
         * Method returns true if there are wavelets with the same type as lambda that intersect with wavelets on patch mu_p!
         * E.g., true is returned in the following setting:
         * lambda = inner wavelet, mu = right boundary wavelet that is extended to lambdas patch. 
         * Observe that if lambda is far away from the boundary over which mu is extended they will not intersect!
         * Therefore this is checked in get_onedim_intersections()
         */
        bool
        get_LMR_info(const unsigned int lambdanum, 
                //const MultiIndex<int,DIM> mu_j,
                const unsigned int mu_p,
                MultiIndex<unsigned int, DIM>& intinfo) const;
        
        /*
         * get_intersection_geometry
         * 
         * This method realizes: lambda_p & LMR_info => something simple to realize the sum over all patches.
         * 
         * Attention: needs specialized_intinfo (needs to include info from both wavelets! not only one wavelet and (j,p) of the other! get_LMR_info is NOT sufficient)
         * 
         * We need to iterate over the patches where lambda and mu intersect.
         * given lambda_patch and LMR_info -> compute:
         * - details of patch intersection: how many, which geometry?
         * - geometry == type (see below)
         * - for each type we pick a "centerpatch". It is the leftmost in each direction.
         *   In the global quadrangulation it has the number centerpatchnumber
         * 
         * We assume that the support intersection is given by a tensor product of intervals,
         * i.e., we DO NOT allow L-shaped support geometry. Note that the theory of tensor product bases allows this!
         * This means that not all possibilities of extension directions and boundary conditions (on \partial\Omega) 
         * can be realized with this implementation!
         * 
         * type = extrended 1st dim? 0:1 + extended in 2nd dim?0:2 + extended in 3rd dim?0:4
         * 
         * Example:
         * types of intersection geometry in 2d:
         * type 0:  (0,1)^2
         * type 1,2 :  (0,2)x(0,1); (0,1)x(0,2)
         * type 3: (0,2)^2
         * 
         * Not implemented: L-shaped:
         * 
         * types of intersection geometry in 3d: (like binary numbers 0,..,7)
         * type 0: (0,1)^3
         * type 1: (0,2)x(0,1)^2; 
         * type 2: (0,1)x(0,2)x(0,1); 
         * type 3: (0,2)^2x(0,1)
         * type 4: (0,1)^2x(0,2)
         * type 5: (0,2)x(0,1)x(0,2)
         * type 6: (0,1)x(0,2)^2
         * type 7: (0,2)^3
         * 
         * In the cache we store the integrals of lambda_i against mu_i, where the lambda_i is reflected if necessary. 
         * I.e., the meaning of the (one dimensional) integrals is w.r.t. the patch of mu.
         * in "orientation" we denote if this is the same orientation as for the "centerpatch", i.e.,
         * orientation[i] == true => nothing to do, false: reverse the ordering
         * 
         * intinfo[i] orientation centerpatch
         * 0 LL refl left_of_lambda
         * 1 LM ok left_of_lambda
         * 2 LR ok left_of_lambda
         * 3 ML refl same_as_lambda
         * 4 MM ok same_as_lambda
         * 5 MR refl same_as_lambda
         * 6 RL refl same_as_lambda
         * 7 RM ok right_of_lambda
         * 8 RR ok same_as_lambda
         * 
         * This method does not check whether the boundary conditions make sense.
         * It is possible to specify illegal BC in the constructor!
         * 
         * Stores results in
         * intersection_geometry_XXXX, therefore not "const"
         */
        void
        get_intersection_geometry(const unsigned int lambda_p, 
                const MultiIndex<unsigned int, DIM> intinfo, 
                unsigned int& type, 
                int& centerpatchnumber, 
                MultiIndex<bool, DIM> & orientation);
        
        /*
         * given the information from get_intersection_info() this method checks whether there is actually an overlap of lambda with some mu.
         * returns gen/wav_intersection == true if there is an intersection with a generator/wavelet
         * false otherwise.
         * The value of gen/wav_inersection determines whether kmingen, kminwav, kmaxgen, kmaxwav
         * is meaningful
         * 
         * calls get_cached_onedim_intersections, therefore not "const"
         */
        
        void get_onedim_intersections(const unsigned int intinfo_i,
            const unsigned int lami_j, 
            const unsigned int lami_e, 
            const int lami_k,
            const unsigned int lami_basisnum,
            const bool min_type_i,
            const unsigned int mui_j,
            const unsigned int mui_basisnum,
            int& kmingen_i,
            int& kmaxgen_i,
            int& kminwav_i,
            int& kmaxwav_i,
            bool& gen_intersection,
            bool& wav_intersection);
        
        /*
         * For the given lambda: return all wavelets mu on level levelnum that have intersecting support.
         * 
         * For each fixed type e_, the set of indices mu such that supp psi_mu intersects supp psi_nu 
         * is a cartesian product of sets corresponding to one dimensional support intersections.
         * Depending on e[i] these one dimensional intersecting sets are fully characterized by the first and last gen/wav in it.
         * If levelnum corresponds to a level (p,j) with j[i] > j0[p][i], then there are no generators for this level and dimension.
         * In this case kmingen[i],kmaxgen[i] do not have a meaning.
         * 
         * intinfo gives all infos about the type of intersection of the (fixed) wavelet lambda and all wavelets from patch p
         * intinfo[i] \in {0,1,...,8} == { L,M,R }^2
         * intinfo[i] = 
         * (L,x) => lambda is a left boundary wavelet and is extended left
         * (M,x) => lambda is a not extended
         * (R,x) => lambda is a right boundary wavelet and is extended right
         * (x,L) => patch levelnum.p is extended left (in this dimension)
         * (x,M) => patch levelnum.p is not extended (in this dimension)
         * (x,R) => patch levelnum.p is extended right (in this dimension)
         * The intersection may hit 1 or 2 onedimensional "segments" (is contained in 1 or 2 shifted unit intervals, but never in 3)
         * Values where 2 patches are involved in this direction: LL,LR, RL, RR
         * 
         * Example: Omega is divided into 4 patches (DIM=2):
         * 1 2
         * 3 4
         * Omega_1 is extended right, Omega_4 is extended left and upwards, Omega_2 and Omega_3 are not extended
         * Let lambda be a wavelet from patch 4 that is extended left and up, i.e., left in x and right in y direction
         * Then intersecting_wavelets(lambda.num, (arbitrary_j,patch_1).num,x,x,x,x,intinfo) yields
         * intinfo = (LR, RM) : both Lx, Rx stem from lambda, the remaining xR and xM describe that patch 1 is extended right in the first dimension and nowhere in the other
         * 
         * usage: intinfo makes it possible to easily compute which Haar wavelets hit the support intersection of lambda and all mus
         * 
         * output: 
         * {kmingen[i], ..., kmaxgen[i]} is the set of possible values that the ith entry of the translation index k of an intersecting mu, where mu[i] is a generator
         * {kminwav[i], ..., kmaxwav[i]} same for mu[i] of wavelet type
         * intinfo
         * 
         * returns false if there are no intersecting wavelets. Other output variables are meaningless in this case.
         * If levelnum corresponds to a level (p,j) with j[i] > j0[p][i], then there are no generators for this level and dimension.
         * In this case kmingen[i],kmaxgen[i] do not have a meaning.
         * 
         * Calls get_onedim_intersections, therefore not "const"
         */
        bool intersecting_wavelets(const unsigned int lambdanum, 
                const unsigned int levelnum,
                MultiIndex<int, DIM>& kmingen,
                MultiIndex<int, DIM>& kmaxgen,
                MultiIndex<int, DIM>& kminwav,
                MultiIndex<int, DIM>& kmaxwav,
                MultiIndex<bool, DIM>& gen_intersection,
                MultiIndex<bool, DIM>& wav_intersection,
                MultiIndex<unsigned int, DIM>& intinfo);
        
        inline
        Array1D<FixedArray1D<bool,2*DIM> > get_bc() const
        {
            return bc_;
        }
        
        inline
        FixedArray1D<IBASIS*, 4> get_bases_infact() const
        {
            return bases_infact_;
        }
        
         
        // read access to neighbours_
        inline
        int get_neighbours(const unsigned int patch, const unsigned int direction) const
        {
            return neighbours_[patch][direction];
        }
        
        inline
        Point<DIM, int> get_corners(const unsigned int patch) const
        {
            return corners_[patch];
        }
        // for faster access
        inline
        int get_corners(const unsigned int patch, const unsigned int direction) const
        {
            return corners_[patch][direction];
        }
        
        inline
        int get_j0(const unsigned int patch, const unsigned int direction) const
        {
            return j0_[patch][direction];
        }
        
        inline
        unsigned int get_numofbw() const
        {
            return numofbw_;
        }
               
        /*!
         * Evaluate a single primal/dual generator or wavelet \psi_\lambda
         * on a dyadic subgrid with 1<<resolution points in each direction and patch
         */
        Array1D<SampledMapping<DIM> > sampled_output(const unsigned int lambdanum,
                const double alpha,
                const bool primal,
                const int resolution);
        /*!
         * Evaluate a single primal/dual generator or wavelet \psi_\lambda
         * Add result to the previous sampling
         */
        void sampled_output(Array1D<SampledMapping<DIM> > & previous_sampling,
                const unsigned int lambdanum,
                const double alpha,
                const bool primal,
                const int resolution);
        /*!
         * Evaluate an arbitrary linear combination of primal/dual wavelets
         */
        Array1D<SampledMapping<DIM> > sampled_output(const InfiniteVector<double, int>& coeffs,
                const bool primal,
                const int resolution);
        
       
        
        
    //protected:
        
        
        //! Coarsest possible level j0
    	Array1D<MultiIndex<int,DIM> > j0_;

        /*
         *  Collection of first and last wavelet indices on all levels up to jmax
         *  Precomputed for speedup
         */
        Array1D<Index> first_wavelets_, last_wavelets_;
        
        /*
         * The ordering of the indices induces an ordering for the pairs (j,p)
         * For speedup, the mapping \N -> Pairs is stored in num_to_level and
         * the mapping Pairs -> \N is stored in level_to_num
         */
        Array1D<std::pair<MultiIndex<int,DIM>, int > > num_to_level_;
        Array1D<Array1D<int> > level_to_num_; // level_to_num[patchnum][num_of_level_wrt_this_patch], i.e. num_of_level[p][jnum] observe that p comes before jnum
        
    	//! Collection of all wavelets between coarsest and finest level
    	Array1D<Index> full_collection_;

        //! The instances of the 1D bases
        // sorted by order of BC
        // (ordered by norm, then lexicographically, i.e., BC = {00,01,10,11})
        // We assume that only BC up the the first order are used!!
// CLEANUP        
        // old (and wrong?) comment: // 0 = homegeneous, 1 = inhomogeneous. Entries correspond to {00,01,10,11}
    	FixedArray1D<IBASIS*, 4> bases_infact_;
        // for boundary conditions beyond first order: (currently not implemented!)
        // ordering of the boundary con ditions by 1-norm, then lexicographically:
        // 00, 01, 10, 02, 11, 20, 03, 12, 21, 30, 04, ...
        //map<IBASIS*> bases_infact;

        Array1D<FixedArray1D<bool,2*DIM> > bc_;
        
        // Redundant:
    	//! For faster access, all relevant pointers to the 1D bases
    	Array1D<FixedArray1D<IBASIS*,DIM> > bases_;

        // offset of the patches relative to the orign
        Array1D<Point<DIM, int> > corners_;

        // for each patch for each face: store number of neighbor the patch is extending to along this facet.
        // value = -1 if no extension along this face
        // Dont forget neighbors in diagonal directions! They aren't stored directly in this array (in the current basis version).
        // This may get complicated (and/or slow) for higher dimensions
        Array1D<FixedArray1D<int,2*DIM> > neighbours_;
        
        // Redundant:
        // neighbours_ together with the one dimensional boundary conditions imply, whether a patch is extended along a facet.
        // However, for convenience this info is stored here.
        Array1D<FixedArray1D<bool,2*DIM> > extinfo_;

    	//! Finest possible level jmax
    	unsigned int jmax_;

        
        // We assume PBasis is used for the one dimensional bases. 
        // Consequently the number of wavelets that has to be extended is (d+dT-2)/2 per side.
        // The number of wavelets with nontrivial values at the boundary may be smaller.
        // We have to extend all boundary adapted wavelets to ensure local support of the dual wavelets.
        // Only 1 generator is extended per side
        // We assume that D,DT are the same for all interval bases
        const unsigned int numofbw_; // number of wavelets that are extended on each side of the interval
        
        // we store which wavelets (or generators) of the Nth one dimensional basis (N=0,1,2,3) intersect 
        // the (reflected) support of the one dimensional wavelet with number num w.r.t. the Mth one dimensional basis
        
        // nonreflected_wavs[M*4+N][num][level] = (kmin,kmax) == all wavelets on level "level" from basis "N" that intersect wavelet "num" from basis "M". 
        // support of "num" is not reflected. (that case is stored in reflected_wavs)
        // intersections with generators from the basis "N" are stored in ..._gens. Those generators are always on the lowest level.
        typedef std::map<int, pair<int, int> > intToPair;
        typedef std::map<int, intToPair> intToIntToPair;
        FixedArray1D< intToIntToPair, 16>  nonreflected_wavs, reflected_wavs;
        FixedArray1D< intToPair, 16>  nonreflected_gens, reflected_gens;
        
        // store output of get_intersection_geometry
        // for each patch (row) and int_info  store the resulting centerpatch_number and geometry type
        // int_info is converted to a number: in 2d: 9*intinfo[0]+intinfo[1]
        Matrix<int> intersection_geometry_centerpatchnumber_;
        Matrix<int> intersection_geometry_type_;
        MultiIndex< Matrix<bool> ,DIM>  intersection_geometry_orientation_;
        
    };

}


#include "qtbasis.cpp"
#endif	/* _WAVELETTL_QTBasis_H */


/*
 * store Array1D<SampledMapping<DIM> > in matlab readable format
 */
void matlab_output (Array1D<SampledMapping<2> > & previous_sampling, std::ostream& os)
{
    for (unsigned p=0; p< previous_sampling.size(); ++p)
    {
        os << "x{" << p+1 << "} = ";
        print_matrix(previous_sampling[p].gridx(), os);
        os << ";" << std::endl;
        os << "y{" << p+1 << "} = ";
        print_matrix(previous_sampling[p].gridy(), os);
        os << ";" << std::endl;
        os << "z{" << p+1 << "} = ";
        print_matrix(previous_sampling[p].values(), os);
        os << ";"
            << std::endl;
    }
};

#if 0
  // dump the internal data to a file for analysis with Matlab
  std::ofstream dump;
  ostringstream filename;
  filename << "DSBasis" << "_"
 	   << "d" << d << "_"
 	   << "dt" << dT << "_"
	   << "bio" << bio << "_"
	   << "bc" << basis.get_s0() << basis.get_s1() << basis.get_sT0() << basis.get_sT1() << "_"
 	   << "data"
 	   << ".m";
  dump.open(filename.str().c_str());
  basis.dump_data(dump);
  dump.close();
#endif
