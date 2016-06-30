// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Ulrich Friedrich                                                   |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_TBASIS_H_
#define _WAVELETTL_TBASIS_H_
#define _TBASIS_DEBUGLEVEL_  1 // more output, more tests

/*
 * replace online computation of first_ and last_ wavelet indices by precomputation
 * This shifts work from the CPU to the Memory. 
 * Its faster in simple tests. Not sure about Programs that use a lot of memory.
 * 1 == precomputation is active
 */
#define _PRECOMPUTE_FIRSTLAST_WAVELETS 1 

#include <list>

#include <algebra/infinite_vector.h>
#include <utils/fixed_array1d.h>
#include <utils/multiindex.h>
#include <cube/tbasis_index.h>
#include <utils/function.h>
#include <geometry/point.h>
#include <utils/array1d.h>

#include <numerics/gauss_data.h>

// for convenience, include also some functionality
#include <cube/tbasis_support.h>

using std::list;
using MathTL::Point;
using MathTL::Function;
using MathTL::FixedArray1D;
using MathTL::MultiIndex;
using MathTL::InfiniteVector;
using MathTL::Array1D;

namespace WaveletTL
{
    /*
     * Template class for tensor product wavelet bases on the d-dimensional unit cube [0,1]^d (or mapped versions thereof) of the type
     * (V_0 \oplus W_0 \oplus W_1 \oplus W_2 \oplus ...)\times (V_0 \oplus W_0 \oplus W_1 \oplus W_2 \oplus ...)
     * For each of the 2*d facets, you can specify the orders s_i, sT_i of the Dirichlet boundary conditions of
     * the primal and dual basis in the normal direction (this is tailored for the DSBasis constructor...).
     */
    template <class IBASIS, unsigned int DIM = 2>
    class TensorBasis
    {
    public:
    	//! Interval basis
    	typedef IBASIS IntervalBasis;
        
        //! Default constructor (no b.c.'s)
        TensorBasis();

        /* Constructor with specified boundary condition orders
         * i-th direction at x=0 <-> index 2*i
         * i-th direction at x=1 <-> index 2*i+1
         */
        /* works only with DS Basis, not with PBasis
    	TensorBasis(const FixedArray1D<int,2*DIM>& s, const FixedArray1D<int,2*DIM>& sT);
         */

    	/*
    	 * Constructor with specified boundary condition orders
    	 * i-th direction at x=0 <-> index 2*i
    	 * i-th direction at x=1 <-> index 2*i+1
    	 */
    	TensorBasis(const FixedArray1D<int,2*DIM>& s);

    	/*
    	 * Constructor with specified Dirichlet boundary conditions for
    	 * the primal functions, the dual functions will be constructed to
    	 * fulfill free b.c.'s
    	 */
    	TensorBasis(const FixedArray1D<bool,2*DIM>& bc);

    	/*
    	 * Constructor with precomputed instances of the 1D bases;
    	 * in this case, the pointers are reused and
    	 * not deleted at destruction time.
    	 */
        TensorBasis(const FixedArray1D<IBASIS*,DIM> bases);

    	//! Destructor
    	~TensorBasis();

    	//! Coarsest possible level j0
    	inline const MultiIndex<int,DIM> j0() const { return j0_; }

    	inline void set_jmax(const MultiIndex<int,DIM> jmax) {
      		jmax_ = multi_degree(jmax);
      		setup_full_collection();
#if _PRECOMPUTE_FIRSTLAST_WAVELETS
                precompute_firstlast_wavelets();
#endif
    	}

        inline const unsigned int get_jmax() const {
            return jmax_;
        }

        inline void set_jmax(const int jmax) {
      		jmax_ = jmax;
      		setup_full_collection();
#if _PRECOMPUTE_FIRSTLAST_WAVELETS
                precompute_firstlast_wavelets();
#endif
    	}
    	//! Wavelet index class
    	typedef TensorIndex<IBASIS,DIM,TensorBasis<IBASIS,DIM> > Index;

    	//! Read access to the bases
    	inline const FixedArray1D<IBASIS*,DIM>& bases() const { return bases_; }

    	/*
    	 * Geometric type of the support sets
         * (j,a,b) <-> 2^{-j_1}[a_1,b_1]x...x2^{-j_DIM}[a_n,b_DIM]
    	 */
    	typedef struct {
      		int j[DIM];
      		int a[DIM];
      		int b[DIM];
    	} Support;


        /*
         * For a given interval basis IBASIS, compute a cube
         * 2^{-j_}<a_,b_> = 2^{-j_1}[a_1,b_1]x...x2^{-j_n}[a_n,b_n]
         * which contains the support of a single primal cube generator
         * or wavelet psi_lambda.
         */
    	void support(const Index& lambda, Support& supp) const;

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

    	//! Index of first generator 
    	Index first_generator() const;

        /*! 
         * For compatibility reasons, returns the same (!) as first_generator()
         * parameter j is ignored !
         */
        Index first_generator(const unsigned int j) const;
        Index first_generator(const MultiIndex<int,DIM> j) const;

        //! Index of last generator on level j0
    	Index last_generator() const;

    	/*! 
         * Index of first wavelet on level j >= j0.
         * Method does not check, whether j is valid, i.e., if j-j0 is nonnegative
         */
    	Index first_wavelet(const MultiIndex<int,DIM> j) const;
        
        //! Index of first wavelet on level j >= ||j0||_1
        Index first_wavelet(const int levelsum) const;

    	/*! 
         * Index of last wavelet on sublevel j >= j0.
         * Method does not check, whether j is valid, i.e., if j-j0 is nonnegative
         */
    	Index last_wavelet(const MultiIndex<int,DIM> j) const;

        //! Index of last wavelet on level j >= ||j0||_1
    	Index last_wavelet(const int levelsum) const;

#if _PRECOMPUTE_FIRSTLAST_WAVELETS
        /*
         * Compute the indices of the first and last wavelet for all wavelet levels up to
         * \|\lambda\| == jmax. Indices are stored. 
         * first/last_wavelet routines do not compute indices any more
         */
        void precompute_firstlast_wavelets();
#endif
        
    	/*!
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
         * 
         * For the special case that f is a tensor a faster routine can be developed.
    	 */
    	void expand(const Function<DIM>* f,
                    const bool primal,
                    const MultiIndex<int,DIM> jmax,
                    InfiniteVector<double,Index>& coeffs) const;

        /*!
         * As above, but with integer max level
         */

        void expand(const Function<DIM>* f,
                    const bool primal,
                    const unsigned int jmax,
                    InfiniteVector<double,Index>& coeffs) const;
    	/*
    	 * Helper function, integrate a smooth function f against a
    	 * primal generator or wavelet
    	 */
    	double integrate(const Function<DIM>* f,
                         const Index& lambda) const;

    	//! Point evaluation of (derivatives of) primal generators or wavelets \psi_\lambda
    	double evaluate(const unsigned int derivative,
                        const Index& lambda,
                        const Point<DIM> x) const;

    	/*!
         * Compute all wavelet indices between beginning at j0_ and the last_wavelet 
         * with \|level\|\leq jmax_, i.e.,
         * degrees_of_freedom = last_wavelet_num<IBASIS,DIM,TensorBasis<IBASIS,DIM> >(this, jmax_) +1
         * many wavelet indices are computed. They are stored in full_collection
         * This codes the mapping \N -> wavelet_indices
         */
    	void setup_full_collection();


    	//! Number of wavelets between coarsest and finest level
    	const int degrees_of_freedom() const { return full_collection.size(); };

    	//! Get the wavelet index corresponding to a specified number
    	const inline Index* get_wavelet (const int number) const {
      		return &full_collection[number];
    	}

    protected:
    	//! Collection of all wavelets between coarsest and finest level
    	Array1D<Index> full_collection;
        
#if _PRECOMPUTE_FIRSTLAST_WAVELETS
        /*
         *  Collection of first and last wavelet indices on all levels up to jmax
         *  Precomputed for speedup
         */
        Array1D<Index> first_wavelets, last_wavelets;
#endif
        
        //! Coarsest possible level j0
    	MultiIndex<int,DIM> j0_;
    	//int j0_[DIM];

    	//! Finest possible level j0
    	//MultiIndex<int,DIM> jmax_;
        // wavelet indices with \|level\|\leq jmax_ are stored in full_collection
    	unsigned int jmax_;

    	/*
    	 * The instances of the 1D bases
    	 */
    	list<IBASIS*> bases_infact;

    	//! For faster access, all relevant pointers to the 1D bases
    	FixedArray1D<IBASIS*,DIM> bases_;

    	//! Flag for deletion of the pointers
    	bool delete_pointers;
  	};
}

#include <cube/tbasis.cpp>

#endif /*_WAVELETTL_TBASIS_H_*/
