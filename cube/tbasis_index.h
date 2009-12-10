// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner, Ulrich Friedrich                   |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_TBASIS_INDEX_H_
#define _WAVELETTL_TBASIS_INDEX_H_

#include <iostream>
using std::cout;
using std::endl;

#include "utils/multiindex.h"

#include <utils/multiindex.h>
#include <utils/fixed_array1d.h>

using MathTL::MultiIndex;

namespace WaveletTL
{
	template <class IBASIS, unsigned int DIM> class TensorBasis;
	/*
	 * An index class for tensor product wavelet bases over the d-dimensional unit cube [0,1]^d (or mapped versions thereof) of the type
	 * (V_0 \oplus W_0 \oplus W_1 \oplus W_2 \oplus ...)\times (V_0 \oplus W_0 \oplus W_1 \oplus W_2 \oplus ...)
	 * as modeled by TensorBasis
  	*/
  	template <class IBASIS, unsigned int DIM, class TENSORBASIS = TensorBasis<IBASIS,DIM> >
  	class TensorIndex
	{
		public:
		// level index type
		typedef MultiIndex<int,DIM> level_type;
		// type index type
		typedef MultiIndex<int,DIM> type_type;
		// translation index type
		typedef MultiIndex<int,DIM> translation_type;

		/*
		 * Constructor with a given tensor basis
		 * (also serves as a default constructor, but yields an invalid index
		 * in this case, because the underlying bases must be specified to work correctly)
    	*/
    	TensorIndex(const TENSORBASIS* basis = 0);

	    /*
	     * Constructor with given j,e,k
    	*/
    	TensorIndex(const level_type& j, const type_type& e, const translation_type& k, const TENSORBASIS* basis);

    	// Copy constructor
    	TensorIndex(const TensorIndex& lambda);

	    // Copy index from const pointer
    	TensorIndex(const TensorIndex* lambda);

	    /*
	     * Constructor.
	     * We always assume to have a lexicographical ordering for the elements in the wavelet bases,
	     * w.r.t. the tuple ( j , k ). According to this ordering there exists a mapping from the
	     * totality of wavelet indices into the nonnegative integers. This routine creates an index
	     * from the given number corresponding to this ordering.
	     * This constructor cannot be used for generators on levels j > j0 (multi index)
    	*/
    	TensorIndex(const int number, const TENSORBASIS* basis);

   		// Assignment
    	TensorIndex& operator = (const TensorIndex& lambda);

    	// Check equality
    	bool operator == (const TensorIndex& lambda) const;

    	// Check non-equality
    	inline bool operator != (const TensorIndex& lambda) const
    	{ return !(*this == lambda); }

    	// Preincrement
    	TensorIndex& operator ++ ();

    	/* Ordering <
         * First by level, i.e. the 1-norm of j,
         * then lexicographically w.r.t. j,e,k
         */
    	bool operator < (const TensorIndex& lambda) const;

    	// Lexicographic order <=
    	bool operator <= (const TensorIndex& lambda) const
    	{ return (*this < lambda || *this == lambda); }

    	// Scale j
    	const level_type& j() const { return j_; }

    	// Type e
    	const type_type& e() const { return e_; }

    	// Translation index k
    	const translation_type& k() const { return k_; }

    	// Underlying basis
    	const TENSORBASIS* basis() const { return basis_; }

    	const int number() const { return num_; }

	protected:

    	// Pointer to the underlying basis
    	const TENSORBASIS* basis_;

    	// Scale
    	MultiIndex<int,DIM> j_;

    	// Type
    	MultiIndex<int,DIM> e_;

    	// Translation
    	MultiIndex<int,DIM> k_;

    	// Number of the index, only for the elements of a wavelet bases
    	int num_;

	};

	//! stream output
	template <class IBASIS, unsigned int DIM, class TENSORBASIS>
	inline std::ostream& operator << (std::ostream& os, const TensorIndex<IBASIS,DIM,TENSORBASIS>& lambda)
	{
		using namespace std;
		os << "("
		<< lambda.j()
		<< ","
		<< lambda.e()
		<< ","
		<< lambda.k()
		<< ")";
	return os;
	}

	/*
	 * index of first generator
	 */
	template <class IBASIS, unsigned int DIM, class TENSORBASIS>
	TensorIndex<IBASIS,DIM,TENSORBASIS>
	first_generator(const TENSORBASIS* basis);

	/*
	 * index of last generator
	 */
	template <class IBASIS, unsigned int DIM, class TENSORBASIS>
	TensorIndex<IBASIS,DIM,TENSORBASIS>
	last_generator(const TENSORBASIS* basis);

	/*
     * index of first wavelet on level j >= j0
     */
	template <class IBASIS, unsigned int DIM, class TENSORBASIS>
	TensorIndex<IBASIS,DIM,TENSORBASIS>
	first_wavelet(const TENSORBASIS* basis, const typename TensorIndex<IBASIS,DIM,TENSORBASIS>::level_type j);

	/*
  	 * index of last wavelet on level j >= j0
	 */
	template <class IBASIS, unsigned int DIM, class TENSORBASIS>
	TensorIndex<IBASIS,DIM,TENSORBASIS>
	last_wavelet(const TENSORBASIS* basis, const typename TensorIndex<IBASIS,DIM,TENSORBASIS>::level_type j);

	/*
     * number of first generator on level j0
     */
	template <class IBASIS, unsigned int DIM, class TENSORBASIS>
	const int
	first_generator_num(const TENSORBASIS* basis);

	/*
	 * number of last generator on level j0
	 */
	template <class IBASIS, unsigned int DIM, class TENSORBASIS>
	const int
	last_generator_num(const TENSORBASIS* basis);

	/*
	 * number of first wavelet on level j >= j0
	 */
	template <class IBASIS, unsigned int DIM, class TENSORBASIS>
	const int
	first_wavelet_num(const TENSORBASIS* basis, const typename TensorIndex<IBASIS,DIM,TENSORBASIS>::level_type j);

	/*
	 * number of last wavelet on level j >= j0
	 */
	template <class IBASIS, unsigned int DIM, class TENSORBASIS>
	const int
	last_wavelet_num(const TENSORBASIS* basis, const typename TensorIndex<IBASIS,DIM,TENSORBASIS>::level_type j);
}

#include "tbasis_index.cpp"

#endif /*TBASIS_INDEX_H_*/

