// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2008                                            |
// | Thorsten Raasch, Manuel Werner, Andreas Schneider                  |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_I_MULTI_INDEX_H
#define _WAVELETTL_I_MULTI_INDEX_H

#include <iostream>
using std::cout;
using std::endl;

namespace WaveletTL
{
  // provide some constants for better readability
  #define E_GENERATOR 0  // index type generator
  #define E_WAVELET 1  // index type wavelet

  //#ifndef N_UNDEFINED
   #define N_UNDEFINED -1
  //#endif
  //#ifndef N_UNKNOWN
   #define N_UNKNOWN -2
  //#endif

  /*!
    A (template) multi-wavelet index class for multiwavelet basis on [0,1]
    modelling a multiindex of the form lambda=(j,e,c,k).

    r is the number of multi-wavelets in the basis.

    We require the template parameter class IBASIS to provide the routines
    j0(), DeltaLmin(), DeltaRmax(j), Nablamin(), Nablamax(j) and the data
    field number_of_components.
  */
  template <class IBASIS>
  class IntervalMultiIndex
  {
  public:
    /* Typedefs ************************************************************/
    //! type index type
    typedef int type_type;

    //! translation index type
    typedef int translation_type;

    //! component type
    typedef unsigned int component_type;

    //! number type
    typedef int number_type;


    /* constructors ********************************************************/
    /*!
      constructor with given interval basis
      (also serves as a default constructor, but yields an invalid index then
      because the underlying interval basis must be specified to work correctly)
    */
    IntervalMultiIndex(const IBASIS* basis = 0);

    //! copy constructor
    IntervalMultiIndex(const IntervalMultiIndex<IBASIS>& lambda);

    //! constructor with specified parameters
    IntervalMultiIndex(const int j, const type_type e, const translation_type k, const component_type c, const IBASIS* basis);

    //! constructor with specified number
    IntervalMultiIndex(const number_type n, const IBASIS* basis);


    /* member functions ****************************************************/
    //! check if it is a valid index
    bool is_valid() const;

    //! assignment
    IntervalMultiIndex<IBASIS>& operator = (const IntervalMultiIndex<IBASIS>& lambda);

    //! check equality
    bool operator == (const IntervalMultiIndex<IBASIS>& lambda) const;

    //! check non-equality
    bool operator != (const IntervalMultiIndex<IBASIS>& lambda) const
    { return !(*this == lambda); }
    
    //! preincrement
    IntervalMultiIndex<IBASIS>& operator ++ ();
    
    //! predecrement
    IntervalMultiIndex<IBASIS>& operator -- ();
    
    //! lexicographic order <
    bool operator < (const IntervalMultiIndex<IBASIS>& lambda) const;

    //! lexicographic order <=
    bool operator <= (const IntervalMultiIndex<IBASIS>& lambda) const
    { return (*this < lambda || *this == lambda); }
    
    //! scale j
    inline const int j() const { return j_; }
    
    //! type e
    inline const type_type& e() const { return e_; }
    
    //! translation index k
    inline const translation_type& k() const { return k_; }

    //! component index c
    inline const component_type& c() const { return c_; }

    //! number of the multi-index
    number_type number();

    //! number of the multi-index, if already computed
    const number_type& number_if_computed() const { return num_; }

    //! underlying basis
    inline const IBASIS* basis() const { return basis_; }

    //! set the values of j, e, k, c directly
    inline void set(int j, type_type e, translation_type k, component_type c);

    /*! encode c and k index of this object into a single integer
        (for usage with an adaptor class which wraps a multiwavelet basis in a classical wavelet basis)
    */
    inline translation_type ck_encode() const { return ck_encode(k_, c_); }

    /*! encode given c and k index into a single integer
        (for usage with an adaptor class which wraps a multiwavelet basis in a classical wavelet basis)
    */
    static translation_type ck_encode(const translation_type& k, const component_type& c);

    /*! decode a single (encoded) integer k into translation index k and component index c
        (for usage with an adaptor class which wraps a multiwavelet basis in a classical wavelet basis)
    */
    static void ck_decode(const translation_type& k_enc, translation_type& k, component_type& c);

  protected:  
    
    //! scale, type, component, translation
    int j_;
    type_type e_;
    translation_type k_;
    component_type c_;

    //! number of the multi-index
    number_type num_;
    
    //! pointer to the underlying interval basis
    const IBASIS* basis_;
  };

  //! stream output
  template <class IBASIS>
  inline std::ostream& operator << (std::ostream& os, const IntervalMultiIndex<IBASIS>& lambda)
  {
    using namespace std;
    os << "("
       << lambda.j()
       << ","
       << lambda.e()
       << ","
       << lambda.k()
       << ","
       << lambda.c()
       << ")";
    return os;
  }

#if 0
  /*!
    index of first (leftmost) generator on level j >= j0
  */
  template <class IBASIS>
  inline
  IntervalMultiIndex<IBASIS> first_generator(const IBASIS* basis, const int j);

  /*!
    index of last (rightmost) generator on level j >= j0
  */
  template <class IBASIS>
  inline
  IntervalMultiIndex<IBASIS> last_generator(const IBASIS* basis, const int j);

  /*!
    index of first (leftmost) wavelet on level j >= j0
  */
  template <class IBASIS>
  inline
  IntervalMultiIndex<IBASIS> first_wavelet(const IBASIS* basis, const int j);

  /*!
    index of last (rightmost) wavelet on level j >= j0
  */
  template <class IBASIS>
  inline
  IntervalMultiIndex<IBASIS> last_wavelet(const IBASIS* basis, const int j);

  /*!
    index of first function with type e
    (mainly for TensorProductBasis)
  */
  template <class IBASIS>
  inline
  IntervalMultiIndex<IBASIS> first_index(const IBASIS* basis, const int j, const typename IntervalMultiIndex<IBASIS>::type_type e);

  /*!
    index of last function with type e
    (mainly for TensorProductBasis)
  */
  template <class IBASIS>
  inline
  IntervalMultiIndex<IBASIS> last_index(const IBASIS* basis, const int j, const typename IntervalMultiIndex<IBASIS>::type_type e);
#endif
}

// include implementation
#include <interval/i_multi_index.cpp>

#endif  // _WAVELETTL_I_MULTI_INDEX_H
