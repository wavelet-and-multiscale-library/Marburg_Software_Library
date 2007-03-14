// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2007                                            |
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

  /*!
    A (template) multi-wavelet index class for multiwavelet basis on [0,1]
    modelling a multiindex of the form lambda=(j,e,c,k).

    r is the number of multi-wavelets in the basis.

    We require the template parameter class IBASIS to provide the routines
    j0(), DeltaLmin(), DeltaRmax(j), Nablamin() and Nablamax(j).
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

    //! underlying basis
    inline const IBASIS* basis() const { return basis_; }

  protected:  
    
    //! scale, type, component, translation
    int j_;
    type_type e_;
    translation_type k_;
    component_type c_;
    
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
}

// include implementation
#include <interval/i_multi_index.cpp>

#endif  // _WAVELETTL_I_MULTI_INDEX_H
