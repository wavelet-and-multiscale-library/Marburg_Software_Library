// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner, Philipp Keding                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_I_Q_INDEX_H
#define _WAVELETTL_I_Q_INDEX_H

#include <iostream>
#include <Rd/r_q_index.h>
#include <interval/i_index.h>

using std::cout;
using std::endl;

namespace WaveletTL
{
  /*!
    A (template) quarklet index class for frames on [0,1] like periodic quarklet bases
    or those from [DKU],[DS] (see template class DSBasis).
    We require the template parameter class IFRAME to provide the routines
    j0(), DeltaLmin(), DeltaRmax(j), Nablamin() and Nablamax(j).
  */
  template <class IFRAME>
  class IntervalQIndex
  {
  public:
    /*!
      constructor with given interval frame
      (also serves as a default constructor, but yields an invalid index then
      because the underlying interval frame must be specified to work correctly)
    */
    IntervalQIndex(const IFRAME* frame = 0);

    //! copy constructor
    IntervalQIndex(const IntervalQIndex& lambda);

    //! copy index from const pointer
    IntervalQIndex(const IntervalQIndex* lambda);
  
    //! constructor with specified parameters
    IntervalQIndex(const int p, const int j, const int e, const int k, const IFRAME* frame);

    //! constructor with specified number
    IntervalQIndex(const int num, const IFRAME* frame);
    
    //! assignment
    IntervalQIndex& operator = (const IntervalQIndex& lambda);

    //! check equality
    bool operator == (const IntervalQIndex& lambda) const;

    //! check non-equality
    bool operator != (const IntervalQIndex& lambda) const
    { return !(*this == lambda); }

    //! preincrement
    IntervalQIndex& operator ++ ();
    
    //! lexicographic order <
    bool operator < (const IntervalQIndex& lambda) const;

    //! lexicographic order <=
    bool operator <= (const IntervalQIndex& lambda) const
    { return (*this < lambda || *this == lambda); }

    //! polynomial degree p
    const int p() const { return p_; }


    //! scale j
    const int j() const { return j_; }

    /*
     * for convenience: level_type (=int)
    */
    typedef int level_type;

    //! type index type
    typedef int type_type;

    //! type e
    const type_type& e() const { return e_; }

    //! translation index type
    typedef int translation_type;

    //! translation index k
    const translation_type& k() const { return k_; }

    //! number
    const int number() const {
//       return num_;
      return (e() == 0
  	      ? k() - frame_->DeltaLmin(p())
  	      : frame_->DeltaNablasize(j(),p()) + k() - frame_->Nablamin(p()));
    }

    /*!
      inverse of constructor
      'IntervalQIndex(const unsigned int number,
                     const IFRAME* frame)'
    */
    void set_number();


    //! underlying basis
    const IFRAME* frame() const { return frame_; }

  protected:  
    
    //! polynomial degree, scale, type, translation
    int p_, j_, e_, k_;
    
//     //! number of the index (not for generators on level j > j0)
//     int num_;
    

    //! pointer to the underlying interval frame
    const IFRAME* frame_;
  };

  //! stream output
  template <class IFRAME>
  std::ostream& operator << (std::ostream& os, const IntervalQIndex<IFRAME>& lambda)
  {
    using namespace std;
    os << "("
       << lambda.p()
       << ","
       << lambda.j()
       << ","
       << lambda.e()
       << ","
       << lambda.k()
       << ")" << " number = " << lambda.number();
    return os;
  }

  /*!
    index of first (leftmost) generator on level j >= j0
  */
  template <class IFRAME>
  IntervalQIndex<IFRAME> first_q_generator(const IFRAME* frame, const int j, const int p = 0);

  /*!
    index of last (rightmost) generator on level j >= j0
  */
  template <class IFRAME>
  IntervalQIndex<IFRAME> last_q_generator(const IFRAME* frame, const int j, const int p = 0);

  /*!
    index of first (leftmost) wavelet on level j >= j0
  */
  template <class IFRAME>
  IntervalQIndex<IFRAME> first_quarklet(const IFRAME* frame, const int j, const int p = 0);

  /*!
    index of last (rightmost) wavelet on level j >= j0
  */
  template <class IFRAME>
  IntervalQIndex<IFRAME> last_quarklet(const IFRAME* frame, const int j, const int p = 0);

  /*!
    index of first function with type e
    (mainly for TensorProductBasis)
   * from this point not changed to quarklet setting @PHK
  */
  template <class IFRAME>
  IntervalQIndex<IFRAME> first_q_index(const IFRAME* frame, const int j, const int e);

  /*!
    index of last function with type e
    (mainly for TensorProductBasis)
  */
  template <class IFRAME>
  IntervalQIndex<IFRAME> last_q_index(const IFRAME* frame, const int j, const int e);

  /*!
    number of first (leftmost) generator on level j0
  */
  template <class IFRAME>
  int first_q_generator_numb(const IFRAME* frame);

  /*!
    number of last (rightmost) generator on level j0
  */
  template <class IFRAME>
  int last_q_generator_numb(const IFRAME* frame);

  /*!
    number of first (leftmost) wavelet on level j >= j0
  */
  template <class IFRAME>
  int first_quarklet_numb(const IFRAME* frame, const int j);

  /*!
    number of last (rightmost) wavelet on level j >= j0
  */
  template <class IFRAME>
  int last_quarklet_numb(const IFRAME* frame, const int j);


  
  
  //
  //
  // from here on new version of IntervalIndex, without references to an instance of the basis


  template <class IFRAME>
  class IntervalQIndex2
    : public RQIndex
  {
  public:
    using RQIndex::p;
    using RQIndex::j;
    using RQIndex::type_type;
    using RQIndex::e;
    using RQIndex::translation_type;
    using RQIndex::k;

    //! default constructor
    IntervalQIndex2();

    //! constructor with specified parameters
    IntervalQIndex2(const int p, const int j, const int e, const int k);

    //! copy index from const pointer
    IntervalQIndex2(const IntervalQIndex2* lambda);
  
    //! constructor from an RQIndex
    IntervalQIndex2(const RQIndex& lambda);

    //! assignment
    IntervalQIndex2& operator = (const IntervalQIndex2& lambda);

    //! check equality
    bool operator == (const IntervalQIndex2& lambda) const;

    //! check non-equality
    bool operator != (const IntervalQIndex2& lambda) const
    { return !(*this == lambda); }

    //! preincrement
    IntervalQIndex2& operator ++ ();
    
    //! lexicographic order <
    bool operator < (const IntervalQIndex2& lambda) const;

    //! lexicographic order <=
    bool operator <= (const IntervalQIndex2& lambda) const
    { return (*this < lambda || *this == lambda); }
    
    //! number
    const int number() const {
     return (e() == 0
  	      ? k() - IFRAME::DeltaLmin()
  	      : IFRAME::Deltasize(j()) + k() - IFRAME::Nablamin());
     }
  };

  //! stream output
  template <class IFRAME>
  std::ostream& operator << (std::ostream& os, const IntervalQIndex2<IFRAME>& lambda)
  {
    using namespace std;
    os << "("
       << lambda.p()
       << ","
       << lambda.j()
       << ","
       << lambda.e()
       << ","
       << lambda.k()
       << ")";
    return os;
  }

  /*!
    index of first (leftmost) generator on level j >= j0, p>=0
  */
  template <class IFRAME>
  IntervalQIndex2<IFRAME> first_q_generator(const int j, const int p = 0);

  /*!
    index of last (rightmost) generator on level j >= j0, p>=0
  */
  template <class IFRAME>
  IntervalQIndex2<IFRAME> last_q_generator(const int j, const int p = 0);

  /*!
    index of first (leftmost) quarklet on level j >= j0, p>=0
  */
  template <class IFRAME>
  IntervalQIndex2<IFRAME> first_quarklet(const int j, const int p = 0);

  /*!
    index of last (rightmost) quarklet on level j >= j0, p>=0
  */
  template <class IFRAME>
  IntervalQIndex2<IFRAME> last_quarklet(const int j, const int p = 0);

}

// include implementation
#include <interval/i_q_index.cpp>

#endif
