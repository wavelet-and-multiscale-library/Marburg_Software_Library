// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2018                                            |
// | Philipp Keding, Alexander Sieber                                    |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_SLITDOMAIN_FRAME_INDEX_H
#define _WAVELETTL_SLITDOMAIN_FRAME_INDEX_H

#include <iostream>
using std::cout;
using std::endl;

#include <utils/multiindex.h>

using MathTL::MultiIndex;

namespace WaveletTL
{
  template <class IFRAME> class SlitDomainFrame;
  
  /*!
    An index class for composite quarklet frames of the slitdomain in R^2.
    In fact, this is essentially a lexicographically ordered multiindex (j,e,p,k), with
    * scale j
    * type e ((0,0),(0,1),(1,0) or (1,1))
    * (logical) patch p (0<=p<=4)
    * translation index k

    The logical patches are:
      p=0 (-1,0)x( 0,1)
      p=1 (-1,0)x(-1,0)
      p=2 ( 0,1)x(-1,0)
      p=3 (0,1)x (0,1)
      p=4 (-1,0) x {0}
      p=5 {0} x (-1,0)
      p=6 (0,1) x {0} 
   * 
   *   */
  template <class IFRAME>
  class SlitDomainFrameIndex
  {
  public:
    
    
     //! type polynomial type
    typedef MultiIndex<int,2> polynomial_type; 
      
    //! type level type
    typedef MultiIndex<int,2> level_type;
      
    //! type index type
    typedef MultiIndex<int,2> type_type;
    
    //! translation index type
    typedef MultiIndex<int,2> translation_type;

    /*!
      constructor with a given slitdomain frame
      (also serves as a default constructor, but yields an invalid index
      in this case, because the underlying frames must be specified to work correctly)
    */
    SlitDomainFrameIndex(const SlitDomainFrame<IFRAME>* frame = 0);

    //! constructor with given p,j,e,patch,k
    SlitDomainFrameIndex(const polynomial_type& p,
                 const level_type& j,
		 const type_type& e,
		 const int patch,
		 const translation_type& k,
                 const unsigned int number,
		 const SlitDomainFrame<IFRAME>* frame);

    //! copy constructor
    SlitDomainFrameIndex(const SlitDomainFrameIndex& lambda);

    //! copy index from const pointer
    SlitDomainFrameIndex(const SlitDomainFrameIndex* lambda);

  
    //! check equality
    bool operator == (const SlitDomainFrameIndex& lambda) const;
    
    //! check non-equality
    inline bool operator != (const SlitDomainFrameIndex& lambda) const
    { return !(*this == lambda); }
    
    //! preincrement
    SlitDomainFrameIndex& operator ++ ();

    //! lexicographic order <
    bool operator < (const SlitDomainFrameIndex& lambda) const;

    //! lexicographic order <=
    bool operator <= (const SlitDomainFrameIndex& lambda) const
    { return (*this < lambda || *this == lambda); }

    //! polynomial p
    const polynomial_type& p() const { return p_; }
    
    //! scale j
    const level_type& j() const { return j_; }

    //! type e
    const type_type& e() const { return e_; }

    //! patch
    const int patch() const { return patch_; }

    //! translation index k
    const translation_type& k() const { return k_; }

    /*!
      By construction, the overall quarklet index set is ordered, so that
      there exists a bijective mapping into the positive integers.
      This routine returns the "number" of the current index, starting with 0
      for the first generator on the coarsest level. If the current index is a
      generator on a higher level j, the number 0 corresponds to the first generator
      on the level j.
    */
//    const int number() const;
    
    const unsigned int& number() const { return num_; }
    
  protected:
    //! pointer to the underlying frame
    const SlitDomainFrame<IFRAME>* frame_;

    //! polynomial
    polynomial_type p_; 
    
    //! scale
    level_type j_;
    
    //! type
    type_type e_;

    //! patch
    int patch_;
    
    //! translation
    translation_type k_;
    
    // Number of the index, only for the elements of a quarklet frame
    unsigned int num_;

  };

  //! stream output
  template <class IFRAME>
  inline std::ostream& operator << (std::ostream& os,
				    const SlitDomainFrameIndex<IFRAME>& lambda)
  {
    using namespace std;
    os << "("
       << lambda.p()
       << ","
       << lambda.j()
       << ","
       << lambda.e()
       << ","
       << lambda.patch()
       << ","
       << lambda.k()
       << ")" << " number = " << lambda.number();
    return os;
    
    
  }

  /*!
    index of first generator on level j >= j0
  */
  template <class IFRAME>
  SlitDomainFrameIndex<IFRAME>
  first_generator(const SlitDomainFrame<IFRAME>* frame, const typename SlitDomainFrameIndex<IFRAME>::level_type& j, const typename SlitDomainFrameIndex<IFRAME>::polynomial_type& p);

  /*!
    index of last generator on level j >= j0
  */
  template <class IFRAME>
  SlitDomainFrameIndex<IFRAME>
  last_generator(const SlitDomainFrame<IFRAME>* frame, const typename SlitDomainFrameIndex<IFRAME>::level_type& j, const typename SlitDomainFrameIndex<IFRAME>::polynomial_type& p);

  /*!
    index of first quarklet on level j >= j0
  */
  template <class IFRAME>
  SlitDomainFrameIndex<IFRAME>
  first_quarklet(const SlitDomainFrame<IFRAME>* frame, const typename SlitDomainFrameIndex<IFRAME>::level_type& j, const typename SlitDomainFrameIndex<IFRAME>::polynomial_type& p);
  
  /*!
    index of first quarklet with a given type on level j >= j0
  */
//  template <class IFRAME>
//  SlitDomainFrameIndex<IFRAME>
//  first_quarklet(const SlitDomainFrame<IFRAME>* frame,
//		const typename SlitDomainFrameIndex<IFRAME>::level_type& j,
//		const typename SlitDomainFrameIndex<IFRAME>::type_type& e, const typename SlitDomainFrameIndex<IFRAME>::polynomial_type& p);
  
  /*!
    index of last quarklet on level j >= j0
  */
  template <class IFRAME>
  SlitDomainFrameIndex<IFRAME>
  last_quarklet(const SlitDomainFrame<IFRAME>* frame, const typename SlitDomainFrameIndex<IFRAME>::level_type& j, const typename SlitDomainFrameIndex<IFRAME>::polynomial_type& p);
}

#include <slitdomain/slitdomain_frame_index.cpp>

#endif
