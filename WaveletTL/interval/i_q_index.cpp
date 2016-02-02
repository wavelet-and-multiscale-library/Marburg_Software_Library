// implementation for i_q_index.h

#include <cassert>
#include <cmath>

namespace WaveletTL
{
  template <class IFRAME>
  IntervalQIndex<IFRAME>::IntervalQIndex(const IFRAME* frame)
    : frame_(frame)
  {
    if (frame_ == 0) {
      p_ = j_ = e_ = k_ = 0; // invalid!
    } else {
      p_ = 0;                   // leftmost
      k_ = frame_->DeltaLmin(); // generator
      e_ = 0;                   // on the coarsest level
      j_ = frame_->j0();        // with polynomial degree 0.

//       num_ = -1;
    }
  }
  
  template <class IFRAME>
  IntervalQIndex<IFRAME>::IntervalQIndex(const IntervalQIndex<IFRAME>& lambda)
  {
    p_ = lambda.p();
    j_ = lambda.j();
    e_ = lambda.e();
    k_ = lambda.k();
    frame_ = lambda.frame();
    //num_ = lambda.number();
  }

  template <class IFRAME>
  IntervalQIndex<IFRAME>::IntervalQIndex(const IntervalQIndex<IFRAME>* lambda)
  {
    p_ = lambda->p();
    j_ = lambda->j();
    e_ = lambda->e();
    k_ = lambda->k();
    frame_ = lambda->frame();
    //num_ = lambda->number();
  }



  template <class IFRAME>
  IntervalQIndex<IFRAME>::IntervalQIndex(const int p, const int j, const int e, const int k, const IFRAME* frame)
  {
    p_ = p;
    j_ = j;
    e_ = e;
    k_ = k;
    frame_ = frame;
//     num_ = -1;
  }


  template <class IFRAME>
  IntervalQIndex<IFRAME>::IntervalQIndex(const int num,
				       const IFRAME* frame)
  :  frame_(frame)
  {
      if(true)
      {
          int j0 = frame_->j0();
          int nabla = frame_->Nablasize(j0);
          int tmp;
          int delta = frame_->Deltasize(j0);
          // using that delta(j0) + sum_{l=0}^{j} 2^l *nabla = Deltasize(j+j0) 
          if (num < delta) 
          {
	      p_ = 0;
              j_ = j0;
              e_ = 0;
          }
          else 
          {
              tmp = num - delta;
              tmp = tmp/nabla +1; // (double)tmp/nabla +1.;
	      p_ = 0;
              j_ = log2((unsigned int) tmp) + j0; // floor(log2(tmp)) + j0;
              e_ = 1;
          }
          if (e_ == 0)
              k_ = frame_->DeltaLmin() + num;
          else
              k_ = frame_->Nablamin() +num - delta + nabla - (1<<j_); // frame_->Nablamin()  + num -delta +nabla -(1<<(j_-j0))*nabla ;
      }
    else
    {
    // num_ = num; 
    int num_ = num; 
    
    // to be decreased successively
    int act_num = num_;
    int j = frame_->j0();
    int tmp2 = 0;

    // determine the level of the index
    while (true) {
      int tmp = frame_->Deltasize(j);
      if (tmp > num)
	break;
      tmp2 = tmp;
      j++;
    }
    if (j == frame_->j0()) {
      p_ = 0;
      j_ = j;
      e_ = 0;
    }
    else {
      p_ = 0;
      j_ = j-1;
      e_ = 1;
    }

    act_num -= tmp2;

    if (e_ == 0)
      k_ = frame_->DeltaLmin() + act_num;
    else
      k_ = frame_->Nablamin()  + act_num;
    }
  }



  template <class IFRAME>
  IntervalQIndex<IFRAME>&
  IntervalQIndex<IFRAME>::operator = (const IntervalQIndex<IFRAME>& lambda)
  {
    p_ = lambda.p();
    j_ = lambda.j();
    e_ = lambda.e();
    k_ = lambda.k();
    frame_ = lambda.frame();
//     num_ = lambda.number();

    return *this;
  }

  template <class IFRAME>
  bool
  IntervalQIndex<IFRAME>::operator == (const IntervalQIndex<IFRAME>& lambda) const
  {
    return (p_ == lambda.p() &&
	    j_ == lambda.j() &&
	    e_ == lambda.e() &&
	    k_ == lambda.k());
  }

  template <class IFRAME>
  bool
  IntervalQIndex<IFRAME>::operator < (const IntervalQIndex<IFRAME>& lambda) const
  {
    return (p_ < lambda.p() ||
	    p_ == lambda.p() && j_ < lambda.j() ||
	    p_ == lambda.p() && j_ == lambda.j() && e_ < lambda.e() ||
	    p_ == lambda.p() && j_ == lambda.j() && e_ == lambda.e() && k_ < lambda.k());
  }
  
  template <class IFRAME>
  IntervalQIndex<IFRAME>&
  IntervalQIndex<IFRAME>::operator ++ ()
  {
//     if (num_ > -1)
//       num_++;
    
    switch (e_) {
    case 0:
      if (k_ == frame_->DeltaRmax(j_)) {
	e_ = 1;
	k_ = frame_->Nablamin();
      }
      else
	k_++;
      break;
    case 1:
      if (k_ == frame_->Nablamax(j_)) {
	j_++;
	k_ = frame_->Nablamin();
      }
      else
	k_++;
      break;
    default:
      break;
    }
    
    return *this;
  }
  
  template <class IFRAME>
  void
  IntervalQIndex<IFRAME>::set_number()
  {
    assert( e_ != 0 || j_ == frame_->j0() );

    int result = 0;

    // determine how many wavelets there are on all the levels
    // below the level of this index
    if (e_ == 1)
      result = frame_->Deltasize(j_);

    // count the wavelets that belong to the same level,
    // whose indices are smaller than this index,
    if (e_ == 0)
      result += k_ - frame_->DeltaLmin();
    else
      result += k_ - frame_->Nablamin();
//     num_ = result;
  }


  template <class IFRAME>
  inline
  IntervalQIndex<IFRAME> first_q_generator(const IFRAME* frame, const int j, const int p)
  {
    assert(j >= frame->j0());
    return IntervalQIndex<IFRAME>(p, j, 0, frame->DeltaLmin(), frame);
  }
  
  template <class IFRAME>
  inline
  IntervalQIndex<IFRAME> last_q_generator(const IFRAME* frame, const int j, const int p)
  {
    assert(j >= frame->j0());
    return IntervalQIndex<IFRAME>(p, j, 0, frame->DeltaRmax(j), frame);
  }

  template <class IFRAME>
  inline
  IntervalQIndex<IFRAME> first_quarklet(const IFRAME* frame, const int j, const int p)
  {
    assert(j >= frame->j0());
    return IntervalQIndex<IFRAME>(p, j, 1, frame->Nablamin(), frame);
  }
  
  template <class IFRAME>
  IntervalQIndex<IFRAME> last_quarklet(const IFRAME* frame, const int j, const int p)
  {
    assert(j >= frame->j0());
    return IntervalQIndex<IFRAME>(p, j, 1, frame->Nablamax(j), frame);
  }

  template <class IFRAME>
  inline
  IntervalQIndex<IFRAME> first_q_index(const IFRAME* frame, const int j, const int e)
  {
    return (e == 0 ? first_q_generator(frame, j) : first_quarklet(frame, j));
  }
  
  template <class IFRAME>
  inline
  IntervalQIndex<IFRAME> last_q_index(const IFRAME* frame, const int j, const int e)
  {
    return (e == 0 ? last_q_generator(frame, j) : last_quarklet(frame, j));
  }

  template <class IFRAME>
  inline
  int first_q_generator_numb(const IFRAME* frame)
  {
    IntervalQIndex<IFRAME> ind(first_q_generator<IFRAME>(frame, frame->j0()));
    ind.set_number();
    return ind.number();
    }
  
  template <class IFRAME>
  inline
  int last_q_generator_numb(const IFRAME* frame)
  {
    IntervalQIndex<IFRAME> ind(last_q_generator<IFRAME>(frame, frame->j0()));
    ind.set_number();
    return ind.number();

  }

  template <class IFRAME>
  inline
  int first_quarklet_numb(const IFRAME* frame, const int j)
  {
    assert(j >= frame->j0());
    IntervalQIndex<IFRAME> ind(first_quarklet<IFRAME>(frame, j));
    ind.set_number();
    return ind.number();

  }
  
  template <class IFRAME>
  inline
  int last_quarklet_numb(const IFRAME* frame, const int j)
  {
    assert(j >= frame->j0());
    IntervalQIndex<IFRAME> ind(last_quarklet<IFRAME>(frame, j));
    ind.set_number();
    return ind.number();
    }




  //
  //
  // from here on new version of IntervalQIndex, without references to an instance of the basis

  template <class IFRAME>
  IntervalQIndex2<IFRAME>::IntervalQIndex2()
  {
    p_ = 0;
    k_ = IFRAME::DeltaLmin(); // leftmost
    e_ = 0;                   // generator
    j_ = IFRAME::j0();        // on the coarsest level with polynomial degree 0
  }

  template <class IFRAME>
  IntervalQIndex2<IFRAME>::IntervalQIndex2(const int p, const int j, const int e, const int k)
    : RQIndex(p, j, e, k)
  {
  }

  template <class IFRAME>
  IntervalQIndex2<IFRAME>::IntervalQIndex2(const IntervalQIndex2<IFRAME>* lambda)
    : RQIndex(*lambda)
  {
  }

  template <class IFRAME>
  IntervalQIndex2<IFRAME>::IntervalQIndex2(const RQIndex& lambda)
    : RQIndex(lambda)
  {
  }
  
  template <class IFRAME>
  inline
  IntervalQIndex2<IFRAME>&
  IntervalQIndex2<IFRAME>::operator = (const IntervalQIndex2<IFRAME>& lambda)
  {
    RQIndex::operator = (lambda);
    return *this;
  }
  
  template <class IFRAME>
  bool
  IntervalQIndex2<IFRAME>::operator == (const IntervalQIndex2<IFRAME>& lambda) const
  {
    return (p_ == lambda.p() &&
	    j_ == lambda.j() &&
	    e_ == lambda.e() &&
	    k_ == lambda.k());
  }
  
  template <class IFRAME>
  bool
  IntervalQIndex2<IFRAME>::operator < (const IntervalQIndex2<IFRAME>& lambda) const
  {
    return (p_ < lambda.p() ||
	    p_ == lambda.p() && j_ < lambda.j() ||
	    p_ == lambda.p() && j_ == lambda.j() && e_ < lambda.e() ||
	    p_ == lambda.p() && j_ == lambda.j() && e_ == lambda.e() && k_ < lambda.k());
  }
  
  template <class IFRAME>
  IntervalQIndex2<IFRAME>&
  IntervalQIndex2<IFRAME>::operator ++ ()
  {
    switch (e_) {
    case 0:
      if (k_ == IFRAME::DeltaRmax(j_)) {
	e_ = 1;
	k_ = IFRAME::Nablamin();
      }
      else
	k_++;
      break;
    case 1:
      if (k_ == IFRAME::Nablamax(j_)) {
	j_++;
        k_ = IFRAME::Nablamin();
      }
      else
	k_++;
      break;
    default:
      break;
    }
    
    
      
    
    return *this;
  }
  
  template <class IFRAME>
  inline
  IntervalQIndex2<IFRAME> first_q_generator(const int j, const int p)
  {
    assert(j >= IFRAME::j0());
    return IntervalQIndex2<IFRAME>(p, j, 0, IFRAME::DeltaLmin());
  }
  
  template <class IFRAME>
  inline
  IntervalQIndex2<IFRAME> last_q_generator(const int j, const int p)
  {
    assert(j >= IFRAME::j0());
    return IntervalQIndex2<IFRAME>(p, j, 0, IFRAME::DeltaRmax(j));
  }

  template <class IFRAME>
  inline
  IntervalQIndex2<IFRAME> first_quarklet(const int j, const int p)
  {
    assert(j >= IFRAME::j0());
    return IntervalQIndex2<IFRAME>(p, j, 1, IFRAME::Nablamin());
  }
  
  template <class IFRAME>
  inline
  IntervalQIndex2<IFRAME> last_quarklet(const int j, const int p)
  {
    assert(j >= IFRAME::j0());
    return IntervalQIndex2<IFRAME>(p, j, 1, IFRAME::Nablamax(j));
  }

}
