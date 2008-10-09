// implementation for i_index.h

#include <cassert>

namespace WaveletTL
{
  template <class IBASIS>
  IntervalIndex<IBASIS>::IntervalIndex(const IBASIS* basis)
    : basis_(basis)
  {
    if (basis_ == 0) {
      j_ = e_ = k_ = 0; // invalid!
    } else {
      k_ = basis_->DeltaLmin(); // leftmost
      e_ = 0;                   // generator
      j_ = basis_->j0();        // on the coarsest level

//       num_ = -1;
    }
  }
  
  template <class IBASIS>
  IntervalIndex<IBASIS>::IntervalIndex(const IntervalIndex<IBASIS>& lambda)
  {
    j_ = lambda.j();
    e_ = lambda.e();
    k_ = lambda.k();
    basis_ = lambda.basis();
    //num_ = lambda.number();
  }

  template <class IBASIS>
  IntervalIndex<IBASIS>::IntervalIndex(const IntervalIndex<IBASIS>* lambda)
  {
    j_ = lambda->j();
    e_ = lambda->e();
    k_ = lambda->k();
    basis_ = lambda->basis();
    //num_ = lambda->number();
  }



  template <class IBASIS>
  IntervalIndex<IBASIS>::IntervalIndex(const int j, const int e, const int k,
				       const IBASIS* basis)
  {
    j_ = j;
    e_ = e;
    k_ = k;
    basis_ = basis;
//     num_ = -1;
  }


  template <class IBASIS>
  IntervalIndex<IBASIS>::IntervalIndex(const int num,
				       const IBASIS* basis)
  :  basis_(basis)
  {
    // num_ = num; 
    int num_ = num; 
    
    // to be decreased successively
    int act_num = num_;
    int j = basis_->j0();
    int tmp2 = 0;

    // determine the level of the index
    while (true) {
      int tmp = basis_->Deltasize(j);
      if (tmp > num)
	break;
      tmp2 = tmp;
      j++;
    }
    if (j == basis_->j0()) {
      j_ = j;
      e_ = 0;
    }
    else {
      j_ = j-1;
      e_ = 1;
    }

    act_num -= tmp2;

    if (e_ == 0)
      k_ = basis_->DeltaLmin() + act_num;
    else
      k_ = basis_->Nablamin()  + act_num;
 
  }

  template <class IBASIS>
  IntervalIndex<IBASIS>&
  IntervalIndex<IBASIS>::operator = (const IntervalIndex<IBASIS>& lambda)
  {
    j_ = lambda.j();
    e_ = lambda.e();
    k_ = lambda.k();
    basis_ = lambda.basis();
//     num_ = lambda.number();

    return *this;
  }

  template <class IBASIS>
  bool
  IntervalIndex<IBASIS>::operator == (const IntervalIndex<IBASIS>& lambda) const
  {
    return (j_ == lambda.j() &&
	    e_ == lambda.e() &&
	    k_ == lambda.k());
  }

  template <class IBASIS>
  bool
  IntervalIndex<IBASIS>::operator < (const IntervalIndex<IBASIS>& lambda) const
  {
    return (j_ < lambda.j() ||
	    (j_ == lambda.j() && (e_ < lambda.e() ||
				  (e_ == lambda.e() && k_ < lambda.k()))));
  }
  
  template <class IBASIS>
  IntervalIndex<IBASIS>&
  IntervalIndex<IBASIS>::operator ++ ()
  {
//     if (num_ > -1)
//       num_++;
    
    switch (e_) {
    case 0:
      if (k_ == basis_->DeltaRmax(j_)) {
	e_ = 1;
	k_ = basis_->Nablamin();
      }
      else
	k_++;
      break;
    case 1:
      if (k_ == basis_->Nablamax(j_)) {
	j_++;
	k_ = basis_->Nablamin();
      }
      else
	k_++;
      break;
    default:
      break;
    }
    
    return *this;
  }
  
  template <class IBASIS>
  void
  IntervalIndex<IBASIS>::set_number()
  {
    assert( e_ != 0 || j_ == basis_->j0() );

    int result = 0;

    // determine how many wavelets there are on all the levels
    // below the level of this index
    if (e_ == 1)
      result = basis_->Deltasize(j_);

    // count the wavelets that belong to the same level,
    // whose indices are smaller than this index,
    if (e_ == 0)
      result += k_ - basis_->DeltaLmin();
    else
      result += k_ - basis_->Nablamin();
//     num_ = result;
  }


  template <class IBASIS>
  inline
  IntervalIndex<IBASIS> first_generator(const IBASIS* basis, const int j)
  {
    assert(j >= basis->j0());
    return IntervalIndex<IBASIS>(j, 0, basis->DeltaLmin(), basis);
  }
  
  template <class IBASIS>
  inline
  IntervalIndex<IBASIS> last_generator(const IBASIS* basis, const int j)
  {
    assert(j >= basis->j0());
    return IntervalIndex<IBASIS>(j, 0, basis->DeltaRmax(j), basis);
  }

  template <class IBASIS>
  inline
  IntervalIndex<IBASIS> first_wavelet(const IBASIS* basis, const int j)
  {
    assert(j >= basis->j0());
    return IntervalIndex<IBASIS>(j, 1, basis->Nablamin(), basis);
  }
  
  template <class IBASIS>
  IntervalIndex<IBASIS> last_wavelet(const IBASIS* basis, const int j)
  {
    assert(j >= basis->j0());
    return IntervalIndex<IBASIS>(j, 1, basis->Nablamax(j), basis);
  }

  template <class IBASIS>
  inline
  IntervalIndex<IBASIS> first_index(const IBASIS* basis, const int j, const int e)
  {
    return (e == 0 ? first_generator(basis, j) : first_wavelet(basis, j));
  }
  
  template <class IBASIS>
  inline
  IntervalIndex<IBASIS> last_index(const IBASIS* basis, const int j, const int e)
  {
    return (e == 0 ? last_generator(basis, j) : last_wavelet(basis, j));
  }

  template <class IBASIS>
  inline
  int first_generator_num(const IBASIS* basis)
  {
    IntervalIndex<IBASIS> ind(first_generator<IBASIS>(basis, basis->j0()));
    ind.set_number();
    return ind.number();
  }
  
  template <class IBASIS>
  inline
  int last_generator_num(const IBASIS* basis)
  {
    IntervalIndex<IBASIS> ind(last_generator<IBASIS>(basis, basis->j0()));
    ind.set_number();
    return ind.number();

  }

  template <class IBASIS>
  inline
  int first_wavelet_num(const IBASIS* basis, const int j)
  {
    assert(j >= basis->j0());
    IntervalIndex<IBASIS> ind(first_wavelet<IBASIS>(basis, j));
    ind.set_number();
    return ind.number();

  }
  
  template <class IBASIS>
  inline
  int last_wavelet_num(const IBASIS* basis, const int j)
  {
    assert(j >= basis->j0());
    IntervalIndex<IBASIS> ind(last_wavelet<IBASIS>(basis, j));
    ind.set_number();
    return ind.number();
  }




  //
  //
  // from here on new version of IntervalIndex, without references to an instance of the basis

  template <class IBASIS>
  IntervalIndex2<IBASIS>::IntervalIndex2()
  {
    k_ = IBASIS::DeltaLmin(); // leftmost
    e_ = 0;                   // generator
    j_ = IBASIS::j0();        // on the coarsest level
  }

  template <class IBASIS>
  IntervalIndex2<IBASIS>::IntervalIndex2(const int j, const int e, const int k)
    : RIndex(j, e, k)
  {
  }

  template <class IBASIS>
  IntervalIndex2<IBASIS>::IntervalIndex2(const IntervalIndex2<IBASIS>* lambda)
    : RIndex(*lambda)
  {
  }

  template <class IBASIS>
  IntervalIndex2<IBASIS>::IntervalIndex2(const RIndex& lambda)
    : RIndex(lambda)
  {
  }
  
  template <class IBASIS>
  inline
  IntervalIndex2<IBASIS>&
  IntervalIndex2<IBASIS>::operator = (const IntervalIndex2<IBASIS>& lambda)
  {
    RIndex::operator = (lambda);
    return *this;
  }
  
  template <class IBASIS>
  bool
  IntervalIndex2<IBASIS>::operator == (const IntervalIndex2<IBASIS>& lambda) const
  {
    return (j_ == lambda.j() &&
	    e_ == lambda.e() &&
	    k_ == lambda.k());
  }
  
  template <class IBASIS>
  bool
  IntervalIndex2<IBASIS>::operator < (const IntervalIndex2<IBASIS>& lambda) const
  {
    return (j_ < lambda.j() ||
	    (j_ == lambda.j() && (e_ < lambda.e() ||
				  (e_ == lambda.e() && k_ < lambda.k()))));
  }
  
  template <class IBASIS>
  IntervalIndex2<IBASIS>&
  IntervalIndex2<IBASIS>::operator ++ ()
  {
    switch (e_) {
    case 0:
      if (k_ == IBASIS::DeltaRmax(j_)) {
	e_ = 1;
	k_ = IBASIS::Nablamin();
      }
      else
	k_++;
      break;
    case 1:
      if (k_ == IBASIS::Nablamax(j_)) {
	j_++;
	k_ = IBASIS::Nablamin();
      }
      else
	k_++;
      break;
    default:
      break;
    }
    
    return *this;
  }
  
  template <class IBASIS>
  inline
  IntervalIndex2<IBASIS> first_generator(const int j)
  {
    assert(j >= IBASIS::j0());
    return IntervalIndex2<IBASIS>(j, 0, IBASIS::DeltaLmin());
  }
  
  template <class IBASIS>
  inline
  IntervalIndex2<IBASIS> last_generator(const int j)
  {
    assert(j >= IBASIS::j0());
    return IntervalIndex2<IBASIS>(j, 0, IBASIS::DeltaRmax(j));
  }

  template <class IBASIS>
  inline
  IntervalIndex2<IBASIS> first_wavelet(const int j)
  {
    assert(j >= IBASIS::j0());
    return IntervalIndex2<IBASIS>(j, 1, IBASIS::Nablamin());
  }
  
  template <class IBASIS>
  inline
  IntervalIndex2<IBASIS> last_wavelet(const int j)
  {
    assert(j >= IBASIS::j0());
    return IntervalIndex2<IBASIS>(j, 1, IBASIS::Nablamax(j));
  }

}
