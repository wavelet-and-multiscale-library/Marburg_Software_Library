// implementation for i_index.h

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
    }
  }
  
  template <class IBASIS>
  IntervalIndex<IBASIS>::IntervalIndex(const IntervalIndex<IBASIS>& lambda)
  {
    j_ = lambda.j();
    e_ = lambda.e();
    k_ = lambda.k();
    basis_ = lambda.basis();
  }

  template <class IBASIS>
  IntervalIndex<IBASIS>::IntervalIndex(const int j, const int e, const int k,
				       const IBASIS* basis)
  {
    j_ = j;
    e_ = e;
    k_ = k;
    basis_ = basis;
  }

  template <class IBASIS>
  IntervalIndex<IBASIS>&
  IntervalIndex<IBASIS>::operator = (const IntervalIndex<IBASIS>& lambda)
  {
    j_ = lambda.j();
    e_ = lambda.e();
    k_ = lambda.k();
    basis_ = lambda.basis();

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
    switch (e_) {
    case 0:
      if (k_ == basis_->DeltaRmax(j_)) {
	e_ = 1;
	k_ = 0;
      }
      else
	k_++;
      break;
    case 1:
      if (k_ == (1<<j_)-1) {
	j_++;
	k_ = 0;
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
  inline
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
}
