// implementation for i_multi_index.h

namespace WaveletTL
{
  /* constructors **********************************************************/
  template <class IBASIS>
  IntervalMultiIndex<IBASIS>::IntervalMultiIndex(const IBASIS* basis)
    : basis_(basis)
  {
    if (basis_ == 0) {
      j_ = e_ = k_ = 0; // invalid!
    } else {
      k_ = basis_->DeltaLmin(); // leftmost
      c_ = 0;                   // first wavelet
      e_ = 0;                   // generator
      j_ = basis_->j0();        // on the coarsest level
    }
  }

  template <class IBASIS>
  IntervalMultiIndex<IBASIS>::IntervalMultiIndex(const IntervalMultiIndex<IBASIS>& lambda)
  {
    operator = (lambda);
  }

  template <class IBASIS>
  IntervalMultiIndex<IBASIS>::IntervalMultiIndex(const int j, const type_type e, const component_type c, const translation_type k, const IBASIS* basis)
    : j_(j), e_(e), c_(c), k_(k), basis_(basis)
  {
  }
  
  /* member functions ******************************************************/
  template <class IBASIS>
  IntervalMultiIndex<IBASIS>&
  IntervalMultiIndex<IBASIS>::operator = (const IntervalMultiIndex<IBASIS>& lambda)
  {
    j_ = lambda.j();
    e_ = lambda.e();
    c_ = lambda.c();
    k_ = lambda.k();
    basis_ = lambda.basis();

    return *this;
  }

  template <class IBASIS>
  bool
  IntervalMultiIndex<IBASIS>::operator == (const IntervalMultiIndex<IBASIS>& lambda) const
  {
    return (j_ == lambda.j() &&
	    e_ == lambda.e() &&
	    c_ == lambda.c() &&
	    k_ == lambda.k());
  }

  template <class IBASIS>
  bool
  IntervalMultiIndex<IBASIS>::operator < (const IntervalMultiIndex<IBASIS>& lambda) const
  {
    return (j_ < lambda.j() ||
	    (j_ == lambda.j() && (e_ < lambda.e() ||
				  (e_ == lambda.e() && (c_ < lambda.c() ||
							(c_ == lambda.c() && k_< lambda.k()))))));
  }
  
  template <class IBASIS>
  IntervalMultiIndex<IBASIS>&
  IntervalMultiIndex<IBASIS>::operator ++ ()
  {
    switch (e_) {
    case 0: // generator
      if (c_ < basis_->number_of_components-1) {
        if (k_ == basis_->DeltaRmax(j_)) {
          c_++;
          k_ = basis_->DeltaLmin();
        }
        else
          k_++;
      }
      else { // c_ == basis_->number_of_components-1
        if (k_ == basis_->DeltaRmax(j_)) {
          e_ = 1;
          c_ = 0;
          k_ = basis_->Nablamin();
        }
        else
          k_++;
      }
      break;
    case 1: // wavelet
      if (c_ < basis_->number_of_components-1) {
        if (k_ == basis_->Nablamax(j_)) {
          c_++;
          k_ = basis_->Nablamin();
        }
        else
          k_++;
      }
      else { // c_ == basis_->number_of_components-1
        if (k_ == basis_->Nablamax(j_)) {
          j_++;
          c_ = 0;
          k_ = basis_->Nablamin();
        }
        else
          k_++;
      }
      break;
    default:
      break;
    }
    
    return *this;
  }
  

  /* non member functions **************************************************/
  template <class IBASIS>
  inline
  IntervalMultiIndex<IBASIS>
  first_generator(const IBASIS* basis, const int j)
  {
    assert(j >= basis->j0());
    return IntervalMultiIndex<IBASIS>(j, 0, basis->DeltaLmin(), basis);
  }
  
  template <class IBASIS>
  inline
  IntervalMultiIndex<IBASIS>
  last_generator(const IBASIS* basis, const int j)
  {
    assert(j >= basis->j0());
    return IntervalMultiIndex<IBASIS>(j, 0, basis->DeltaRmax(j), basis);
  }

  template <class IBASIS>
  inline
  IntervalMultiIndex<IBASIS>
  first_wavelet(const IBASIS* basis, const int j)
  {
    assert(j >= basis->j0());
    return IntervalMultiIndex<IBASIS>(j, 1, basis->Nablamin(), basis);
  }
  
  template <class IBASIS>
  inline
  IntervalMultiIndex<IBASIS>
  last_wavelet(const IBASIS* basis, const int j)
  {
    assert(j >= basis->j0());
    return IntervalMultiIndex<IBASIS>(j, 1, basis->Nablamax(j), basis);
  }

  template <class IBASIS>
  inline
  IntervalMultiIndex<IBASIS>
  first_index(const IBASIS* basis, const int j, const typename IntervalMultiIndex<IBASIS>::type_type e)
  {
    return (e == 0 ? first_generator(basis, j) : first_wavelet(basis, j));
  }
  
  template <class IBASIS>
  inline
  IntervalMultiIndex<IBASIS>
  last_index(const IBASIS* basis, const int j, const typename IntervalMultiIndex<IBASIS>::type_type e)
  {
    return (e == 0 ? last_generator(basis, j) : last_wavelet(basis, j));
  }

}
