// implementation for i_multi_index.h

#include <interval/i_multi_index.h>

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
      c_ = 0;                   // first wavelet
      k_ = basis_->DeltaLmin(); // leftmost
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
  IntervalMultiIndex<IBASIS>::IntervalMultiIndex(const int j, const type_type e, const translation_type k, const component_type c, const IBASIS* basis)
    : j_(j), e_(e), k_(k), c_(c), basis_(basis)
  {
  }
  
  /* member functions ******************************************************/
  template <class IBASIS>
  bool
  IntervalMultiIndex<IBASIS>::is_valid()
  {
    return ( (j_ >= basis_->j0()) && ((e_ == E_GENERATOR) || (e_ == E_WAVELET))
             && (k_ >= basis_->DeltaLmin()) && (k_ <= basis_->DeltaRmax(j_))
             && (c_ >= 0) && (c_ < basis_->number_of_components) );
  }

  template <class IBASIS>
  IntervalMultiIndex<IBASIS>&
  IntervalMultiIndex<IBASIS>::operator = (const IntervalMultiIndex<IBASIS>& lambda)
  {
    j_ = lambda.j();
    e_ = lambda.e();
    k_ = lambda.k();
    c_ = lambda.c();
    basis_ = lambda.basis();

    return *this;
  }

  template <class IBASIS>
  bool
  IntervalMultiIndex<IBASIS>::operator == (const IntervalMultiIndex<IBASIS>& lambda) const
  {
    return (j_ == lambda.j() &&
            e_ == lambda.e() &&
            k_ == lambda.k() &&
            c_ == lambda.c());
  }

  template <class IBASIS>
  bool
  IntervalMultiIndex<IBASIS>::operator < (const IntervalMultiIndex<IBASIS>& lambda) const
  {
    return (j_ < lambda.j() ||
           (j_ == lambda.j() && (e_ < lambda.e() ||
                                (e_ == lambda.e() && (k_ < lambda.k() ||
                                                     (k_ == lambda.k() && c_< lambda.c()))))));
  }
  
  template <class IBASIS>
  IntervalMultiIndex<IBASIS>&
  IntervalMultiIndex<IBASIS>::operator ++ ()
  {
    int kmax;

    if (c_ < basis_->number_of_components-1)
      c_++;
    else { // c_ == basis_->number_of_components-1
      c_ = 0;
      kmax = (e_ == E_GENERATOR) ? basis_->DeltaRmax(j_) : basis_->Nablamax(j_);
      if (k_ < kmax)
        k_++;
      else { // k_ == kmax
        if (e_ == E_GENERATOR) {
          k_ = basis_->Nablamin();
          e_ = E_WAVELET;
        }
        else { // e_ == E_WAVELET
          k_ = basis_->Nablamin();
          j_++;
        }
      }
    }
    
    return *this;
  }
  
  template <class IBASIS>
  IntervalMultiIndex<IBASIS>&
  IntervalMultiIndex<IBASIS>::operator -- ()
  {
    int kmin;

    if (c_ > 0)
      c_--;
    else { // c_ == 0
      c_ = basis_->number_of_components-1;
      kmin = (e_ == E_GENERATOR) ? basis_->DeltaLmin() : basis_->Nablamin();
      if (k_ > kmin)
        k_--;
      else { // k_ == kmin
        k_ = (e_ == E_GENERATOR) ? basis_->DeltaRmax(j_) : basis_->Nablamax(j_);
        if (j_ == basis_->j0()) {
          assert((j_ > basis_->j0()) || (e_ == E_WAVELET));
          e_ = E_GENERATOR;
        }
        else {
          j_--;
        }
      }
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
    return IntervalMultiIndex<IBASIS>(j, E_GENERATOR, basis->DeltaLmin(), 0, basis);
  }
  
  template <class IBASIS>
  inline
  IntervalMultiIndex<IBASIS>
  last_generator(const IBASIS* basis, const int j)
  {
    assert(j >= basis->j0());
    return IntervalMultiIndex<IBASIS>(j, E_GENERATOR, basis->DeltaRmax(j), basis->number_of_components-1, basis);
  }

  template <class IBASIS>
  inline
  IntervalMultiIndex<IBASIS>
  first_wavelet(const IBASIS* basis, const int j)
  {
    assert(j >= basis->j0());
    return IntervalMultiIndex<IBASIS>(j, E_WAVELET, basis->Nablamin(), 0, basis);
  }
  
  template <class IBASIS>
  inline
  IntervalMultiIndex<IBASIS>
  last_wavelet(const IBASIS* basis, const int j)
  {
    assert(j >= basis->j0());
    return IntervalMultiIndex<IBASIS>(j, E_WAVELET, basis->Nablamax(j), basis->number_of_components-1, basis);
  }

  template <class IBASIS>
  inline
  IntervalMultiIndex<IBASIS>
  first_index(const IBASIS* basis, const int j, const typename IntervalMultiIndex<IBASIS>::type_type e)
  {
    return (e == E_GENERATOR ? first_generator(basis, j) : first_wavelet(basis, j));
  }
  
  template <class IBASIS>
  inline
  IntervalMultiIndex<IBASIS>
  last_index(const IBASIS* basis, const int j, const typename IntervalMultiIndex<IBASIS>::type_type e)
  {
    return (e == E_GENERATOR ? last_generator(basis, j) : last_wavelet(basis, j));
  }

}
