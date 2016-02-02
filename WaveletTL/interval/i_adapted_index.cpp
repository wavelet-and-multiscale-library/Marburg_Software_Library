// implementation for i_adapted_index.h

#include <cassert>

namespace WaveletTL
{
  template <class IAdaptedBasis>
  IntervalAdaptedIndex<IAdaptedBasis>::IntervalAdaptedIndex(const IAdaptedBasis* basis)
    : basis_(basis)
  {
    if (basis_ == 0) {
      index_ = 0;
    } else {
      index_ = new IMultiIndex(basis->multi_basis());
    }
  }

  template <class IAdaptedBasis>
  IntervalAdaptedIndex<IAdaptedBasis>::IntervalAdaptedIndex(const IMultiIndex& mu, const IAdaptedBasis* basis)
    : basis_(basis)
  {
    index_ = new IMultiIndex(mu); // make a deep copy
  }
  
  template <class IAdaptedBasis>
  IntervalAdaptedIndex<IAdaptedBasis>::IntervalAdaptedIndex(const IntervalAdaptedIndex<IAdaptedBasis>& lambda)
    : basis_(lambda.basis())
  {
    index_ = new IMultiIndex(*lambda.multi_index()); // make a deep copy of the encapsulated multi-index
  }

  template <class IAdaptedBasis>
  IntervalAdaptedIndex<IAdaptedBasis>::IntervalAdaptedIndex(const IntervalAdaptedIndex<IAdaptedBasis>* lambda)
    : basis_(lambda->basis())
  {
    index_ = new IMultiIndex(*lambda->multi_index()); // make a deep copy of the encapsulated multi-index
  }

  template <class IAdaptedBasis>
  IntervalAdaptedIndex<IAdaptedBasis>::IntervalAdaptedIndex(const int j, const int e, const int k,
                                       const IAdaptedBasis* basis)
  {
    typename IMultiIndex::translation_type k_multi;
    typename IMultiIndex::component_type c_multi;

    IMultiIndex::ck_decode(k, k_multi, c_multi);
    index_ = new IMultiIndex(j, e, k_multi, c_multi, basis->multi_basis());
    basis_ = basis;
  }

  template <class IAdaptedBasis>
  IntervalAdaptedIndex<IAdaptedBasis>::IntervalAdaptedIndex(const int num,
                                       const IAdaptedBasis* basis)
  :  basis_(basis)
  {
    index_ = new IMultiIndex(num, basis_->multi_basis());
  }

  template <class IAdaptedBasis>
  IntervalAdaptedIndex<IAdaptedBasis>::~IntervalAdaptedIndex()
  {
    if (index_ != 0)
      delete index_; // free memory of capsulated multi-index object
  }


  template <class IAdaptedBasis>
  IntervalAdaptedIndex<IAdaptedBasis>&
  IntervalAdaptedIndex<IAdaptedBasis>::operator = (const IntervalAdaptedIndex<IAdaptedBasis>& lambda)
  {
    basis_ = lambda.basis();
    if (index_ == 0)
      index_ = new IMultiIndex(*lambda.multi_index());
    else
      index_->operator= (*lambda.multi_index());

    return *this;
  }

  template <class IAdaptedBasis>
  bool
  IntervalAdaptedIndex<IAdaptedBasis>::operator == (const IntervalAdaptedIndex<IAdaptedBasis>& lambda) const
  {
    return index_->operator== (*lambda.multi_index());
  }

  template <class IAdaptedBasis>
  bool
  IntervalAdaptedIndex<IAdaptedBasis>::operator < (const IntervalAdaptedIndex<IAdaptedBasis>& lambda) const
  {
    return index_->operator< (*lambda.multi_index());
  }
  
  template <class IAdaptedBasis>
  IntervalAdaptedIndex<IAdaptedBasis>&
  IntervalAdaptedIndex<IAdaptedBasis>::operator ++ ()
  {
    index_->operator++ ();
    
    return *this;
  }


  template <class IAdaptedBasis>
  inline
  IntervalAdaptedIndex<IAdaptedBasis> first_generator(const IAdaptedBasis* basis, const int j)
  {
    assert(j >= basis->j0());
    return IntervalAdaptedIndex<IAdaptedBasis>(j, 0, basis->DeltaLmin(), basis);
  }
  
  template <class IAdaptedBasis>
  inline
  IntervalAdaptedIndex<IAdaptedBasis> last_generator(const IAdaptedBasis* basis, const int j)
  {
    assert(j >= basis->j0());
    return IntervalAdaptedIndex<IAdaptedBasis>(j, 0, basis->DeltaRmax(j), basis);
  }

  template <class IAdaptedBasis>
  inline
  IntervalAdaptedIndex<IAdaptedBasis> first_wavelet(const IAdaptedBasis* basis, const int j)
  {
    assert(j >= basis->j0());
    return IntervalAdaptedIndex<IAdaptedBasis>(j, 1, basis->Nablamin(), basis);
  }
  
  template <class IAdaptedBasis>
  inline
  IntervalAdaptedIndex<IAdaptedBasis> last_wavelet(const IAdaptedBasis* basis, const int j)
  {
    assert(j >= basis->j0());
    return IntervalAdaptedIndex<IAdaptedBasis>(j, 1, basis->Nablamax(j), basis);
  }

  template <class IAdaptedBasis>
  inline
  IntervalAdaptedIndex<IAdaptedBasis> first_index(const IAdaptedBasis* basis, const int j, const int e)
  {
    return (e == 0 ? first_generator(basis, j) : first_wavelet(basis, j));
  }
  
  template <class IAdaptedBasis>
  inline
  IntervalAdaptedIndex<IAdaptedBasis> last_index(const IAdaptedBasis* basis, const int j, const int e)
  {
    return (e == 0 ? last_generator(basis, j) : last_wavelet(basis, j));
  }

  template <class IAdaptedBasis>
  inline
  int first_generator_num(const IAdaptedBasis* basis)
  {
    IntervalAdaptedIndex<IAdaptedBasis> ind(first_generator<IAdaptedBasis>(basis, basis->j0()));
    ind.set_number();
    return ind.number();
  }
  
  template <class IAdaptedBasis>
  inline
  int last_generator_num(const IAdaptedBasis* basis)
  {
    IntervalAdaptedIndex<IAdaptedBasis> ind(last_generator<IAdaptedBasis>(basis, basis->j0()));
    ind.set_number();
    return ind.number();

  }

  template <class IAdaptedBasis>
  inline
  int first_wavelet_num(const IAdaptedBasis* basis, const int j)
  {
    assert(j >= basis->j0());
    IntervalAdaptedIndex<IAdaptedBasis> ind(first_wavelet<IAdaptedBasis>(basis, j));
    ind.set_number();
    return ind.number();

  }
  
  template <class IAdaptedBasis>
  inline
  int last_wavelet_num(const IAdaptedBasis* basis, const int j)
  {
    assert(j >= basis->j0());
    IntervalAdaptedIndex<IAdaptedBasis> ind(last_wavelet<IAdaptedBasis>(basis, j));
    ind.set_number();
    return ind.number();
  }

}
