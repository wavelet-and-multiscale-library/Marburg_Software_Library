// implementation of some (inline) Tensor<RANK, DIM>:: methods

#include <cassert>

namespace MathTL
{
  template <unsigned int RANK, unsigned int DIM>
  const unsigned int Tensor<RANK, DIM>::dimension;

  template <unsigned int DIM>
  const unsigned int Tensor<1, DIM>::dimension;

  template <unsigned int RANK, unsigned int DIM>
  const unsigned int Tensor<RANK, DIM>::rank;

  template <unsigned int DIM>
  const unsigned int Tensor<1, DIM>::rank;

  template <unsigned int RANK, unsigned int DIM>
  inline
  Tensor<RANK, DIM>::Tensor()
  {
    // due to the recursive structure and the specialization to Tensor<1, DIM>,
    // we don't need to initialize anything here
  }

  template <unsigned int RANK, unsigned int DIM>
  inline
  void Tensor<RANK, DIM>::clear()
  {
    for (unsigned int i(0); i < DIM; ++i)
      subtensor[i].clear();
  }

  template <unsigned int RANK, unsigned int DIM>
  inline
  const typename Tensor<RANK, DIM>::value_type&
  Tensor<RANK, DIM>::operator [] (const size_type i) const
  {
    assert(i < DIM);
    return subtensor[i];
  }

  template <unsigned int RANK, unsigned int DIM>
  inline
  typename Tensor<RANK, DIM>::value_type&
  Tensor<RANK, DIM>::operator [] (const size_type i)
  {
    assert(i < DIM);
    return subtensor[i];
  }

  template <unsigned int RANK, unsigned int DIM>
  inline
  Tensor<RANK, DIM> &
  Tensor<RANK, DIM>::operator = (const Tensor<RANK, DIM>& T)
  {
    for (unsigned int i(0); i < DIM; ++i)
      subtensor[i] = T.subtensor[i];
    return *this;
  }

  template <unsigned int RANK, unsigned int DIM>
  inline
  bool Tensor<RANK, DIM>::operator == (const Tensor<RANK, DIM>& T) const
  {
    for (unsigned int i(0); i < DIM; ++i)
      if (subtensor[i] != T.subtensor[i]) return false;
    return true;
  }

  template <unsigned int RANK, unsigned int DIM>
  inline
  bool Tensor<RANK, DIM>::operator != (const Tensor<RANK, DIM>& T) const
  {
    return !(*this == T);
  }

  template <unsigned int RANK, unsigned int DIM>
  inline
  Tensor<RANK, DIM>&
  Tensor<RANK, DIM>::operator += (const Tensor<RANK, DIM>& T)
  {
    for (unsigned int i(0); i < DIM; ++i)
      subtensor[i] += T.subtensor[i];
    return *this;
  }
  
  template <unsigned int RANK, unsigned int DIM>
  inline
  Tensor<RANK, DIM>&
  Tensor<RANK, DIM>::operator -= (const Tensor<RANK, DIM>& T)
  {
    for (unsigned int i(0); i < DIM; ++i)
      subtensor[i] -= T.subtensor[i];
    return *this;
  }

  template <unsigned int RANK, unsigned int DIM>
  inline
  Tensor<RANK, DIM>&
  Tensor<RANK, DIM>::operator *= (const double s)
  {
    for (unsigned int i(0); i < DIM; ++i)
      subtensor[i] *= s;
    return *this;
  }

  template <unsigned int RANK, unsigned int DIM>
  inline
  Tensor<RANK, DIM>&
  Tensor<RANK, DIM>::operator /= (const double s)
  {
    for (unsigned int i(0); i < DIM; ++i)
      subtensor[i] /= s;
    return *this;
  }
  
  template <unsigned int RANK, unsigned int DIM>
  inline
  Tensor<RANK, DIM>
  Tensor<RANK, DIM>::operator + (const Tensor<RANK, DIM>& T) const
  {
    Tensor<RANK, DIM> r(*this); // uses default copy constructor
    for (unsigned int i(0); i < DIM; ++i)
      r.subtensor[i] += T.subtensor[i];
    return r;
  }

  template <unsigned int RANK, unsigned int DIM>
  inline
  Tensor<RANK, DIM>
  Tensor<RANK, DIM>::operator - (const Tensor<RANK, DIM>& T) const
  {
    Tensor<RANK, DIM> r(*this); // uses default copy constructor
    for (unsigned int i(0); i < DIM; ++i)
      r.subtensor[i] -= T.subtensor[i];
    return r;
  }
  
  template <unsigned int RANK, unsigned int DIM>
  inline
  Tensor<RANK, DIM>
  Tensor<RANK, DIM>::operator - () const
  {
    Tensor<RANK, DIM> r;
    for (unsigned int i(0); i < DIM; ++i)
      r.subtensor[i] = -subtensor[i];
    return r;
  }

  template <unsigned int RANK, unsigned int DIM>
  inline
  const typename Tensor<RANK, DIM>::size_type
  Tensor<RANK, DIM>::memory_consumption()
  {
    return sizeof(Tensor<RANK, DIM>);
  }

  
  //
  //
  // implementation for external Tensor<RANK, DIM> functionality

  template <unsigned int RANK, unsigned int DIM>
  inline
  std::ostream& operator << (std::ostream& os, const Tensor<RANK, DIM>& T)
  {
    os << "[";
    for (unsigned int i(0); i < DIM-1; ++i)
      os << T[i] << ' ';
    os << T[DIM-1] << "]";

    return os;
  }

  template <unsigned int RANK>
  inline
  std::ostream& operator << (std::ostream& os, const Tensor<RANK, 1>& T)
  {
    os << "[" << T[0] << "]";
    
    return os;
  }

}
