// implementation of some (inline) Tensor<RANK, DIM, VALUE>:: methods

#include <cassert>

namespace MathTL
{
  template <unsigned int RANK, unsigned int DIM, class VALUE>
  const unsigned int Tensor<RANK, DIM, VALUE>::dimension;

  template <unsigned int DIM, class VALUE>
  const unsigned int Tensor<1, DIM, VALUE>::dimension;

  template <unsigned int RANK, unsigned int DIM, class VALUE>
  const unsigned int Tensor<RANK, DIM, VALUE>::rank;

  template <unsigned int DIM, class VALUE>
  const unsigned int Tensor<1, DIM, VALUE>::rank;

  template <unsigned int RANK, unsigned int DIM, class VALUE>
  inline
  Tensor<RANK, DIM, VALUE>::Tensor()
  {
    // due to the recursive structure and the specialization to
    // Tensor<1, DIM, VALUE>,
    // we don't need to initialize anything here
  }

  template <unsigned int RANK, unsigned int DIM, class VALUE>
  inline
  void Tensor<RANK, DIM, VALUE>::clear()
  {
    for (unsigned int i(0); i < DIM; ++i)
      subtensor[i].clear();
  }

  template <unsigned int RANK, unsigned int DIM, class VALUE>
  inline
  const typename Tensor<RANK, DIM, VALUE>::value_type&
  Tensor<RANK, DIM, VALUE>::operator [] (const size_type i) const
  {
    assert(i < DIM);
    return subtensor[i];
  }

  template <unsigned int RANK, unsigned int DIM, class VALUE>
  inline
  typename Tensor<RANK, DIM, VALUE>::value_type&
  Tensor<RANK, DIM, VALUE>::operator [] (const size_type i)
  {
    assert(i < DIM);
    return subtensor[i];
  }

  template <unsigned int RANK, unsigned int DIM, class VALUE>
  inline
  Tensor<RANK, DIM, VALUE> &
  Tensor<RANK, DIM, VALUE>::operator = (const Tensor<RANK, DIM, VALUE>& T)
  {
    for (unsigned int i(0); i < DIM; ++i)
      subtensor[i] = T.subtensor[i];
    return *this;
  }

  template <unsigned int RANK, unsigned int DIM, class VALUE>
  inline
  bool Tensor<RANK, DIM, VALUE>::operator == (const Tensor<RANK, DIM, VALUE>& T) const
  {
    for (unsigned int i(0); i < DIM; ++i)
      if (subtensor[i] != T.subtensor[i]) return false;
    return true;
  }

  template <unsigned int RANK, unsigned int DIM, class VALUE>
  inline
  bool Tensor<RANK, DIM, VALUE>::operator != (const Tensor<RANK, DIM, VALUE>& T) const
  {
    return !(*this == T);
  }

  template <unsigned int RANK, unsigned int DIM, class VALUE>
  inline
  Tensor<RANK, DIM, VALUE>&
  Tensor<RANK, DIM, VALUE>::operator += (const Tensor<RANK, DIM, VALUE>& T)
  {
    for (unsigned int i(0); i < DIM; ++i)
      subtensor[i] += T.subtensor[i];
    return *this;
  }
  
  template <unsigned int RANK, unsigned int DIM, class VALUE>
  inline
  Tensor<RANK, DIM, VALUE>&
  Tensor<RANK, DIM, VALUE>::operator -= (const Tensor<RANK, DIM, VALUE>& T)
  {
    for (unsigned int i(0); i < DIM; ++i)
      subtensor[i] -= T.subtensor[i];
    return *this;
  }

  template <unsigned int RANK, unsigned int DIM, class VALUE>
  inline
  Tensor<RANK, DIM, VALUE>&
  Tensor<RANK, DIM, VALUE>::operator *= (const double s)
  {
    for (unsigned int i(0); i < DIM; ++i)
      subtensor[i] *= s;
    return *this;
  }

  template <unsigned int RANK, unsigned int DIM, class VALUE>
  inline
  Tensor<RANK, DIM, VALUE>&
  Tensor<RANK, DIM, VALUE>::operator /= (const double s)
  {
    for (unsigned int i(0); i < DIM; ++i)
      subtensor[i] /= s;
    return *this;
  }
  
  template <unsigned int RANK, unsigned int DIM, class VALUE>
  inline
  Tensor<RANK, DIM, VALUE>
  Tensor<RANK, DIM, VALUE>::operator + (const Tensor<RANK, DIM, VALUE>& T) const
  {
    Tensor<RANK, DIM, VALUE> r(*this); // uses default copy constructor
    for (unsigned int i(0); i < DIM; ++i)
      r.subtensor[i] += T.subtensor[i];
    return r;
  }

  template <unsigned int RANK, unsigned int DIM, class VALUE>
  inline
  Tensor<RANK, DIM, VALUE>
  Tensor<RANK, DIM, VALUE>::operator - (const Tensor<RANK, DIM, VALUE>& T) const
  {
    Tensor<RANK, DIM, VALUE> r(*this); // uses default copy constructor
    for (unsigned int i(0); i < DIM; ++i)
      r.subtensor[i] -= T.subtensor[i];
    return r;
  }
  
  template <unsigned int RANK, unsigned int DIM, class VALUE>
  inline
  Tensor<RANK, DIM, VALUE>
  Tensor<RANK, DIM, VALUE>::operator - () const
  {
    Tensor<RANK, DIM, VALUE> r;
    for (unsigned int i(0); i < DIM; ++i)
      r.subtensor[i] = -subtensor[i];
    return r;
  }

  template <unsigned int RANK, unsigned int DIM, class VALUE>
  inline
  const typename Tensor<RANK, DIM, VALUE>::size_type
  Tensor<RANK, DIM, VALUE>::memory_consumption()
  {
    return sizeof(Tensor<RANK, DIM, VALUE>);
  }

  
  //
  //
  // implementation for external Tensor<RANK, DIM, VALUE> functionality

  template <unsigned int RANK, unsigned int DIM, class VALUE>
  inline
  std::ostream& operator << (std::ostream& os,
			     const Tensor<RANK, DIM, VALUE>& T)
  {
    os << "[";
    for (unsigned int i(0); i < DIM-1; ++i)
      os << T[i] << ' ';
    os << T[DIM-1] << "]";

    return os;
  }

  template <unsigned int RANK, class VALUE>
  inline
  std::ostream& operator << (std::ostream& os,
			     const Tensor<RANK, 1, VALUE>& T)
  {
    os << "[" << T[0] << "]";
    
    return os;
  }

}
