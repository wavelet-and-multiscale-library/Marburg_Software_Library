// implementation of some (inline) Tensor<1, DIM, VALUE>:: methods

#include <cassert>
#include <io/vector_io.h>

namespace MathTL
{
  template <unsigned int DIM, class VALUE>
  inline
  Tensor<1, DIM, VALUE>::Tensor(const bool initialize)
    : values()
  {
    if (initialize)
      {
	for (unsigned int i(0); i < DIM; ++i)
	  values[i] = 0;
      }
  }

  template <unsigned int DIM, class VALUE>
  inline
  void Tensor<1, DIM, VALUE>::clear()
  {
    for (unsigned int i(0); i < DIM; ++i)
      values[i] = 0;
  }

  template <unsigned int DIM, class VALUE>
  inline
  const typename Tensor<1, DIM, VALUE>::size_type
  Tensor<1, DIM, VALUE>::size() const
  {
    return dimension;
  }

  template <unsigned int DIM, class VALUE>
  inline
  const VALUE Tensor<1, DIM, VALUE>::operator [] (const size_type i) const
  {
    assert(i >= 0 && i < DIM);
    return values[i];
  }

  template <unsigned int DIM, class VALUE>
  inline
  VALUE& Tensor<1, DIM, VALUE>::operator [] (const size_type i)
  {
    assert(i >= 0 && i < DIM);
    return values[i];
  }

  template <unsigned int DIM, class VALUE>
  inline
  bool Tensor<1, DIM, VALUE>::operator == (const Tensor<1, DIM, VALUE>& T) const
  {
    for (unsigned int i(0); i < DIM; ++i)
      if (values[i] != T.values[i]) return false;
    return true;
  }

  template <unsigned int DIM, class VALUE>
  inline
  bool Tensor<1, DIM, VALUE>::operator != (const Tensor<1, DIM, VALUE>& T) const
  {
    return !((*this) == T);
  }

  template <unsigned int DIM, class VALUE>
  inline
  Tensor<1, DIM, VALUE>&
  Tensor<1, DIM, VALUE>::operator += (const Tensor<1, DIM, VALUE>& T)
  {
    for (unsigned int i(0); i < DIM; ++i)
      values[i] += T.values[i];
    return *this;
  }
  
  template <unsigned int DIM, class VALUE>
  inline
  void Tensor<1, DIM, VALUE>::add(const VALUE s, const Tensor<1, DIM, VALUE>& T)
  {
    for (unsigned int i(0); i < DIM; ++i)
      values[i] += s * T.values[i];
  }

  template <unsigned int DIM, class VALUE>
  inline
  Tensor<1, DIM, VALUE>&
  Tensor<1, DIM, VALUE>::operator -= (const Tensor<1, DIM, VALUE>& T)
  {
    for (unsigned int i(0); i < DIM; ++i)
      values[i] -= T.values[i];
    return *this;
  }

  template <unsigned int DIM, class VALUE>
  inline
  Tensor<1, DIM, VALUE>&
  Tensor<1, DIM, VALUE>::operator *= (const VALUE s)
  {
    for (unsigned int i(0); i < DIM; ++i)
      values[i] *= s;
    return *this;
  }

  template <unsigned int DIM, class VALUE>
  inline
  Tensor<1, DIM, VALUE>&
  Tensor<1, DIM, VALUE>::operator /= (const VALUE s)
  {
    assert(s != 0);
    for (unsigned int i(0); i < DIM; ++i)
      values[i] /= s;
    return *this;
  }

  template <unsigned int DIM, class VALUE>
  inline
  Tensor<1, DIM, VALUE>
  Tensor<1, DIM, VALUE>::operator + (const Tensor<1, DIM, VALUE>& T) const
  {
    return (Tensor<1, DIM, VALUE>(*this) += T);
  }

  template <unsigned int DIM, class VALUE>
  inline
  Tensor<1, DIM, VALUE>
  Tensor<1, DIM, VALUE>::operator - (const Tensor<1, DIM, VALUE>& T) const
  {
    return (Tensor<1, DIM, VALUE>(*this) -= T);
  }
  
  template <unsigned int DIM, class VALUE>
  inline
  Tensor<1, DIM, VALUE>
  Tensor<1, DIM, VALUE>::operator - () const
  {
    return (Tensor<1, DIM, VALUE>() -= (*this));
  }

  template <unsigned int DIM, class VALUE>
  inline
  VALUE Tensor<1, DIM, VALUE>::operator * (const Tensor<1, DIM, VALUE>& T) const
  {
    VALUE r(0);
    for (unsigned int i(0); i < DIM; ++i)
      r += values[i] * T.values[i];
    return r;
  }

  template <unsigned int DIM, class VALUE>
  inline
  const typename Tensor<1, DIM, VALUE>::size_type
  Tensor<1, DIM, VALUE>::memory_consumption()
  {
    return sizeof(Tensor<1, DIM, VALUE>);
  }

  //
  //
  // implementation for external Tensor<RANK, DIM, VALUE> functionality

  template <unsigned int DIM, class VALUE>
  inline
  std::ostream& operator << (std::ostream& os,
			     const Tensor<1, DIM, VALUE>& T)
  {
    print_vector(T, os);
    return os;
  }

  template <class VALUE>
  inline
  std::ostream& operator << (std::ostream& os,
			     const Tensor<1, 1, VALUE>& T)
  {
    print_vector(T, os);
    return os;
  }
}
