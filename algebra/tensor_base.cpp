// implementation of some (inline) Tensor<1, DIM>:: methods

#include <cassert>
#include <io/vector_io.h>

namespace MathTL
{
  template <unsigned int DIM>
  inline
  Tensor<1, DIM>::Tensor(const bool initialize)
    : values(DIM > 0 ? DIM : 1)
  {
    if (initialize)
      {
	for (unsigned int i(0); i < DIM; ++i)
	  values[i] = 0;
      }
  }

  template <unsigned int DIM>
  inline
  void Tensor<1, DIM>::clear()
  {
    for (unsigned int i(0); i < DIM; ++i)
      values[i] = 0;
  }

  template <unsigned int DIM>
  inline
  const typename Tensor<1, DIM>::size_type
  Tensor<1, DIM>::size() const
  {
    return dimension;
  }

  template <unsigned int DIM>
  inline
  const double Tensor<1, DIM>::operator [] (const size_type i) const
  {
    assert(i >= 0 && i < DIM);
    return values[i];
  }

  template <unsigned int DIM>
  inline
  double& Tensor<1, DIM>::operator [] (const size_type i)
  {
    assert(i >= 0 && i < DIM);
    return values[i];
  }

  template <unsigned int DIM>
  inline
  bool Tensor<1, DIM>::operator == (const Tensor<1, DIM>& T) const
  {
    for (unsigned int i(0); i < DIM; ++i)
      if (values[i] != T.values[i]) return false;
    return true;
  }

  template <unsigned int DIM>
  inline
  bool Tensor<1, DIM>::operator != (const Tensor<1, DIM>& T) const
  {
    return !((*this) == T);
  }

  template <unsigned int DIM>
  inline
  Tensor<1, DIM>&
  Tensor<1, DIM>::operator += (const Tensor<1, DIM>& T)
  {
    for (unsigned int i(0); i < DIM; ++i)
      values[i] += T.values[i];
    return *this;
  }
  
  template <unsigned int DIM>
  inline
  Tensor<1, DIM>&
  Tensor<1, DIM>::operator -= (const Tensor<1, DIM>& T)
  {
    for (unsigned int i(0); i < DIM; ++i)
      values[i] -= T.values[i];
    return *this;
  }

  template <unsigned int DIM>
  inline
  Tensor<1, DIM>&
  Tensor<1, DIM>::operator *= (const double s)
  {
    for (unsigned int i(0); i < DIM; ++i)
      values[i] *= s;
    return *this;
  }

  template <unsigned int DIM>
  inline
  Tensor<1, DIM>&
  Tensor<1, DIM>::operator /= (const double s)
  {
    assert(s != 0);
    for (unsigned int i(0); i < DIM; ++i)
      values[i] /= s;
    return *this;
  }

  template <unsigned int DIM>
  inline
  Tensor<1, DIM>
  Tensor<1, DIM>::operator + (const Tensor<1, DIM>& T) const
  {
    return (Tensor<1, DIM>(*this) += T);
  }

  template <unsigned int DIM>
  inline
  Tensor<1, DIM>
  Tensor<1, DIM>::operator - (const Tensor<1, DIM>& T) const
  {
    return (Tensor<1, DIM>(*this) -= T);
  }
  
  template <unsigned int DIM>
  inline
  Tensor<1, DIM>
  Tensor<1, DIM>::operator - () const
  {
    return (Tensor<1, DIM>() -= (*this));
  }

  template <unsigned int DIM>
  inline
  double Tensor<1, DIM>::operator * (const Tensor<1, DIM>& T) const
  {
    double r(0);
    for (unsigned int i(0); i < DIM; ++i)
      r += values[i] * T.values[i];
    return r;
  }

  template <unsigned int DIM>
  inline
  const typename Tensor<1, DIM>::size_type
  Tensor<1, DIM>::memory_consumption()
  {
    return sizeof(Tensor<1, DIM>);
  }

  //
  //
  // implementation for external Tensor<RANK, DIM> functionality

  template <unsigned int DIM>
  inline
  std::ostream& operator << (std::ostream& os, const Tensor<1, DIM>& T)
  {
    print_vector(T, os);
    return os;
  }

  inline
  std::ostream& operator << (std::ostream& os, const Tensor<1, 1>& T)
  {
    print_vector(T, os);
    return os;
  }
}
