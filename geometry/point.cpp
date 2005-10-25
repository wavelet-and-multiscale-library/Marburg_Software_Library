// implementation of the Point<DIM>:: inline methods

#include <io/vector_io.h>

namespace MathTL
{
  template <unsigned int DIM, class VALUE>
  inline
  Point<DIM, VALUE>::Point()
    : Tensor<1, DIM, VALUE>(true)
  {
  }

  template <unsigned int DIM, class VALUE>
  inline
  Point<DIM, VALUE>::Point(const Tensor<1, DIM, VALUE>& T)
    : Tensor<1, DIM, VALUE>(T)
  {
  }

  template <unsigned int DIM, class VALUE>
  inline
  Point<DIM, VALUE>::Point(const VALUE x)
  {
    for (unsigned int i = 0; i < DIM; i++)
      values[i] = x;
  }

  template <unsigned int DIM, class VALUE>
  inline
  Point<DIM, VALUE>::Point(const VALUE x, const VALUE y)
  {
    assert(DIM == 2);
    this->values[0] = x;
    this->values[1] = y;
  }

  template <unsigned int DIM, class VALUE>
  inline
  Point<DIM, VALUE>::Point(const VALUE x, const VALUE y, const VALUE z)
  {
    assert(DIM == 3);
    this->values[0] = x;
    this->values[1] = y;
    this->values[2] = z;
  }

  template <unsigned int DIM, class VALUE>
  Point<DIM, VALUE>& Point<DIM, VALUE>::operator = (const VALUE x)
  {
    for (unsigned int i = 0; i < DIM; i++)
      values[i] = x;
    return *this;
  }

//   template <unsigned int DIM, class VALUE>
//   Point<DIM, VALUE>& Point<DIM, VALUE>::operator = (const Point<DIM, VALUE>& p)
//   {
//     return (Tensor<1, DIM, VALUE>::operator = (p));
//   }

  template <unsigned int DIM, class VALUE>
  inline
  const typename Point<DIM, VALUE>::size_type
  Point<DIM, VALUE>::size() const
  {
    return this->dimension;
  }

  template <unsigned int DIM, class VALUE>
  inline
  const VALUE Point<DIM, VALUE>::operator () (const size_type i) const
  {
    assert(i < DIM);
    return this->values[i];
  }

  template <unsigned int DIM, class VALUE>
  inline
  VALUE& Point<DIM, VALUE>::operator () (const size_type i)
  {
    assert(i < DIM);
    return this->values[i];
  }

  template <unsigned int DIM, class VALUE>
  inline
  std::ostream& operator << (std::ostream& os, const Point<DIM, VALUE>& p)
  {
    print_vector(p, os);
    return os;
  }
}
