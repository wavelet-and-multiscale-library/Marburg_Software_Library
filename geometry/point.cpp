// implementation of the Point<DIM>:: inline methods

#include <io/vector_io.h>

namespace MathTL
{
  template <unsigned int DIM>
  inline
  Point<DIM>::Point()
    : Tensor<1, DIM>(true)
  {
  }

  template <unsigned int DIM>
  inline
  Point<DIM>::Point(const Tensor<1, DIM>& T)
    : Tensor<1, DIM>(T)
  {
  }

  template <unsigned int DIM>
  inline
  Point<DIM>::Point(const double x)
  {
    assert(DIM == 1);
    this->values[0] = x;
  }

  template <unsigned int DIM>
  inline
  Point<DIM>::Point(const double x, const double y)
  {
    assert(DIM == 2);
    this->values[0] = x;
    this->values[1] = y;
  }

  template <unsigned int DIM>
  inline
  Point<DIM>::Point(const double x, const double y, const double z)
  {
    assert(DIM == 3);
    this->values[0] = x;
    this->values[1] = y;
    this->values[2] = z;
  }

  template <unsigned int DIM>
  inline
  const typename Point<DIM>::size_type
  Point<DIM>::size() const
  {
    return dimension;
  }

  template <unsigned int DIM>
  inline
  const double Point<DIM>::operator () (const size_type i) const
  {
    assert(i < DIM);
    return this->values[i];
  }

  template <unsigned int DIM>
  inline
  double& Point<DIM>::operator () (const size_type i)
  {
    assert(i < DIM);
    return this->values[i];
  }

  template <unsigned int DIM>
  inline
  std::ostream& operator << (std::ostream& os, const Point<DIM>& p)
  {
    print_vector(p, os);
    return os;
  }
}
