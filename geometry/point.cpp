// implementation of the Point<DIM>:: inline methods

namespace MathTL
{
  template <int DIM>
  inline
  Point<DIM>::Point()
    : Tensor<1, DIM>(true)
  {
  }

  template <int DIM>
  inline
  Point<DIM>::Point(const Tensor<1, DIM>& T)
    : Tensor<1, DIM>(T)
  {
  }

  template <int DIM>
  inline
  Point<DIM>::Point(const double x)
  {
    assert(DIM == 1);
    this->values[0] = x;
  }

  template <int DIM>
  inline
  Point<DIM>::Point(const double x, const double y)
  {
    assert(DIM == 2);
    this->values[0] = x;
    this->values[1] = y;
  }

  template <int DIM>
  inline
  Point<DIM>::Point(const double x, const double y, const double z)
  {
    assert(DIM == 3);
    this->values[0] = x;
    this->values[1] = y;
    this->values[2] = z;
  }

  template <int DIM>
  inline
  const double Point<DIM>::operator () (const size_type i) const
  {
    assert(i < DIM);
    return this->values[i];
  }

  template <int DIM>
  inline
  double& Point<DIM>::operator () (const size_type i)
  {
    assert(i < DIM);
    return this->values[i];
  }
}
