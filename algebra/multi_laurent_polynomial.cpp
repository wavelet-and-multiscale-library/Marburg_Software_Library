// implementation for multi_laurent_polynomial.h

namespace MathTL
{

  template <class R, unsigned int DIMENSION>
  R MultiLaurentPolynomial<R, DIMENSION>::value
  (const Point<DIMENSION>& p,
   const unsigned int component) const
  {
    // TODO: implement Horner scheme! (multivariate?)
    return R(42);
//     return value(p[0]);
  }

  template <class R, unsigned int DIMENSION>
  inline
  void MultiLaurentPolynomial<R, DIMENSION>::vector_value
  (const Point<DIMENSION> &p,
   Vector<R>& values) const
  {
    values.resize(1, false);
    values[0] = value(p);
  }
}
