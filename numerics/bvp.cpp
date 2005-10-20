// implementation for bvp.h

namespace MathTL
{
  template <unsigned int DIM>
  inline
  TwoPointBVP<DIM>::~TwoPointBVP()
  {
  }

  template <unsigned int DIM>
  EllipticBVP<DIM>::EllipticBVP(const Function<DIM>* a,
				const Function<DIM>* q,
				const Function<DIM>* f)
    : a_(a), q_(q), f_(f)
  {
  }

  template <unsigned int DIM>
  PoissonBVP<DIM>::PoissonBVP(const Function<DIM>* f)
    : EllipticBVP<DIM>(f, f, f)
  {
  }
}
