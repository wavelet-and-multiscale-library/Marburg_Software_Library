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
  void
  EllipticBVP<DIM>::set_f(const Function<DIM>* f)
  {
    f_ = f;
  }


  template <unsigned int DIM>
  PoissonBVP<DIM>::PoissonBVP(const Function<DIM>* f)
    : EllipticBVP<DIM>(f, f, f)
  {
  }

  template <unsigned int DIM>
  PoissonBVP_Coeff<DIM>::PoissonBVP_Coeff(const Function<DIM>* a, const Function<DIM>* f)
    : EllipticBVP<DIM>(a, f, f)
  {
  }

  template <unsigned int DIM>
  PoissonBVP_Coeff<DIM>::PoissonBVP_Coeff(const Function<DIM>* a, const Function<DIM>* f, const Function<DIM>* f_ar)
    : EllipticBVP<DIM>(a, f, f), f_ar_(f_ar)
  {
  }

  template <unsigned int DIM>
  void
  PoissonBVP_Coeff<DIM>::set_a(const Function<DIM>* a)
  {
    this->a_ = a;
  }



  template <unsigned int DIM>
  IdentityBVP<DIM>::IdentityBVP(const Function<DIM>* f)
    : EllipticBVP<DIM>(f, f, f)
  {
  }

  template <unsigned int DIM>
  BiharmonicBVP<DIM>::BiharmonicBVP(const Function<DIM>* f)
    :f_(f)
  {
  }

  template <unsigned int DIM>
  void
  BiharmonicBVP<DIM>::set_f(const Function<DIM>* f)
  {
    f_ = f;
  }

}
