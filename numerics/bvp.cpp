// implementation for bvp.h

namespace MathTL
{
  template <unsigned int DIM>
  inline
  TwoPointBVP<DIM>::~TwoPointBVP()
  {
  }

  template <unsigned int DIM, class ATLAS>
  EllipticBVP<DIM,ATLAS>::EllipticBVP(const ATLAS* atlas,
				      const Array1D<FixedArray1D<int,2*DIM> >& bc,
				      const Function<DIM>* a,
				      const Function<DIM>* q,
				      const Function<DIM>* f)
    : atlas_(atlas), delete_atlas(false), bc_(bc), delete_functions(false), a_(a), q_(q), f_(f)
  {
  }
  
  template <unsigned int DIM, class ATLAS>
  EllipticBVP<DIM,ATLAS>::~EllipticBVP()
  {
    if (delete_atlas) delete atlas_;
    if (delete_functions) {
      delete a_;
      delete q_;
      delete f_;
    }
  }
}
