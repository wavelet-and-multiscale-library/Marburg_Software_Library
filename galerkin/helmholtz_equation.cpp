// implementation for helmholtz_equation.h

namespace WaveletTL
{
  template <int d, int dT>
  HelmholtzEquation1D<d,dT>::HelmholtzEquation1D(const Function<1>* f,
						 const double alpha,
						 const bool precompute_rhs)
    : f_(f), alpha_(alpha),
      basis_("P","",1,1,0,0) // PBasis, complementary b.c.'s
  {
  }
}
