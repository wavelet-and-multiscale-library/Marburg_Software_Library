// implementation for ldomain_equation.h

namespace WaveletTL
{
  template <class IBASIS>
  LDomainEquation<IBASIS>::LDomainEquation(const EllipticBVP<2>* bvp)
    : bvp_(bvp), basis_(), normA(0.0), normAinv(0.0)
  {
    // do something
//     compute_rhs();
  }
}
