// implementation for periodic_gramian.h

namespace WaveletTL
{
  template <class RBASIS>
  PeriodicIntervalGramian<RBASIS>::PeriodicIntervalGramian
  (const PeriodicBasis<RBASIS>& basis,
   const InfiniteVector<double, typename PeriodicBasis<RBASIS>::Index>& y)
    : basis_(basis), y_(y),
      normA(1.0),
      normAinv(ldexp(1.0, 2*(RBASIS::primal_polynomial_degree()-1))) // lower bound from [Bittner]
  {
  }
  
  template <class RBASIS>
  inline
  double
  PeriodicIntervalGramian<RBASIS>::a(const typename WaveletBasis::Index& lambda,
				     const typename WaveletBasis::Index& mu) const
  {
    return basis_.integrate(lambda, mu);
  }
  
}
