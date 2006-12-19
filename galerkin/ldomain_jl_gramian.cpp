// implementation for ldomain_jl_gramian.h

namespace WaveletTL
{
  LDomainJLGramian::LDomainJLGramian(const WaveletBasis& basis,
				     const InfiniteVector<double,Index>& y)
    : basis_(basis), y_(y), normA(0.0), normAinv(0.0)
  {
  }
  
}
