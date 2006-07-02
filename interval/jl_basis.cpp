// implementation for jl_basis.h

namespace WaveletTL
{
  JLBasis::JLBasis(const int s0, const int s1)
  {
    this->s0 = s0;
    this->s1 = s1;
    
    setup();
  }
  
  void
  JLBasis::setup() {
    j0_ = 1;
  }
  
  void
  JLBasis::reconstruct(const InfiniteVector<double, Index>& c,
		       const int j,
		       InfiniteVector<double, Index>& v) const {
    v.clear();
    for (InfiniteVector<double, Index>::const_iterator it(c.begin()), itend(c.end());
	 it != itend; ++it) {
      InfiniteVector<double, Index> help;
      reconstruct_1(it.index(), j, help);
      v.add(*it, help);
    }
  }
  
  void
  JLBasis::reconstruct_1(const Index& lambda,
			 const int j,
			 InfiniteVector<double, Index>& c) const {
    c.clear();
    
    if (lambda.j() >= j)
      c.set_coefficient(lambda, 1.0); 
    else {
      // For the reconstruction of psi_lambda, we have to compute
      // the corresponding column of the transformation matrix Mj=(Mj0, Mj1).
      
    }
  }
  
}
