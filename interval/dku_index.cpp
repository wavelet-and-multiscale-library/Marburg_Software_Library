// implementation for dku_index.h

namespace WaveletTL
{
  template <int d, int dT>
  DKUIndex<d, dT>::DKUIndex(const DKUBasis<d, dT>* basis)
    : basis_(basis)
  {
    if (basis_ == 0) {
      j_ = e_ = k_ = 0; // invalid!
    } else {
      k_ = basis_->DeltaLmin(); // leftmost
      e_ = 0;                   // generator
      j_ = basis_->j0();        // on the coarsest level
    }
  }
  
  template <int d, int dT>
  DKUIndex<d, dT>::DKUIndex(const DKUIndex<d, dT>& lambda)
  {
    j_ = lambda.j();
    e_ = lambda.e();
    k_ = lambda.k();
  }

  template <int d, int dT>
  DKUIndex<d, dT>::DKUIndex(const int j, const int e, const int k)
  {
    j_ = j;
    e_ = e;
    k_ = k;
  }

  template <int d, int dT>
  DKUIndex<d, dT>& DKUIndex<d, dT>::operator = (const DKUIndex<d, dT>& lambda)
  {
    j_ = lambda.j();
    e_ = lambda.e();
    k_ = lambda.k();

    return *this;
  }

  template <int d, int dT>
  bool DKUIndex<d, dT>::operator == (const DKUIndex<d, dT>& lambda) const
  {
    return (j_ == lambda.j() &&
	    e_ == lambda.e() &&
	    k_ == lambda.k());
  }

  template <int d, int dT>
  bool DKUIndex<d, dT>::operator < (const DKUIndex<d, dT>& lambda) const
  {
    return (j_ < lambda.j() ||
	    (j_ == lambda.j() && e_ < lambda.e()) ||
	    (j_ == lambda.j() && e_ == lambda.e() && k_ < lambda.k()));
  }
}
