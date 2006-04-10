// implementation for ldomain_index.h

namespace WaveletTL
{
  template <class IBASIS>
  LDomainIndex<IBASIS>::LDomainIndex(const LDomainBasis<IBASIS>* basis)
    : basis_(basis), p_(0)
  {
    if (basis_ == 0) {
      j_ = 0; // invalid (e and k are initialized by zero automatically)
    } else {
//       j_ = basis_->j0(); // coarsest level;
//       // e_ is zero by default: generator
//       for (unsigned int i = 0; i < DIM; i++)
// 	k_[i] = basis_->bases()[i]->DeltaLmin();
    }
  }

}
