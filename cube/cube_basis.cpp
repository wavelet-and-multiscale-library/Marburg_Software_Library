// implementation for cube_basis.h

namespace WaveletTL
{
  template <class IBASIS, unsigned int DIM>
  CubeBasis<IBASIS,DIM>::CubeBasis()
    : bases_()
  {
    for (unsigned int i = 0; i < DIM; i++)
      bases_[i] = new IBASIS();
    j0_ = bases_[0]->j0();
  }

  template <class IBASIS, unsigned int DIM>
  CubeBasis<IBASIS,DIM>::~CubeBasis()
  {
    for (unsigned int i = 0; i < DIM; i++)
      delete bases_[i];
  }

}
