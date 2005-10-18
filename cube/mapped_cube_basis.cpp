// implementation for mapped_cube_basis.h

using MathTL::AffineLinearMapping;

namespace WaveletTL
{
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  MappedCubeBasis<IBASIS,DIM_d,DIM_m>::MappedCubeBasis()
    : CubeBasis<IBASIS,DIM_d>(), delete_kappa(true)
  {
    kappa_ = new AffineLinearMapping<DIM_d>();
  }

  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  MappedCubeBasis<IBASIS,DIM_d,DIM_m>::MappedCubeBasis(const Chart<DIM_d,DIM_m>* kappa,
						       const FixedArray1D<int,2*DIM_d>& s,
						       const FixedArray1D<int,2*DIM_d>& sT)
    : CubeBasis<IBASIS,DIM_d>(s, sT), kappa_(kappa), delete_kappa(false)
  {
  }

  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  MappedCubeBasis<IBASIS,DIM_d,DIM_m>::~MappedCubeBasis()
  {
    if (delete_kappa)
      delete kappa_;
  }
}
