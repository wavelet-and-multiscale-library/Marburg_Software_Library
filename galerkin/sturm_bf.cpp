// implementation for sturm_bf.h

namespace WaveletTL
{
  template <class WBASIS>
  SturmBilinearForm<WBASIS>::SturmBilinearForm(const simpleSturmBVP& bvp)
    : bvp_(bvp), wbasis_(bvp.bc_left(), bvp.bc_right())
  {
  }
}
