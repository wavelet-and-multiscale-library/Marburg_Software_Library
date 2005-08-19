// implementation for sturm_bf.h

namespace WaveletTL
{
  template <class WBASIS>
  SturmBilinearForm<WBASIS>::SturmBilinearForm(const simpleSturmBVP& bvp)
    : bvp_(bvp), wbasis_(bvp.bc_left(), bvp.bc_right())
  {
  }

  template <class WBASIS>
  double
  SturmBilinearForm<WBASIS>::operator () (const typename WBASIS::Index& lambda,
					  const typename WBASIS::Index& nu) const
  {
    // a(u,v) = \int_0^1 [p(t)u'(t)v'(t)+q(t)u(t)v(t)] dt

    double r = 0;
    
    return r;
  }
}
