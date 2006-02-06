// implementation for precond.h

namespace WaveletTL
{
  template <class PROBLEM>
  inline
  void
  DiagonalPreconditioner<PROBLEM>::apply_preconditioner
  (const InfiniteVector<double, typename PROBLEM::Index>& y,
   InfiniteVector<double, typename PROBLEM::Index>& x) const
  {
    x = y;
    x.scale(*this, -1);
  };

}
