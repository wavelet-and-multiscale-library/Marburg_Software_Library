// implementation for precond.h

namespace WaveletTL
{
  template <class PROBLEM>
  void
  DiagonalPreconditioner<PROBLEM>::apply_preconditioner
  (const InfiniteVector<double, typename PROBLEM::Index>& y,
   InfiniteVector<double, typename PROBLEM::Index>& x) const
  {
    x = y;
    for (typename InfiniteVector<double, typename PROBLEM::Index>::iterator it(x.begin());
	 it != x.end(); ++it)
      {
      }
  };

}
