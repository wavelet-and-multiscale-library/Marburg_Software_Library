// implementation for precond.h

#include <cmath>

namespace WaveletTL
{
  template <class INDEX>
  InfinitePreconditioner<INDEX>::~InfinitePreconditioner() {}

  template <class INDEX>
  inline
  void
  InfiniteSymmetricPreconditioner<INDEX>::apply_left_preconditioner
  (const InfiniteVector<double,INDEX>& y,
   InfiniteVector<double,INDEX>& x) const
  {
    apply_preconditioner(y, x);
  }
    
  template <class INDEX>
  inline
  void
  InfiniteSymmetricPreconditioner<INDEX>::apply_right_preconditioner
  (const InfiniteVector<double,INDEX>& y,
   InfiniteVector<double,INDEX>& x) const
  {
    apply_preconditioner(y, x);
  }
  
  template <class INDEX>
  inline
  void
  FullyDiagonalPreconditioner<INDEX>::apply_preconditioner
  (const InfiniteVector<double,INDEX>& y,
   InfiniteVector<double,INDEX>& x) const
  {
    x = y;
    x.scale(this, -1);
  };
  
  template <class INDEX>
  inline
  double
  WaveletNEPreconditioner<INDEX>::diag(const INDEX& lambda) const
  {
    return pow(ldexp(1.0, lambda.j()), operator_order());
  }

  template <class INDEX>
  inline
  double
  EnergyNormPreconditioner<INDEX>::diag(const INDEX& lambda) const
  {
    return sqrt(a(lambda, lambda));
  };
}
