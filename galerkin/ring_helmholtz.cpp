// implementation for ring_helmholtz.h

namespace WaveletTL
{
  template <int d, int dt, int s0, int s1>
  RingHelmholtzEquation<d,dt,s0,s1>::RingHelmholtzEquation
  (const WaveletBasis& basis,
   const char* G_file,
   const char* A_file,
   const int jmax,
   const double alpha,
   const InfiniteVector<double,Index>& y)
    : basis_(basis), alpha_(alpha), y_(y),
      G_(basis, InfiniteVector<double,Index>()),
      GC_(&G_, G_file, jmax, 1.0, 1.0), // dirty
      A_(basis, InfiniteVector<double,Index>()),
      AC_(&A_, A_file, jmax, 1.0, 1.0), // dirty
      normA(1.0), normAinv(1.0) // dirty
  {
    y_precond = y_;
    y_precond.scale(this, -1);
  }

  template <int d, int dt, int s0, int s1>
  inline
  double
  RingHelmholtzEquation<d,dt,s0,s1>::a
  (const Index& lambda,
   const Index& mu) const
  {
    return alpha_ * GC_.a(lambda, mu) + A_.a(lambda, mu);
  }
  
  template <int d, int dt, int s0, int s1>
  void
  RingHelmholtzEquation<d,dt,s0,s1>::set_alpha(const double alpha) const
  {
    assert(alpha >= 0);
    alpha_ = alpha;
    y_precond = y_;
    y_precond.scale(this, -1);
  }

}
