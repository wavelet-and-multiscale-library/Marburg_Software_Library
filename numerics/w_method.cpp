// implementation for w_method.h

namespace MathTL
{
  template <class VECTOR>
  void
  WMethod<VECTOR>::increment(const AbstractIVP<VECTOR>& ivp,
			     const double t_m, const VECTOR& u_m,
			     const double tau,
			     VECTOR& u_mplus1,
			     VECTOR& error_estimate,
			     const double tolerance) const
  {
  }
}
