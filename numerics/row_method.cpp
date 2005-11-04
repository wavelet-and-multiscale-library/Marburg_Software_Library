// implementation for row_method.h

namespace MathTL
{
  template <class VECTOR>
  ROWMethod<VECTOR>::ROWMethod(const typename WMethod<VECTOR>::Method method,
			       const WMethodStageEquationSolver<VECTOR>& stage_equation_solver)
    : WMethod<VECTOR>(method, stage_equation_solver)
  {
  }
  
  template <class VECTOR>
  void
  ROWMethod<VECTOR>::increment(const AbstractIVP<VECTOR>& ivp,
			       const double t_m,
			       const VECTOR& u_m,
			       const double tau,
			       VECTOR& u_mplus1,
			       VECTOR& error_estimate,
			       const double tolerance) const
  {
  }

  template <class VECTOR>
  void
  ROWMethod<VECTOR>::solve_stage_equation(const AbstractIVP<VECTOR>& ivp,
					  const double alpha,
					  const VECTOR& y,
					  const double tolerance,
					  VECTOR& x) const
  {
  }
   
  template <class VECTOR>
  void
  ROWMethod<VECTOR>::g(const AbstractIVP<VECTOR>& ivp,
		       const double t_m,
		       const VECTOR& u_m,
		       const double tolerance,
		       VECTOR& result) const
  {
  }

}
