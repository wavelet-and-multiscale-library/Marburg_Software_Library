// implementation for lin_par_equation.h

#include <adaptive/apply.h>
#include <adaptive/cdd1.h>

namespace WaveletTL
{
  template <class ELLIPTIC_EQ>
  LinParEqROWStageEquationHelper<ELLIPTIC_EQ>
  ::LinParEqROWStageEquationHelper
  (const double a,
   const ELLIPTIC_EQ* elliptic,
   const InfiniteVector<double, typename ELLIPTIC_EQ::Index>& z)
    : alpha(a), T(elliptic), y(z)
  {
  }
  
  template <class ELLIPTIC_EQ>
  void
  LinParEqROWStageEquationHelper<ELLIPTIC_EQ>
  ::add_level (const Index& lambda,
	       InfiniteVector<double, Index>& w, const int j,
	       const double factor,
	       const int J,
	       const CompressionStrategy strategy) const
  {
    T->add_level(lambda, w, j, factor, J, strategy);
    if (lambda.j() == j)
      w.add_coefficient(lambda, alpha*factor/(T->D(lambda)*T->D(lambda)));
  }

  template <class ELLIPTIC_EQ>
  LinearParabolicEquation<ELLIPTIC_EQ>
  ::LinearParabolicEquation(const ELLIPTIC_EQ* helper,
			    const InfiniteVector<double,Index>& initial)
    : elliptic(helper)
  {
    AbstractIVP<InfiniteVector<double,Index> >::u0 = initial;
  }

  template <class ELLIPTIC_EQ>
  void
  LinearParabolicEquation<ELLIPTIC_EQ>
  ::evaluate_f(const double t,
	       const InfiniteVector<double,Index>& v,
	       const double tolerance,
	       InfiniteVector<double,Index>& result) const
  {
    result.clear();
    InfiniteVector<double,Index> w(v);
    elliptic->rescale(w, 1); // w = Dv
    APPLY(*elliptic, w, tolerance, result, 8, St04a); // yields -D^{-1}AD^{-1}w
    elliptic->rescale(result, 1);
    result.scale(-1.0);

    
    // no source term, for the moment
  }
  
  template <class ELLIPTIC_EQ>
  void
  LinearParabolicEquation<ELLIPTIC_EQ>
  ::evaluate_ft(const double t,
		const InfiniteVector<double,Index>& v,
		const double tolerance,
		InfiniteVector<double,Index>& result) const
  {
    result.clear(); // constant source term, for the moment
  }
    
  template <class ELLIPTIC_EQ>
  void
  LinearParabolicEquation<ELLIPTIC_EQ>
  ::solve_ROW_stage_equation(const double t,
			     const InfiniteVector<double,Index>& v,
			     const double alpha,
			     const InfiniteVector<double,Index>& y,
			     const double tolerance,
			     InfiniteVector<double,Index>& result) const
  {
    LinParEqROWStageEquationHelper<ELLIPTIC_EQ> helper(alpha, elliptic, y);
    CDD1_SOLVE(helper, tolerance, result, 8); // D^{-1}(alpha*I-T)D^{-1}*Dx = D^{-1}y    
    elliptic->rescale(result, -1); // Dx -> x
  }
}
