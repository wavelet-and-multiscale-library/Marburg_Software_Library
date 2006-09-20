// implementation for row_stage_equation.h

// #include <adaptive/apply.h>
// #include <adaptive/cdd1.h>

namespace WaveletTL
{

//   template <class ELLIPTIC_EQ>
//   LinParEqROWStageEquationHelper<ELLIPTIC_EQ>
//   ::LinParEqROWStageEquationHelper
//   (const double a,
//    const ELLIPTIC_EQ* elliptic,
//    const CachedProblem<IntervalGramian<typename ELLIPTIC_EQ::WaveletBasis> >& GC,
//    const InfiniteVector<double, typename ELLIPTIC_EQ::Index>& z)
//     : alpha(a), T(elliptic), G(GC), y(z)
//   {
//   }
  
//   template <class ELLIPTIC_EQ>
//   void
//   LinParEqROWStageEquationHelper<ELLIPTIC_EQ>
//   ::add_level (const Index& lambda,
// 	       InfiniteVector<double, Index>& w, const int j,
// 	       const double factor,
// 	       const int J,
// 	       const CompressionStrategy strategy) const
//   {
//     // Gramian part, we have to take care of the preconditioning factors
//     InfiniteVector<double,Index> g;
//     G.add_level(lambda, g, j, factor * alpha/T->D(lambda), J, strategy);
//     g.scale(this, -1);
//     w.add(g);
    
//     T->add_level(lambda, w, j, factor, J, strategy);

// //     G.add_level(lambda, w, j, factor * alpha/(T->D(lambda)*T->D(lambda))
// //     if (lambda.j() == j)
// //       w.add_coefficient(lambda, alpha*factor/(T->D(lambda)*T->D(lambda)));
//   }
  
//   template <class ELLIPTIC_EQ>
//   LinearParabolicEquation<ELLIPTIC_EQ>
//   ::LinearParabolicEquation(const ELLIPTIC_EQ* helper,
//  			    const InfiniteVector<double,typename ELLIPTIC_EQ::Index>& initial,
//  			    const InfiniteVector<double,typename ELLIPTIC_EQ::Index>& f,
//  			    const int jmax)
//     : elliptic(helper), G(helper->basis(), InfiniteVector<double,typename ELLIPTIC_EQ::Index>()),
//       GC(&G), constant_f_(f), f_(0), jmax_(jmax)
//   {
//     AbstractIVP<InfiniteVector<double,typename ELLIPTIC_EQ::Index> >::u0 = initial;
//   }

//   template <class ELLIPTIC_EQ>
//   void
//   LinearParabolicEquation<ELLIPTIC_EQ>
//   ::evaluate_f(const double t,
// 	       const InfiniteVector<double,Index>& v,
// 	       const double tolerance,
// 	       InfiniteVector<double,Index>& result) const
//   {
//     result.clear();
//     InfiniteVector<double,Index> w(v), temp;
//     w.scale(elliptic, 1); // w = Dv
//     APPLY(*elliptic, w, tolerance, temp, jmax_, St04a); // yields -D^{-1}AD^{-1}w
//     temp.scale(elliptic, 1);
//     temp.scale(-1.0); // result = -D(-D^{-1}AD^{-1}Dv) = Av

// //     // multiply with inverse primal gramian (i.e., switch from dual to primal basis)
// //     G.set_rhs(temp);
// //     CDD1_SOLVE(GC, tolerance, result, jmax_);

//     result = temp;
    
//     // add constant driving term (if present)
//     if (!constant_f_.empty())
//       result.add(constant_f_);

//     // add time-dependent driving term (if present)
//     if (f_ != 0) {
//       f_->set_time(t);
//       w.clear();
// //       expand(f_, elliptic->basis(), false, jmax_, w); // expand in the primal basis
//       expand(f_, elliptic->basis(), true, jmax_, w); // expand in the dual (!) basis
//       result.add(w);
//     }
//   }
  
//   template <class ELLIPTIC_EQ>
//   void
//   LinearParabolicEquation<ELLIPTIC_EQ>
//   ::evaluate_ft(const double t,
// 		const InfiniteVector<double,Index>& v,
// 		const double tolerance,
// 		InfiniteVector<double,Index>& result) const
//   {
//     result.clear(); // from the constant driving term

//     // approximate derivative of time-dependent driving term (if present)
//     if (f_ != 0) {
//       const double h = 1e-6;
//       InfiniteVector<double,Index> fhelp;
//       f_->set_time(t);
// //       expand(f_, elliptic->basis(), false, jmax_, fhelp); // expand in the primal basis
//       expand(f_, elliptic->basis(), true, jmax_, fhelp); // expand in the dual (!) basis
//       f_->set_time(t+h);
// //       expand(f_, elliptic->basis(), false, jmax_, result);
//       expand(f_, elliptic->basis(), true, jmax_, result);
//       result.add(-1., fhelp);
//       result.scale(1./h);
//     }
//   }
    
//   template <class ELLIPTIC_EQ>
//   void
//   LinearParabolicEquation<ELLIPTIC_EQ>
//   ::solve_ROW_stage_equation(const double t,
// 			     const InfiniteVector<double,Index>& v,
// 			     const double alpha,
// 			     const InfiniteVector<double,Index>& y,
// 			     const double tolerance,
// 			     InfiniteVector<double,Index>& result) const
//   {
// //     // multiply everything with the Gramian
// //     InfiniteVector<double,Index> Gy;
// //     APPLY(GC, y, tolerance, Gy, jmax_, St04a);
// //     LinParEqROWStageEquationHelper<ELLIPTIC_EQ> helper(alpha, elliptic, GC, Gy);

//     LinParEqROWStageEquationHelper<ELLIPTIC_EQ> helper(alpha, elliptic, GC, y);
//     CDD1_SOLVE(helper, tolerance, result, jmax_); // D^{-1}(alpha*I-T)D^{-1}*Dx = D^{-1}y    
//     result.scale(elliptic, -1); // Dx -> x
//   }
}
