// implementation for cdd1.h

#include <cmath>
#include <set>

#include <adaptive/apply.h>

using std::set;

namespace WaveletTL
{
  template <class PROBLEM>
  void CDD1_SOLVE(const PROBLEM& P, const double epsilon,
		  InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& u_epsilon)
  {
    // INIT, cf. [BB+] 

    CDD1Parameters params;

    params.c1 = 1.0/P.norm_Ainv();
    params.c2 = P.norm_A();
    params.kappa = params.c2/params.c1;
    params.gamma = 0.5;
    params.F = P.F_norm();

    // determination of q=q1=q2=q3,q4 according to [CDD1, (7.23)ff]
    params.q4 = 1 / (20 * params.kappa);
    const double A = params.c1 / (20 * (3 + params.c1 / params.c2));
    const double B = params.c2 * (0.1 - params.q4 * sqrt(params.kappa));
    const double C = params.q4 / (1 / params.c2 + 6 * (params.gamma + 1) / (params.gamma * params.c1));
    params.q1 = params.q2 = params.q3 = std::min(A, std::min(B, C));

    params.q0 = sqrt(params.kappa) + params.q3/params.c2;

    params.theta = sqrt(1 - params.c1 * params.gamma * params.gamma / (4 * params.c2));
    params.thetabar = 1 - 1 / (6 * params.kappa);

    params.K = (unsigned int) floor(log(20 * params.kappa) / fabs(log(params.theta))) + 1;

    typedef typename PROBLEM::WaveletBasis::Index Index;
    set<Index> Lambda, Lambdahat;
    u_epsilon.clear();
    double delta = params.F;

    InfiniteVector<double,Index> vhat, rhat, ubar, f;
    while (delta > sqrt(params.c1)*epsilon) { // check the additional factor c1^{1/2} in [BB+]
      cout << "CDD1_SOLVE: delta=" << delta << endl;
      NPROG(P, params, Lambda, u_epsilon, delta, vhat, Lambdahat, rhat, ubar);
      if (l2_norm(rhat)+(params.q1+params.q2+(1+1/params.kappa)*params.q3)*delta <= params.c1*epsilon)
	{
	  u_epsilon = ubar;
	  break;
	}
      else
	{
	  u_epsilon = vhat;
	}
      delta *= 0.5;
    }
  }

  template <class PROBLEM>
  void NPROG(const PROBLEM& P, const CDD1Parameters& params,
	     const set<typename PROBLEM::WaveletBasis::Index>& Lambda,
	     const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& v,
	     const double delta,
	     InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& v_hat,
	     set<typename PROBLEM::WaveletBasis::Index>& Lambda_hat,
	     InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& r_hat,
	     InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& u_Lambda_k)
  {
    typedef typename PROBLEM::WaveletBasis::Index Index;
    set<Index> Lambda_k(Lambda), Lambda_kplus1;
    unsigned int k = 0;
    GALERKIN(P, params, Lambda_k, v, delta, params.q3*delta/params.c2, u_Lambda_k);
    while (true) {
      NGROW(P, params, Lambda_k, u_Lambda_k, params.q1*delta, params.q2*delta, Lambda_kplus1, r_hat);
      if (l2_norm(r_hat) <= params.c1*delta/20 || k == params.K) {
	u_Lambda_k.COARSE(2*delta/5, v_hat);
	v_hat.support(Lambda_hat);
	break;
      }
      GALERKIN(P, params, Lambda_kplus1, u_Lambda_k, params.q0*delta, params.q3*delta/params.c2, v_hat);
      u_Lambda_k = v_hat;
      Lambda_k = Lambda_kplus1;
      k++;
    }
  }

  template <class PROBLEM>
  void GALERKIN(const PROBLEM& P, const CDD1Parameters& params,
		const set<typename PROBLEM::WaveletBasis::Index>& Lambda,
 		const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& v,
 		const double delta,
		const double eta,
 		InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& ubar)
  {
    typedef typename PROBLEM::WaveletBasis::Index Index;
    InfiniteVector<double,Index> r;
    INRESIDUAL(P, params, Lambda, v, params.c1*eta/6, params.c1*eta/6, r);
  }

  template <class PROBLEM>
  void NGROW(const PROBLEM& P, const CDD1Parameters& params,
	     const set<typename PROBLEM::WaveletBasis::Index>& Lambda,
 	     const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& ubar,
 	     const double xi1,
 	     const double xi2,
 	     set<typename PROBLEM::WaveletBasis::Index>& LambdaTilde,
 	     InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& r)
  {
  }

  template <class PROBLEM>
  void INRESIDUAL(const PROBLEM& P, const CDD1Parameters& params,
		  const set<typename PROBLEM::WaveletBasis::Index>& Lambda,
		  const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& v,
		  const double eta1,
		  const double eta2,
		  InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& r)
  {
  }
}
