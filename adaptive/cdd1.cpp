// implementation for cdd1.h

#include <cmath>
#include <set>
#include <algorithm>

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
    params.gamma = 0.9;
    params.F = P.F_norm();

    // determination of q=q1=q2=q3,q4 according to [CDD1, (7.23)ff]
    params.q4 = 1 / (20 * params.kappa);
    const double A = params.c1 / (20 * (3 + params.c1 / params.c2));
    const double B = params.c2 * (0.1 - params.q4 * sqrt(params.kappa));
    const double C = params.q4 / (1 / params.c2 + 6 * (params.gamma + 1) / (params.gamma * params.c1));
    params.q1 = params.q2 = params.q3 = std::min(A, std::min(B, C));

    params.q0 = sqrt(params.kappa) + params.q3/params.c2;

    params.theta = sqrt(1 - params.c1 * params.gamma * params.gamma / (4 * params.c2));
    params.theta_bar = 1 - 1 / (6 * params.kappa);

    params.K = (unsigned int) floor(log(20 * params.kappa) / fabs(log(params.theta))) + 1;

    cout << "CDD1_SOLVE parameters:" << endl;
    cout << "c1=" << params.c1 << ", c2=" << params.c2 << ", kappa=" << params.kappa << endl;
    cout << "gamma=" << params.gamma << endl;
    cout << "F=" << params.F << endl;
    cout << "q0=" << params.q0 << ", q1=" << params.q1 << ", q2=" << params.q2
	 << ", q3=" << params.q3 << ", q4=" << params.q4 << endl;
    cout << "theta=" << params.theta << ", theta_bar=" << params.theta_bar << endl;
    cout << "K=" << params.K << endl;

    typedef typename PROBLEM::WaveletBasis::Index Index;
    set<Index> Lambda, Lambda_hat;
    u_epsilon.clear();
    double delta = params.F;

    InfiniteVector<double,Index> v_hat, r_hat, u_bar, F;
    P.RHS(2*params.q2*epsilon, F);
    while (delta > sqrt(params.c1)*epsilon) { // check the additional factor c1^{1/2} in [BB+] !?
      cout << "CDD1_SOLVE: delta=" << delta << endl;
      NPROG(P, params, F, Lambda, u_epsilon, delta, v_hat, Lambda_hat, r_hat, u_bar);
      if (l2_norm(r_hat)+(params.q1+params.q2+(1+1/params.kappa)*params.q3)*delta <= params.c1*epsilon)
	{
	  u_epsilon = u_bar;
	  break;
	}
      else
	{
	  u_epsilon = v_hat;
	}
      delta *= 0.5;
    }
  }

  template <class PROBLEM>
  void NPROG(const PROBLEM& P, const CDD1Parameters& params,
	     const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& F,
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
    GALERKIN(P, params, F, Lambda_k, v, delta, params.q3*delta/params.c2, u_Lambda_k);
    while (true) {
      NGROW(P, params, F, Lambda_k, u_Lambda_k, params.q1*delta, params.q2*delta, Lambda_kplus1, r_hat);
      if (l2_norm(r_hat) <= params.c1*delta/20 || k == params.K) {
	u_Lambda_k.COARSE(2*delta/5, v_hat);
	v_hat.support(Lambda_hat);
	break;
      }
      GALERKIN(P, params, F, Lambda_kplus1, u_Lambda_k, params.q0*delta, params.q3*delta/params.c2, v_hat);
      u_Lambda_k.swap(v_hat);
      Lambda_k.swap(Lambda_kplus1);
      k++;
    }
  }

  template <class PROBLEM>
  void GALERKIN(const PROBLEM& P, const CDD1Parameters& params,
		const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& F,
		const set<typename PROBLEM::WaveletBasis::Index>& Lambda,
 		const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& v,
 		const double delta,
		const double eta,
 		InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& u_bar)
  {
    cout << "GALERKIN called..." << endl;
    typedef typename PROBLEM::WaveletBasis::Index Index;
    InfiniteVector<double,Index> r;
    u_bar = v;
    double mydelta = delta;
    while (true) {
      INRESIDUAL(P, params, F, Lambda, u_bar, params.c1*eta/6, params.c1*eta/6, r);
      const double inresidual_norm = l2_norm(r);
      if (eta >= std::min(params.theta_bar*mydelta, inresidual_norm/params.c1+eta/3)) {
	cout << "... GALERKIN done, norm of internal residual: " << inresidual_norm << endl;
	break;
      }
      u_bar += 1/params.c2 * r;
      mydelta *= params.theta_bar;
    }
  }

  template <class PROBLEM>
  void NGROW(const PROBLEM& P, const CDD1Parameters& params,
	     const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& F,
	     const set<typename PROBLEM::WaveletBasis::Index>& Lambda,
 	     const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& u_bar,
 	     const double xi1,
 	     const double xi2,
 	     set<typename PROBLEM::WaveletBasis::Index>& Lambda_tilde,
 	     InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& r)
  {
    cout << "NGROW called..." << endl;
    typedef typename PROBLEM::WaveletBasis::Index Index;
    set<Index> Lambda_c;
    NRESIDUAL(P, params, F, Lambda, u_bar, xi1, xi2, r, Lambda_c);
    const double residual_norm = l2_norm(r);
    cout << "* in NGROW, current residual norm is " << residual_norm << endl;
    InfiniteVector<double,Index> pr;
    r.COARSE(sqrt(1-params.gamma*params.gamma)*residual_norm, pr);
    pr.support(Lambda_c);
    Lambda_tilde.clear();
    cout << "* in NGROW, size of Lambda is " << Lambda.size() << ", size of Lambda_c is " << Lambda_c.size() << endl;
    std::set_union(Lambda.begin(), Lambda.end(),
		   Lambda_c.begin(), Lambda_c.end(),
		   inserter(Lambda_tilde, Lambda_tilde.begin()));
    cout << "... NGROW done, size of new index set: " << Lambda_tilde.size() << endl;
  }

  template <class PROBLEM>
  void INRESIDUAL(const PROBLEM& P, const CDD1Parameters& params,
		  const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& F,
		  const set<typename PROBLEM::WaveletBasis::Index>& Lambda,
		  const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& v,
		  const double eta1,
		  const double eta2,
		  InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& r)
  {
    typedef typename PROBLEM::WaveletBasis::Index Index;
    InfiniteVector<double,Index> w, g(F);

    // TODO: speed up the following two lines
    APPLY(P, v, eta1, w);
    w.clip(Lambda);

    g.clip(Lambda);
    g.COARSE(eta2, r);
    r -= w;
  }

  template <class PROBLEM>
  void NRESIDUAL(const PROBLEM& P, const CDD1Parameters& params,
		 const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& F,
		 const set<typename PROBLEM::WaveletBasis::Index>& Lambda,
		 const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& v,
		 const double eta1,
		 const double eta2,
		 InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& r,
		 set<typename PROBLEM::WaveletBasis::Index>& Lambda_tilde)
  {
    typedef typename PROBLEM::WaveletBasis::Index Index;
    InfiniteVector<double,Index> w;
    APPLY(P, v, eta1, w);
    F.COARSE(eta2, r);
    r -= w;
    r.support(Lambda_tilde);
  }
}
