// implementation for cdd2.h

#include <cmath>
#include <set>

#include <adaptive/apply.h>

using std::set;

namespace WaveletTL
{
  template <class PROBLEM>
  void CDD2_SOLVE(const PROBLEM& P, const double nu, const double epsilon,
		  InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& u_epsilon)
  {
    typedef typename PROBLEM::WaveletBasis::Index Index;

    // compute optimal relaxation parameter omega
    const double omega = 2.0 / (P.norm_A() + 1.0/P.norm_Ainv());
    cout << "CDD2_SOLVE: omega=" << omega << endl;

    // compute spectral norm rho
    const double cond_A = P.norm_A() * P.norm_Ainv();
    const double rho = (cond_A - 1.0) / (cond_A + 1.0);
    cout << "CDD2_SOLVE: rho=" << rho << endl;
    
    // desired error reduction factor theta < 1/3
    const double theta = 2.0/7.0;
    cout << "CDD2_SOLVE: theta=" << theta << endl;

    // compute minimal K such that 3*rho^K < theta
    const int K = (int) ceil(log10(theta/3.0) / log10(rho));
    cout << "CDD2_SOLVE: K=" << K << endl;
    
    u_epsilon.clear();

    double epsilon_k = nu;
    InfiniteVector<double,Index> f, v;
    P.RHS(1e-6, f);
    while (epsilon_k > epsilon) {
      epsilon_k *= 3*pow(rho, K) / theta;
      cout << "CDD2_SOLVE: epsilon_k=" << epsilon_k << endl;
      const double eta = theta * epsilon_k / (6*omega*K);
//       P.RHS(eta, f);
      v = u_epsilon;
      for (int j = 1; j <= K; j++) {
	APPLY(P, u_epsilon, eta, v);
	v += omega * (f - v);
      }
      v.COARSE((1-theta)*epsilon_k, u_epsilon);
    } 
    
    set<Index> Lambda;
  }
}
