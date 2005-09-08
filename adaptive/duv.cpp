// implementation for duv.h

#include <adaptive/apply.h>

namespace WaveletTL
{
  template <class PROBLEM>
  void DUV_SOLVE_SD(const PROBLEM& P, const double nu, const double epsilon,
		    InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& u_epsilon)
  {
    typedef typename PROBLEM::WaveletBasis::Index Index;

    // compute spectral norm rho
    const double cond_A = P.norm_A() * P.norm_Ainv();
    const double rho = (cond_A - 1.0) / (cond_A + 1.0);
    cout << "DUV_SOLVE_SD: rho=" << rho << endl;
    
    // desired error reduction factor theta < 1/3
    const double theta = 2.0/7.0;
    cout << "DUV_SOLVE: theta=" << theta << endl;

    // compute minimal K such that 3*rho^K < theta
    const int K = (int) ceil(log10(theta/3.0) / log10(rho));
    cout << "DUV_SOLVE: K=" << K << endl;
    
    u_epsilon.clear();

    double epsilon_k = nu;
    InfiniteVector<double,Index> f, rj, v, Av, Arj;
    P.RHS(1e-6, f);
    while(epsilon_k > epsilon) {
      cout << "DUV_SOLVE: epsilon_k=" << epsilon_k << endl;
      v = u_epsilon;
      for (int j = 0; j < K; j++) {
	const double eta = pow(rho, (double)j) * epsilon_k;
	APPLY(P, v, eta, Av);
	rj = f - Av;
 	cout << "current residual error: " << l2_norm(rj) << endl;
	APPLY(P, rj, eta/5.0, Arj);
// 	cout << "rj * rj = " << rj * rj << endl;
// 	cout << "rj * Arj = " << rj * Arj << endl;
	const double alphaj = (rj * rj) / (rj * Arj);
// 	cout << "descent parameter alphaj=" << alphaj << endl;
	v += alphaj * rj;
      }
      v.COARSE(2.0*epsilon_k/5.0, u_epsilon);
      epsilon_k /= 2.0;
    }
  }
}
