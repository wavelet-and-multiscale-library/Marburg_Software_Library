// implementation for steepest_descent.h

#include <cmath>
#include <set>

#include <adaptive/apply.h>

using std::set;

namespace FrameTL
{
  template <class PROBLEM>
  void steepest_descent_SOLVE(const PROBLEM& P,  const double epsilon,
			      InfiniteVector<double, typename PROBLEM::Index>& u_epsilon)
  {
    unsigned int loops = 0;
    const int jmax = 3;
    typedef typename PROBLEM::Index Index;

    double a_inv     = P.norm_Ainv();
    double kappa     = P.norm_A()/a_inv;
    double omega_i   = a_inv*P.F_norm();
    cout << "omega_i = " << omega_i << endl;
    double delta     = 1./(5.*kappa+a_inv);
    cout << "delta = " << delta << endl;
    const double A = 1 + delta;
    const double C = 1.0 / ((1 - ((kappa*(delta*delta+2.*delta)+a_inv*delta)/((1-delta)*(1-delta))))
      * (((1-delta)*(1-delta))/(a_inv)));
    cout << "C = " << C << endl;
    const double B = C * (A*A);
    cout << "B = " << B << endl;
    //double lambda = (kappa-1)/(kappa+1) + P.norm_A()*std::max(3.*A*A*B,C*(1./(1-delta)))*delta;
    double lambda = (kappa-1)/(kappa+1)+0.0001;
    cout << "lambda = " << lambda << endl;
    const double C3 = B;
    cout << "C3 = " << C3 << endl;

    double mu        = 2; //shall be > 1
    //beta in (0,1)
    double beta      = 0.5;
    //let K be such that beta^K * omega <= epsilon
    unsigned int K   = (int) (log(epsilon/omega_i) / log(beta) + 1);
    //let M be such that lambda^M <= ((1-delta) / (1+delta)) * (beta / ((1+3*mu)*kappa))
    int M            = (int) ((log( ((1-delta)/(1+delta)) * (beta / ((1+3.0*mu)*kappa)) )
			       / log(lambda)) + 1);
    

    InfiniteVector<double, Index> w, tilde_r;
    
    cout << "K = " << K << endl;
    cout << "M = " << M << endl;

    for (unsigned int i = 1; i < K; i++) {
      //omega_i *= 0.000001;/*beta;*/
      omega_i *= beta;
      double xi_i = omega_i / ((1+3.0*mu)*C3*M);
      cout << "xi = " << xi_i << endl;
      double nu_i = 0.;
      RES(P, w, xi_i, delta, omega_i/((1+3.*mu)*a_inv), jmax,
	  tilde_r, nu_i, CDD1);
      //cout << tilde_r << endl;
      while ( nu_i > omega_i/((1+3.*mu)*a_inv)) {
	InfiniteVector<double, Index> z_i;
	cout << "apply tol = " << delta*l2_norm(tilde_r) << endl;
	APPLY(P, tilde_r, delta*l2_norm(tilde_r), z_i, jmax, CDD1);
	cout << "after apply " << endl; 
	//APPLY_COARSE(P, tilde_r, delta*l2_norm(tilde_r), z_i, 0.000001, jmax, CDD1);
	w += (l2_norm_sqr(tilde_r)/(z_i*tilde_r))*tilde_r;
	cout << "xi = " << xi_i << endl;
	RES(P, w, xi_i, delta, omega_i/((1+3.*mu)*a_inv), jmax,
	    tilde_r, nu_i, CDD1);
	cout << "loop: " << ++loops << " nu = " << nu_i << endl;
	cout << "active indices: " << w.size() << endl;
	if (loops==20) {
	  u_epsilon = w;
	  return;
	}
      }
      InfiniteVector<double, Index> tmp;
      w.COARSE(((3.*mu*omega_i)/(1+3.*mu)),tmp);
      w = tmp;
    }
    
    //condition number

  }
}
