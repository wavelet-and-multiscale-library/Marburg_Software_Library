#include <cmath>
#include <vector>

namespace MathTL
{
  template <class VECTOR, class IVP>
  Rosenbrock<VECTOR, IVP>::Rosenbrock(const Method method)
    : method_(method)
  {
    switch(method_)
      {
      case ROS2:
	alpha_ = LowerTriangularMatrix<double>(2);
	alpha_(1,0) = 1.0;
	gamma_ = LowerTriangularMatrix<double>(2);
	gamma_(0,0) = gamma_(1,1) = 1 + M_SQRT1_2;
	gamma_(1,0) = -2*(1 + M_SQRT1_2);
	b_.resize(2);
	b_[0] = b_[1] = 0.5;
	break;
      case Euler:
      default:
	alpha_ = LowerTriangularMatrix<double>(1, 1, "0");
	gamma_ = LowerTriangularMatrix<double>(1, 1, "1");
	b_.resize(1);
	b_[0] = 1.0;
	break;
      }
  }

  template <class VECTOR, class IVP>
  void Rosenbrock<VECTOR, IVP>::increment(const IVP& ivp,
					  const double t_m, const VECTOR& u_m,
					  const double tau,
					  VECTOR& u_mplus1) const
  {
    unsigned int s(b_.size()); // number of stages, for readability
    std::vector<VECTOR> k(s);
    VECTOR rhs1(u_m), help(u_m), uhelp(u_m); // ensure the correct size

    // precompute the coefficients alpha_i = sum_{j=0}^{i-1}alpha_{i,j}
    Vector<double> alphai(s);
    for (unsigned int i(0); i < s; i++)
      {
	for (unsigned int j(0); j < i; j++)
	  alphai[i] += alpha_(i,j);
      }

    // precompute the coefficients gamma_i = sum_{j=0}^i\gamma_{i,j}
    Vector<double> gammai(s);
    for (unsigned int i(0); i < s; i++)
      {
	for (unsigned int j(0); j <= i; j++)
	  gammai[i] += gamma_(i,j);
      }
    
    // compute the stage solutions k_i
    for (unsigned int i(0); i < s; i++)
      {
	uhelp = u_m;
	for (unsigned int j(0); j < i; j++)
	  uhelp.add(tau*alpha_(i,j), k[j]);
  	ivp.apply_f(t_m+tau*alphai[i], uhelp, rhs1);

	ivp.apply_ft(t_m, u_m, help);
	rhs1.add(tau*gammai[i], help);

	uhelp = 0;
	for (unsigned int j(0); j < i; j++)
	  uhelp.add(gamma_(i,j)/gamma_(i,i), k[j]);
	rhs1 += uhelp;
	
	ivp.solve_jacobian(t_m, rhs1, tau*gamma_(i, i), help);

	k[i] = help - uhelp;
      }

    // update u^{(m)} -> u^{(m+1)} by the k_i
    u_mplus1 = u_m;
    for (unsigned int i(0); i < s; i++)
      {
	u_mplus1.add(tau*b_[i], k[i]);
      }
  }
}
