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
      case ROS3:
	alpha_ = LowerTriangularMatrix<double>(3);
	gamma_ = LowerTriangularMatrix<double>(3);
	alpha_(1,0) = alpha_(2,0) = gamma_(0,0) = gamma_(1,1) = gamma_(2,2)
	  = 0.43586652150845899941601945119356;
	gamma_(1,0) = -0.19294655696029095575009695436041;
	gamma_(2,1) = 1.74927148125794685173529749738960;
	b_.resize(3);
	b_(0) = -0.75457412385404315829818998646589;
	b_(1) = 1.94100407061964420292840123379419;
	b_(2) = -0.18642994676560104463021124732829;
	break;
      case ROWDA3:
	alpha_ = LowerTriangularMatrix<double>(3);
	alpha_(1,0) = alpha_(2,0) = 0.7;
	gamma_ = LowerTriangularMatrix<double>(3);
	gamma_(0,0) = gamma_(1,1) = gamma_(2,2) = 0.435866521508459;
	gamma_(1,0) = 0.1685887625570998;
	gamma_(2,0) = 4.943922277836421;
	gamma_(2,1) = 1.0;
	b_.resize(3);
	b_[0] = 0.3197278911564624;
	b_[1] = 0.7714777906171382;
	b_[2] = -0.09120568177360061;
	break;
      case RODAS3:
	alpha_ = LowerTriangularMatrix<double>(4);
	alpha_(2,0) = 1.0;
	alpha_(3,0) = 0.75; alpha_(3,1) = -0.25; alpha_(3,2) = 0.5;
	gamma_ = LowerTriangularMatrix<double>(4);
	gamma_(0,0) = gamma_(1,1) = gamma_(2,2) = gamma_(3,3) = 0.5;
	gamma_(1,0) = 1.0;
	gamma_(2,0) = gamma_(2,1) = -0.25;
	gamma_(3,0) = gamma_(3,1) = 1.0/12.0;
	gamma_(3,2) = -2.0/3.0;
	b_.resize(4);
	b_[0] = 5.0/6.0;
	b_[1] = b_[2] = -1.0/6.0;
	b_[3] = 0.5;
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
  void Rosenbrock<VECTOR, IVP>::increment(IVP& ivp,
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
