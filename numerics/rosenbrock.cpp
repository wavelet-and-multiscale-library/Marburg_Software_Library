#include <vector>

namespace MathTL
{
  template <class VECTOR, class IVP>
  Rosenbrock<VECTOR, IVP>::Rosenbrock()
    : alpha_(1), gamma_(1), b_(1)
  {
    gamma_(0, 0) = 1;
    b_[0] = 1;
  }

  template <class VECTOR, class IVP>
  Rosenbrock<VECTOR, IVP>::Rosenbrock(const LowerTriangularMatrix<double>& alpha,
				 const LowerTriangularMatrix<double>& gamma,
				 const Vector<double>& b)
    : alpha_(alpha), gamma_(gamma), b_(b)
  {
  }

  template <class VECTOR, class IVP>
  void Rosenbrock<VECTOR, IVP>::increment(const IVP& ivp,
					  const double t_m, const VECTOR& u_m,
					  const double tau,
					  VECTOR& umplus1) const
  {
#if 0
    // for a first test, we only solve the first stage equation
    VECTOR k1(u_m), rhs1(u_m), rhs2(u_m); // just to ensure the same size

    ivp.apply_f(t_m+tau*alpha_(0,0), u_m, rhs1);
    ivp.apply_ft(t_m, u_m, rhs2);
    rhs2 *= tau*gamma_(0,0);
    rhs1 += rhs2;

    ivp.solve_jacobian(t_m, u_m, tau*gamma_(0,0), k1);

    k1 *= tau;
    umplus1 = u_m + k1;
#else
    unsigned int s(b_.size()); // number of stages, for simplicity
    std::vector<VECTOR> k(s);
    VECTOR rhs1(u_m), rhs2(u_m); // ensure the correct size
    for (unsigned int i(0); i < s; i++)
      {
      }
    
    
#endif
  }
}
