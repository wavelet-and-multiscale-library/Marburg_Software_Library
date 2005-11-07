// implementation for onestep.h

#include <cmath>
#include <iostream>

using std::cout;
using std::endl;

namespace MathTL
{
  template <class VECTOR>
  OneStepScheme<VECTOR>::~OneStepScheme()
  {
  }

  template <class VECTOR>
  void solve_IVP(const AbstractIVP<VECTOR>* ivp,
		 const OneStepScheme<VECTOR>* scheme,
		 const double T,
		 const double tolerance,
		 const double q,
		 const double tau_max,
		 IVPSolution<VECTOR>& result)
  {
    // IVP solver version from Deuflhard/Bornemann, Numerik II, Alg. 5.2 (p. 208)

    result.t.clear();
    result.u.clear();

    double t_m = 0;
    VECTOR u_m(ivp->u0);

    unsigned int m = 0;
    result.t.push_back(t_m);
    result.u.push_back(u_m);
    cout << "t_{" << m << "}=" << t_m << " accepted!" << endl;
    m++;

    const double rho = 0.9; // safety factor

    double tau_m = (T - t_m)/10.0; // first guess for the time stepsize
    
    VECTOR u_mplus1, error_estimate;
    while (t_m < T)
      {
	// try to advance one step
	scheme->increment(ivp, t_m, u_m, tau_m, u_mplus1, error_estimate);
	
	double tau = std::min(pow(rho*tolerance/l2_norm(error_estimate),1./(scheme->order()+1))*tau_m,
			      std::min(q*tau_m,tau_max));
	
	if (l2_norm(error_estimate) <= tolerance)
	  {
	    // accept the time step
	    t_m += tau_m;
	    result.t.push_back(t_m);
	    u_m = u_mplus1;
	    result.u.push_back(u_m);

	    cout << "t_{" << m << "}=" << t_m << " accepted!" << endl;
	    
	    tau_m = std::min(tau, T-t_m);

	    m++;
	  }
	else
	  {
	    // reject the time step
	    cout << "t_{" << m << "}=" << t_m+tau_m << " rejected!" << endl;
	    tau_m = tau;
	  }
      }
  }
}
