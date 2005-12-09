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
		 const double atol,
		 const double rtol,
		 const double q,
		 const double tau_max,
		 IVPSolution<VECTOR>& result)
  {
#if 0
    // IVP solver a la Hairer/Wanner
    result.t.clear();
    result.u.clear();
    
    double t_m = 0;
    VECTOR u_m(ivp->u0);
    
    unsigned int m = 0;
    result.t.push_back(t_m);
    result.u.push_back(u_m);
    cout << "t_{" << m << "}=" << t_m << " accepted!" << endl;
    m++;
    
    const double rho = 0.5; // safety factor
    
    double tau_m = (T - t_m)/10.0; // first crude guess for the time stepsize (TODO...)
    
    VECTOR u_mplus1, error_estimate;
    while (t_m < T)
      {
	// try to advance one step
	scheme->increment(ivp, t_m, u_m, tau_m, u_mplus1, error_estimate);
	
	// compute an error estimator
	const double u_m_norm = linfty_norm(u_m);
	const double u_mplus1_norm = linfty_norm(u_mplus1);
	const double maxnorm = std::max(u_m_norm, u_mplus1_norm);
	double errest = 0;
	for (typename VECTOR::const_iterator it(error_estimate.begin());
	     it != error_estimate.end(); ++it)
	  errest += ((*it * *it)/((atol+maxnorm*rtol)*(atol+maxnorm*rtol)));
	errest = sqrt(errest / error_estimate.size());
	
	// estimate new stepsize
	double tau = std::min(tau_max, tau_m * std::min(q, rho*pow(1./errest, 1./(scheme->order()+1))));
	
	if (errest <= 1)
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
#else
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
    
    const double rho = 0.8; // safety factor
    
    double tau_m = (T - t_m)/10.0; // first guess for the time stepsize
    
    VECTOR u_mplus1, error_estimate;
    while (t_m < T)
      {
	// try to advance one step
	scheme->increment(ivp, t_m, u_m, tau_m, u_mplus1, error_estimate);
	
	double tau = std::min(pow(rho*atol/l2_norm(error_estimate),1./(scheme->order()+1))*tau_m,
			      std::min(q*tau_m,tau_max));
	
	if (l2_norm(error_estimate) <= atol) // + std::max(linfty_norm(u_m), linfty_norm(u_mplus1)) * rtol)
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
#endif
  }
}
