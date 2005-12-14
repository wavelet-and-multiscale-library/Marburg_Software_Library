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
#if 1
    // IVP solver a la Hairer/Wanner
    result.t.clear();
    result.u.clear();
    
    double t_m = 0;
    VECTOR u_m(ivp->u0);
    
    unsigned int m = 0;
    result.t.push_back(t_m);
    result.u.push_back(u_m);

#if _MATHTL_ONESTEPSCHEME_VERBOSITY >= 1
    cout << "t_{" << m << "}=" << t_m << " accepted!" << endl;
#endif

    m++;
    
    const double rho = 0.8; // overall safety factor

    // guess the initial time stepsize (cf. Hairer/Wanner, p. 169)
    double u0_norm = linfty_norm(ivp->u0);
    VECTOR yp0(ivp->u0);
    ivp->evaluate_f(0, ivp->u0, atol, yp0); // initial slope
    double ft0u0_norm = linfty_norm(yp0);
    double tau0 = (u0_norm < 1e-5 || ft0u0_norm < 1e-5)
      ? 1e-6 : 1e-2*u0_norm/ft0u0_norm;

    // TODO: approximate ypp
    
    double tau_m = std::min(tau_max, tau0); // 100*tau0);
    
    VECTOR u_mplus1, error_estimate;
    bool done = false;

    // parameters for step size selection;
    const double fac1 = 5.0;
    const double fac2 = 1./q;
    double fac = 0, facgus = 0, hacc = 0, erracc = 0; // will be set in the loop

    while (!done)
      {
 	// jump to T if within 10% of T-t_m
 	if (1.1*tau_m >= fabs(T-t_m) || t_m+tau_m >= T)
 	  {
 	    tau_m = T-t_m;
 	    done = true; // we would be "done" if the step is accepted
 	  }

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
//  	errest = linfty_norm(error_estimate) / (atol+maxnorm*rtol); // not good

 	// estimate new stepsize
	fac = std::max(fac2, std::min(fac1, pow(errest, 1./scheme->order()) / rho));
	double tau_new = tau_m / fac;
	
	if (errest <= 1)
	  {
	    // accept the time step

#if _MATHTL_ONESTEPSCHEME_VERBOSITY >= 1
	    cout << "t_{" << m << "}=" << t_m+tau_m << " accepted!" << endl;
#endif

	    t_m += tau_m;
	    result.t.push_back(t_m);
	    u_m = u_mplus1;
	    result.u.push_back(u_m);

	    // predictive controller of Gustafsson
	    if (m >= 2) {
	      facgus = (hacc/tau_m) * pow(errest*errest/erracc, 1./scheme->order()) / rho;
	      facgus = std::max(fac2, std::min(fac1, facgus));
	      fac = std::max(fac, facgus);
	    }
	    hacc = tau_m;
	    erracc = std::max(1e-2, errest);

	    tau_m = tau_new;
	    
	    m++;
	  }
	else
	  {
	    // reject the time step

#if _MATHTL_ONESTEPSCHEME_VERBOSITY >= 1
	    cout << "t_{" << m << "}=" << t_m+tau_m << " rejected!"
		 << " (errest=" << errest << ")" << endl;
#endif
	    
	    tau_m = tau_new;
	    
	    done = false;
	  }
      }

    



//     // IVP solver a la Hairer/Wanner
//     result.t.clear();
//     result.u.clear();
    
//     double t_m = 0;
//     VECTOR u_m(ivp->u0);
    
//     unsigned int m = 0;
//     result.t.push_back(t_m);
//     result.u.push_back(u_m);

// #if _MATHTL_ONESTEPSCHEME_VERBOSITY >= 1
//     cout << "t_{" << m << "}=" << t_m << " accepted!" << endl;
// #endif

//     m++;
    
//     const double rho = 0.8; // safety factor

//     // guess the initial time stepsize (cf. Hairer/Wanner, p. 169)
//     double u0_norm = linfty_norm(ivp->u0);
//     VECTOR yp0(ivp->u0);
//     ivp->evaluate_f(0, ivp->u0, atol, yp0); // initial slope
//     double ft0u0_norm = linfty_norm(yp0);
//     double tau0 = (u0_norm < 1e-5 || ft0u0_norm < 1e-5)
//       ? 1e-6 : 1e-2*u0_norm/ft0u0_norm;

//     // TODO: approximate ypp
    
//     double tau_m = 100*tau0;

//     m++;
    
//     bool done = false;
//     VECTOR u_mplus1, error_estimate;
//     while (!done)
//       {
// 	double tau_min = 1e-10 * t_m;
// 	tau_m = std::min(tau_max, std::max(tau_min, tau_m)); // this will be the trial stepsize

// 	// jump to T if within 10% of T-t_m
// 	if (1.1*tau_m >= fabs(T-t_m))
// 	  {
// 	    tau_m = T-t_m;
// 	    done = true; // we would be "done" if the step is accepted
// 	  }

// 	// try to advance one step
// 	bool nofailed = true;
// 	while(true) {
// 	  scheme->increment(ivp, t_m, u_m, tau_m, u_mplus1, error_estimate);

// 	  // compute an error estimator
// 	  const double u_m_norm = linfty_norm(u_m);
// 	  const double u_mplus1_norm = linfty_norm(u_mplus1);
// 	  const double maxnorm = std::max(u_m_norm, u_mplus1_norm);
// 	  double errest = 0;
// 	  for (typename VECTOR::const_iterator it(error_estimate.begin());
// 	       it != error_estimate.end(); ++it)
// 	    errest += ((*it * *it)/((atol+maxnorm*rtol)*(atol+maxnorm*rtol)));
// 	  errest = sqrt(errest / error_estimate.size());
// //  	errest = linfty_norm(error_estimate) / (atol+maxnorm*rtol);
	  
// 	  if (errest > 1)
// 	    {
// 	      // failed step
	      
// #if _MATHTL_ONESTEPSCHEME_VERBOSITY >= 1
// 	      cout << "t_{" << m << "}=" << t_m+tau_m << " rejected!" << endl;
// #endif
	      
// 	      nofailed = false;
// 	      tau_m = std::max(tau_min, tau_m * std::max(0.1, rho * pow(1./errest, 1./scheme->order())));
// 	      done = false;
// 	    }
// 	  else
// 	    {
// 	      // accept the time step
	      
// 	      m++;

// #if _MATHTL_ONESTEPSCHEME_VERBOSITY >= 1
// 	      cout << "t_{" << m << "}=" << t_m+tau_m << " accepted!" << endl;
// #endif
	      
// 	      t_m += tau_m;
// 	      result.t.push_back(t_m);
// 	      u_m = u_mplus1;
// 	      result.u.push_back(u_m);

// 	      if (done)
// 		{
// 		  t_m = T; // adjust t_m
// 		}
// 	      else
// 		{
// 		  // guess tau_m for the next step
// 		  tau_m = std::max(tau_min, tau_m * std::max(0.1, rho * pow(1./errest, 1./scheme->order())));
// 		}
	      
// 	      break;
// 	    }
// 	}
//       }




#else
    // IVP solver version from Deuflhard/Bornemann, Numerik II, Alg. 5.2 (p. 208)
    
    result.t.clear();
    result.u.clear();
    
    double t_m = 0;
    VECTOR u_m(ivp->u0);
    
    unsigned int m = 0;
    result.t.push_back(t_m);
    result.u.push_back(u_m);

#if _MATHTL_ONESTEPSCHEME_VERBOSITY >= 1
    cout << "t_{" << m << "}=" << t_m << " accepted!" << endl;
#endif

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
	
	if (linfty_norm(error_estimate) <= atol) // + std::max(linfty_norm(u_m), linfty_norm(u_mplus1)) * rtol)
	  {
	    // accept the time step
	    t_m += tau_m;
	    result.t.push_back(t_m);
	    u_m = u_mplus1;
	    result.u.push_back(u_m);

#if _MATHTL_ONESTEPSCHEME_VERBOSITY >= 1
	    cout << "t_{" << m << "}=" << t_m << " accepted!" << endl;
#endif
	    
	    tau_m = std::min(tau, T-t_m);

	    m++;
	  }
	else
	  {
	    // reject the time step

#if _MATHTL_ONESTEPSCHEME_VERBOSITY >= 1
	    cout << "t_{" << m << "}=" << t_m+tau_m << " rejected!" << endl;
#endif

	    tau_m = std::min(tau, T-t_m);
	    // tau_m = tau;
	  }
      }
#endif
  }
}
