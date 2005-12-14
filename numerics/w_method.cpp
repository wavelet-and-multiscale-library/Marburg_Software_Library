// implementation for w_method.h

#include <cmath>

namespace MathTL
{
  template <class VECTOR>
  WMethodStageEquationHelper<VECTOR>::~WMethodStageEquationHelper()
  {
  }

  template <class VECTOR>
  WMethod<VECTOR>::WMethod(const Method method,
			   const WMethodStageEquationHelper<VECTOR>* s)
    : stage_equation_helper(s), preprocessor(0)
  {
    switch(method)
      {
      case GRK4T:
	A.resize(4,4); // corresponds to the [HW] notation
	A(1,0) = 0.2000000000000000e+01;
        A(2,0) = A(3,0) = 0.4524708207373116e+01;
        A(2,1) = A(3,1) = 0.4163528788597648e+01;
	
	C.resize(4,4); // corresponds to the [HW] notation
	C(0,0) = C(1,1) = C(2,2) = C(3,3) = 0.231; // gamma
	C(1,0) = -0.5071675338776316e+01;
	C(2,0) = 0.6020152728650786e+01;
 	C(2,1) = 0.1597506846727117e+00;
 	C(3,0) = -0.1856343618686113e+01;
 	C(3,1) = -0.8505380858179826e+01;
 	C(3,2) = -0.2084075136023187e+01;
	
	gamma_vector.resize(4); // corresponds to Di in [HW]
	gamma_vector[0] = 0.2310000000000000e+00;
        gamma_vector[1] = -0.3962966775244303e-01;
	gamma_vector[2] = 0.5507789395789127e+00;
        gamma_vector[3] = -0.5535098457052764e-01;

	alpha_vector.resize(4); // corresponds to Ci in [HW]
	alpha_vector[1] = 0.4620000000000000e+00;
	alpha_vector[2] = alpha_vector[3]
	  = 0.8802083333333334e+00;

 	m.resize(4); // corresponds to Bi in [HW] (the m_i in the [HW II] book)
 	m[0] = 0.3957503746640777e+01;
	m[1] = 0.4624892388363313e+01;
	m[2] = 0.6174772638750108e+00;
	m[3] = 0.1282612945269037e+01;

 	e.resize(4);
	e[0] = 0.2302155402932996e+01;
	e[1] = 0.3073634485392623e+01;
	e[2] = -0.8732808018045032e+00;
	e[3] = -0.1282612945269037e+01;

	p = 4;
	break;
      case RODAS3:
	A.resize(4,4);
	A(2,0) = A(3,0) = 2.0;
	A(3,2) = 1.0;

	C.resize(4,4);
	C(0,0) = C(1,1) = C(2,2) = C(3,3) = 0.5; // = gamma
	C(1,0) = 4.0;
	C(2,0) = C(3,0) = 1.0;
	C(2,1) = C(3,1) = -1.0;
	C(3,2) = -8.0/3.0;
	
	gamma_vector.resize(4);
	gamma_vector[0] = 0.5;
	gamma_vector[1] = 1.5;

	alpha_vector.resize(4);
	alpha_vector[2] = alpha_vector[3] = 1.0;

 	m.resize(4);
	m[0] = 2.;
	m[2] = m[3] = 1.0;

 	e.resize(4);
	e[3] = 1.0;

 	p = 3;
 	break;
      case ROS2:
	A.resize(2,2);
	A(1,0) = 2 / (2 + M_SQRT2);
	
	C.resize(2,2);
	C(0,0) = C(1,1) = 1 + M_SQRT1_2; // = gamma
	C(1,0) = -4 / (2 + M_SQRT2);
	
	gamma_vector.resize(2);
	gamma_vector[0] = 1 + M_SQRT1_2;
	gamma_vector[1] = -1 - M_SQRT1_2;

	alpha_vector.resize(2);
	alpha_vector[1] = 1.0;

 	m.resize(2);
 	m[0] = 3 / (2 + M_SQRT2);
	m[1] = 1 / (2 + M_SQRT2);

 	e.resize(2);
	e[0] = e[1] = 1 / (2 + M_SQRT2);

 	p = 2;
 	break;
      case ROS3:
	A.resize(3,3);
	A(1,0) = A(2,0) = 1.0;

	C.resize(3,3);
	C(0,0) = C(1,1) = C(2,2) = 0.43586652150845899941601945119356; // = gamma
	C(1,0) = -1.0156171083877702091975600115544914;
	C(2,0) = 4.0759956452537699824805835358065385;
	C(2,1) = 9.2076794298330791242156818474001883;

	gamma_vector.resize(3);
	gamma_vector[0] = .43586652150845899941601945119356;
	gamma_vector[1] = .24291996454816804366592249683315;
	gamma_vector[2] = 2.18513800276640585115131694858316;

	alpha_vector.resize(3);
	alpha_vector[1] = alpha_vector[2] = .43586652150845899941601945119356;

	m.resize(3);
	m[0] = 1.0;
	m[1] = 6.1697947043828245592553615689728652;
	m[2] = -.42772256543218573326238373806511718;

	e.resize(3);
	e[0] = 0.5;
	e[1] = -2.9079558716805469821718236208016638;
	e[2] = .22354069897811569627360909276198390;

 	p = 3;
 	break;
      case ROWDA3:
	A.resize(3,3);
	A(1,0) = A(2,0) = 1.6059962522264972810030329434237128;

	C.resize(3,3);
	C(0,0) = C(1,1) = C(2,2) = 0.435866521508459; // = gamma
	C(1,0) = .88740444110022637192938576568447141;
	C(2,0) = 23.987479717241899735266914145517441;
	C(2,1) = 5.2637223717664389241008661915179357;

	gamma_vector.resize(3);
	gamma_vector[0] = .43586652150;
	gamma_vector[1] = .6044552840570998;
	gamma_vector[2] = 6.379788799336421;

	alpha_vector.resize(3);
	alpha_vector[1] = alpha_vector[2] = 0.7;

	m.resize(3);
	m[0] = 2.2367270453655677993222954788155933;
	m[1] = 2.2500677310226296985367412611749105;
	m[2] = -.20925140444309304448380305345383128;

	e.resize(3);
	e[0] = .1773708776809320419636338468173019;
	e[1] = 2.0806662990846894958960573746516369;
	e[2] = -.20925140444309304448380305345383128;

 	p = 3;
 	break;
      default:
	break;
      }
  }
  
  template <class VECTOR>
  void
  WMethod<VECTOR>::increment(const AbstractIVP<VECTOR>* ivp,
			     const double t_m, const VECTOR& u_m,
			     const double tau,
			     VECTOR& u_mplus1,
			     VECTOR& error_estimate,
			     const double tolerance) const
  {
    const unsigned int stages = A.row_dimension(); // for readability

    Array1D<VECTOR> u(stages);
    u[0] = u_m; u[0].scale(0.0); // ensures correct size
    for (unsigned int i = 1; i < stages; i++)
      u[i] = u[0];

    VECTOR rhs(u[0]), help(u[0]); // ensures correct size
    
    // solve stage equations (TODO: adjust the tolerances appropriately)
    for (unsigned int i(0); i < stages; i++) {
      // setup i-th right-hand side
      //   f(t_m + \tau * \alpha_i, u^{(m)} + \sum_{j=1}^{i-1} a_{i,j} * u_j)
      //   + \sum_{j=1}^{i-1} \frac{c_{i,j}}{\tau} * u_j
      //   + \tau * \gamma_i * g

      help = u_m;
      for (unsigned int j(0); j < i; j++)
  	help.add(A(i,j), u[j]);
      ivp->evaluate_f(t_m+tau*alpha_vector[i], help, tolerance/(3*stages), rhs);

      if (preprocessor == 0) {
	for (unsigned int j(0); j < i; j++)
	  rhs.add(C(i,j)/tau, u[j]);
      } else {
	// TODO: insert action of the preprocessor here!
      }
	
      stage_equation_helper->approximate_ft(ivp, t_m, u_m, tolerance/(3*stages), help);
      rhs.add(tau*gamma_vector[i], help);
      
      // solve i-th stage equation
      // (\tau*\gamma_{i,i})^{-1}I - T) u_i = rhs
      stage_equation_helper->solve_W_stage_equation(ivp, t_m, u_m, 1./(tau*C(i,i)), rhs, tolerance/(3*stages), u[i]);
    }
    
    // update u^{(m)} -> u^{(m+1)} by the k_i
    u_mplus1 = u_m;
    for (unsigned int i(0); i < stages; i++)
      u_mplus1.add(m[i], u[i]);
    
    // error estimate
    error_estimate = u_m; error_estimate.scale(0.0); // ensures correct size
    for (unsigned int i = 0; i < stages; i++)
      error_estimate.add(e[i], u[i]);
  }
}
