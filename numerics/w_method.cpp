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
    : stage_equation_helper(s)
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
// 	alpha_matrix.resize(4,4);
// 	alpha_matrix(2,0) = 1.0;
// 	alpha_matrix(3,0) = 0.75;
// 	alpha_matrix(3,1) = -0.25;
// 	alpha_matrix(3,2) = 0.5;
// 	gamma_matrix.resize(4,4);
// 	gamma_matrix(0,0) = gamma_matrix(1,1)
// 	  = gamma_matrix(2,2) = gamma_matrix(3,3) = 0.5;
// 	gamma_matrix(1,0) = 1.0;
// 	gamma_matrix(2,0) = gamma_matrix(2,1) = -0.25;
// 	gamma_matrix(3,0) = gamma_matrix(3,1) = 1.0/12.0;
// 	gamma_matrix(3,2) = -2.0/3.0;
// 	b.resize(4);
// 	b[0] = 5.0/6.0;
// 	b[1] = b[2] = -1.0/6.0;
// 	b[3] = 0.5;
// 	bhat.resize(4);
// 	bhat[0] = 0.75;
// 	bhat[1] = -0.25;
// 	bhat[2] = 0.5;
// 	p = 3;
 	break;
      case ROS2:
	A.resize(2,2); // corresponds to the [HW] notation
	A(1,0) = 2 / (2 + M_SQRT2);
	
	C.resize(2,2); // corresponds to the [HW] notation
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
//       case ROS3:
// 	alpha_matrix.resize(3,3);
// 	gamma_matrix.resize(3,3);
// 	alpha_matrix(1,0) = alpha_matrix(2,0)
// 	  = gamma_matrix(0,0) = gamma_matrix(1,1) = gamma_matrix(2,2)
// 	  = 0.43586652150845899941601945119356;
// 	gamma_matrix(1,0) = -0.19294655696029095575009695436041;
// 	gamma_matrix(2,1) = 1.74927148125794685173529749738960;
// 	b.resize(3);
// 	b[0] = -0.75457412385404315829818998646589;
// 	b[1] = 1.94100407061964420292840123379419;
// 	b[2] = -0.18642994676560104463021124732829;
// 	bhat.resize(3);
// 	bhat[0] = -1.53358745784149585370766523913002;
// 	bhat[1] = 2.81745131148625772213931745457622;
// 	bhat[2] = -0.28386385364476186843165221544619;
// 	p = 3;
// 	break;
//       case ROWDA3:
// 	alpha_matrix.resize(3,3);
// 	alpha_matrix(1,0) = alpha_matrix(2,0) = 0.7;
// 	gamma_matrix.resize(3,3);
// 	gamma_matrix(0,0) = gamma_matrix(1,1)
// 	  = gamma_matrix(2,2) = 0.435866521508459;
// 	gamma_matrix(1,0) = 0.1685887625570998;
// 	gamma_matrix(2,0) = 4.943922277836421;
// 	gamma_matrix(2,1) = 1.0;
// 	b.resize(3);
// 	b[0] = 0.3197278911564624;
// 	b[1] = 0.7714777906171382;
// 	b[2] = -0.09120568177360061;
// 	bhat.resize(3);
// 	bhat[0] = 0.926163587124091;
// 	bhat[1] = 0.073836412875909;
// 	p = 3;
// 	break;
      default:
	break;
      }

//     const unsigned int stages = alpha_matrix.row_dimension();
    
//     alpha_vector.resize(stages);
//     for (unsigned int i = 0; i < stages; i++) {
//       double help = 0;
//       for (unsigned int j = 0; j < i; j++)
// 	help += alpha_matrix(i,j);
//       alpha_vector[i] = help; // alpha_i = sum_{j=1}^{i-1} alpha_{i,j}
//     }
    
//     gamma_vector.resize(stages);
//     for (unsigned int i = 0; i < stages; i++) {
//       double help = 0;
//       for (unsigned int j = 0; j <= i; j++)
// 	help += gamma_matrix(i,j);
//       gamma_vector[i] = help; // gamma_i = sum_{j=1}^i gamma_{i,j}
//     }
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
    u[0] = u_m; u[0] = 0; // ensures correct size
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
      ivp->evaluate_f(t_m+tau*alpha_vector[i], help, tolerance/stages, rhs);
      
      for (unsigned int j(0); j < i; j++)
	rhs.add(C(i,j)/tau, u[j]);
      
      stage_equation_helper->approximate_ft(ivp, t_m, u_m, tolerance/stages, help);
      rhs.add(tau*gamma_vector[i], help);
      
      // solve i-th stage equation
      // (\tau*\gamma_{i,i})^{-1}I - T) u_i = rhs
      stage_equation_helper->solve_W_stage_equation(ivp, t_m, u_m, 1./(tau*C(i,i)), rhs, tolerance/stages, u[i]);
    }
    
    // update u^{(m)} -> u^{(m+1)} by the k_i
    u_mplus1 = u_m;
    for (unsigned int i(0); i < stages; i++)
      u_mplus1.add(m[i], u[i]);
    
    // error estimate
    error_estimate = u_m; error_estimate = 0; // ensures correct size
    for (unsigned int i = 0; i < stages; i++)
      error_estimate.add(e[i], u[i]);
  }
}
