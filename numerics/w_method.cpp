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
      case ROS2:
	alpha_matrix.resize(2,2);
	alpha_matrix(1,0) = 1.0;
	gamma_matrix.resize(2,2);
	gamma_matrix(0,0) = gamma_matrix(1,1) = 1 + M_SQRT1_2;
	gamma_matrix(1,0) = -2*(1 + M_SQRT1_2);
	b.resize(2);
	b[0] = b[1] = 0.5;
	bhat.resize(2);
	bhat[0] = 1.0;
	p = 2;
	break;
      case ROWDA3:
	alpha_matrix.resize(3,3);
	alpha_matrix(1,0) = alpha_matrix(2,0) = 0.7;
	gamma_matrix.resize(3,3);
	gamma_matrix(0,0) = gamma_matrix(1,1)
	  = gamma_matrix(2,2) = 0.435866521508459;
	gamma_matrix(1,0) = 0.1685887625570998;
	gamma_matrix(2,0) = 4.943922277836421;
	gamma_matrix(2,1) = 1.0;
	b.resize(3);
	b[0] = 0.3197278911564624;
	b[1] = 0.7714777906171382;
	b[2] = -0.09120568177360061;
	bhat.resize(3);
	bhat[0] = 0.926163587124091;
	bhat[1] = 0.073836412875909;
	p = 3;
	break;
      case ROS3:
	alpha_matrix.resize(3,3);
	gamma_matrix.resize(3,3);
	alpha_matrix(1,0) = alpha_matrix(2,0)
	  = gamma_matrix(0,0) = gamma_matrix(1,1) = gamma_matrix(2,2)
	  = 0.43586652150845899941601945119356;
	gamma_matrix(1,0) = -0.19294655696029095575009695436041;
	gamma_matrix(2,1) = 1.74927148125794685173529749738960;
	b.resize(3);
	b[0] = -0.75457412385404315829818998646589;
	b[1] = 1.94100407061964420292840123379419;
	b[2] = -0.18642994676560104463021124732829;
	bhat.resize(3);
	bhat[0] = -1.53358745784149585370766523913002;
	bhat[1] = 2.81745131148625772213931745457622;
	bhat[2] = -0.28386385364476186843165221544619;
	p = 3;
	break;
      case RODAS3:
	alpha_matrix.resize(4,4);
	alpha_matrix(2,0) = 1.0;
	alpha_matrix(3,0) = 0.75;
	alpha_matrix(3,1) = -0.25;
	alpha_matrix(3,2) = 0.5;
	gamma_matrix.resize(4,4);
	gamma_matrix(0,0) = gamma_matrix(1,1)
	  = gamma_matrix(2,2) = gamma_matrix(3,3) = 0.5;
	gamma_matrix(1,0) = 1.0;
	gamma_matrix(2,0) = gamma_matrix(2,1) = -0.25;
	gamma_matrix(3,0) = gamma_matrix(3,1) = 1.0/12.0;
	gamma_matrix(3,2) = -2.0/3.0;
	b.resize(4);
	b[0] = 5.0/6.0;
	b[1] = b[2] = -1.0/6.0;
	b[3] = 0.5;
	bhat.resize(4);
	bhat[0] = 0.75;
	bhat[1] = -0.25;
	bhat[2] = 0.5;
	p = 3;
	break;
      default:
	break;
      }

    const unsigned int stages = alpha_matrix.row_dimension();
    
    alpha_vector.resize(stages);
    for (unsigned int i = 0; i < stages; i++) {
      double help = 0;
      for (unsigned int j = 0; j < i; j++)
	help += alpha_matrix(i,j);
      alpha_vector[i] = help; // alpha_i = sum_{j=1}^{i-1} alpha_{i,j}
    }
    
    gamma_vector.resize(stages);
    for (unsigned int i = 0; i < stages; i++) {
      double help = 0;
      for (unsigned int j = 0; j <= i; j++)
	help += gamma_matrix(i,j);
      gamma_vector[i] = help; // gamma_i = sum_{j=1}^i gamma_{i,j}
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
    const unsigned int stages = alpha_matrix.row_dimension(); // for readability

    Array1D<VECTOR> k(stages);
    k[0] = u_m; k[0] = 0; // ensures correct size
    for (unsigned int i = 1; i < stages; i++)
      k[i] = k[0];

    VECTOR rhs1(k[0]), help(k[0]), uhelp(k[0]); // ensures correct size
    
    // solve stage equations (TODO: adjust the tolerances appropriately)
    for (unsigned int i(0); i < stages; i++) {
      uhelp = u_m;
      for (unsigned int j(0); j < i; j++)
	uhelp.add(tau*alpha_matrix(i,j), k[j]);
      ivp->evaluate_f(t_m+tau*alpha_vector[i], uhelp, tolerance/stages, rhs1);
      
      ivp->evaluate_ft(t_m, u_m, tolerance/stages, help);
      rhs1.add(tau*gamma_vector[i], help);
      
      uhelp = 0;
      for (unsigned int j(0); j < i; j++)
	uhelp.add(gamma_matrix(i,j)/gamma_matrix(i,i), k[j]);
      rhs1.add(uhelp);
      
      ivp->solve_ROW_stage_equation(t_m, u_m, tau*gamma_matrix(i,i), rhs1, tolerance/stages, help);
      k[i] = help - uhelp;
    }
    
    // update u^{(m)} -> u^{(m+1)} by the k_i
    u_mplus1 = u_m;
    for (unsigned int i(0); i < stages; i++)
      u_mplus1.add(tau*b[i], k[i]);
    
    // error estimate
    error_estimate = u_m; error_estimate = 0; // ensures correct size
    for (unsigned int i = 0; i < stages; i++)
      error_estimate.add(tau*(b[i]-bhat[i]), k[i]);
  }
}
