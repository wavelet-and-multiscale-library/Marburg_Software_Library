// implementation for w_method.h

#include <cmath>

namespace MathTL
{
  template <class VECTOR>
  WMethod<VECTOR>::WMethod(const Method method,
			   const WMethodStageEquationSolver<VECTOR>& s)
    : stage_equation_solver(s)
  {
    switch(method)
      {
      case ROS2:
	alpha_matrix = LowerTriangularMatrix<double>(2);
	alpha_matrix(1,0) = 1.0;
	gamma_matrix = LowerTriangularMatrix<double>(2);
	gamma_matrix(0,0) = gamma_matrix(1,1) = 1 + M_SQRT1_2;
	gamma_matrix(1,0) = -2*(1 + M_SQRT1_2);
	b.resize(2);
	b[0] = b[1] = 0.5;
	bhat.resize(2);
	bhat[0] = 1.0;
	break;
      default:
	break;
      }

    const unsigned int stages = alpha_matrix.row_dimension();
    
    alpha_vector.resize(s);
    for (unsigned int i = 0; i < stages; i++) {
      double help = 0;
      for (unsigned int j = 0; j < i; j++)
	help += alpha_matrix(i, j);
      alpha_vector[i] = help; // alpha_i = sum_{j=1}^{i-1} alpha_{i,j}
    }
    
    gamma_vector.resize(stages);
    for (unsigned int i = 0; i < stages; i++) {
      double help = 0;
      for (unsigned int j = 0; j <= i; j++)
	help += gamma_matrix(i, j);
      gamma_vector[i] = help; // gamma_i = sum_{j=1}^i gamma_{i,j}
    }
  }
  
  template <class VECTOR>
  void
  WMethod<VECTOR>::increment(const AbstractIVP<VECTOR>& ivp,
			     const double t_m, const VECTOR& u_m,
			     const double tau,
			     VECTOR& u_mplus1,
			     VECTOR& error_estimate,
			     const double tolerance) const
  {
  }
}
