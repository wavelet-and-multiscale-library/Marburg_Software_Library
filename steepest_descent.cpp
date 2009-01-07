// implementation for steepest_descent.h

#include <cmath>
#include <set>
#include <utils/plot_tools.h>
#include <adaptive/apply.h>
#include <numerics/corner_singularity.h>
#include <frame_evaluate.h>

using std::set;

namespace FrameTL{

/*!
*/
template<class VALUE = double>
class Singularity1D_RHS_2
  : public Function<1, VALUE>
{
public:
  Singularity1D_RHS_2() {};
  virtual ~Singularity1D_RHS_2() {};
  VALUE value(const Point<1>& p,
	      const unsigned int component = 0) const
  {
    return -sin(3.*M_PI*p[0])*9.*M_PI*M_PI - 4.;
  }
  
  void vector_value(const Point<1> &p,
		    Vector<VALUE>& values) const { ; }
  
};

/*!
  special function with steep gradients
  near the right end of the interval
*/
template<class VALUE = double>
class Singularity1D_2
  : public Function<1, VALUE>
{
public:
  Singularity1D_2() {};
  virtual ~Singularity1D_2() {};
  VALUE value(const Point<1>& p,
	      const unsigned int component = 0) const
  {
    if (0. <= p[0] && p[0] < 0.5)
      return -sin(3.*M_PI*p[0]) + 2.*p[0]*p[0];

    if (0.5 <= p[0] && p[0] <= 1.0)
      return -sin(3.*M_PI*p[0]) + 2.*(1-p[0])*(1-p[0]);

    return 0.;

  }
  
  void vector_value(const Point<1> &p,
		    Vector<VALUE>& values) const { ; }
  
};



  template <class PROBLEM>
  void steepest_descent_SOLVE(const PROBLEM& P,  const double epsilon,
			      InfiniteVector<double, typename PROBLEM::Index>& u_epsilon,
			      Array1D<InfiniteVector<double, typename PROBLEM::Index> >& approximations)
  {
    //typedef DSBasis<3,3> Basis1D;
    typedef PBasis<3,3> Basis1D;

//     Point<2> origin;
//     origin[0] = 0.0;
//     origin[1] = 0.0;
    
//     CornerSingularity sing2D(origin, 0.5, 1.5);
//     CornerSingularityRHS singRhs(origin, 0.5, 1.5);

//    Singularity1D_RHS_2<double> sing1D;
//    Singularity1D_2<double> exactSolution1D;

    unsigned int loops = 0;
    unsigned int niter = 0;

    const int jmax = JMAX;
    typedef typename PROBLEM::Index Index;

    double a_inv     = P.norm_Ainv();
    double kappa     = P.norm_A()*a_inv;
    double omega_i   = a_inv*P.F_norm();
    cout << "a_inv = " << a_inv << endl;
    cout << "omega_i = " << omega_i << endl;

    //double delta     = 1./(5.*kappa+a_inv);
    double delta = 1.;
    cout << "delta = " << delta << endl;

    //const double A = 1 + delta;
    const double A = 1.;
    //const double C = 1.0 / ((1 - ((kappa*(delta*delta+2.*delta)+a_inv*delta)/((1-delta)*(1-delta))))
    //			    * (((1-delta)*(1-delta))/(a_inv)));
    const double C = 1.0;
    cout << "C = " << C << endl;
    const double B = C * (A*A);
    cout << "B = " << B << endl;

    //double lambda = (kappa-1)/(kappa+1) + P.norm_A()*std::max(3.*A*A*B,C*(1./(1-delta)))*delta;
    //double lambda = ((kappa-1)/(kappa+1)+1.)/2.;
    double lambda = 0.95;
    cout << "lambda = " << lambda << endl;

    const double C3 = B;
    cout << "C3 = " << C3 << endl;

    double mu        = 1.0001; //shall be > 1

    //beta in (0,1)
    //double beta      = 0.98;
    //double beta      = 0.85;
    double beta      = 0.85;
    
    cout << "beta = " << beta << endl;

    //let K be such that beta^K * omega <= epsilon
    unsigned int K   = (int) (log(epsilon/omega_i) / log(beta) + 1);
    //let M be such that lambda^M <= ((1-delta) / (1+delta)) * (beta / ((1+3*mu)*kappa))
    int M            = std::max((int) ((log( ((1-delta)/(1+delta)) * (beta / ((1+3.0*mu)*kappa)) )
    			       / log(lambda)) + 1),1);

    cout << "K = " << K << endl;
    cout << "M = " << M << endl;


    InfiniteVector<double, Index> w, tilde_r, help, f, Av;

    map<double,double> log_10_residual_norms;
    map<double,double> degrees_of_freedom;
    map<double,double> asymptotic;
    map<double,double> time_asymptotic;
    map<double,double> descent_params;


    
    bool exit = 0;
    double time = 0.;
    double residual_norm = 5.0;
    clock_t tstart, tend;
    tstart = clock();

    //EvaluateFrame<Basis1D,2,2> evalObj;

    double d = 0.5;

    for (unsigned int i = 1; i <= K; i++) {
      omega_i *= beta;
      double xi_i = omega_i / ((1+3.0*mu)*C3*M);
      double nu_i = 0.;

      RES(P, w, xi_i, delta, omega_i/((1+3.*mu)*a_inv), jmax,
	  tilde_r, nu_i, niter, CDD1); 

      while ( nu_i > omega_i/((1+3.*mu)*a_inv)) {

	InfiniteVector<double, Index> z_i;
	APPLY_COARSE(P, tilde_r, delta*l2_norm(tilde_r), z_i, 1.0e-6, jmax, CDD1);
	//APPLY_COARSE(P, tilde_r, delta*l2_norm(tilde_r), z_i, 0.5, jmax, CDD1);
 	double g = z_i*tilde_r;
	if  (g != 0.)
	  d = (tilde_r*tilde_r)/g;
	
	w += d*tilde_r;
 	InfiniteVector<double, Index> tmp;
 	w.COARSE(1.0/100.0*residual_norm, tmp);
 	w = tmp;


	cout << "descent param = " << d << endl;
	++loops;
	++niter;

	RES(P, w, xi_i, delta, omega_i/((1+3.*mu)*a_inv), jmax,
	    tilde_r, nu_i, niter, CDD1);

	cout << "loop: " << loops << " nu = " 
	     << nu_i << " epsilon = " << omega_i/((1+3.*mu)*a_inv) << endl;
	cout << "xi: " << xi_i << endl; 

	// ############# output #############
	tend = clock();
	time += ((double) (tend-tstart))/((double) CLOCKS_PER_SEC);

  	P.RHS(1.0e-6,f);
 	APPLY(P, w, 1.0e-6, Av, jmax, CDD1);
 	help = f-Av;
	  
	//residual_norm = l2_norm(help);
	double tmp1 = log10(residual_norm);
	cout << "residual norm = " << residual_norm << endl;

	cout << "active indices: " << w.size() << endl;
	asymptotic[log10( (double)w.size() )] = tmp1;
	time_asymptotic[log10(time)] = tmp1;
	
	descent_params[loops] = d;
	


// 	if ((loops <= 10) || ((loops % 10 == 0) && (loops <= 100))
// 	    || ((loops % 50 == 0) && (loops <= 600))
// 	    || ((loops % 200 == 0) && (loops <= 1000))
// 	    || (loops % 500 == 0)
// 	    ) {
// 	  u_epsilon = w;
// 	  u_epsilon.scale(&P,-1);  
//  	  char filename1[50];
//  	  char filename2[50];
//  	  sprintf(filename1, "%s%d%s%d%s", "approx_sol_steep35_2D_2509out_", loops, "_nactive_", w.size(),".m");
//  	  sprintf(filename2, "%s%d%s%d%s", "error_steep35_2D_2509out_", loops, "_nactive_", w.size(),".m");
// 	  Array1D<SampledMapping<2> > U = evalObj.evaluate(P.basis(), u_epsilon, true, 6);
// 	  cout << "done plotting approximate solution" << endl;
// 	  Array1D<SampledMapping<2> > Error = evalObj.evaluate_difference(P.basis(), u_epsilon, sing2D, 6);
// 	  cout << "done plotting pointwise error" << endl;
// 	  std::ofstream ofs5(filename1);
// 	  matlab_output(ofs5,U);
// 	  ofs5.close();
// 	  std::ofstream ofs6(filename2);
// 	  matlab_output(ofs6,Error);
// 	  ofs6.close();
// 	}

	int d  = Basis1D::primal_polynomial_degree();
	int dT = Basis1D::primal_vanishing_moments();
	char name1[60];
	char name2[60];
	char name3[60];
	
	sprintf(name1, "%s%d%s%d%s", "./steepest_descent_results/steep1D_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name2, "%s%d%s%d%s", "./steepest_descent_results/steep1D_time_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name3, "%s%d%s%d%s", "./steepest_descent_results/steep1D_descent_params_P_jmax18_", d, "_dT", dT, ".m");
	
	
	std::ofstream os1(name1);
	matlab_output(asymptotic,os1);
	os1.close();
	
	std::ofstream os2(name2);
	matlab_output(time_asymptotic,os2);
	os2.close();
	
	std::ofstream os3(name3);
	matlab_output(descent_params,os3);
	os3.close();
	
	
	tstart = clock();
	// ############# end output #############

	//if (residual_norm < 3.1623e-04 || loops == 1000) {
	if (residual_norm < 1.0e-3 || loops == 1000) {

	u_epsilon = w;
	  exit = true;
	  break;
	}

      }//end while
      
      cout << "#######################" << endl;
      cout << "exiting inner loop" << endl;
      cout << "#######################" << endl;
      cout << "number of calls of APPLY = " << niter << endl;

      if (exit)
	break;
//       cout << "tolerance for COARSE = " << ((3.*mu*omega_i)/(1+3.*mu)) << endl;
//       InfiniteVector<double, Index> tmp;
//       w.COARSE(residual_norm, tmp);
//       w = tmp;
    }// end for 
    
    
    // collect final approximation and its local parts
    approximations[P.basis().n_p()] = u_epsilon;
    
    for (int i = 0; i < P.basis().n_p(); i++) {
      approximations[i].clear();
      for (typename InfiniteVector<double, Index>::const_iterator it = u_epsilon.begin(), itend = u_epsilon.end();
	   it != itend; ++it)
	if (it.index().p() == i)
	  approximations[i].set_coefficient(it.index(),*it);
    }
  }
}
