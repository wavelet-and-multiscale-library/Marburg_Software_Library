// implementation for steepest_descent_ks.h

#include <cmath>
#include <set>
#include <utils/plot_tools.h>
#include <adaptive/apply.h>
#include <numerics/corner_singularity.h>
#include <interval/p_basis.h>


using std::set;

namespace WaveletTL{

  template <class PROBLEM>
  void steepest_descent_ks_SOLVE(const PROBLEM& P,  const double epsilon,
			      InfiniteVector<double, typename PROBLEM::Index>& approximations)
  {

    // promal and dual spline orders of the wavelets
    const int d  = PRIMALORDER;
    const int dt = DUALORDER;
    typedef PBasis<d,dt> Basis1D;
  
    //     Point<2> origin;
    //     origin[0] = 0.0;
    //     origin[1] = 0.0;
    
    //     CornerSingularity sing2D(origin, 0.5, 1.5);
    //     CornerSingularityRHS singRhs(origin, 0.5, 1.5);

    //    Singularity1D_RHS_2<double> sing1D;
    //    Singularity1D_2<double> exactSolution1D;

    unsigned int loops = 0;
    unsigned int niter = 0;

    // the maximal level
    const int jmax = JMAX;
    typedef typename PROBLEM::Index Index;

    // #####################################################################################
    // We have to set up the various constants steering the decay of the accuracy.
    // However, most theoretical estimates turn out to be too pessimistic. To obtain
    // a good performance, we have to manually choose them more optimistic.
    // #####################################################################################

    // norm of the pseudo inverse
    double a_inv     = P.norm_Ainv();
    
    // spectral condition number
    double kappa     = P.norm_A()*a_inv;

    // upper bound for the \ell_2-norm of the exact discrete solution in the range of the
    // stiffness matrix
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

    //double beta      = 0.85; // this was chosen for d=2,3
    //double beta      = 0.75; // this was chosen for d=4
    double beta      = 0.5; // this was chosen for d=4
    
    cout << "beta = " << beta << endl;

    //let K be such that beta^K * omega <= epsilon
    unsigned int K   = (int) (log(epsilon/omega_i) / log(beta) + 1);
    //let M be such that lambda^M <= ((1-delta) / (1+delta)) * (beta / ((1+3*mu)*kappa))
    int M            = std::max((int) ((log( ((1-delta)/(1+delta)) * (beta / ((1+3.0*mu)*kappa)) )
					/ log(lambda)) + 1),1);

    cout << "K = " << K << endl;
    cout << "M = " << M << endl;
    // #####################################################################################
    // End setting up constants.
    // #####################################################################################

    // InfiniteVector's used in the iterative algorithm
    InfiniteVector<double, Index> w, tilde_r, help, f, Av;

    // map's used for generating output
    map<double,double> log_10_residual_norms;
    map<double,double> degrees_of_freedom;
    map<double,double> asymptotic;
    map<double,double> time_asymptotic;
    map<double,double> descent_params;
    map<double,double> weak_ell_tau_norms;

    
    bool exit = 0;
    double time = 0.;
    double residual_norm = 5.0;

    // variables for runtime measurement
    clock_t tstart, tend;
    // get the current time
    tstart = clock();

    //EvaluateFrame<Basis1D,2,2> evalObj;

    double dd = 0.5;

    // the adaptive algorithm
    for (unsigned int i = 1; i <= K; i++) {
      omega_i *= beta;
      double xi_i = omega_i / ((1+3.0*mu)*C3*M);
      double nu_i = 0.;

      RES(P, w, xi_i, delta, omega_i/((1+3.*mu)*a_inv), jmax,
	  tilde_r, nu_i, niter, CDD1); 

      while ( nu_i > omega_i/((1+3.*mu)*a_inv)) {

	InfiniteVector<double, Index> z_i;

	// Instead of using APPLY only, we use a call of apply followed by an 
	// immeadiate call of COARSE. This is done to prevent that the number
	// of non-zeros in the iterates grow very quickly.
	APPLY_COARSE(P, tilde_r, delta*l2_norm(tilde_r), z_i, 1.0e-6, jmax, CDD1);
	//APPLY_COARSE(P, tilde_r, delta*l2_norm(tilde_r), z_i, 0.5, jmax, CDD1);
 	double g = z_i*tilde_r;
	if  (g != 0.)
	  dd = (tilde_r*tilde_r)/g;
	
	w += dd*tilde_r;
	//  	InfiniteVector<double, Index> tmp;
	//  	w.COARSE(1.0/100.0*residual_norm, tmp);
	//  	w = tmp;


	cout << "descent param = " << dd << endl;
	++loops;
	++niter;

	RES(P, w, xi_i, delta, omega_i/((1+3.*mu)*a_inv), jmax,
	    tilde_r, nu_i, niter, CDD1);

	cout << "loop: " << loops << " nu = " 
	     << nu_i << " epsilon = " << omega_i/((1+3.*mu)*a_inv) << endl;
	cout << "xi: " << xi_i << endl; 


	// take the elapsed time
	tend = clock();
	time += ((double) (tend-tstart))/((double) CLOCKS_PER_SEC);

	// #####################################################################################
	// Approximate the EXACT residual using a sufficiently small precision
	// and perform output.
	// #####################################################################################
  	P.RHS(1.0e-6,f);
 	APPLY(P, w, 1.0e-6, Av, jmax, CDD1);
 	help = f-Av;

	residual_norm = l2_norm(help);
	double tmp1 = log10(residual_norm);
	cout << "residual norm = " << residual_norm << endl;

	cout << "active indices: " << w.size() << endl;
	asymptotic[log10( (double)w.size() )] = tmp1;
	time_asymptotic[log10(time)] = tmp1;
	
	descent_params[loops] = dd;

	int d  = Basis1D::primal_polynomial_degree();
	int dT = Basis1D::primal_vanishing_moments();

	char name1[128];
	char name2[128];
	char name3[128];
	char name4[128];

	// setup filenames for output files for the one-dimensional cases
#ifdef ONE_D
	switch (d) {
	case 2: {
	  weak_ell_tau_norms[loops] = w.weak_norm(1./1.5);// (d,dT)=(2,2)=1./1.5 (d,dT)=(3,3)=1./2.5, (d,dT)=(4,6)=1./3.5
	  sprintf(name1, "%s%d%s%d%s", "./sd_results22/steep1D_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name2, "%s%d%s%d%s", "./sd_results22/steep1D_time_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name3, "%s%d%s%d%s", "./sd_results22/steep1D_descent_params_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name4, "%s%d%s%d%s", "./sd_results22/steep1D_weak_ell_tau_norms_P_jmax18_", d, "_dT", dT, ".m");
	  break;
	}
	case 3: {
          weak_ell_tau_norms[loops] = w.weak_norm(1./2.5);
	  sprintf(name1, "%s%d%s%d%s", "./sd_results33_basis/steep1D_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name2, "%s%d%s%d%s", "./sd_results33_basis/steep1D_time_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name3, "%s%d%s%d%s", "./sd_results33_basis/steep1D_descent_params_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name4, "%s%d%s%d%s", "./sd_results33_basis/steep1D_weak_ell_tau_norms_P_jmax18_", d, "_dT", dT, ".m");
	  break;
	}
	case 4: {
	  weak_ell_tau_norms[loops] = w.weak_norm(1./3.5);
	  sprintf(name1, "%s%d%s%d%s", "./sd_results46/steep1D_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name2, "%s%d%s%d%s", "./sd_results46/steep1D_time_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name3, "%s%d%s%d%s", "./sd_results46/steep1D_descent_params_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name4, "%s%d%s%d%s", "./sd_results46/steep1D_weak_ell_tau_norms_P_jmax18_", d, "_dT", dT, ".m");  
	  break;
	}
	};
#endif
	// setup filenames for output files for the two-dimensional cases
#ifdef TWO_D
	switch (d) {
	case 2: {
	  weak_ell_tau_norms[loops] = w.weak_norm(1./1.0);
	  sprintf(name1, "%s%d%s%d%s", "./sd_results2D_22/steep2D_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name2, "%s%d%s%d%s", "./sd_results2D_22/steep2D_time_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name3, "%s%d%s%d%s", "./sd_results2D_22/steep2D_descent_params_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name4, "%s%d%s%d%s", "./sd_results2D_22/steep2D_weak_ell_tau_norms_P_jmax18_", d, "_dT", dT, ".m");
	  break;
	}
	case 3: {
	  weak_ell_tau_norms[loops] = w.weak_norm(1./1.5);
	  sprintf(name1, "%s%d%s%d%s", "./sd_results2D_33/steep2D_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name2, "%s%d%s%d%s", "./sd_results2D_33/steep2D_time_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name3, "%s%d%s%d%s", "./sd_results2D_33/steep2D_descent_params_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name4, "%s%d%s%d%s", "./sd_results2D_33/steep2D_weak_ell_tau_norms_P_jmax18_", d, "_dT", dT, ".m");
	  break;
	}
	case 4: {
	  weak_ell_tau_norms[loops] = w.weak_norm(1./2.0);
	  sprintf(name1, "%s%d%s%d%s", "./sd_results2D_46/steep2D_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name2, "%s%d%s%d%s", "./sd_results2D_46/steep2D_time_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name3, "%s%d%s%d%s", "./sd_results2D_46/steep2D_descent_params_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name4, "%s%d%s%d%s", "./sd_results2D_46/steep2D_weak_ell_tau_norms_P_jmax18_", d, "_dT", dT, ".m");  
	  break;
	}
	};

#endif
	
	// perform matlab output
	std::ofstream os1(name1);
	matlab_output(asymptotic,os1);
	os1.close();
	
	std::ofstream os2(name2);
	matlab_output(time_asymptotic,os2);
	os2.close();
	
	std::ofstream os3(name3);
	matlab_output(descent_params,os3);
	os3.close();
	
	std::ofstream os4(name4);
	matlab_output(weak_ell_tau_norms,os4);
	os3.close();

	// #####################################################################################
	// End performing output.
	// #####################################################################################

	// restart the stopwatch.
	tstart = clock();


	// create MANUAL stopping criterion for the one-dimensional case
#ifdef ONE_D
	double tol=1.0;
	switch (d) {
	case 2: {
	  tol=0.015;
	  break;
	}
	case 3: {
	  tol=3.1623e-04;
	  break;
	}
	case 4: {
	  tol=3.1623e-04;
	  break;
	}
	}
	if (residual_norm < tol || loops == 5000) {
#endif

	  // create MANUAL stopping criterion for the two-dimensional case
#ifdef TWO_D
	  if (residual_norm < 0.01 || loops == 5000) {
#endif
	    approximations = w;
	    //u_epsilon = w;
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
      
      
	// #####################################################################################
	// Perform the coarsening. This can usually be dropped when COARSE is always applied
	// directly after the APPLY. On the theoretical side, however, it is still mandatory
	// for this algorithm.
	// #####################################################################################
	//       cout << "tolerance for COARSE = " << ((3.*mu*omega_i)/(1+3.*mu)) << endl;
	//       InfiniteVector<double, Index> tmp;
	//       w.COARSE(residual_norm, tmp);
	//       w = tmp;
      }// end for 
      
      // #####################################################################################
      // The algorithm is finished.
      // Collect final approximation and its local parts
      // #####################################################################################
      
      
      
//    approximations.clear();
//    //for (typename InfiniteVector<double, Index>::const_iterator it = u_epsilon.begin(), itend = u_epsilon.end();
//    for (typename InfiniteVector<double, Index>::const_iterator it = approximations.begin(),
//           itend = approximations.end();
//         it != itend; ++it)
//      
//        approximations.set_coefficient(it.index(),*it);
//      
    }
  
    
  
  
  template <class PROBLEM>
  void steepest_descent_ks_QUARKLET_SOLVE(const PROBLEM& P,  const double epsilon,
			      InfiniteVector<double, typename PROBLEM::Index>& approximations, const CompressionStrategy strategy,
                              const double a, const double b)
  {

    // promal and dual spline orders of the wavelets
    const int d  = PRIMALORDER;
    const int dt = DUALORDER;
    typedef PQFrame<d,dt> Basis1D;
  
    //     Point<2> origin;
    //     origin[0] = 0.0;
    //     origin[1] = 0.0;
    
    //     CornerSingularity sing2D(origin, 0.5, 1.5);
    //     CornerSingularityRHS singRhs(origin, 0.5, 1.5);

    //    Singularity1D_RHS_2<double> sing1D;
    //    Singularity1D_2<double> exactSolution1D;

    unsigned int loops = 0;
    unsigned int niter = 0;

    // the maximal level
    const int jmax = JMAX;
    const int pmax = PMAX;
    typedef typename PROBLEM::Index Index;

    // #####################################################################################
    // We have to set up the various constants steering the decay of the accuracy.
    // However, most theoretical estimates turn out to be too pessimistic. To obtain
    // a good performance, we have to manually choose them more optimistic.
    // #####################################################################################

    // norm of the pseudo inverse
    double a_inv     = P.norm_Ainv();
    
    // spectral condition number
    double kappa     = P.norm_A()*a_inv;

    // upper bound for the \ell_2-norm of the exact discrete solution in the range of the
    // stiffness matrix
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

    double beta      = 0.85; // this was chosen for d=2,3
    //double beta      = 0.75; // this was chosen for d=4
//    double beta      = 0.5; // this was chosen for d=4
    
    cout << "beta = " << beta << endl;

    //let K be such that beta^K * omega <= epsilon
    unsigned int K   = (int) (log(epsilon/omega_i) / log(beta) + 1);
    //let M be such that lambda^M <= ((1-delta) / (1+delta)) * (beta / ((1+3*mu)*kappa))
    int M            = std::max((int) ((log( ((1-delta)/(1+delta)) * (beta / ((1+3.0*mu)*kappa)) )
					/ log(lambda)) + 1),1);

    cout << "K = " << K << endl;
    cout << "M = " << M << endl;
//    abort();
    // #####################################################################################
    // End setting up constants.
    // #####################################################################################

    // InfiniteVector's used in the iterative algorithm
    InfiniteVector<double, Index> w, tilde_r, help, f, Av;

    // map's used for generating output
    map<double,double> log_10_residual_norms;
    map<double,double> degrees_of_freedom;
    map<double,double> asymptotic;
    map<double,double> time_asymptotic;
    map<double,double> descent_params;
    map<double,double> weak_ell_tau_norms;

    
    bool exit = 0;
    double time = 0.;
    double residual_norm = 5.0;

    // variables for runtime measurement
    clock_t tstart, tend;
    // get the current time
    tstart = clock();

    //EvaluateFrame<Basis1D,2,2> evalObj;

    double dd = 0.5;

    // the adaptive algorithm
    for (unsigned int i = 1; i <= K; i++) {
      omega_i *= beta;
      double xi_i = omega_i / ((1+3.0*mu)*C3*M);
      double nu_i = 0.;

      RES_QUARKLET(P, w, xi_i, delta, omega_i/((1+3.*mu)*a_inv), jmax,
	  tilde_r, nu_i, niter, strategy, pmax, a, b); 

      while ( nu_i > omega_i/((1+3.*mu)*a_inv)) {

	InfiniteVector<double, Index> z_i;

	// Instead of using APPLY only, we use a call of apply followed by an 
	// immeadiate call of COARSE. This is done to prevent that the number
	// of non-zeros in the iterates grow very quickly.
	APPLY_QUARKLET_COARSE(P, tilde_r, delta*l2_norm(tilde_r), z_i, jmax, strategy, pmax, a, b, 1e-6);
	//APPLY_COARSE(P, tilde_r, delta*l2_norm(tilde_r), z_i, 0.5, jmax, CDD1);
 	double g = z_i*tilde_r;
	if  (g != 0.)
	  dd = (tilde_r*tilde_r)/g;
	
	w += dd*tilde_r;
	//  	InfiniteVector<double, Index> tmp;
	//  	w.COARSE(1.0/100.0*residual_norm, tmp);
	//  	w = tmp;


	cout << "descent param = " << dd << endl;
	++loops;
	++niter;

	RES_QUARKLET(P, w, xi_i, delta, omega_i/((1+3.*mu)*a_inv), jmax,
	    tilde_r, nu_i, niter, strategy, pmax, a, b);

	cout << "loop: " << loops << " nu = " 
	     << nu_i << " epsilon = " << omega_i/((1+3.*mu)*a_inv) << endl;
	cout << "xi: " << xi_i << endl; 


	// take the elapsed time
	tend = clock();
	time += ((double) (tend-tstart))/((double) CLOCKS_PER_SEC);

	// #####################################################################################
	// Approximate the EXACT residual using a sufficiently small precision
	// and perform output.
	// #####################################################################################
  	P.RHS(1.0e-6,f);
 	APPLY_QUARKLET(P, w, 1.0e-6, Av, jmax, strategy, pmax, a, b);
 	help = f-Av;

	residual_norm = l2_norm(help);
	double tmp1 = log10(residual_norm);
	cout << "residual norm = " << residual_norm << endl;

	cout << "active indices: " << w.size() << endl;
	asymptotic[log10( (double)w.size() )] = tmp1;
	time_asymptotic[log10(time)] = tmp1;
	
	descent_params[loops] = dd;

	int d  = Basis1D::primal_polynomial_degree();
	int dT = Basis1D::primal_vanishing_moments();

	char name1[128];
	char name2[128];
	char name3[128];
	char name4[128];

	// setup filenames for output files for the one-dimensional cases
#ifdef ONE_D
	switch (d) {
	case 2: {
	  weak_ell_tau_norms[loops] = w.weak_norm(1./1.5);// (d,dT)=(2,2)=1./1.5 (d,dT)=(3,3)=1./2.5, (d,dT)=(4,6)=1./3.5
	  sprintf(name1, "%s%d%s%d%s", "./sd_results22/steep1D_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name2, "%s%d%s%d%s", "./sd_results22/steep1D_time_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name3, "%s%d%s%d%s", "./sd_results22/steep1D_descent_params_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name4, "%s%d%s%d%s", "./sd_results22/steep1D_weak_ell_tau_norms_P_jmax18_", d, "_dT", dT, ".m");
	  break;
	}
	case 3: {
          weak_ell_tau_norms[loops] = w.weak_norm(1./2.5);
	  sprintf(name1, "%s%d%s%d%s", "./sd_results33/steep1D_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name2, "%s%d%s%d%s", "./sd_results33/steep1D_time_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name3, "%s%d%s%d%s", "./sd_results33/steep1D_descent_params_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name4, "%s%d%s%d%s", "./sd_results3/steep1D_weak_ell_tau_norms_P_jmax18_", d, "_dT", dT, ".m");
	  break;
	}
	case 4: {
	  weak_ell_tau_norms[loops] = w.weak_norm(1./3.5);
	  sprintf(name1, "%s%d%s%d%s", "./sd_results46/steep1D_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name2, "%s%d%s%d%s", "./sd_results46/steep1D_time_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name3, "%s%d%s%d%s", "./sd_results46/steep1D_descent_params_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name4, "%s%d%s%d%s", "./sd_results46/steep1D_weak_ell_tau_norms_P_jmax18_", d, "_dT", dT, ".m");  
	  break;
	}
	};
#endif
	// setup filenames for output files for the two-dimensional cases
#ifdef TWO_D
	switch (d) {
	case 2: {
	  weak_ell_tau_norms[loops] = w.weak_norm(1./1.0);
	  sprintf(name1, "%s%d%s%d%s", "./sd_results2D_22/steep2D_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name2, "%s%d%s%d%s", "./sd_results2D_22/steep2D_time_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name3, "%s%d%s%d%s", "./sd_results2D_22/steep2D_descent_params_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name4, "%s%d%s%d%s", "./sd_results2D_22/steep2D_weak_ell_tau_norms_P_jmax18_", d, "_dT", dT, ".m");
	  break;
	}
	case 3: {
	  weak_ell_tau_norms[loops] = w.weak_norm(1./1.5);
	  sprintf(name1, "%s%d%s%d%s", "./sd_quarklets_results/steep2D_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name2, "%s%d%s%d%s", "./sd_quarklets_results/steep2D_time_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name3, "%s%d%s%d%s", "./sd_quarklets_results/steep2D_descent_params_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name4, "%s%d%s%d%s", "./sd_quarklets_results/steep2D_weak_ell_tau_norms_P_jmax18_", d, "_dT", dT, ".m");
	  break;
	}
	case 4: {
	  weak_ell_tau_norms[loops] = w.weak_norm(1./2.0);
	  sprintf(name1, "%s%d%s%d%s", "./sd_results2D_46/steep2D_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name2, "%s%d%s%d%s", "./sd_results2D_46/steep2D_time_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name3, "%s%d%s%d%s", "./sd_results2D_46/steep2D_descent_params_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name4, "%s%d%s%d%s", "./sd_results2D_46/steep2D_weak_ell_tau_norms_P_jmax18_", d, "_dT", dT, ".m");  
	  break;
	}
	};

#endif
	
	// perform matlab output
	std::ofstream os1(name1);
	matlab_output(asymptotic,os1);
	os1.close();
	
	std::ofstream os2(name2);
	matlab_output(time_asymptotic,os2);
	os2.close();
	
	std::ofstream os3(name3);
	matlab_output(descent_params,os3);
	os3.close();
	
	std::ofstream os4(name4);
	matlab_output(weak_ell_tau_norms,os4);
	os3.close();

	// #####################################################################################
	// End performing output.
	// #####################################################################################

	// restart the stopwatch.
	tstart = clock();


	// create MANUAL stopping criterion for the one-dimensional case
#ifdef ONE_D
	double tol=1.0;
	switch (d) {
	case 2: {
//	  tol=0.015;
          tol=1e-5;
	  break;
	}
	case 3: {
	  tol=3.1623e-04;
	  break;
	}
	case 4: {
	  tol=3.1623e-04;
	  break;
	}
        case 5: {
	  tol=3.1623e-04;
	  break;
	}
	}
	if (residual_norm < tol || loops == 5000) {
#endif

	  // create MANUAL stopping criterion for the two-dimensional case
#ifdef TWO_D
	  if (residual_norm < 0.01 || loops == 5000) {
#endif
	    approximations = w;
	    //u_epsilon = w;
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
      
      
	// #####################################################################################
	// Perform the coarsening. This can usually be dropped when COARSE is always applied
	// directly after the APPLY. On the theoretical side, however, it is still mandatory
	// for this algorithm.
	// #####################################################################################
	//       cout << "tolerance for COARSE = " << ((3.*mu*omega_i)/(1+3.*mu)) << endl;
	//       InfiniteVector<double, Index> tmp;
	//       w.COARSE(residual_norm, tmp);
	//       w = tmp;
      }// end for 
      
      // #####################################################################################
      // The algorithm is finished.
      // Collect final approximation and its local parts
      // #####################################################################################
      
      
      
//    approximations.clear();
//    //for (typename InfiniteVector<double, Index>::const_iterator it = u_epsilon.begin(), itend = u_epsilon.end();
//    for (typename InfiniteVector<double, Index>::const_iterator it = approximations.begin(),
//           itend = approximations.end();
//         it != itend; ++it)
//      
//        approximations.set_coefficient(it.index(),*it);
//      
    }
  
    
  
  
  template <class PROBLEM>
  void steepest_descent_ks_QUARKLET_SOLVE(const PROBLEM& P,  const double epsilon,
			      InfiniteVector<double, int>& approximations, const CompressionStrategy strategy,
                              const double a, const double b)
  {

    // promal and dual spline orders of the wavelets
    const int d  = PRIMALORDER;
    const int dt = DUALORDER;
    typedef PQFrame<d,dt> Basis1D;
  
    //     Point<2> origin;
    //     origin[0] = 0.0;
    //     origin[1] = 0.0;
    
    //     CornerSingularity sing2D(origin, 0.5, 1.5);
    //     CornerSingularityRHS singRhs(origin, 0.5, 1.5);

    //    Singularity1D_RHS_2<double> sing1D;
    //    Singularity1D_2<double> exactSolution1D;

    unsigned int loops = 0;
    unsigned int niter = 0;

    // the maximal level
    const int jmax = JMAX;
    const int pmax = PMAX;
//    typedef typename PROBLEM::Index Index;

    // #####################################################################################
    // We have to set up the various constants steering the decay of the accuracy.
    // However, most theoretical estimates turn out to be too pessimistic. To obtain
    // a good performance, we have to manually choose them more optimistic.
    // #####################################################################################

    // norm of the pseudo inverse
    double a_inv     = P.norm_Ainv();
    
    // spectral condition number
    double kappa     = P.norm_A()*a_inv;

    // upper bound for the \ell_2-norm of the exact discrete solution in the range of the
    // stiffness matrix
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

    //double beta      = 0.85; // this was chosen for d=2,3
    //double beta      = 0.75; // this was chosen for d=4
    double beta      = 0.5; // this was chosen for d=4
    
    cout << "beta = " << beta << endl;

    //let K be such that beta^K * omega <= epsilon
    unsigned int K   = (int) (log(epsilon/omega_i) / log(beta) + 1);
    //let M be such that lambda^M <= ((1-delta) / (1+delta)) * (beta / ((1+3*mu)*kappa))
    int M            = std::max((int) ((log( ((1-delta)/(1+delta)) * (beta / ((1+3.0*mu)*kappa)) )
					/ log(lambda)) + 1),1);

    cout << "K = " << K << endl;
    cout << "M = " << M << endl;
//    abort();
    // #####################################################################################
    // End setting up constants.
    // #####################################################################################

    // InfiniteVector's used in the iterative algorithm
    InfiniteVector<double, int> w, tilde_r, help, f, Av;

    // map's used for generating output
    map<double,double> log_10_residual_norms;
    map<double,double> degrees_of_freedom;
    map<double,double> asymptotic;
    map<double,double> time_asymptotic;
    map<double,double> descent_params;
    map<double,double> weak_ell_tau_norms;

    
    bool exit = 0;
    double time = 0.;
    double residual_norm = 5.0;

    // variables for runtime measurement
    clock_t tstart, tend;
    // get the current time
    tstart = clock();

    //EvaluateFrame<Basis1D,2,2> evalObj;

    double dd = 0.5;

    // the adaptive algorithm
    for (unsigned int i = 1; i <= K; i++) {
      omega_i *= beta;
      double xi_i = omega_i / ((1+3.0*mu)*C3*M);
      double nu_i = 0.;

      RES_QUARKLET(P, w, xi_i, delta, omega_i/((1+3.*mu)*a_inv), jmax,
	  tilde_r, nu_i, niter, strategy, pmax, a, b); 

      while ( nu_i > omega_i/((1+3.*mu)*a_inv)) {

	InfiniteVector<double, int> z_i;

	// Instead of using APPLY only, we use a call of apply followed by an 
	// immeadiate call of COARSE. This is done to prevent that the number
	// of non-zeros in the iterates grow very quickly.
	APPLY_QUARKLET_COARSE(P, tilde_r, delta*l2_norm(tilde_r), z_i, jmax, strategy, pmax, a, b, 1e-6);
	//APPLY_COARSE(P, tilde_r, delta*l2_norm(tilde_r), z_i, 0.5, jmax, CDD1);
 	double g = z_i*tilde_r;
	if  (g != 0.)
	  dd = (tilde_r*tilde_r)/g;
	
	w += dd*tilde_r;
	//  	InfiniteVector<double, Index> tmp;
	//  	w.COARSE(1.0/100.0*residual_norm, tmp);
	//  	w = tmp;


	cout << "descent param = " << dd << endl;
	++loops;
	++niter;

	RES_QUARKLET(P, w, xi_i, delta, omega_i/((1+3.*mu)*a_inv), jmax,
	    tilde_r, nu_i, niter, strategy, pmax, a, b);

	cout << "loop: " << loops << " nu = " 
	     << nu_i << " epsilon = " << omega_i/((1+3.*mu)*a_inv) << endl;
	cout << "xi: " << xi_i << endl; 


	// take the elapsed time
	tend = clock();
	time += ((double) (tend-tstart))/((double) CLOCKS_PER_SEC);

	// #####################################################################################
	// Approximate the EXACT residual using a sufficiently small precision
	// and perform output.
	// #####################################################################################
  	P.RHS(1.0e-6,f);
 	APPLY_QUARKLET(P, w, 1.0e-6, Av, jmax, strategy, pmax, a, b);
 	help = f-Av;

	residual_norm = l2_norm(help);
	double tmp1 = log10(residual_norm);
	cout << "residual norm = " << residual_norm << endl;

	cout << "active indices: " << w.size() << endl;
	asymptotic[log10( (double)w.size() )] = tmp1;
	time_asymptotic[log10(time)] = tmp1;
	
	descent_params[loops] = dd;

	int d  = Basis1D::primal_polynomial_degree();
	int dT = Basis1D::primal_vanishing_moments();

	char name1[128];
	char name2[128];
	char name3[128];
	char name4[128];

	// setup filenames for output files for the one-dimensional cases
#ifdef ONE_D
	switch (d) {
	case 2: {
	  weak_ell_tau_norms[loops] = w.weak_norm(1./1.5);// (d,dT)=(2,2)=1./1.5 (d,dT)=(3,3)=1./2.5, (d,dT)=(4,6)=1./3.5
	  sprintf(name1, "%s%d%s%d%s", "./sd_results22/steep1D_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name2, "%s%d%s%d%s", "./sd_results22/steep1D_time_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name3, "%s%d%s%d%s", "./sd_results22/steep1D_descent_params_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name4, "%s%d%s%d%s", "./sd_results22/steep1D_weak_ell_tau_norms_P_jmax18_", d, "_dT", dT, ".m");
	  break;
	}
	case 3: {
          weak_ell_tau_norms[loops] = w.weak_norm(1./2.5);
	  sprintf(name1, "%s%d%s%d%s", "./sd_results33_basis/steep1D_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name2, "%s%d%s%d%s", "./sd_results33_basis/steep1D_time_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name3, "%s%d%s%d%s", "./sd_results33_basis/steep1D_descent_params_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name4, "%s%d%s%d%s", "./sd_results33_basis/steep1D_weak_ell_tau_norms_P_jmax18_", d, "_dT", dT, ".m");
	  break;
	}
	case 4: {
	  weak_ell_tau_norms[loops] = w.weak_norm(1./3.5);
	  sprintf(name1, "%s%d%s%d%s", "./sd_results46/steep1D_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name2, "%s%d%s%d%s", "./sd_results46/steep1D_time_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name3, "%s%d%s%d%s", "./sd_results46/steep1D_descent_params_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name4, "%s%d%s%d%s", "./sd_results46/steep1D_weak_ell_tau_norms_P_jmax18_", d, "_dT", dT, ".m");  
	  break;
	}
	};
#endif
	// setup filenames for output files for the two-dimensional cases
#ifdef TWO_D
	switch (d) {
	case 2: {
	  weak_ell_tau_norms[loops] = w.weak_norm(1./1.0);
	  sprintf(name1, "%s%d%s%d%s", "./sd_results2D_22/steep2D_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name2, "%s%d%s%d%s", "./sd_results2D_22/steep2D_time_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name3, "%s%d%s%d%s", "./sd_results2D_22/steep2D_descent_params_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name4, "%s%d%s%d%s", "./sd_results2D_22/steep2D_weak_ell_tau_norms_P_jmax18_", d, "_dT", dT, ".m");
	  break;
	}
	case 3: {
	  weak_ell_tau_norms[loops] = w.weak_norm(1./1.5);
	  sprintf(name1, "%s%d%s%d%s", "./sd_quarklets_results/steep2D_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name2, "%s%d%s%d%s", "./sd_quarklets_results/steep2D_time_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name3, "%s%d%s%d%s", "./sd_quarklets_results/steep2D_descent_params_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name4, "%s%d%s%d%s", "./sd_quarklets_results/steep2D_weak_ell_tau_norms_P_jmax18_", d, "_dT", dT, ".m");
	  break;
	}
	case 4: {
	  weak_ell_tau_norms[loops] = w.weak_norm(1./2.0);
	  sprintf(name1, "%s%d%s%d%s", "./sd_results2D_46/steep2D_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name2, "%s%d%s%d%s", "./sd_results2D_46/steep2D_time_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name3, "%s%d%s%d%s", "./sd_results2D_46/steep2D_descent_params_P_jmax18_", d, "_dT", dT, ".m");
	  sprintf(name4, "%s%d%s%d%s", "./sd_results2D_46/steep2D_weak_ell_tau_norms_P_jmax18_", d, "_dT", dT, ".m");  
	  break;
	}
	};

#endif
	
	// perform matlab output
	std::ofstream os1(name1);
	matlab_output(asymptotic,os1);
	os1.close();
	
	std::ofstream os2(name2);
	matlab_output(time_asymptotic,os2);
	os2.close();
	
	std::ofstream os3(name3);
	matlab_output(descent_params,os3);
	os3.close();
	
	std::ofstream os4(name4);
	matlab_output(weak_ell_tau_norms,os4);
	os3.close();

	// #####################################################################################
	// End performing output.
	// #####################################################################################

	// restart the stopwatch.
	tstart = clock();


	// create MANUAL stopping criterion for the one-dimensional case
#ifdef ONE_D
	double tol=1.0;
	switch (d) {
	case 2: {
	  tol=0.015;
	  break;
	}
	case 3: {
	  tol=3.1623e-04;
	  break;
	}
	case 4: {
	  tol=3.1623e-04;
	  break;
	}
        case 5: {
	  tol=3.1623e-04;
	  break;
	}
	}
	if (residual_norm < tol || loops == 5000) {
#endif

	  // create MANUAL stopping criterion for the two-dimensional case
#ifdef TWO_D
	  if (residual_norm < 0.01 || loops == 5000) {
#endif
	    approximations = w;
	    //u_epsilon = w;
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
      
      
	// #####################################################################################
	// Perform the coarsening. This can usually be dropped when COARSE is always applied
	// directly after the APPLY. On the theoretical side, however, it is still mandatory
	// for this algorithm.
	// #####################################################################################
	//       cout << "tolerance for COARSE = " << ((3.*mu*omega_i)/(1+3.*mu)) << endl;
	//       InfiniteVector<double, Index> tmp;
	//       w.COARSE(residual_norm, tmp);
	//       w = tmp;
      }// end for 
      
      // #####################################################################################
      // The algorithm is finished.
      // Collect final approximation and its local parts
      // #####################################################################################
      
      
      
//    approximations.clear();
//    //for (typename InfiniteVector<double, Index>::const_iterator it = u_epsilon.begin(), itend = u_epsilon.end();
//    for (typename InfiniteVector<double, Index>::const_iterator it = approximations.begin(),
//           itend = approximations.end();
//         it != itend; ++it)
//      
//        approximations.set_coefficient(it.index(),*it);
//      
    }
    
  }
  

  
  


