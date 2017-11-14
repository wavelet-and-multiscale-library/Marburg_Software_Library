// implementation for richardson.h

#include <cmath>
#include <set>
#include <utils/plot_tools.h>
#include <adaptive/apply.h>
#include <numerics/corner_singularity.h>
#include <frame_evaluate.h>
#include <poisson_1d_testcase.h>

using std::set;

namespace FrameTL
{
  template <class PROBLEM>
  void richardson_SOLVE(const PROBLEM& P, const double epsilon,
			InfiniteVector<double, typename PROBLEM::Index>& u_epsilon,
			Array1D<InfiniteVector<double, typename PROBLEM::Index> >& approximations)
  {


    const int d  = PRIMALORDER;
    const int dt = DUALORDER;
    typedef PBasis<d,dt> Basis1D;

//     Singularity1D_RHS_2<double> sing1D;
//     Singularity1D_2<double> exactSolution1D;
     
//      Point<2> origin;
//      origin[0] = 0.0;
//      origin[1] = 0.0;

//      CornerSingularity sing2D(origin, 0.5, 1.5);
//      CornerSingularityRHS singRhs(origin, 0.5, 1.5);

    const unsigned int jmax = JMAX;

    const double nu = P.norm_Ainv()*P.F_norm();
    typedef typename PROBLEM::Index Index;

    cout << "Rich_SOLVE: nu=" << nu << endl;

    // compute optimal relaxation parameter omega
    //const double omega = 2.0 / (P.norm_A() + 1.0/P.norm_Ainv());
    
    // ####### 1D #######
//    const double omega = 0.4;
    const double omega = 0.15;

    // ####### 2D #######
    // bad one
//    const double omega = 0.05;
    // good one
    //const double omega = 0.25;

    //const double omega = 0.3;


    //const double omega = 2.0/2.47-0.5;
    cout << "Rich_SOLVE: omega=" << omega << endl;

    // compute spectral norm rho
    //const double cond_A = P.norm_A() * P.norm_Ainv();
    //const double rho = (cond_A - 1.0) / (cond_A + 1.0);
    const double rho = 0.8;
    cout << "Rich_SOLVE: rho=" << rho << endl;
    
    // desired error reduction factor theta < 1/3
    //const double theta = 2.0/7.0;
    const double theta = 0.333;
    cout << "Rich_SOLVE: theta=" << theta << endl;

    // compute minimal K such that 3*rho^K < theta
    const int K = (int) ceil(log(theta/3.0) / log(rho));
    cout << "Rich_SOLVE: K=" << K << endl;
    
    u_epsilon.clear();

    map<double,double> log_10_residual_norms;
    map<double,double> degrees_of_freedom;
    map<double,double> asymptotic;
    map<double,double> time_asymptotic;
    map<double,double> log_10_L2_error;
    map<double,double> weak_ell_tau_norms;
    
    bool exit = 0;
    unsigned int loops = 0;

    double epsilon_k = nu;
    InfiniteVector<double,Index> f, v, Av;

    double time = 0;
    clock_t tstart, tend;
    tstart = clock();

//     EvaluateFrame<Basis1D,2,2> evalObj;
    
    
     while (epsilon_k > epsilon) {
      epsilon_k *= 3*pow(rho, K) / theta;
      cout << "Rich_SOLVE: epsilon_k=" << epsilon_k << endl;
      double eta = theta * epsilon_k / (6*omega*K);
      cout << "eta= " << eta << endl;
      cout << "Rich_SOLVE: eta=" << eta << endl;
      P.RHS(eta, f);
      for (int j = 1; j <= 1/*K*/; j++) {
	APPLY_COARSE(P, v, eta, Av, 1.0e-6, jmax, CDD1);

	v += omega * (f - Av);

	++loops;
	tend = clock();
	time += (double)(tend-tstart)/CLOCKS_PER_SEC;
	
	// ############ output #############
	P.RHS(1.0e-6, f);
 	APPLY(P, v, 1.0e-6, Av, jmax, CDD1);
  	double residual_norm = l2_norm(f - Av);
 	double tmp1 = log10(residual_norm);
	cout << "current residual error ||f-Av||=" << residual_norm << endl;
	
	asymptotic[log10( (double)v.size() )] = tmp1;
	time_asymptotic[log10(time)] = tmp1;


#ifdef ONE_D
	weak_ell_tau_norms[loops] = v.weak_norm(1./2.5); // d=3, dT=3
#endif
#ifdef TWO_D
	weak_ell_tau_norms[loops] = v.weak_norm(1./1.5); // d=3, dT=3
#endif


	cout << "active indices: " << v.size() << endl;
 	cout << "loop: " << loops << endl;


	int d  = Basis1D::primal_polynomial_degree();
	int dT = Basis1D::primal_vanishing_moments();
	char name1[128];
	char name2[128];
	char name3[128];

#ifdef ONE_D
	sprintf(name1, "%s%d%s%d%s", "./Richardson_results33_alpha_0p15/rich1D_asymptotic_P_jmax18_d", d, "_dT", dT, ".m");
	sprintf(name2, "%s%d%s%d%s", "./Richardson_results33_alpha_0p15/rich1D_time_asymptotic_P_jmax18_d", d, "_dT", dT, ".m");
	sprintf(name3, "%s%d%s%d%s", "./Richardson_results33_alpha_0p15/rich1D_weak_ell_tau_norms_P_jmax18_d", d, "_dT", dT, ".m");
#endif	

#ifdef TWO_D
	sprintf(name1, "%s%d%s%d%s", "./Richardson_results33_2D_alpha_0p05/rich2D_asymptotic_P_jmax18_d", d, "_dT", dT, ".m");
	sprintf(name2, "%s%d%s%d%s", "./Richardson_results33_2D_alpha_0p05/rich2D_time_asymptotic_P_jmax18_d", d, "_dT", dT, ".m");
	sprintf(name3, "%s%d%s%d%s", "./Richardson_results33_2D_alpha_0p05/rich2D_weak_ell_tau_norms_P_jmax18_d", d, "_dT", dT, ".m");
#endif

	std::ofstream os1(name1);
	matlab_output(asymptotic,os1);
	os1.close();
	
	std::ofstream os2(name2);
	matlab_output(time_asymptotic,os2);
	os2.close();

	std::ofstream os3(name3);
	matlab_output(weak_ell_tau_norms,os3);
	os3.close();

	// ############ end output #############	
	tstart = clock();

#ifdef ONE_D
	if (residual_norm < 3.1623e-04 || loops == 5000) {
#endif
#ifdef TWO_D
	if (residual_norm < 0.01 || loops == 5000) {
#endif
	  u_epsilon = v;
	  exit = true;
	  break;
	}
      }
      //v.COARSE((1-theta)*epsilon_k, u_epsilon);
      
      if (exit)
	break;
      
    } 
    
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
