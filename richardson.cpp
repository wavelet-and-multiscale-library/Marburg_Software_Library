// implementation for richardson.h

#include <cmath>
#include <set>
#include <utils/plot_tools.h>
#include <adaptive/apply.h>
#include <numerics/corner_singularity.h>
#include <frame_evaluate.h>

using std::set;

namespace FrameTL
{
  template <class PROBLEM>
  void richardson_SOLVE(const PROBLEM& P, const double epsilon,
			InfiniteVector<double, typename PROBLEM::Index>& u_epsilon)
  {

    typedef DSBasis<2,2> Basis1D;
//     Singularity1D_RHS_2<double> sing1D;
//     Singularity1D_2<double> exactSolution1D;
     
//      Point<2> origin;
//      origin[0] = 0.0;
//      origin[1] = 0.0;

//      CornerSingularity sing2D(origin, 0.5, 1.5);
//      CornerSingularityRHS singRhs(origin, 0.5, 1.5);

    const unsigned int jmax = 4;

    const double nu = P.norm_Ainv()*P.F_norm();
    typedef typename PROBLEM::Index Index;

    cout << "CDD2_SOLVE: nu=" << nu << endl;

    // compute optimal relaxation parameter omega
    //const double omega = 2.0 / (P.norm_A() + 1.0/P.norm_Ainv());
    const double omega = 2.0/2.47-0.5;
    cout << "CDD2_SOLVE: omega=" << omega << endl;

    // compute spectral norm rho
    const double cond_A = P.norm_A() * P.norm_Ainv();
    const double rho = (cond_A - 1.0) / (cond_A + 1.0);
    cout << "CDD2_SOLVE: rho=" << rho << endl;
    
    // desired error reduction factor theta < 1/3
    //const double theta = 2.0/7.0;
    const double theta = 0.333;
    cout << "CDD2_SOLVE: theta=" << theta << endl;

    // compute minimal K such that 3*rho^K < theta
    const int K = (int) ceil(log(theta/3.0) / log(rho));
    cout << "CDD2_SOLVE: K=" << K << endl;
    
    u_epsilon.clear();

    map<double,double> log_10_residual_norms;
    map<double,double> degrees_of_freedom;
    map<double,double> asymptotic;
    map<double,double> time_asymptotic;
    map<double,double> log_10_L2_error;
    
    bool exit = 0;
    unsigned int loops = 0;

    double epsilon_k = nu;
    InfiniteVector<double,Index> f, v, Av;

    double time = 0;
    clock_t tstart, tend;
    tstart = clock();

    EvaluateFrame<Basis1D,2,2> evalObj;
    
    double eta = theta * epsilon_k / (6*omega*K);
    while (epsilon_k > epsilon) {
      epsilon_k *= 3*pow(rho, K) / theta;
      cout << "CDD2_SOLVE: epsilon_k=" << epsilon_k << endl;
      double eta = theta * epsilon_k / (6*omega*K);
      P.RHS(eta, f);
      //cout << f << endl;
      //v = u_epsilon;
      cout << "eta= " << eta << endl;
      for (int j = 1; j <= 1; j++) {
        //eta = pow(rho,j) * epsilon_k / (6*omega*K);
	cout << "CDD2_SOLVE: eta=" << eta << endl;
	//APPLY(P, v, 0.000001/*eta*/, Av, jmax,  CDD1);
	//P.RHS(eta, f);
	APPLY_COARSE(P, v, theta * epsilon_k / (6*omega*K), Av, 0.00000001, jmax, CDD1);
	double tmp = l2_norm(f - Av);
	cout << "current residual error ||f-Av||=" << tmp << endl;
	v += omega * (f - Av);

	++loops;

// 	if ( loops == 1 || loops == 2 || loops == 4 || loops == 6 || loops == 8
// 	     || loops == 10 || (loops % 20 == 0)) {
// 	  u_epsilon = v;
// 	  P.rescale(u_epsilon,-1);
// 	  //cout << "computing L_2 error..." << endl;
// 	  //double L2err = evalObj.L_2_error(P.basis(), u_epsilon, sing2D, 5, 0.0, 1.0);
// 	  //log_10_L2_error[loops] = log10(L2err);
// 	  //cout << "L_2 error = " << L2err << endl;
// 	}

	double tmp1 = log10(tmp);
	asymptotic[log10( (double)v.size() )] = tmp1;
	log_10_residual_norms[loops] = tmp1;
	tend = clock();
	time += (double)(tend-tstart)/CLOCKS_PER_SEC;
	time_asymptotic[log10(time)] = tmp1;
	degrees_of_freedom[loops] = v.size();
	cout << "active indices: " << v.size() << endl;
	cout << "loop: " << loops << endl;

	std::ofstream os3("rich_orig_2D_asymptotic_3101.m");
	matlab_output(asymptotic,os3);
	os3.close();
	
// 	std::ofstream os4("rich_orig_2D_L2_errors.m");
// 	matlab_output(log_10_L2_error,os4);
// 	os4.close();
	
	std::ofstream os5("rich_orig_2D_time_asymptotic3101.m");
	matlab_output(asymptotic,os5);
	os5.close();

	if (tmp < 0.004 || loops == 400) {
	  u_epsilon = v;
	  exit = true;
	  break;
	}
      }
      if (exit)
	break;
      //v.COARSE((1-theta)*epsilon_k, u_epsilon);
    } 
    
    set<Index> Lambda;
    
    std::ofstream os1("residual_norms.m");
    matlab_output(log_10_residual_norms,os1);
    os1.close();
    
    std::ofstream os2("degrees_of_freedom.m");
    matlab_output(degrees_of_freedom,os2);
    os2.close();

    std::ofstream os3("rich_asymptotic.m");
    matlab_output(asymptotic,os3);
    os3.close();

    std::ofstream os4("time_asymptotic.m");
    matlab_output(time_asymptotic,os4);
    os4.close();


  }


    

}
