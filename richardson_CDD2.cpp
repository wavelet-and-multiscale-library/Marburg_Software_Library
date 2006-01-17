// implementation for richardson.h

#include <cmath>
#include <set>
#include <utils/plot_tools.h>
#include <adaptive/apply.h>
#include <frame_evaluate.h>
#include <numerics/corner_singularity.h>

using std::set;

namespace FrameTL
{


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
  void richardson_SOLVE_CDD2(const PROBLEM& P, const double epsilon,
			InfiniteVector<double, typename PROBLEM::Index>& u_epsilon)
  {
    typedef DSBasis<2,2> Basis1D;
    //     Singularity1D_RHS_2<double> sing1D;
    //     Singularity1D_2<double> exactSolution1D;
     
//     Point<2> origin;
//     origin[0] = 0.0;
//     origin[1] = 0.0;

//     CornerSingularity sing2D(origin, 0.5, 1.5);
//     CornerSingularityRHS singRhs(origin, 0.5, 1.5);
     
    const int jmax = 10;

    double nu = P.norm_Ainv()*P.F_norm();
    cout << "nu = " << nu << endl;
    typedef typename PROBLEM::Index Index;

    // compute optimal relaxation parameter omega
    //const double omega = 2.0 / (P.norm_A() + 1.0/P.norm_Ainv());
    //const double omega = 2.0/5.0048-0.01;//0.337721577225;//0.476369621181
    const double omega = 0.52;
    //const double omega = 2.0/10;
    //const double omega = 2.0/100;
    cout << "CDD2_SOLVE: omega=" << omega << endl;

    // compute spectral norm rho
    const double cond_A = P.norm_A() * P.norm_Ainv();
    const double rho = (cond_A-1.) / (cond_A+1.);
    cout << "CDD2_SOLVE: rho=" << rho << endl;
    
    // desired error reduction factor theta < 1/3
    //     const double theta = 2.0/7.0;
    const double theta = 2.0/7.0;
    cout << "CDD2_SOLVE: theta=" << theta << endl;

    // compute minimal K such that 3*rho^K < theta
    const int K = (int) ceil(log(theta/2.0) / log(rho))+1;
    //const int K = 100;
    cout << "CDD2_SOLVE: K=" << K << endl;
    
    u_epsilon.clear();

    map<double,double> log_10_residual_norms;
    map<double,double> log_10_L2_error;
    map<double,double> degrees_of_freedom;
    map<double,double> asymptotic;
    map<double,double> time_asymptotic;
    
    double min_error = 100.0;

    bool exit = 0;
    unsigned int loops = 0;

    InfiniteVector<double,Index> f, v, Av, res, help, help2;

    unsigned int innerloop = 0;

    //EvaluateFrame<Basis1D,2,2> evalObj;


    double time = 0;
    clock_t tstart, tend;
    tstart = clock();

    double acctime = 0;

    double d = pow(rho,K)*nu/(2.0*omega*K);


    double tmp = 0;
    double tmp1 = 0;
    while (nu > epsilon) {
 
      for (int j = 1; j <= K; j++) {


	//APPLY(P, v, pow(rho,j)*nu/(2.0*omega*K), Av, jmax,  CDD1);
	double eta = 0;
	if (innerloop == 0)
	  eta = pow(rho,j)*nu/(2.0*omega*K);
	else
	  eta = d;
	P.RHS(eta, f);
	cout << "CDD2_SOLVE: eta = " << pow(rho,j)*nu/(2.0*omega*K) << endl;
	tend = clock();
	time = (double)(tend-tstart)/CLOCKS_PER_SEC;
	APPLY_COARSE(P, v, eta, Av, 0.00000001, jmax, CDD1);
	tend = clock();
	acctime += ((double)(tend-tstart)/CLOCKS_PER_SEC - time);
	cout << "time = " << acctime  << endl;

	//APPLY(P, v, eta, Av, jmax, CDD1);

	res = f - Av;
	tmp = l2_norm(res);
	tmp1 = log10(tmp);
	cout << "current residual error ||f-Av||=" << tmp << endl;
	loops++;
	//########### aswymptotics ###############
	// 	if (j == K) {
	// 	  if (tmp1 < min_error)
	// 	    min_error = tmp1;
	  
	// 	  log_10_residual_norms[loops] = min_error;

	// 	  tend = clock();
	// 	  time = (double)(tend-tstart)/CLOCKS_PER_SEC;
	// 	  time_asymptotic[log10(time)] = tmp1;
	//  	  if (tmp < 0.5e-1) {
	//  	    u_epsilon = v;
	//  	    exit = true;
	//  	    break;
	//	}  
	//	}
	// 	if (v.size() != 0)
	// 	  asymptotic[log10( (double)v.size() )] = tmp1;
	//########################################

	v += omega * res;

	 u_epsilon = v;
	 P.rescale(u_epsilon,-1);

	//	if ((loops < 20) || loops % 10 == 0) {
	//	  u_epsilon = v;
	//	  P.rescale(u_epsilon,-1);
	//	  cout << "computing L_2 error..." << endl;
	// 	  double L2err = evalObj.L_2_error(P.basis(), u_epsilon, sing2D, 5, 0.0, 1.0);
	// 	  log_10_L2_error[loops] = log10(L2err);
	//cout << "L_2 error = " << L2err << endl;
	//	}
	//asymptotic[log10( (double)v.size() )] = tmp1;
	log_10_residual_norms[loops] = tmp1;

	cout << "active indices: " << v.size() << endl;
	cout << "loop: " << loops << endl;
	// 	if (loops % 10 == 0) {
	// 	  std::ofstream os3("richCDD2_2D_asymptotic.m");
	// 	  matlab_output(asymptotic,os3);
	// 	  os3.close();

	// 	  std::ofstream os4("richCDD2_2D_L2_errors.m");
	// 	  matlab_output(log_10_L2_error,os4);
	// 	  os4.close();
	// 	}
 	if (tmp < 0.004 || loops == 300) {
 	  u_epsilon = v;
  	  exit = true;
  	  break;
  	}
	// 	if (innerloop >= 1) {
	// 	  //u_epsilon = v;
	//   	  break;
	// 	}


	// 	if (((loops % 10 == 0) && loops <= 100) || loops % 100 == 0){
	// 	  u_epsilon = v;
	// 	  P.rescale(u_epsilon,-1);
	// 	  char filename1[50];
	// 	  char filename2[50];
	  
	// 	  sprintf(filename1, "%s%d%s%d%s", "approx_sol_rich33_1D_out_", loops, "_nactive_", v.size(),".m");
	// 	  sprintf(filename2, "%s%d%s%d%s", "error_rich_1D33_out_", loops, "_nactive_", v.size(),".m");

	// 	  EvaluateFrame<Basis1D,1,1> evalObj;
	  
	// 	  Array1D<SampledMapping<1> > U = evalObj.evaluate(P.basis(), u_epsilon, true, 11);//expand in primal basis
	// 	  cout << "...plotting approximate solution" << endl;
	// 	  Array1D<SampledMapping<1> > Error = evalObj.evaluate_difference(P.basis(), u_epsilon, exactSolution1D, 11);
	// 	  cout << "...plotting error" << endl;
	// 	  std::ofstream ofs5(filename1);
	// 	  matlab_output(ofs5,U);
	// 	  ofs5.close();
	  
	// 	  std::ofstream ofs6(filename2);
	// 	  matlab_output(ofs6,Error);
	// 	  ofs6.close();

	// 	}


      }
  
      cout << "#######################" << endl;
      cout << "exiting inner loop" << endl;
      //cout << "nu = " << 2.0*pow(rho,K)*nu / theta << endl;
      cout << "#######################" << endl;
      
      ++innerloop;
      
      nu = 2.0*pow(rho,K)*nu / theta;
      if (exit)
	break;
      //v.COARSE((1-theta)*nu, v);
      
      

    }
    
    set<Index> Lambda;
    
    std::ofstream os1("residual_norms_rich.m");
    matlab_output(log_10_residual_norms,os1);
    os1.close();
    
    //     std::ofstream os2("degrees_of_freedom.m");
    //     matlab_output(degrees_of_freedom,os2);
    //     os2.close();

    //     std::ofstream os3("richCDD2_asymptotic.m");
    //     matlab_output(asymptotic,os3);
    //     os3.close();

    std::ofstream os4("time_asymptotic_rich.m");
    matlab_output(time_asymptotic,os4);
    os4.close();


  }


    

}
