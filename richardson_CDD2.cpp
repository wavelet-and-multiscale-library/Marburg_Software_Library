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
    typedef DSBasis<3,5> Basis1D;
    //typedef PBasis<2,2> Basis1D;
    //     Singularity1D_RHS_2<double> sing1D;
    //     Singularity1D_2<double> exactSolution1D;
     
//     Point<2> origin;
//     origin[0] = 0.0;
//     origin[1] = 0.0;

//     CornerSingularity sing2D(origin, 0.5, 1.5);
//     CornerSingularityRHS singRhs(origin, 0.5, 1.5);
     
    const int jmax = 13;

    double nu = P.norm_Ainv()*P.F_norm();
    //nu *= 10;
    //double nu = 50;
    cout << "nu = " << nu << endl;
    typedef typename PROBLEM::Index Index;

    // compute optimal relaxation parameter omega
    //const double omega = 2.0 / (P.norm_A() + 1.0/P.norm_Ainv());
    //const double omega = 2.0/5.0048-0.01;//0.337721577225;//0.476369621181
    
    //13.09.06
    //const double omega = 0.15;
    //22.09.06
    const double omega = 0.2;
    //const double omega = 2./4.;//0.337721577225;//0.476369621181
    //const double omega = 0.2;//that was the bad alpha estimate
    //const double omega = 0.1;// 2D good one
    //const double omega = 2.0/10;
    //const double omega = 2.0/100;
    cout << "CDD2_SOLVE: omega=" << omega << endl;



    //const double omega = 0.2;



    // compute spectral norm rho
    double cond_A = P.norm_A() * P.norm_Ainv();
    //double cond_A = 100;
    cout << "cond A = " << cond_A << endl;
    //cond_A = 32.8823;
    const double rho = (cond_A-1.) / (cond_A+1.);
    //const double rho = 0.92;
    cout << "CDD2_SOLVE: rho=" << rho << endl;
    
    // desired error reduction factor theta < 1/3
    //     const double theta = 2.0/7.0;
    const double theta = 2.0/7.0;
    cout << "CDD2_SOLVE: theta=" << theta << endl;

    // compute minimal K such that 3*rho^K < theta
    const int K = (int) ceil(log(theta/2.0) / log(rho))+1;
    //const int K = 488;
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

    InfiniteVector<double,Index> f, v, Av, resid, help, help2, r_exact;

    unsigned int innerloop = 0;

    //EvaluateFrame<Basis1D,2,2> evalObj;

    double time = 0;
    clock_t tstart, tend;
    tstart = clock();

    double acctime = 0;

    double d = pow(rho,K)*nu/(2.0*omega*K);


    double tmp = 0;
    double tmp1 = 0;
    //while (nu > epsilon) {
    while (true) {
      for (int j = 1; j <= K; j++) {
	//APPLY(P, v, pow(rho,j)*nu/(2.0*omega*K), Av, jmax,  CDD1);
	double eta = 0;
	if (innerloop == 0)
	  eta = pow(rho,K)*nu/(2.0*omega*K);
	else
	  eta = pow(rho,K)*nu/(2.0*omega*K);
	P.RHS(eta, f);
	cout << "CDD2_SOLVE: eta = " << pow(rho,K)*nu/(2.0*omega*K) << endl;
	//tend = clock();
	//time = (double)(tend-tstart)/CLOCKS_PER_SEC;
	APPLY_COARSE(P, v, eta, Av, 0.00000001, jmax, CDD1);
	//tend = clock(); 
	//acctime += ((double)(tend-tstart)/CLOCKS_PER_SEC - time);
	//cout << "time = " << acctime  << endl;

	//APPLY(P, v, eta, Av, jmax, CDD1);

	resid = f-Av;

	//approximate residual
	P.RHS(0., f);
	APPLY(P, v, 0., help, jmax, CDD1);
	r_exact = f - help;
	

	tmp = l2_norm(r_exact);
 	tmp1 = log10(tmp);
	
	tend = clock();

	time += (double)(tend-tstart)/CLOCKS_PER_SEC;
	time_asymptotic[log10(time)] = tmp1;

	cout << "current residual error ||f-Av||=" << tmp << endl;

	loops++;
	//########### aswymptotics ###############
	if (v.size() != 0)
	  asymptotic[log10( (double)v.size() )] = tmp1;
	//########################################

	v += omega * resid;

	u_epsilon = v;
	//u_epsilon.scale(&P,-1);

	//	if ((loops < 20) || loops % 10 == 0) {
	//	  u_epsilon = v;
	//	  P.rescale(u_epsilon,-1);
	//	  cout << "computing L_2 error..." << endl;
	// 	  double L2err = evalObj.L_2_error(P.basis(), u_epsilon, sing2D, 5, 0.0, 1.0);
	// 	  log_10_L2_error[loops] = log10(L2err);
	//cout << "L_2 error = " << L2err << endl;
	//	}
	// 	asymptotic[log10( (double)v.size() )] = tmp1;
	// 	log_10_residual_norms[loops] = tmp1;

	cout << "active indices: " << v.size() << endl;
	cout << "loop: " << loops << endl;
	if (loops % 1 == 0) {
	  std::ofstream os3("rich_asymptotic_2D_2209_35_al_0_2_0_98931213.m");
	  matlab_output(asymptotic,os3);
	  os3.close();

	  std::ofstream os4("time_asymptotic_rich_2D_35_2209_al_0_2_0_98931213.m");
	  matlab_output(time_asymptotic,os4);
	  os4.close();
	}
	tstart = clock();

 	if ((tmp < 0.0001 || loops == 2000) && innerloop) {
 	  u_epsilon = v;
  	  exit = true;
  	  break;
  	}
	if (innerloop >= 0) {
	  //u_epsilon = v;
	  break;
	}
      }
  
      cout << "#######################" << endl;
      cout << "exiting inner loop" << endl;
      //cout << "nu = " << 2.0*pow(rho,K)*nu / theta << endl;
      cout << "#######################" << endl;
      
      ++innerloop;
      
      //nu = 2.0*pow(rho,K)*nu / theta;
      //nu *= 0.98779464;
      nu *= 0.98931213;
      //nu *= 0.99083194;
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

    std::ofstream os3("richCDD2_asymptotic.m");
    matlab_output(asymptotic,os3);
    os3.close();

    std::ofstream os4("time_asymptotic_rich.m");
    matlab_output(time_asymptotic,os4);
    os4.close();


  }
    

}
