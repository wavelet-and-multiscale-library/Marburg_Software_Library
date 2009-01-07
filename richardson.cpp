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
  void richardson_SOLVE(const PROBLEM& P, const double epsilon,
			InfiniteVector<double, typename PROBLEM::Index>& u_epsilon,
			Array1D<InfiniteVector<double, typename PROBLEM::Index> >& approximations)
  {

    typedef PBasis<3,3> Basis1D;
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

    cout << "CDD2_SOLVE: nu=" << nu << endl;

    // compute optimal relaxation parameter omega
    //const double omega = 2.0 / (P.norm_A() + 1.0/P.norm_Ainv());
    const double omega = 0.4;
    //const double omega = 0.15;
    //const double omega = 2.0/2.47-0.5;
    cout << "CDD2_SOLVE: omega=" << omega << endl;

    // compute spectral norm rho
    const double cond_A = P.norm_A() * P.norm_Ainv();
    //const double rho = (cond_A - 1.0) / (cond_A + 1.0);
    const double rho = 0.8;
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

//     EvaluateFrame<Basis1D,2,2> evalObj;
    
    
     while (epsilon_k > epsilon) {
      epsilon_k *= 3*pow(rho, K) / theta;
      cout << "CDD2_SOLVE: epsilon_k=" << epsilon_k << endl;
      double eta = theta * epsilon_k / (6*omega*K);
      cout << "eta= " << eta << endl;
      cout << "CDD2_SOLVE: eta=" << eta << endl;
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

	cout << "active indices: " << v.size() << endl;
 	cout << "loop: " << loops << endl;


	int d  = Basis1D::primal_polynomial_degree();
	int dT = Basis1D::primal_vanishing_moments();
	char name1[60];
	char name2[60];
	char name3[60];
	
	sprintf(name1, "%s%d%s%d%s", "./Richardson_results/rich1D_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name2, "%s%d%s%d%s", "./Richardson_results/rich1D_time_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	
	std::ofstream os1(name1);
	matlab_output(asymptotic,os1);
	os1.close();
	
	std::ofstream os2(name2);
	matlab_output(time_asymptotic,os2);
	os2.close();

	// ############ end output #############	
	tstart = clock();

	if (residual_norm < 1.0e-3 || loops == 1000) {
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
