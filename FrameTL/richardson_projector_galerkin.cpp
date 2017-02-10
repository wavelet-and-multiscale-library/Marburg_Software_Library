// implementation for richardson_projector_galerkin.h

#include <cmath>
#include <set>
#include <utils/plot_tools.h>
#include <adaptive/apply.h>
#include <numerics/corner_singularity.h>
#include <frame_evaluate.h>
#include <poisson_1d_testcase.h>
#include <projector_equation.h>
#include <projector_utils.h>
#include <galerkin/cached_problem.h>
#include <algebra/infinite_vector.h>
#include <richardson.h>
#include <frame_index.h>
#include <lhs_equation.h>
#include <cdd1_local.h>

//#undef SINGULARITY_RHS

using std::set;
using namespace std;
using namespace FrameTL;
using namespace MathTL;
using namespace WaveletTL;
using FrameTL::EvaluateFrame;

namespace FrameTL
{
  template <class PROBLEM>
  void richardson_SOLVE_projector_galerkin(const PROBLEM& P, const double epsilon,
			InfiniteVector<double, typename PROBLEM::Index>& u_epsilon,
			Array1D<InfiniteVector<double, typename PROBLEM::Index> >& approximations)
  {

    const int d  = PRIMALORDER;
    const int dt = DUALORDER;

    
#ifdef ONE_D
    const int DIM = 1;
#else
    const int DIM = 2;
#endif

    typedef PBasis<d,dt> Basis1D;


    const unsigned int jmax = JMAX;

#ifdef SINGULARITY_RHS
    double nu = P.norm_Ainv();
#else
    const double nu = P.norm_Ainv()*P.F_norm();
#endif
    typedef typename PROBLEM::Index Index;
    typedef typename InfiniteVector<double, Index>::const_iterator Iter;

    cout << "Rich_SOLVE: nu=" << nu << endl;

    // compute optimal relaxation parameter omega
    //const double omega = 2.0 / (P.norm_A() + 1.0/P.norm_Ainv());

    // ####### 1D #######
    //const double omega = 0.4;
    //const double omega = 0.15;

    // ####### 2D #######
    // bad one
#ifdef SINGULARITY_RHS
    const double omega = 0.25;
#else
    const double omega = 0.05;
#endif

    // good one
    //const double omega = 0.25;

    //const double omega = 0.3;


    //const double omega = 2.0/2.47-0.5;
    cout << "Rich_SOLVE: omega=" << omega << endl;

    // compute spectral norm rho
    // const double cond_A = P.norm_A() * P.norm_Ainv();
    // const double rho = (cond_A - 1.0) / (cond_A + 1.0);
    const double rho = 0.9;
    cout << "Rich_SOLVE: rho=" << rho << endl;

    // desired error reduction factor theta < 1/3
    //const double theta = 2.0/7.0;
    const double theta = 0.333;
    cout << "Rich_SOLVE: theta=" << theta << endl;

    // compute minimal K such that 3*rho^K < theta
    const int K = (int) ceil(log(theta/3.0) / log(rho));
    cout << "Rich_SOLVE: K=" << K << endl;

    // get number of patches - should be 2 anyway
    const int m = P.get_problem()->frame().n_p();
    cout << "m=" << m << endl;

    // prepare indices for patches - needed for clip
    set<Index> patch[m];
    for (Index iterator = P.get_problem() -> frame().first_generator( P.get_problem()->frame().j0() );; ++iterator)
    {
        patch[iterator.p()].insert(iterator);
        if (iterator == P.get_problem()->frame().last_wavelet(JMAX)) break;
    }


    // Setup problem for left-hand-size of (6.2.17)
    LHS_Equation <Basis1D, DIM> lhs_equation_0(&P.get_problem()->get_bvp(), &P.get_problem()->frame(), jmax, 0);
    CachedProblemLocal<LHS_Equation<Basis1D, DIM> >  lhs_problem_0(&lhs_equation_0, 1.0, 1.0);

    // Setup problem for left-hand-size of (6.2.17)
    LHS_Equation <Basis1D, DIM> lhs_equation_1(&P.get_problem()->get_bvp(), &P.get_problem()->frame(), jmax, 1);
    CachedProblemLocal<LHS_Equation<Basis1D, DIM> >  lhs_problem_1(&lhs_equation_1, 1.0, 1.0);

    // Problem for calculation of RHS in (6.2.17)
    ProjectorEquation <Basis1D, DIM> rhs_equation_0(&P.get_problem()->get_bvp(), &P.get_problem()->frame(), jmax, true, 0);
    CachedProblem<ProjectorEquation<Basis1D, DIM> >  rhs_problem_0(&rhs_equation_0, 1.0, 1.0);

    // Problem for calculation of RHS in (6.2.17)
    ProjectorEquation <Basis1D, DIM> rhs_equation_1(&P.get_problem()->get_bvp(), &P.get_problem()->frame(), jmax, true, 1);
    CachedProblem<ProjectorEquation<Basis1D, DIM> >  rhs_problem_1(&rhs_equation_1, 1.0, 1.0);

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
    InfiniteVector<double,Index> f, v, Av, RHS, temp, u_scaled, start, null;
    InfiniteVector<double,Index> vi[m];

// in case we use a singularity function as right-hand-side
#ifdef SINGULARITY_RHS
    InfiniteVector<double,Index> sing_RHS;
    MultiIndex<int, 2> e;
    e[0] = 1; e[1] = 1;
    MultiIndex<int, 2> k;
    k[0] = 28; k[1] = 20;
    Index lambda(&(P.get_problem()->frame()),6,e,1,k);
    //lambda = P.get_problem()->frame().last_wavelet(JMAX);
    cout << "Lambda = " << lambda << endl;
    sing_RHS[lambda] = 100;
    sing_RHS.scale(&P,-1);
    cout << "Singularity function prepared" << endl;
    nu *= l2_norm(sing_RHS);
#endif

    double residual_norm = 0;
    double eta;
    double time = 0;
    clock_t tstart, tend;
    tstart = clock();

//     EvaluateFrame<Basis1D,2,2> evalObj;


     while (epsilon_k > epsilon) {
      epsilon_k *= 3*pow(rho, K) / theta;
      cout << "Rich_SOLVE: epsilon_k=" << epsilon_k << endl;
      eta = theta * epsilon_k / (6*omega*K);
      cout << "eta= " << eta << endl;
      cout << "Rich_SOLVE: eta=" << eta << endl;
// bad hack for testing with constant RHS
#ifdef SINGULARITY_RHS
      f = sing_RHS;
#else
      P.RHS(eta, f);
#endif
      for (int j = 1; j <= 1;/*K;*/ j++) {

        APPLY(P, v, eta, Av, jmax, CDD1);

        v += omega * (f - Av);

        v.COARSE(eta, temp);
        v = temp;

	++loops;
	tend = clock();
	time += (double)(tend-tstart)/CLOCKS_PER_SEC;

	// ############ output #############
#ifdef SINGULARITY_RHS
        f = sing_RHS;
#else
        P.RHS(10e-6, f);
#endif
 	APPLY(P, v, 1.0e-6, Av, jmax, CDD1);
  	residual_norm = l2_norm(f - Av);
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
	sprintf(name1, "%s%d%s%d%s", "./Richardson_Proj_1D/rich1D_asymptotic_P_jmax18_d", d, "_dT", dT, ".m");
	sprintf(name2, "%s%d%s%d%s", "./Richardson_Proj_1D/rich1D_time_asymptotic_P_jmax18_d", d, "_dT", dT, ".m");
	sprintf(name3, "%s%d%s%d%s", "./Richardson_Proj_1D/rich1D_weak_ell_tau_norms_P_jmax18_d", d, "_dT", dT, ".m");
#endif

#ifdef TWO_D
	sprintf(name1, "%s%d%s%d%s", "./Richardson_Proj_2D/rich2D_asymptotic_P_jmax18_d", d, "_dT", dT, ".m");
	sprintf(name2, "%s%d%s%d%s", "./Richardson_Proj_2D/rich2D_time_asymptotic_P_jmax18_d", d, "_dT", dT, ".m");
	sprintf(name3, "%s%d%s%d%s", "./Richardson_Proj_2D/rich2D_weak_ell_tau_norms_P_jmax18_d", d, "_dT", dT, ".m");
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

	if (residual_norm < 3.1623e-04 || loops == 5000 || epsilon_k <= epsilon) {
	  u_epsilon = v;
	  exit = true;
	  break;
	}
      }

#ifdef PROJECTOR
#ifdef CONSTFUN
      if (false && residual_norm >= 0.005 && loops % (3*K) == 0){
#elif defined SINGULARITY_RHS
      if (residual_norm >= 0.005 && loops % (3*K) == 0){
#else
      if (residual_norm >= 0.05 && loops % (3*K) == 0){
#endif
          cout << "START PROJECTION STEP" << endl;

          // tolerance for the projection step on each patch, according to Stevenson
          double epsilon_proj = theta*epsilon_k/(3*m);
#ifdef TWO_D
          // Adjust parameter in 2D-case. Initial guess is too pessimistic
          // and leads to too many degrees of freedom
          // epsilon_proj *= 10e2;
#endif
          cout << "epsilon_proj = " << epsilon_proj << endl;

          for(int i=0; i<m; i++)
          {



              // scale result from solver
              u_scaled = v;
              u_scaled.scale(&P, -1);

              start = u_scaled;
              start.clip(patch[i]);
              if (i==0)
                lhs_equation_0.rescale(start,1);//start.scale(&lhs_problem_0, 1);
              else if (i==1)
                lhs_equation_1.rescale(start,1);//start.scale(&lhs_problem_1, 1);

              if (i==0)
                rhs_equation_0.rescale(u_scaled,1);//u_scaled.scale(&rhs_problem_0, 1);
              else if (i==1)
                rhs_equation_1.rescale(u_scaled,1);//u_scaled.scale(&rhs_problem_1, 1);
#ifdef ONE_D
              double tol_apply = 10e-12;
#else
 #ifdef CONSTFUN
              cout << "l2_norm(u_scaled) = " << l2_norm(u_scaled) << endl;
              double tol_apply = 10e-6;//epsilon_proj;
 #elif defined SINGULARITY_RHS
              double tol_apply = 10e-6;
 #else
              double tol_apply = 10e-6;
 #endif
#endif
              cout << "start APPLY for rhs..." << endl;
              // prepare RHS for patch 0 and 1
              if (i == 0)
                APPLY(rhs_problem_0, u_scaled, tol_apply, RHS, jmax, CDD1);
              else if(i == 1)
                APPLY(rhs_problem_1, u_scaled, tol_apply, RHS, jmax, CDD1);
              cout << "APPLY done!" << endl;

              RHS.clip(patch[i]);

              // scale right-hand-side
              if (i==0){
                rhs_equation_0.rescale(RHS,1);//RHS.scale(&rhs_problem_0, 1);
                lhs_equation_0.rescale(RHS,-1);//RHS.scale(&lhs_problem_0, -1);
              }
              else if (i==1){
                rhs_equation_1.rescale(RHS,1);//RHS.scale(&rhs_problem_1, 1);
                lhs_equation_1.rescale(RHS,-1);//RHS.scale(&lhs_problem_1, -1);
              }

              vi[i] = v;
              vi[i].scale(&P, -1);
              if (i==0)
                lhs_equation_0.rescale(vi[i],1);//vi[i].scale(&lhs_problem_0, 1);
              else if(i==1)
                lhs_equation_1.rescale(vi[i],1);//vi[i].scale(&lhs_problem_1, 1);
              

              cout << "RHS on patch " << i << " done." << endl;
              cout << "RHS.size() = " << RHS.size() << endl;

              // due to zeros on the main diagonal in the matrix needed to set up the RHS, we need to remove nan-entries here
              for(typename InfiniteVector<double, Index>::const_iterator it = RHS.begin(); it != RHS.end(); ++it)
              {
                  if(RHS[it.index()] != RHS[it.index()]) {cout << "nan!" << endl; RHS[it.index()] = 0;}
              }

              // solve local problems using the Galerkin method from CDD
              if (i==0){
                lhs_equation_0.set_rhs(RHS);
                #ifdef SINGULARITY_RHS
                    CDD1_LOCAL_SOLVE(lhs_problem_0, 0, epsilon_proj*10e-4/*10e-8*/, start , vi[i], null, jmax, CDD1);
                #else
                    CDD1_LOCAL_SOLVE(lhs_problem_0, 0, eta/*epsilon_proj*/, start , vi[i], null, jmax, CDD1);
                #endif
                lhs_equation_0.rescale(vi[i],-1);//vi[i].scale(&lhs_problem_0,-1);
              }
              else if (i==1){
                lhs_equation_1.set_rhs(RHS);
                #ifdef SINGULARITY_RHS
                    CDD1_LOCAL_SOLVE(lhs_problem_1, 1, epsilon_proj*10e-4/*10e-8*/, start , vi[i], null, jmax, CDD1);
                #else
                    CDD1_LOCAL_SOLVE(lhs_problem_1, 1, eta/*epsilon_proj*/ , start , vi[i], null, jmax, CDD1);
                #endif
                lhs_equation_1.rescale(vi[i],-1);//vi[i].scale(&lhs_problem_1,-1);
              }
              cout << "Projection on patch " << i << " done." << endl;
              
          }

          // collect local parts
          v.clear();
          for(int i=0; i<m; i++) v += vi[i];

          // scale back
          v.scale(&P,1);

#ifdef ONE_D
          // in 1D: extra COARSE step after scaling to reduce degrees of freedom
          v.COARSE(10*m*epsilon_proj,temp);
          v = temp;
#endif

          cout << "PROJECTION STEP DONE" << endl;
          cout << "v.size() = " << v.size() << endl;
      }
     
#endif

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
#ifdef SINGULARITY_RHS
    sing_RHS.scale(&P,1);
    EvaluateFrame<Basis1D,2,2> evalObj;
    Array1D<SampledMapping<2> > evaluate_RHS = evalObj.evaluate(P.get_problem()->frame(), sing_RHS, true, 6);//expand in primal basis
    std::ofstream ofs5("./Richardson_Proj_2D/singularity_rhs.m");
    matlab_output(ofs5,evaluate_RHS);
    ofs5.close();
#endif

  }




}
