// implementation for richardson.h

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
#include <frame_index.h>

#define PROJECTOR

//using std::set;
using namespace std;
using namespace FrameTL;
using namespace MathTL;
using namespace WaveletTL;

namespace FrameTL
{
  template <class PROBLEM>
  void richardson_SOLVE_projector(const PROBLEM& P, const double epsilon,
			InfiniteVector<double, typename PROBLEM::Index>& u_epsilon,
			Array1D<InfiniteVector<double, typename PROBLEM::Index> >& approximations)
  {


    const int d  = PRIMALORDER;
    const int dt = DUALORDER;

    
    // Maximal number of Iterations within the projection step
    const int maxAnzahlIt = 200;

    const int DIM = 1;
    

    double tol_apply = 10e-6;
    double tol;

    typedef PBasis<d,dt> Basis1D;

    const unsigned int jmax = JMAX;

    const double nu = P.norm_Ainv()*P.F_norm();
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
    const double omega = 0.05;
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

    // Number of iterations after which the projection step will be applied
    // in theory, this should be K
    const int projectionFreq = 3*K;

    // constants for the projection step
    double omega_proj = omega;
    const double rho_proj = rho;
    const double theta_proj = theta;
    double epsilon_proj;


    // get number of patches - should be 2 anyway
    const int m = P.get_problem()->frame().n_p();
    cout << "m=" << m << endl;

    // prepare indices for patches - needed for clip
    std:set<Index> patch[m];
    for (Index iterator = P.get_problem() -> frame().first_generator( P.get_problem()->frame().j0() );; ++iterator)
    {
        patch[iterator.p()].insert(iterator);
        if (iterator == P.get_problem()->frame().last_wavelet(JMAX)) break;
    }


    // Setup problem for left-hand-size of (6.2.17)
    ProjectorEquation <Basis1D, DIM> lhs_equation(&P.get_problem()->get_bvp(), &P.get_problem()->frame(), jmax, false, 0);
    CachedProblem<ProjectorEquation<Basis1D, DIM> >  lhs_problem(&lhs_equation, 1.0, 1.0);

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
    InfiniteVector<double,Index> f, v, Av, RHS, temp, u_scaled, residual;
    InfiniteVector<double,Index> vi[m];

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
      P.RHS(eta, f);
      for (int j = 1; j <= 1;/*K;*/ j++) {

        APPLY(P, v, eta, Av, jmax, CDD1);

        v += omega * (f - Av);

        v.COARSE(eta, temp);
        v = temp;

	++loops;
	tend = clock();
	time += (double)(tend-tstart)/CLOCKS_PER_SEC;

	// ############ output #############
	P.RHS(1.0e-6, f);
 	APPLY(P, v, 1.0e-6, Av, jmax, CDD1);
  	residual_norm = l2_norm(f - Av);
 	double tmp1 = log10(residual_norm);
	cout << "current residual error ||f-Av||=" << residual_norm << endl;

	asymptotic[log10( (double)v.size() )] = tmp1;
	time_asymptotic[log10(time)] = tmp1;

	cout << "active indices: " << v.size() << endl;
 	cout << "loop: " << loops << endl;

	// ############ end output #############
	tstart = clock();

	if (residual_norm < 3.1623e-04 || loops == 5000 || epsilon_k <= epsilon) {
	  u_epsilon = v;
	  exit = true;
	  break;
	}
      }

#ifdef PROJECTOR


      if (loops % K == 0){
          cout << "START PROJECTION STEP" << endl;

          for(int i=0; i<m; i++)
          {

              // tolerance for the projection step on each patch

              epsilon_proj = theta*epsilon_k/(3*m);
              cout << "epsilon_proj = " << epsilon_proj << endl;

              // scale result from solver
              u_scaled = v;
              u_scaled.scale(&P, -1);
              if (i==0)
                u_scaled.scale(&rhs_problem_0, 1);
              else if (i==1)
                u_scaled.scale(&rhs_problem_1, 1);


              // prepare RHS for patch 0 and 1
              if (i == 0)
                APPLY(rhs_problem_0, u_scaled, tol_apply, RHS, jmax, CDD1);
              else if(i == 1)
                APPLY(rhs_problem_1, u_scaled, tol_apply, RHS, jmax, CDD1);

              RHS.clip(patch[i]);

              // scale right-hand-side
              if (i==0)
                RHS.scale(&rhs_problem_0, 1);
              else if (i==1)
                RHS.scale(&rhs_problem_1, 1);

              RHS.scale(&lhs_problem, -1);

              vi[i] = v;
              vi[i].clip(patch[i]);
              double epsilon_proj_k = l2_norm(vi[i]);
              vi[i].scale(&P, -1);
              vi[i].scale(&lhs_problem, 1);


              int K_proj = (int) ceil(log(theta_proj/3.0) / log(rho_proj));
              cout << "Rich_SOLVE: K_proj=" << K_proj << endl;

              tol = tol_apply;
              RHS.clip(patch[i]);

              // due to zeros on the main diagonal in the matrix needed to set up the RHS, we need to remove nan-entries here
              for(typename InfiniteVector<double, Index>::const_iterator it = RHS.begin(); it != RHS.end(); ++it)
              {
                  if(RHS[it.index()] != RHS[it.index()]) RHS[it.index()] = 0;
              }


              int count = 0;

              while(true)
              {
                  
                  count++;
                  double eta_proj = theta_proj * epsilon_proj_k / (6*omega_proj*K_proj);
                  epsilon_proj_k *= 3 * pow(rho_proj,K_proj) / theta_proj;

                  APPLY(lhs_problem, vi[i], 10e-6, Av, jmax, CDD1);
                  Av.clip(patch[i]);
                  residual = RHS-Av;
                  vi[i] += omega_proj*(residual);
                  vi[i].COARSE(10e-6,temp);
                  vi[i] = temp;

                  if (epsilon_proj_k < epsilon_proj) break;
                 
              }

              cout << "Projection step on patch " << i << ": " << count << " iterations" << endl;

              vi[i].scale(&lhs_problem,-1);

          }

          // collect local parts
          v.clear();
          for(int i=0; i<m; i++) v += vi[i];

          // scale back
          v.scale(&P,1);

          // extra COARSE step after scaling to reduce degrees of freedom

          v.COARSE(eta,temp);
          temp = v;

          cout << "PROJECTION STEP DONE" << endl;
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

  }




}
