// implementation for steepest_descent.h

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
  void GALERKIN (PROBLEM& P,
		 const set<typename PROBLEM::Index>& Lambda,
		 const InfiniteVector<double, typename PROBLEM::Index>& rhs,
		 Vector<double>& u)
  {
    typedef typename PROBLEM::Index Index;
    // setup stiffness matrix
    cout << "setting up full stiffness matrix..." << endl;
    SparseMatrix<double> A_Lambda;
    WaveletTL::setup_stiffness_matrix(P, Lambda, A_Lambda);

    cout << "setting up full right hand side..." << endl;
    Vector<double> F(Lambda.size());
    unsigned int id = 0;
    typename set<Index>::const_iterator it = Lambda.begin();
    for (; it != Lambda.end(); ++it, ++id) {
      F[id] = rhs.get_coefficient(*it);
    }

    cout << "size = " << F.size() << endl;
    unsigned int iterations = 0;
    CG(A_Lambda, F, u, 1.0e-15, 500, iterations);
    cout << "CG done!!!!" << " Needed " << iterations << " iterations" << endl;
  }



  template <class PROBLEM>
  void additive_Schwarz_SOLVE(PROBLEM& P,  const double epsilon,
			      InfiniteVector<double, typename PROBLEM::Index>& u_epsilon)
  {
    //typedef DSBasis<2,2> Basis1D;
    typedef PBasis<3,3> Basis1D;	

    Point<2> origin;
    origin[0] = 0.0;
    origin[1] = 0.0;
    
    CornerSingularity sing2D(origin, 0.5, 1.5);
    CornerSingularityRHS singRhs(origin, 0.5, 1.5);

    const int jmax = 7;
    typedef typename PROBLEM::Index Index;

    double a_inv     = P.norm_Ainv();
    double omega_i   = a_inv*P.F_norm();

    InfiniteVector<double, Index> u_k, u_k_1_2, r, help, f, Av, r_exact;
    
    map<double,double> log_10_residual_norms;
    map<double,double> degrees_of_freedom;
    map<double,double> asymptotic;
    map<double,double> time_asymptotic;
    map<double,double> log_10_L2_error;

    EvaluateFrame<Basis1D,2,2> evalObj;

    double eta = 0.5;
    //double eta = 0.00001;

    const int number_patches = P.basis().n_p();

    //const double alpha = 0.35;//pbasis 1D 3 3, 0.7x0.7
    const double alpha = 0.7;
    //double alpha = 0.5;
    //const double alpha = 0.3;
    //const double alpha = 0.2;

    unsigned int global_iterations = 0;
    double tmp = 5.;

    bool exit = 0;
    double time = 0;

    clock_t tstart, tend;
    tstart = clock();

    Array1D<InfiniteVector<double, Index> > old_xk(2);

    

    while (tmp > 1.0e-3) {

      cout << "reentering global loop " << endl;

      //approximate residual
      P.RHS(eta, f);
      APPLY_COARSE(P, u_k, eta, help, 0.00000001, jmax, CDD1);
      r = f - help;

      //################# OUTPUT ####################
      tend = clock();
      //approximate residual
      P.RHS(1.0e-8, f);
      //cout << f << endl;
      //abort();
      APPLY(P, u_k, 1.0e-8, help, jmax, CDD1);
      r_exact = f - help;
      
      //r_exact.COARSE(eta, r);

      tmp = l2_norm(r_exact);
      cout << "residual norm = " << tmp  << endl;
      double tmp1 = log10(tmp);

      //time = (double)((tend-tstart)/CLOCKS_PER_SEC);
      time += ((double) (tend-tstart))/((double) CLOCKS_PER_SEC);
      time_asymptotic[log10(time)] = tmp1;


      if (u_k.size() != 0)
	asymptotic[log10( (double)u_k.size() )] = tmp1;
      std::ofstream os3("add_schwarz_asymptotic_35_2D_al_02_0_9_2709.m");
      //std::ofstream os3("add_schwarz_asymptotic_P_35_2D_1108.m");
      matlab_output(asymptotic,os3);
      os3.close();

      std::ofstream os4("add_schwarz_time_asymptotic_35_2D_al_02_0_9_2709.m");
      //std::ofstream os4("add_schwarz_time_asymptotic_P_35_2D_1108.m");
      matlab_output(time_asymptotic,os4);
      os4.close();
      
      tstart = clock();
      //################# END OUTPUT #################


      // setup local index set
      set<Index> Lambda;
      r.support(Lambda);

#if 0
      //#### computing descent parameter ####
      help.clear();
      InfiniteVector<double, Index> Ahelp;

      APPLY_COARSE(P, r, eta/2., help, 0.00000001, jmax, CDD1);

      set<Index> Lambda2;
      help.support(Lambda2);
      
      // compute relaxation parameter
      for (int i = 0; i < number_patches; i++) {
	cout << "doing patch for relaxation parameter computation " << i << endl;
	set<Index> local_index_set;
	typename set<Index>::const_iterator it = Lambda2.begin();

	for (; it != Lambda2.end(); ++it) {
	  if ((*it).p() == i) {
	    local_index_set.insert(*it);
	  }
	}

	if (local_index_set.size() == 0)
	  continue;
	
	Vector<double> xk(local_index_set.size());
	GALERKIN(P, local_index_set, help, xk);

	int id = 0;
	for (typename set<Index>::const_iterator it = local_index_set.begin(), itend = local_index_set.end();
	     it != itend; ++it, ++id)
	  Ahelp.set_coefficient(*it, xk[id]);		
      }
      
      if (Ahelp.size() != 0 && r.size() != 0)
	alpha = ((r*r)/(r*Ahelp));
      cout << "descent parameter = " << alpha << endl;
      //#### end computing descent parameter ####
#endif

      for (int i = 0; i < number_patches; i++) {
	cout << "doing patch " << i << endl;
	set<Index> local_index_set;
	set<int> local_index_set2;
	typename set<Index>::const_iterator it = Lambda.begin();

	for (; it != Lambda.end(); ++it) {
	  if ((*it).p() == i) {
	    local_index_set.insert(*it);
	  }
	}

	if (local_index_set.size() == 0)
	  continue;
	
	Vector<double> xk(local_index_set.size());
	GALERKIN(P, local_index_set, r, xk);

	help.clear();
	int id = 0;
	for (typename set<Index>::const_iterator it = local_index_set.begin(), itend = local_index_set.end();
	     it != itend; ++it, ++id)
	  help.set_coefficient(*it, xk[id]);
	
	u_k = u_k + alpha*help;
		
      }// end for loop over patches
       
      global_iterations++;
      
      eta *= 0.8;
      
      u_epsilon = u_k;
      
      cout << global_iterations <<" loops completed" << endl;
      
      //       if (tmp < 0.0009) {
      // 	u_epsilon = u_k;
      // 	u_epsilon.scale(&P,-1);
      // 	char filename1[50];
      // 	char filename2[50];
	
      // 	sprintf(filename1, "%s%d%s%d%s", "approx_sol_add_schwarz33_2D_out_", global_iterations, "_nactive_", u_k.size(),".m");
      // 	sprintf(filename2, "%s%d%s%d%s", "error_add_schwarz_2D_out_", global_iterations, "_nactive_", u_k.size(),".m");
      // 	cout << "...plotting approximate solution" << endl;
      // 	Array1D<SampledMapping<2> > U = evalObj.evaluate(P.basis(), u_epsilon, true, 6);//expand in primal basis
	
      // 	std::ofstream ofs5(filename1);
      // 	matlab_output(ofs5,U);
      // 	ofs5.close();



      // 	Array1D<SampledMapping<2> > Error = evalObj.evaluate_difference(P.basis(), u_epsilon, sing2D, 6);
      //  	cout << "...plotting error" << endl;
      // 	std::ofstream ofs6(filename2);
      // 	matlab_output(ofs6, Error);
      // 	ofs6.close();
	
	
      //       }
    }

  }
}
