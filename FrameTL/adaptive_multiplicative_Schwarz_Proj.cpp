
#include <cmath>
#include <set>
#include <utils/plot_tools.h>
#include <adaptive/apply.h>
#include <numerics/corner_singularity.h>
#include <frame_evaluate.h>
#include <cdd1_local.h>
#include <error_H_scale.h>
#include <poisson_1d_testcase.h>
#include <projector_equation.h>
#include <lhs_equation.h>
#include <galerkin/cached_problem.h>

using std::set;

namespace FrameTL
{


  // forward declaration
  template <class IBASIS, int DIM>
  double
  H_1_error_interval(const AggregatedFrame<IBASIS,DIM,DIM>& frame,
		     const InfiniteVector<double, typename AggregatedFrame<IBASIS,DIM,DIM>::Index>& coeffs,
		     const Function<1>& f);
  
  // forward declaration
  template <class IBASIS, int DIM>
  double
  error_H_scale_interval (const int order,
			  const AggregatedFrame<IBASIS,DIM,DIM>& frame,
			  const InfiniteVector<double, typename AggregatedFrame<IBASIS,DIM,DIM>::Index>& coeffs,
			  const Function<1>& f);

  // #####################################################################################
  // Some dummy test cases. (Only for testing purposes.)
  // #####################################################################################

  /*!
   */
  template<class VALUE = double>
  class PolySolBiharmonic
    : public Function<1, VALUE>
  {
  public:
    PolySolBiharmonic() {};
    virtual ~PolySolBiharmonic() {};
    VALUE value(const Point<1>& p,
		const unsigned int component = 0) const
    {
      return 16*(p[0]*p[0])*(1-p[0])*(1-p[0]);
    }
  
    void vector_value(const Point<1> &p,
		      Vector<VALUE>& values) const { ; }
  
  };


  /*!
   */
  template<class VALUE = double>
  class SimpleTest
    : public Function<2, VALUE>
  {
  public:
    SimpleTest() {};
    virtual ~SimpleTest() {};
    VALUE value(const Point<2>& p,
		const unsigned int component = 0) const
    {
      return 10.0*p[0]*(1-p[0])*p[1]*(1-p[1]);
    }
  
    void vector_value(const Point<2> &p,
		      Vector<VALUE>& values) const {
      ; 
    }
  
  };

  /*!
   */
  template<class VALUE = double>
  class SimpleTestRHS
    : public Function<2, VALUE>
  {
  public:
    SimpleTestRHS() {};
    virtual ~SimpleTestRHS() {};
    VALUE value(const Point<2>& p,
		const unsigned int component = 0) const
    {
      return 10.0*2*(p[0]*(1-p[0]) + p[1]*(1-p[1]));
    }
  
    void vector_value(const Point<2> &p,
		      Vector<VALUE>& values) const {
      ; 
    }
  
  };

  /*!
   */
  template<class VALUE = double>
  class SimpleTestGradient
    : public Function<2, VALUE>
  {
  public:
    SimpleTestGradient() {};
    virtual ~SimpleTestGradient() {};
    VALUE value(const Point<2>& p,
		const unsigned int component = 0) const
    {
      double res = 0.;
      if (component == 0) {
	res = 10.0*(1-(2*p[0]))*(p[1]*(1-p[1]));
      }
      else if (component == 1) {
	res = 10.0*(1-(2*p[1]))*(p[0]*(1-p[0]));
      }
      
      return res;
      
    }
  
    void vector_value(const Point<2> &p,
		      Vector<VALUE>& values) const {
      values[0] = value(p,0);
      values[1] = value(p,1);
    }
  
  };
  // #####################################################################################



  // A helper routine for the adaptive algorithm. It works for two types
  // of domain coverings. The first one is the unit interval covered by two
  // subintervals [0,1-OVERLAP] \cup [OVERLAP,1]. A typical choice is OVERLAP=0.7.
  // The second is the L-shaped domain (-1,1)^2\setminus [0,1)^2 covered by
  // two patches [-OVERLAP,1]\times[-1,0] \cup [-1,0]\times[-1,1].
  // We take the sequence u and put into u1 those coefficients of u
  // that belong to patch 1-i and that corrrespond to wavelets the supports of which
  // nontrivially intersect the overlapping region.
  // Moreover, we put into u2 all those coefficients of u that correspond to patch 1-i.
  template <class PROBLEM>
  void split (PROBLEM& P, const int i,
	      const InfiniteVector<double, typename PROBLEM::Index>& u,
	      InfiniteVector<double, typename PROBLEM::Index>& u1,
	      InfiniteVector<double, typename PROBLEM::Index>& u2)
  {
    typedef typename PROBLEM::WaveletBasis::Support SuppType;
    u1.clear();
    u2.clear();
    typename InfiniteVector<double, typename PROBLEM::Index>::const_iterator it = u.begin();
#ifdef ONE_D
    if (i==0) {
      for (; it != u.end(); ++it) {
	const SuppType* supp = &(P.basis().all_patch_supports[it.index().number()]);
	if ((it.index().p() == 1) && (supp->a[0] < OVERLAP))
	  u1.set_coefficient(it.index(), *it);
	if (it.index().p() == 1)
	  u2.set_coefficient(it.index(), *it);
      }
    }
    else if (i==1) {
      for (; it != u.end(); ++it) {
	const SuppType* supp = &(P.basis().all_patch_supports[it.index().number()]);
	if ((it.index().p() == 0) && (supp->b[0] > 1-OVERLAP))
	  u1.set_coefficient(it.index(), *it);
	
	if (it.index().p() == 0)
	  u2.set_coefficient(it.index(), *it);
      }
    }
#endif
  }

  // A helper routine for the adaptive algorithm. It works for the case
  // of the rectangular ring-shaped domain as in Sect. 7.2.3 of Manuel's
  // PhD thesis, i.e. (-1,2)^2\setminus[0,1]^2, covered with 4 congruent
  // rectangles.
  // We put into u_sparse those coefficients of u that correspond to a patch
  // different from i and that correspond to wavelets not being fully
  // supported in patch i.
  // We put into u_very_sparse those coefficients of u that correspond to a patch
  // different from i and that correspond to wavelets which intersect with patch i
  // but which are not fully contained in it.
  template <class PROBLEM>
  void thin_out_ring (PROBLEM& P, const int i,
		      const InfiniteVector<double, typename PROBLEM::Index>& u,
		      InfiniteVector<double, typename PROBLEM::Index>& u_sparse,
		      InfiniteVector<double, typename PROBLEM::Index>& u_very_sparse)
  {

  }


  // A helper routine for the adaptive algorithm. It works for the case
  // of the L-shaped domain (-1,1)^2\setminus[0,1)^2, covered with the two rectangles
  // [-OVERLAP,1]\times[-1,0] \cup [-1,0]\times[-1,1].
  // We put into u_sparse those coefficients of u that correspond to the patch
  // 1-i and that correspond to wavelets not being fully supported in patch i.
  // We put into u_very_sparse those coefficients of u that correspond to the patch
  // 1-i and that correspond to wavelets which intersect with patch i
  // but which are not fully contained in it.
  template <class PROBLEM>
  void thin_out (PROBLEM& P, const int i,
		 const InfiniteVector<double, typename PROBLEM::Index>& u,
		 InfiniteVector<double, typename PROBLEM::Index>& u_sparse,
		 InfiniteVector<double, typename PROBLEM::Index>& u_very_sparse)
  {
 
    typedef typename PROBLEM::WaveletBasis::Support SuppType;
    u_sparse.clear();
    u_very_sparse.clear();
    typename InfiniteVector<double, typename PROBLEM::Index>::const_iterator it = u.begin();
    
#ifdef ONE_D
    if (i==0) {
      Point<1> x(OVERLAP);
      for (; it != u.end(); ++it) {
	if (it.index().p() == 1) {
	  if (in_support(P.basis(), it.index(),x))
	    u_very_sparse.set_coefficient(it.index(), *it);
	  
	  const SuppType* supp = &(P.basis().all_patch_supports[it.index().number()]);
	  if (supp->b[0] > OVERLAP)
	    u_sparse.set_coefficient(it.index(), *it);
	}
      }
    }
    else if (i==1) {
      Point<1> x(1-OVERLAP);
      for (; it != u.end(); ++it) {
	if (it.index().p() == 0) {
	  if (in_support(P.basis(), it.index(),x))
	    u_very_sparse.set_coefficient(it.index(), *it);
	  
	  const SuppType* supp = &(P.basis().all_patch_supports[it.index().number()]);
	  if (supp->a[0] < 1-OVERLAP)
	    u_sparse.set_coefficient(it.index(), *it);
	}
      }
    }
#endif

    
  }

  // Delete all coefficients of u corresponding to patch i.
  template <class PROBLEM>
  void remove_i (const int i, InfiniteVector<double, typename PROBLEM::Index>& u)
  {
    typename InfiniteVector<double, typename PROBLEM::Index>::const_iterator it = u.begin();
    for (; it != u.end(); ++it) {
      if (it.index().p() == i)
	u.set_coefficient(it.index(), 0.0);
    }
    u.compress();
  }

  
  template <class PROBLEM>
  void  MultSchw_Proj(const PROBLEM& P, const double epsilon,
		 Array1D<InfiniteVector<double, typename PROBLEM::Index> >& approximations)
  {
    // Exact solution and its first order derivative for the
    // one dimensional Poisson equation. This shall be used
    // for the computation of L_2- and H^1-errors.
#ifdef ONE_D
    Singularity1D_2<double> exact1D;
    Singularity1D_2_prime<double> exact1D_prime;
#endif

    // Exact solution and its gradient for the
    // two dimensional Poisson equation. This shall be used
    // for the computation of L_2- and H^1-errors.

    typedef typename PROBLEM::WaveletBasis::IntervalBasis Basis1D;
    const int jmax = JMAX;
    typedef typename PROBLEM::Index Index;
    typedef typename PROBLEM::WaveletBasis Frame;

    int d  = Basis1D::primal_polynomial_degree();
    int dT = Basis1D::primal_vanishing_moments();

#ifdef TWO_D
#ifdef RINGDOMAIN
    // Energy norm of the exact solution of the Poisson equation
    // in the ring-shaped domain.
    // 1.6544 is the estimate for the classical corner singularity.
    // Since we have four of them, we simply take 4.0*1.6544.
    const double mu = 4.0*1.6544;
#else
#ifdef BIHARMONIC
    // H^2 semi-norm of exact solution.
    // H_2_semi_norm_singularity is a macro defined in the test file
    // from where this routine is called. The estimate has been
    // computed i matlab using the .m-file estim_H2_seminorm_biharmL.m
    // in the directory diss/numerics/plotscripts. There, from a
    // plot of the solution, using a difference formula, the H^2-seminorm
    // has been estimated.
    const double mu = H_2_semi_norm_singularity;
#else
    // Energy norm of the exact solution of the Poisson equation
    // in the L-shaped domain.
    const double mu = 1.6544;
#endif
#endif
    const double M = 1.0; // we choose the most optimistic case

    // (d,dt) = (2,2) jmin = 3
    //const double M = sqrt(5.01773);
    // (d,dt) = (2,2) jmin = 4
    //const double M = sqrt(4.60975);

    // (d,dt) = (3,3) jmin = 3
    //const double M = sqrt(6.98681);

    // (d,dt) = (3,3) jmin = 4
    //const double M = sqrt(6.98986);

    // (d,dt) = (4,4) jmin = 4
    //const double M = sqrt(12.3335);

    // Next we have to set up the error reduction rate of the iterative solver.
    // Instead of using theoretical upper bounds, we work with manually chosen
    // values gained from some experiments with the full system matrix.
    //const double rho = sqrt(1 - 1.0/4.0); // [Xu92] Theorem 4.4 with \omega_1=1, K_0 = K_1=1
#ifdef RINGDOMAIN
    const double rho = 0.2;
#else
    const double rho = 0.2996;
#endif
    //const double rho = 0.5;
    //const double rho = 0.1;
    //const double rho = 0.1;
#endif
#ifdef ONE_D
    // Energy norm of the exact solution of the Poisson equation
    // from eq. (4.4.2) in Manuel's PhD thesis.
    const double mu = 7.44609; 
    const double M = 1.0; // we choose the most optimistic case

    // (d,dt) = (2,2)
    //const double M = sqrt(3.37459); // spectral radius of stiffness matrix

    // (d,dt) = (3,3), jmin = 4
    //const double M = sqrt(4.17833); // spectral radius of stiffness matrix
        
    // (d,dt) = (3,3), jmin = 3
    //const double M = sqrt(4.87718); // spectral radius of stiffness matrix

    // (d,dt) = (4,4)
    //const double M = sqrt(5.62803); // spectral radius of stiffness matrix

    // Next we have to set up the error reduction rate of the iterative solver.
    // Instead of using theoretical upper bounds, we work with manually chosen
    // values gained from some experiments with the full system matrix.
    //const double rho = sqrt(1 - 1.0/4.0); // [Xu92] Theorem 4.4 with \omega_1=1, K_0 = K_1=1
    const double rho = 0.1837;
#endif

    // #####################################################################################
    // Setup of constants.
    // #####################################################################################
    double rho_1 = pow(rho, 1./P.basis().n_p());
    const double C = 0.5 * rho_1*((1./rho)-1) / (1-rho_1);

    const double sigma = std::max(1./M, C + 1./M) + 0.1;//10.0
    cout << "sigma = " << sigma << endl; 
    const int K = std::max(1,(int)ceil(log(1.0/(2.0 * M * sigma)) / log(rho)));
    // #####################################################################################

    // #####################################################################################
    // We manually set the outer loop index L instead of taking the theoretical one.
    // #####################################################################################
#ifdef ONE_D
    int L;
    switch (d) {
    case 2: {
      L = 10;
      break;
    }
    case 3: {
      L = 8; // for j0=4 use 8, // for j0=3 use 10, for FULL alg. take 8 also
      break;
    }
    case 4: {
      L = 9;
      break;
    }
    };
#endif

    cout << "epsilon = " << epsilon << endl
	 << "mu = " << mu << endl
	 << "M = " << M << endl
	 << "rho = " <<  rho << endl
	 << "sigma = " << sigma << endl
	 << "K = " << K << endl
	 << "L = " << L << endl;
    // #####################################################################################
    // End of constant setup.
    // #####################################################################################


    double epsilon_proj;

    // number of patches
    const int m = P.basis().n_p();
    int k = 0;

    const int DIM = 1;
    
    // prepare indices for patches - needed for clip
    std:set<Index> patch[m];
    for (Index iterator = P.get_problem()->frame().first_generator(P.get_problem()->frame().j0() );; ++iterator)
    {
        patch[iterator.p()].insert(iterator);
        if (iterator == P.get_problem()->frame().last_wavelet(jmax)) break;
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

    // InfiniteVector's used in the adaptive algorithm
    InfiniteVector<double,Index> v, RHS, temp, u_scaled, null, start;
    InfiniteVector<double,Index> vi[m];


    // InfiniteVector's used in the adaptive algorithm
    InfiniteVector<double, Index> f, w, r, tmp, tmp_w;
    InfiniteVector<double, Index> u_k, u_k_sparse,u_k_very_sparse;
    InfiniteVector<double, Index> precond_r_i;

    Array1D<InfiniteVector<double, Index> > xks(P.basis().n_p()); // stores local approximations which are used as starting vectors
                                                                  // for the next cycle of local solves
    // map's used for generating output
    map<double,double> log_10_H1_error;
    map<double,double> log_10_H1_error_time;
    map<double,double> log_10_L2_error;
    map<double,double> log_10_L2_error_time;
    map<double,double> log_10_residual_error;
    map<double,double> log_10_residual_error_time;
    map<double,double> tolerances;
    map<double,double> rho_estimates;
    map<double,double> weak_ell_tau_norms;


    // variables for runtime measurement
    double time = 0.;
    clock_t tstart, tend;
    tstart = clock();
    double local_eps = 1.0;

    // #####################################################################################
    // The adaptive algorithm.
    // #####################################################################################
    for (int l = 1; l < L; l++) {
      //for (int p = 2; p <= K; p++) {// THAT WAS USED FOR THE 2D PLAIN DD CASE!!!
	for (int p = 1; p <= K; p++) {
	  for (int i = 0; i < P.basis().n_p(); i++) {
	  k = (l-1)*m*K+(p-1)*m+i+1;
	  cout << "################################" << endl;
	  cout << "number of iteration = " << k << endl;
	  cout << "################################" << endl;

	  // Setup tolerance for the solution of the local problems.
	  // This needs to be manually tuned not to end up with
	  // slow performance.
#ifdef RINGDOMAIN
 	  local_eps = mu*pow(2.0*pow(rho,K)*M*sigma,l-1)*pow(rho,p)/(m*K)*10;
#else
	  local_eps = mu*pow(2.0*pow(rho,K)*M*sigma,l-1)*pow(rho,p)/(m*K)*100;
#endif
 	  
	  // write the tolerances
	  tend = clock();
	  time += ((double) (tend-tstart))/((double) CLOCKS_PER_SEC);
	  tolerances[k] = local_eps;
	  std::ofstream os2d("tolerances.m");
	  matlab_output(tolerances,os2d);
	  os2d.close();
	  cout << "tolerance for solution of local problem = " << local_eps << endl;
	  tstart = clock();

	  precond_r_i.clear();

	  
	  // #####################################################################################
	  // Now follows the solution of the local problems.
	  // #####################################################################################
	  set<Index> Lambda_i;
	  
#if 1 // in this branch we perform the adaptive algorithm
#ifdef SPARSE

#ifdef RINGDOMAIN
	  // Preparation for the ring-shaped domain case.
	  thin_out_ring(P, i, u_k, u_k_sparse, u_k_very_sparse);
#else
	  // Preparation for the L-shaped domain case.
	  thin_out(P, i, u_k, u_k_sparse, u_k_very_sparse);
#endif
#endif
	  // CDD1
	  cout << "entering CDD solver..." << endl;
#ifdef SPARSE
	  // Solution of the local problem in case we use the sparse version of the algorithm
	  // as proposed in Stevenson, Werner 2009, where it is proposed to throw away
	  // all degrees of freedom that are contained in the current subdomain before the local solve.
	  CDD1_LOCAL_SOLVE(P, i, local_eps, xks[i], precond_r_i, u_k_very_sparse, jmax, CDD1);
	  //CDD1_LOCAL_SOLVE(P, i, local_eps, xks[i], precond_r_i, u_k_very_sparse, jmax, CDD1);
#endif
#ifdef FULL
	  // Solution of the local problem in case we use a plain multiplicative Schwarz adaptive 
	  // method, i.e., without removing degrees of freedom in the overlapping region before 
	  // the local solve as in the SPARSE-branch.
	  //remove_i<PROBLEM>(i, u_k);
	  CDD1_LOCAL_SOLVE(P, i, local_eps, xks[i], precond_r_i, u_k, jmax, CDD1); // THAT WAS USED FOR THE 2D CASE
	  //CDD1_LOCAL_SOLVE(P, i, 10.0*local_eps, xks[i], precond_r_i, u_k, jmax, CDD1); // THAT WAS USED FOR THE 1D CASE
#endif
	  cout << "CDD 1 solve completed, size of output is " << precond_r_i.size() << endl;

	  // Store the calculated local solution. These are always used as the initial guess in the
	  // next call of CDD1_LOCAL_SOLVE.
	  xks[i] = precond_r_i;
	  // ######################################################################################################
#else // In this branch we use a different strategy for the solution of the local problem: Approximate the right-hand side
      // for the local problem, apply COARSE to it, take the support of the right-hand side as Galerkin index set,
      // and approximate the Galerkin solution with the cg solver.

	  // approximate right hand side for local problem
	  P.RHS(local_eps, i, f);
	  cout << "fsize = " << f.size() << endl;
#ifdef FULL
	  APPLY(P, i, u_k, local_eps, w, jmax, CDD1);
#endif
#ifdef SPARSE
	  APPLY(P, i, u_k_very_sparse, local_eps, w, jmax, CDD1);
#endif
	  
	  r = f-w;
	  r.COARSE(local_eps, tmp_w);
	  r = tmp_w;
	  tmp_w.clear();

	  r.support(Lambda_i);
// 	  for (Index lambda = FrameTL::first_generator<Basis1D,2,2,Frame>(&P.basis(), P.basis().j0());
// 	       lambda <= FrameTL::last_wavelet<Basis1D,2,2,Frame>(&P.basis(), jmax); ++lambda) {
// 	    if (lambda.p() == i)
// 	      Lambda_i.insert(lambda);
// 	  }

	  cout << "setting up full stiffness matrix..." << endl;
	  cout << "size of local index set = " << Lambda_i.size() << endl;
	  if (Lambda_i.size() > 0) {
	    SparseMatrix<double> A_Lambda;
	    WaveletTL::setup_stiffness_matrix(P, Lambda_i, A_Lambda);
	    
	    cout << "setting up full right hand side..." << endl;
	    Vector<double> F(Lambda_i.size()), xk(Lambda_i.size());
	    unsigned int id = 0;
	    typename set<Index>::const_iterator it = Lambda_i.begin();
	    for (; it != Lambda_i.end(); ++it, ++id) {
	      F[id] = r.get_coefficient(*it);
	      xk[id] = xks[i].get_coefficient(*it);
	    }
	    unsigned int iterations = 0;
	    CG(A_Lambda, F, xk, 1.0e-15, 500, iterations);
	    cout << "CG done!!!!" << " Needed " << iterations << " iterations" << endl;
	    id = 0;
	    for (typename set<Index>::const_iterator it = Lambda_i.begin(), itend = Lambda_i.end();
		 it != itend; ++it, ++id)
	    precond_r_i.set_coefficient(*it, xk[id]);
	  }
#endif 
	  // ######################################################################################################  

	  // #####################################################################################
	  // setup next global iterate
	  // #####################################################################################
#ifdef SPARSE
	  u_k = precond_r_i + u_k_sparse;
#endif
#ifdef FULL
	  u_k = precond_r_i + u_k;
#endif
	  cout << "degrees of freedom: " << u_k.size() << endl;
	  // #####################################################################################


 	  xks[i].clear();
 	  xks[i] = precond_r_i;
      	}// end loop over patches
	cout << "############## full cycle of local solves completed ##############"<< endl;



	}// end loop p

      	          // projection_step
          if (/*l <= L/2 && i == m-1 && p == 1*/ false && l < L-1){
              cout << "START PROJECTION STEP" << endl;
              double solve_tol = 0.5*(sigma - 1./M)*2.0*pow(rho,K)*mu*pow(2.0*pow(rho,K)*M*sigma,l-1);
              for(int i=0; i<m; i++)
              {



                  // scale result from solver
                  u_scaled.clear();
                  u_scaled = u_k;
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
                  double tol_apply = 10e-6;
    #endif
                  // prepare RHS for patch 0 and 1
                  if (i == 0)
                    APPLY(rhs_problem_0, u_scaled, tol_apply, RHS, jmax, CDD1);
                  else if(i == 1)
                    APPLY(rhs_problem_1, u_scaled, tol_apply, RHS, jmax, CDD1);

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

                  vi[i] = u_k;
                  vi[i].scale(&P, -1);
                  if (i==0)
                    lhs_equation_0.rescale(vi[i],1);//vi[i].scale(&lhs_problem_0, 1);
                  else if(i==1)
                    lhs_equation_1.rescale(vi[i],1);//vi[i].scale(&lhs_problem_1, 1);


                  RHS.clip(patch[i]);

                  cout << "RHS on patch " << i << " done." << endl;

                  // due to zeros on the main diagonal in the matrix needed to set up the RHS, we need to remove nan-entries here
                  for(typename InfiniteVector<double, Index>::const_iterator it = RHS.begin(); it != RHS.end(); ++it)
                  {
                      if(RHS[it.index()] != RHS[it.index()]) {cout << "nan!" << endl; RHS[it.index()] = 0;}
                  }


                  // solve local problems
                  if (i==0){
                    lhs_equation_0.set_rhs(RHS);
                    CDD1_LOCAL_SOLVE(lhs_problem_0, 0, solve_tol, null , vi[i], null, jmax, CDD1);
                    lhs_equation_0.rescale(vi[i],-1);//vi[i].scale(&lhs_problem_0,-1);
                  }
                  else if (i==1){
                    lhs_equation_1.set_rhs(RHS);
                    CDD1_LOCAL_SOLVE(lhs_problem_1, 1, solve_tol, null , vi[i], null, jmax, CDD1);
                    lhs_equation_1.rescale(vi[i],-1);//vi[i].scale(&lhs_problem_1,-1);
                  }
                  cout << "Projection on patch " << i << " done." << endl;

              }

              // collect local parts
              u_k.clear();
              for(int i=0; i<m; i++) u_k += vi[i];

              // scale back
              u_k.scale(&P,1);

              // extra COARSE step after scaling to reduce degrees of freedom
              u_k.COARSE(solve_tol,temp);
              temp = v;

              cout << "PROJECTION STEP DONE" << endl;
          }
	// setup tolerance for coarsening
#ifdef RINGDOMAIN
      double coarse_tol = (sigma - 1./M)*2.0*pow(rho,K)*mu*pow(2.0*pow(rho,K)*M*sigma,l-1)*0.1;
#else
      double coarse_tol = (sigma - 1./M)*2.0*pow(rho,K)*mu*pow(2.0*pow(rho,K)*M*sigma,l-1);
#endif
      cout << "tolerance for coarsening = " << coarse_tol << endl;
      cout << "norm of u_k = " << l2_norm(u_k) << endl;
      u_k.COARSE(coarse_tol, tmp_w);
      u_k = tmp_w;
      cout << "degrees of freedom after coarsening: " << u_k.size() << endl;

      // #####################################################################################
      //  Approximate global EXACT residual and perform output.
      // #####################################################################################
      tend = clock();
      time += ((double) (tend-tstart))/((double) CLOCKS_PER_SEC);
      
      tmp = u_k;
      tmp.scale(&P,-1);

      //compute global residual
      P.RHS(1.0e-8, f);
      cout << "fsize exact res = " << f.size() << endl;
      tmp_w.clear(); 
      for (int i = 0; i < P.basis().n_p(); i++) {
      	APPLY(P, i, u_k, 1.0e-8, w, jmax, CDD1);
      	tmp_w += w;
      }
      double residual_norm = l2_norm(f-tmp_w);
      cout << "norm of global residual = " << residual_norm  << endl;

      char name1[128];
      char name2[128];
      char name3[128];
      char name4[128];
      char name5[128];
      char name6[128];
      char name7[128];

      	// setup filenames for output files for the one-dimensional cases
#ifdef ONE_D
      switch (d) {
      case 2: {
	weak_ell_tau_norms[k] = u_k.weak_norm(1./1.5);// (d,dT)=(2,2)=1./1.5 (d,dT)=(3,3)=1./2.5, (d,dT)=(4,6)=1./3.5
	sprintf(name1, "%s%d%s%d%s", "./ms_results22/ms1D_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name2, "%s%d%s%d%s", "./ms_results22/ms1D_time_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name3, "%s%d%s%d%s", "./ms_results22/ms1D_weak_ell_tau_norms_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name4, "%s%d%s%d%s", "./ms_results22/ms1D_H1err_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name5, "%s%d%s%d%s", "./ms_results22/ms1D_H1err_time_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name6, "%s%d%s%d%s", "./ms_results22/ms1D_L2err_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name7, "%s%d%s%d%s", "./ms_results22/ms1D_L2err_time_P_jmax18_", d, "_dT", dT, ".m");
	break;
      }
      case 3: {
	weak_ell_tau_norms[k] = u_k.weak_norm(1./2.5);
	sprintf(name1, "%s%d%s%d%s", "./ms_results33/ms1D_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name2, "%s%d%s%d%s", "./ms_results33/ms1D_time_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name3, "%s%d%s%d%s", "./ms_results33/ms1D_weak_ell_tau_norms_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name4, "%s%d%s%d%s", "./ms_results33/ms1D_H1err_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name5, "%s%d%s%d%s", "./ms_results33/ms1D_H1err_time_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name6, "%s%d%s%d%s", "./ms_results33/ms1D_L2err_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name7, "%s%d%s%d%s", "./ms_results33/ms1D_L2err_time_P_jmax18_", d, "_dT", dT, ".m");
	break;
      }
      case 4: {
	weak_ell_tau_norms[k] = u_k.weak_norm(1./3.5);
	sprintf(name1, "%s%d%s%d%s", "./ms_results46/ms1D_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name2, "%s%d%s%d%s", "./ms_results46/ms1D_time_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name3, "%s%d%s%d%s", "./ms_results46/ms1D_weak_ell_tau_norms_P_jmax18_", d, "_dT", dT, ".m");  
	sprintf(name4, "%s%d%s%d%s", "./ms_results46/ms1D_H1err_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name5, "%s%d%s%d%s", "./ms_results46/ms1D_H1err_time_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name6, "%s%d%s%d%s", "./ms_results46/ms1D_L2err_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name7, "%s%d%s%d%s", "./ms_results46/ms1D_L2err_time_P_jmax18_", d, "_dT", dT, ".m");
	break;
      }
      };

      // perform matlab output
      log_10_residual_error[log10(u_k.size())] = log10(l2_norm(f-tmp_w));
      std::ofstream os2c(name1);
      matlab_output(log_10_residual_error,os2c);
      os2c.close();
      log_10_residual_error_time[log10(time)] = log10(l2_norm(f-tmp_w));
      std::ofstream os2e(name2);
      matlab_output(log_10_residual_error_time,os2e);
      os2e.close();

      std::ofstream os2g(name3);
      matlab_output(weak_ell_tau_norms,os2g);
      os2g.close();


#endif
	
	// setup filenames for output files for the two-dimensional cases
#ifdef TWO_D
      switch (d) {
      case 2: {
	weak_ell_tau_norms[k] = u_k.weak_norm(1./1.0);
	sprintf(name1, "%s%d%s%d%s", "./ms_results2D_22/ms2D_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name2, "%s%d%s%d%s", "./ms_results2D_22/ms2D_time_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name3, "%s%d%s%d%s", "./ms_results2D_22/ms2D_weak_ell_tau_norms_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name4, "%s%d%s%d%s", "./ms_results2D_22/ms2D_H1err_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name5, "%s%d%s%d%s", "./ms_results2D_22/ms2D_H1err_time_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name6, "%s%d%s%d%s", "./ms_results2D_22/ms2D_L2err_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name7, "%s%d%s%d%s", "./ms_results2D_22/ms2D_L2err_time_P_jmax18_", d, "_dT", dT, ".m");
	break;
      }
      case 3: {
#ifdef BIHARMONIC
	weak_ell_tau_norms[k] = u_k.weak_norm(1./1.0);
#else
	weak_ell_tau_norms[k] = u_k.weak_norm(1./1.5);
#endif
	sprintf(name1, "%s%d%s%d%s", "./ms_results2D_33_biharmL_v2/ms2D_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name2, "%s%d%s%d%s", "./ms_results2D_33_biharmL_v2/ms2D_time_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name3, "%s%d%s%d%s", "./ms_results2D_33_biharmL_v2/ms2D_weak_ell_tau_norms_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name4, "%s%d%s%d%s", "./ms_results2D_33_biharmL_v2/ms2D_H1err_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name5, "%s%d%s%d%s", "./ms_results2D_33_biharmL_v2/ms2D_H1err_time_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name6, "%s%d%s%d%s", "./ms_results2D_33_biharmL_v2/ms2D_L2err_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name7, "%s%d%s%d%s", "./ms_results2D_33_biharmL_v2/ms2D_L2err_time_P_jmax18_", d, "_dT", dT, ".m");
	break;
      }
      case 4: {
#ifdef BIHARMONIC
	weak_ell_tau_norms[k] = u_k.weak_norm(1./1.5);
#else
	weak_ell_tau_norms[k] = u_k.weak_norm(1./2.0);
#endif
	sprintf(name1, "%s%d%s%d%s", "./ms_results2D_46_biharmL/ms2D_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name2, "%s%d%s%d%s", "./ms_results2D_46_biharmL/ms2D_time_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name3, "%s%d%s%d%s", "./ms_results2D_46_biharmL/ms2D_weak_ell_tau_norms_P_jmax18_", d, "_dT", dT, ".m");  
	sprintf(name4, "%s%d%s%d%s", "./ms_results2D_46_biharmL/ms2D_H1err_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name5, "%s%d%s%d%s", "./ms_results2D_46_biharmL/ms2D_H1err_time_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name6, "%s%d%s%d%s", "./ms_results2D_46_biharmL/ms2D_L2err_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name7, "%s%d%s%d%s", "./ms_results2D_46_biharmL/ms2D_L2err_time_P_jmax18_", d, "_dT", dT, ".m");
	break;
      }
      };

      // perform matlab output
      log_10_residual_error[log10(u_k.size())] = log10(l2_norm(f-tmp_w));
      std::ofstream os2c(name1);
      matlab_output(log_10_residual_error,os2c);
      os2c.close();

      log_10_residual_error_time[log10(time)] = log10(l2_norm(f-tmp_w));
      std::ofstream os2e(name2);
      matlab_output(log_10_residual_error_time,os2e);
      os2e.close();

      std::ofstream os2g(name3);
      matlab_output(weak_ell_tau_norms,os2g);
      os2g.close();
    

#endif
      // #####################################################################################
      //  End performing output
      // #####################################################################################

      tstart = clock();


    }// end loop L


    // #####################################################################################
    // The adaptive algorithm is finished here.
    // #####################################################################################    
    
    // collect final approximation and its local parts
    approximations[P.basis().n_p()] = u_k;
    
    for (int i = 0; i < P.basis().n_p(); i++) {
      approximations[i].clear();
      for (typename InfiniteVector<double, Index>::const_iterator it = u_k.begin(), itend = u_k.end();
	   it != itend; ++it)
	if (it.index().p() == i)
	  approximations[i].set_coefficient(it.index(),*it);
    }




  }


   template <class PROBLEM>
   void  adaptive_multiplicative_Schwarz_SOLVE(const PROBLEM& P, const double epsilon,
 					      Array1D<InfiniteVector<double, typename PROBLEM::Index> >& approximations)
   {
    }
}
