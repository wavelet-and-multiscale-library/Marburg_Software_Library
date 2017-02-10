// implementation for steepest_descent.h

#include <cmath>
#include <set>
#include <utils/plot_tools.h>
#include <adaptive/apply.h>
#include <numerics/corner_singularity.h>
#include <frame_evaluate.h>
#include <cdd1_local.h>
#include <error_H_scale.h>
#include <poisson_1d_testcase.h>

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
#ifdef TWO_D
    if (i==0) {
      for (; it != u.end(); ++it) {
	const SuppType* supp = &(P.basis().all_patch_supports[it.index().number()]);
	if ((it.index().p() == 1) && (supp->b[0] > -OVERLAP) && (supp->a[1] < 0.)) {
	  u1.set_coefficient(it.index(), *it);
	}
	if (it.index().p() == 1)
	  u2.set_coefficient(it.index(), *it);
      }
    }
    else if (i==1) {
      for (; it != u.end(); ++it) {
	const SuppType* supp = &(P.basis().all_patch_supports[it.index().number()]);
	if ((it.index().p() == 0) && (supp->a[0] < 0.)) {
	  u1.set_coefficient(it.index(), *it);
	}
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
#ifdef TWO_D 
    typedef typename PROBLEM::WaveletBasis::Support SuppType;
    u_sparse.clear();
    u_very_sparse.clear();
    typename InfiniteVector<double, typename PROBLEM::Index>::const_iterator it = u.begin();

    switch (i) {
    case 0: {
      for (; it != u.end(); ++it) {
	const SuppType* supp = &(P.basis().all_patch_supports[it.index().number()]);
	if ((it.index().p() != i) && (supp->a[1] < 0.) && (0. < supp->b[1]) ) {
	  u_very_sparse.set_coefficient(it.index(), *it);
	}
	if ((it.index().p() != i) && (supp->b[1] > 0.) ) {
	  u_sparse.set_coefficient(it.index(), *it);
	}
      }
      break;
    }
    case 1: {
      for (; it != u.end(); ++it) {
	const SuppType* supp = &(P.basis().all_patch_supports[it.index().number()]);
	if ((it.index().p() != i) && (supp->a[0] < 1.0) && (1.0 < supp->b[0]) ) {
	  u_very_sparse.set_coefficient(it.index(), *it);
	}
	if ((it.index().p() != i) && (supp->a[0] < 1.0) ) {
	  u_sparse.set_coefficient(it.index(), *it);
	}
      }
      break;
    }
    case 2: {
      for (; it != u.end(); ++it) {
	const SuppType* supp = &(P.basis().all_patch_supports[it.index().number()]);
	if ((it.index().p() != i) && (supp->a[1] < 1.0) && (1.0 < supp->b[1]) ) {
	  u_very_sparse.set_coefficient(it.index(), *it);
	}
	if ((it.index().p() != i) && (supp->a[1] < 1.0) ) {
	  u_sparse.set_coefficient(it.index(), *it);
	}
      }
      break;
    }
    case 3: {
      for (; it != u.end(); ++it) {
	const SuppType* supp = &(P.basis().all_patch_supports[it.index().number()]);
	if ((it.index().p() != i) && (supp->a[0] < 0.0) && (0.0 < supp->b[0]) ) {
	  u_very_sparse.set_coefficient(it.index(), *it);
	}
	if ((it.index().p() != i) && (supp->b[0] > 0.0) ) {
	  u_sparse.set_coefficient(it.index(), *it);
	}
      }
      break;
    }
    }
#endif
  

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
#ifdef TWO_D
    if (i==0) {
      for (; it != u.end(); ++it) {
	const SuppType* supp = &(P.basis().all_patch_supports[it.index().number()]);
	
	// check whether first line is intersected
	if ((it.index().p() == 1) && (supp->b[0] > -OVERLAP) && (supp->a[1] < 0.) && (0. < supp->b[1]) ) {
	  u_very_sparse.set_coefficient(it.index(), *it);
	}
	// check whether second line is intersected
	if ((it.index().p() == 1) && (supp->a[1] < 0.) && (supp->a[0] < -OVERLAP) && (-OVERLAP < supp->b[0])) {
	  u_very_sparse.set_coefficient(it.index(), *it);
	}
	if ((it.index().p() == 1) && ((supp->a[0] < -OVERLAP) || (supp->b[1] > 0.)) ) {
	  u_sparse.set_coefficient(it.index(), *it);
	}
// 	if ((it.index().p() == 1) && !((supp->a[0] > -OVERLAP) && (supp->b[1] < 0.)) ) {
// 	  u_sparse.set_coefficient(it.index(), *it);
//	}
      }
    }
    else if (i==1) {
      for (; it != u.end(); ++it) {
	const SuppType* supp = &(P.basis().all_patch_supports[it.index().number()]);
	if ((it.index().p() == 0) && (supp->a[0] < 0.) && (0. < supp->b[0])) {
	  u_very_sparse.set_coefficient(it.index(), *it);
	}
	if ( (it.index().p() == 0) && ( supp->b[0] > 0. ) ) {
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
  void  MultSchw(const PROBLEM& P, const double epsilon,
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
#ifdef TWO_D
    SimpleTestGradient<double> simple_sol_grad;
    SimpleTest<double> simple_sol;

    Point<2> origin;
    origin[0] = 0.0;
    origin[1] = 0.0;
    
    CornerSingularity exact2D(origin, 0.5, 1.5);
    CornerSingularityGradient singGrad(origin, 0.5, 1.5);
#endif

    typedef typename PROBLEM::WaveletBasis::IntervalBasis Basis1D;
    const int jmax = JMAX;
    typedef typename PROBLEM::Index Index;
    //typedef typename PROBLEM::WaveletBasis Frame;

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
#ifdef TWO_D
    int L;
    switch (d) {
    case 2: {
      L = 8; // formerly 7
      break;
    }
    case 3: {
#ifdef RINGDOMAIN
      L = 7;
#else
      L = 11;
#endif
      break;
    }
    case 4: {
      L = 11;
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

    // number of patches
    const int m = P.basis().n_p();
    int k = 0;

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

//       double H1err = error_H_scale_interval<Basis1D>(1,P.basis(), tmp, exact1D_prime);
//       log_10_H1_error[log10(u_k.size())] = log10(H1err);
//       std::ofstream os2h(name4);
//       matlab_output(log_10_H1_error,os2h);
//       os2h.close();

//       log_10_H1_error_time[log10(time)] = log10(H1err);
//       std::ofstream os2i(name5);
//       matlab_output(log_10_H1_error_time,os2i);
//       os2i.close();

//       double L2err = error_H_scale_interval<Basis1D>(0, P.basis(), tmp, exact1D);
//       log_10_L2_error[log10(u_k.size())] = log10(L2err);
//       std::ofstream os2j(name6);
//       matlab_output(log_10_L2_error,os2j);
//       os2j.close();

//       log_10_L2_error_time[log10(time)] = log10(L2err);
//       std::ofstream os2k(name7);
//       matlab_output(log_10_L2_error_time,os2k);
//       os2k.close();


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
      
//       double H1err = error_H_scale_Lshaped<Basis1D>(1, P.basis(), tmp, singGrad);
//       log_10_H1_error[log10(u_k.size())] = log10(H1err);
//       std::ofstream os2h(name4);
//       matlab_output(log_10_H1_error,os2h);
//       os2h.close();

//       log_10_H1_error_time[log10(time)] = log10(H1err);
//       std::ofstream os2i(name5);
//       matlab_output(log_10_H1_error_time,os2i);
//       os2i.close();

//       double L2err = error_H_scale_Lshaped<Basis1D>(0, P.basis(), tmp, exact2D);
//       log_10_L2_error[log10(u_k.size())] = log10(L2err);
//       std::ofstream os2j(name6);
//       matlab_output(log_10_L2_error,os2j);
//       os2j.close();

//       log_10_L2_error_time[log10(time)] = log10(L2err);
//       std::ofstream os2k(name7);
//       matlab_output(log_10_L2_error_time,os2k);
//       os2k.close();


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

//   template <class PROBLEM>
//   double compute_exact_residual_norm (const PROBLEM& P,
// 				      InfiniteVector<double, typename PROBLEM::Index>& u_k) {
//     typedef typename PROBLEM::Index Index;
//     typedef typename PROBLEM::WaveletBasis::IntervalBasis Basis1D;
//     typedef typename PROBLEM::WaveletBasis Frame;
//     const int jmax = JMAX;
//     set<Index> Lambda1;
//     for (Index lambda = FrameTL::first_generator<Basis1D,1,1,Frame>(&P.basis(), P.basis().j0());
// 	 lambda <= FrameTL::last_wavelet<Basis1D,1,1,Frame>(&P.basis(), jmax); ++lambda) {
//       Lambda1.insert(lambda);
//     }
//     set<Index> Lambda2;
//     u_k.support(Lambda2);
//     SparseMatrix<double> A(Lambda1.size(), Lambda2.size());
//     WaveletTL::setup_stiffness_matrix(P, Lambda1, Lambda2, A);

//     // copy u_k in vector
//     Vector<double> U_k(u_k.size());
//     Vector<double> A_U_k(Lambda1.size());
//     unsigned int id = 0;
//     typename set<Index>::const_iterator it = Lambda2.begin();
//     for (; it != Lambda2.end(); ++it, ++id) {
//       U_k[id] = u_k.get_coefficient(*it);
//     }

//     Vector<double> F(Lambda1.size());    
//     WaveletTL::setup_righthand_side(P, Lambda1, F);

//     A.apply(U_k, A_U_k);
//     return l2_norm(F-A_U_k);

//   }

   template <class PROBLEM>
   void  adaptive_multiplicative_Schwarz_SOLVE(const PROBLEM& P, const double epsilon,
 					      Array1D<InfiniteVector<double, typename PROBLEM::Index> >& approximations)
   {
// #ifdef ONE_D
//     Singularity1D_2<double> exact1D;
//     Singularity1D_2_prime<double> exact1D_prime;
// #endif

// #ifdef TWO_D
//     SimpleTestGradient<double> simple_sol_grad;
//     SimpleTest<double> simple_sol;

//     Point<2> origin;
//     origin[0] = 0.0;
//     origin[1] = 0.0;
    
//     CornerSingularity exact2D(origin, 0.5, 1.5);
//     CornerSingularityGradient singGrad(origin, 0.5, 1.5);
// #endif

//     //typedef DSBasis<2,2> Basis1D;
//     typedef PBasis<3,3> Basis1D;

//     const int jmax = JMAX;
//     typedef typename PROBLEM::Index Index;
//     typedef typename PROBLEM::WaveletBasis Frame;


//     double a_inv     = P.norm_Ainv();
//     double omega_i   = a_inv*P.F_norm();

//     InfiniteVector<double, Index> f, w, r, tmp_w, tmp;
//     InfiniteVector<double, Index> u_k, u_k_sparse,u_k_very_sparse, u1, u2;
//     InfiniteVector<double, Index> precond_r_i, r_i;
//     //Array1D<InfiniteVector<double, Index> > u_k_i(P.basis().n_p());


//     Array1D<InfiniteVector<double, Index> > global_residual_parts(P.basis().n_p());
      
//     map<double,double> log_10_residual_norms;
//     map<double,double> weak_ell_tau_norms;
//     map<double,double> degrees_of_freedom;
//     map<double,double> asymptotic;
//     map<double,double> time_asymptotic;
//     map<double,double> log_10_H1_error;
//     map<double,double> log_10_L2_error;

//     EvaluateFrame<Basis1D,2,2> evalObj;

//     //double eta = 0.5;
//     //double eta = 0.000001;

//     const int number_patches = P.basis().n_p();

//     const double alpha = 0.1;

//     unsigned int global_iterations = 1;
//     double eps = 1.;

//     Array1D<InfiniteVector<double, Index> > xks(P.basis().n_p());
// //     while (global_iterations < 250) {
// //       P.RHS(eps, f);
// //       APPLY(P, u_k, eps, w, jmax, CDD1);
// //       r = f-w;
// //       cout << l2_norm(r) << endl;
// //       u_k = u_k + alpha * r;
// //       u_k.COARSE(eps/2,w);
// //       u_k = w;
// //       eps *= 0.95;
// //       cout << "eps = " << eps << endl;
// //       cout << "degrees of freedom = " << u_k.size() << endl;
// //       cout << "loop = " << global_iterations << endl;
// //       global_iterations++;
// //     }    

//     double time = 0.;
//     clock_t tstart, tend;
//     tstart = clock();

//     while (global_iterations < 20) {
//       cout << eps << endl;
//       // loop over patches
//       for (int i = 0; i < P.basis().n_p(); i++) {
// 	//P.RHS(1.0e-5, i, f);
// 	P.RHS(eps, i, f);
// #ifdef SPARSE
//  	thin_out(P, i, u_k, u_k_sparse, u_k_very_sparse);
// 	//APPLY(P, i, u_k_very_sparse, 1.0e-10, w, jmax, CDD1);
// 	APPLY(P, i, u_k_very_sparse, eps, w, jmax, CDD1);
// 	cout << "wsize = " << w.size() << endl;
// #endif

// #ifdef FULL
//   	split(P, i, u_k, u1, u2);
//   	APPLY(P, i, u1, eps, w, jmax, CDD1);
// // 	APPLY(P, i, u_k, eps, w, jmax, CDD1);
// //  	w.COARSE(eps/2,tmp_w);
// //  	w = tmp_w;
	
// #endif

//  	r = f - w;
//     	r.COARSE(2*eps,tmp_w);
// 	r = tmp_w;

// 	cout << "fsize = " << f.size() << endl;
// 	//cout << f << endl;
// 	// local solve
// 	precond_r_i.clear();
// #if 0
// 	// CDD1
// 	CDD1_LOCAL_SOLVE(P, i, 100*eps, xks[i], precond_r_i, r, jmax, CDD1);
// #endif	

	

// #if 0
// 	// adaptive Richardson
// 	for (int k = 0; k < 200; k++) {
// 	  APPLY(P, i, precond_r_i, 1.0e-8, w, jmax, CDD1);
// 	  r_i = r - w;
// 	  //cout << r_i << endl;
// 	  precond_r_i = precond_r_i + alpha * r_i;
// 	  //cout << precond_r_i << endl;
// 	  cout << "norm of local residal = " << l2_norm(r_i) << endl;
// 	}
// 	break;
// #endif
// #if 1
// 	// single step of CDD1
// 	set<Index> Lambda_i;
// 	r.support(Lambda_i);

// // 	for (Index lambda = FrameTL::first_generator<Basis1D,2,2,Frame>(&P.basis(), P.basis().j0());
// // 	     lambda <= FrameTL::last_wavelet<Basis1D,2,2,Frame>(&P.basis(), jmax); ++lambda) {
// // 	  if (lambda.p() == i)
// // 	    Lambda_i.insert(lambda);
// // 	  cout << lambda << endl;
// // 	}
// 	cout << "setting up full stiffness matrix..." << endl;
// 	cout << "size of local index set = " << Lambda_i.size() << endl;
// 	if (Lambda_i.size() > 0) {
// 	  SparseMatrix<double> A_Lambda;
// 	  WaveletTL::setup_stiffness_matrix(P, Lambda_i, A_Lambda);
	  
// 	  cout << "setting up full right hand side..." << endl;
// 	  Vector<double> F(Lambda_i.size()), xk(Lambda_i.size());
// 	  unsigned int id = 0;
// 	  typename set<Index>::const_iterator it = Lambda_i.begin();
// 	  for (; it != Lambda_i.end(); ++it, ++id) {
// 	    F[id] = r.get_coefficient(*it);
// 	    xk[id] = xks[i].get_coefficient(*it);
// 	  }
	  
// 	  //cout << "size = " << F.size() << endl;
// 	  unsigned int iterations = 0;
// 	  CG(A_Lambda, F, xk, 1.0e-15, 500, iterations);
// 	  cout << "CG done!!!!" << " Needed " << iterations << " iterations" << endl;
// 	  id = 0;
// 	  for (typename set<Index>::const_iterator it = Lambda_i.begin(), itend = Lambda_i.end();
// 	       it != itend; ++it, ++id)
// 	    precond_r_i.set_coefficient(*it, xk[id]);
// 	}
// #endif
// 	// put together next iterate
// #ifdef SPARSE
// 	//precond_r_i.COARSE(eps,tmp_w);
// 	//precond_r_i = tmp_w;
// 	u_k = precond_r_i + u_k_sparse;
// 	//u_k = precond_r_i + u_k;
// 	cout << "degrees of freedom: " << u_k.size() << endl;
// 	xks[i].clear();
// 	xks[i] = precond_r_i;
//       	u_k.COARSE(2*eps,tmp_w);
//        	u_k = tmp_w;
// 	cout << "degrees of freedom after coarsening: " << u_k.size() << endl;
// #endif
// #ifdef FULL
// 	//u_k = precond_r_i + u_k;
// 	u_k = precond_r_i + u2;
//    	u_k.COARSE(2*eps,tmp_w);
//    	u_k = tmp_w;

// 	cout << "degrees of freedom: " << u_k.size() << endl;

// #endif
// 	cout << "norm of local residal = " << l2_norm(r) << endl;
// 	cout << "weak ell tau norm = " << u_k.weak_norm(1./1.5) << endl;
	
// 	tend = clock();
// 	time += ((double) (tend-tstart))/((double) CLOCKS_PER_SEC);
      
// 	if (global_iterations % P.basis().n_p() == 0) {
	  
// 	  weak_ell_tau_norms[global_iterations] = u_k.weak_norm(1./1.5);
// 	  std::ofstream os3("weak_ell_tau_norms.m");
// 	  matlab_output(weak_ell_tau_norms,os3);
// 	  os3.close();

// // 	  asymptotic[log10(u_k.size())] = log10(l2_norm(r));
// // 	  std::ofstream os4("asymptotic_test.m");
// // 	  matlab_output(asymptotic,os4);
// // 	  os4.close();
// 	}
// 	tstart = clock();
//       }

      
//       tend = clock();
//       time += ((double) (tend-tstart))/((double) CLOCKS_PER_SEC);

//       cout << "############### loop " << global_iterations << " ###############" << endl;
//       global_iterations++;
//        if (global_iterations % 1 == 0 || (global_iterations == 1)) {
//          tmp = u_k;
//          tmp.scale(&P,-1);
// 	 //         double H1err = error_H_scale_Lshaped<Basis1D>(1, P.basis(), tmp, singGrad);
// 	 double H1err = error_H_scale_interval<Basis1D>(1,P.basis(), tmp, exact1D_prime);
// 	 cout << "H_1 error = " <<  H1err << endl;
// // 	//       //double L_2_error = error_H_scale_interval<Basis1D,2>(0, P.basis(), tmp, exact1D);
// // 	// //       double L_2_error = error_H_scale_Lshaped<Basis1D>(0, P.basis(), tmp, exact2D);
// // 	// //       cout << "L_2 error = " << L_2_error << endl;
// // 	// //       log_10_L2_error[log10(u_k.size())] = log10(L_2_error);
	
//          log_10_H1_error[log10(u_k.size())] = log10(H1err);
//          std::ofstream os2a("log_10_H1_error1D.m");
//          matlab_output(log_10_H1_error,os2a);
//          os2a.close();
// // 	//        std::ofstream os2("log_10_L_2_error2D.m");
// // 	//        matlab_output(log_10_L2_error,os2);
// // 	//        os2.close();
	
// 	 time_asymptotic[log10(time)] = log10(H1err);
// 	 std::ofstream os1("time_asymptotic1D.m");
// 	 matlab_output(time_asymptotic,os1);
// 	 os1.close();
	
	
//        }
	
//      tstart = clock();
     
//      //eps *= 0.3;
//      eps *= 0.3;
//  }
    
//     //tmp = u_k;
//     //tmp.scale(&P,-1);
//      //double H1err = error_H_scale_Lshaped<Basis1D>(1, P.basis(), tmp, singGrad);
//      //double H1err = error_H_scale_interval<Basis1D>(1, P.basis(), tmp, exact1D_prime);
//      //cout << "H_1 error = " <<  H1err << endl;
//      //double L_2_error = error_H_scale_interval<Basis1D>(0, P.basis(), tmp, exact1D);
//      //double L_2_error = error_H_scale_Lshaped<Basis1D>(0, P.basis(), tmp, simple_sol);
//      //double L_2_error = error_H_scale_Lshaped<Basis1D>(0, P.basis(), tmp, exact2D);
//      //cout << "L_2 error = " << L_2_error << endl;
//     //log_10_L2_error[log10(u_k.size())] = log10(L_2_error);



//     approximations[P.basis().n_p()] = u_k;
    
//     for (int i = 0; i < P.basis().n_p(); i++) {
//       for (typename InfiniteVector<double, Index>::const_iterator it = u_k.begin(), itend = u_k.end();
// 	   it != itend; ++it)
// 	if (it.index().p() == i)
// 	  approximations[i].set_coefficient(it.index(),*it);
//     }
 


    }
}
