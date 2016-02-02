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
  // 1-i and that correspond to wavelets which intersect with patch i,
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
#if PATCHES == 2
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
#if PATCHES == 3
    if (i==0) {
      for (; it != u.end(); ++it) {
	const SuppType* supp = &(P.basis().all_patch_supports[it.index().number()]);
	// check whether first line is intersected
	if ((it.index().p() != 0) && (supp->b[0] > -0.5) && (supp->a[1] < 0.) && (0. < supp->b[1]) ) {
	  u_very_sparse.set_coefficient(it.index(), *it);
	}
	// check whether second line is intersected
	if ((it.index().p() != 0) && (supp->a[1] < 0.) && (supp->a[0] < -0.5) && (-0.5 < supp->b[0])) {
	  u_very_sparse.set_coefficient(it.index(), *it);
	}
	if ((it.index().p() != 0) && ((supp->a[0] < -0.5) || (supp->b[1] > 0.)) ) {
	  u_sparse.set_coefficient(it.index(), *it);
	}
      }
    }
    if (i==1) {
      for (; it != u.end(); ++it) {
	const SuppType* supp = &(P.basis().all_patch_supports[it.index().number()]);
	// check whether first line is intersected
	if ((it.index().p() != 1) && (supp->b[1] > -0.5) && (supp->a[0] < 0.) && (0. < supp->b[0]) ) {
	  u_very_sparse.set_coefficient(it.index(), *it);
	}
	// check whether second line is intersected
	if ((it.index().p() != 1) && (supp->a[0] < 0.) && (supp->a[1] < -0.5) && (-0.5 < supp->b[1])) {
	  u_very_sparse.set_coefficient(it.index(), *it);
	}
	if ((it.index().p() != 1) && ((supp->a[1] < -0.5) || (supp->b[0] > 0.)) ) {
	  u_sparse.set_coefficient(it.index(), *it);
	}
      }
    }
    if (i==2) {
      for (; it != u.end(); ++it) {
	const SuppType* supp = &(P.basis().all_patch_supports[it.index().number()]);
	// check whether first line is intersected
	if ((it.index().p() != 2) &&
	    (
	    ((supp->a[0] < 0.) && (0. < supp->b[0]))
	    ||
	    ((supp->a[1] < 0.) && (0. < supp->b[1]))
	    )
	    ) {
	  u_very_sparse.set_coefficient(it.index(), *it);
	}
	if ((it.index().p() != 2) && (((supp->b[1] > 0.) || (supp->b[0] > 0.))) ) {
	  u_sparse.set_coefficient(it.index(), *it);
	}
      }
    }
#endif
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
  void  AddSchw(const PROBLEM& P, const double epsilon,
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
    Point<2> origin;
    origin[0] = 0.0;
    origin[1] = 0.0;
    
    CornerSingularity exact2D(origin, 0.5, 1.5);
    CornerSingularityGradient singGrad(origin, 0.5, 1.5);
#endif

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
    const double mu = 4.0*1.6544; // energy norm of the exact solution
#else
    // Energy norm of the exact solution of the Poisson equation
    // in the L-shaped domain.
    const double mu = 1.6544; // energy norm of the exact solution
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


    //const double rho = sqrt(1 - 1.0/4.0); // [Xu92] Theorem 4.4 with \omega_1=1, K_0 = K_1=1

    // Next we have to set up the error reduction rate of the iterative solver.
    // Instead of using theoretical upper bounds, we work with manually chosen
    // values gained from some experiments with the full system matrix.
    //const double rho = sqrt(1 - 1.0/4.0); // [Xu92] Theorem 4.4 with \omega_1=1, K_0 = K_1=1
#ifdef RINGDOMAIN
    const double rho = 0.2;
#else
    // We try the same value as for the multiplicative case.
    const double rho = 0.2996;
#endif

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

    // number of patches
    const int m = P.basis().n_p();

    // #####################################################################################
    // Setup of constants.
    // #####################################################################################
    const double C = m;

    const double sigma = std::max(1./M, C + 1./M) + 0.1;//10.0
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
      L = 8;
#else
      L = 23; // j0=4, then take L=16, j0=3, then take L=14
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
    InfiniteVector<double, Index> f, w, r, tmp, tmp2;
    InfiniteVector<double, Index> u_k, u_k_sparse,u_k_very_sparse;
    InfiniteVector<double, Index> precond_r_i;

    Array1D<InfiniteVector<double, Index> > xks(m); // stores local approximations which are used as starting vectors
                                                                  // for the next cycle of local solves
    Array1D<InfiniteVector<double, Index> > uks(m);

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

    // relaxation parameter
#ifdef RINGDOMAIN
    const double alpha = 0.5;
#else
    // if the relaxation parameter is chosen equal to 1/m, then
    // the right-hand sides for the local auxiliary problems
    // are as sparse as for the multiplicative algorithm, cf. Manuel's
    // PhD thesis page 165/166 and Remark 6.4.
    const double alpha = 1.0/m;
#endif

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
      for (int p = 1; p <= K; p++) {
	for (int i = 0; i < m; i++) {
	  k = (l-1)*m*K+(p-1)*m+i+1;
	  cout << "################################" << endl;
	  cout << "number of iteration = " << k << endl;
	  cout << "################################" << endl;

	  // Setup tolerance for the solution of the local problems.
	  // This needs to be manually tuned not to end up with
	  // slow performance.
#ifdef RINGDOMAIN
	  local_eps = 10.0*mu*pow(2.0*pow(rho,K)*M*sigma,l-1)*pow(rho,p)/(alpha*m*K);
#else
	  local_eps = 10.0*mu*pow(2.0*pow(rho,K)*M*sigma,l-1)*pow(rho,p)/(alpha*m*K); // Faktor 10 bei 3 Patches
	  //local_eps = 100.0*mu*pow(2.0*pow(rho,K)*M*sigma,l-1)*pow(rho,p)/(alpha*m*K); // Faktor 100 bei Test mit 2 Patches
#endif
	  cout << "tolerance for solution of local problem = " << local_eps << endl;
	  precond_r_i.clear();


#ifdef RINGDOMAIN
	  // Preparation for the ring-shaped domain case.
	  thin_out_ring(P, i, u_k, u_k_sparse, u_k_very_sparse);
#else
	  // Preparation for the L-shaped domain case.
	  thin_out(P, i, u_k, u_k_sparse, u_k_very_sparse);
#endif
	  tmp = u_k-u_k_sparse;
	  tmp.compress(1.0e-15);
	  tmp = u_k-((1./(m*alpha))*tmp);
	  tmp.compress(1.0e-15);

// 	  tmp = u_k;
// 	  tmp=u_k-((1./(m*alpha))*(tmp-u_k_sparse));
	  
	  cout << "entering CDD solver..." << endl;

	  // Solution of the local problem in case we use the sparse version of the algorithm
	  // as proposed Manuel's PhD thesis, where it is proposed to throw away
	  // all degrees of freedom that are contained in the current subdomain before the local solve.
#ifdef RINGDOMAIN
	  CDD1_LOCAL_SOLVE(P, i, local_eps, xks[i], precond_r_i, tmp, jmax, CDD1);
#else
	  CDD1_LOCAL_SOLVE(P, i, local_eps, xks[i], precond_r_i, tmp, jmax, CDD1);
#endif
	  // Store the calculated local solution. These are always used as the initial guess in the
	  // next call of CDD1_LOCAL_SOLVE.
	  xks[i] = precond_r_i;

	  // setup global(!) intermediate approximation
	  uks[i] = u_k_sparse + (m*alpha*precond_r_i);
	  //xks[i] = precond_r_i;	  
	} // end loop over the patches

	// compute the average of the global intermediate iterates
	tmp.clear();
	for (int j = 0; j < m; j++) {
	  tmp = tmp + uks[j];
	}
	u_k = (1./m)*tmp;
	cout << "degrees of freedom: " << u_k.size() << endl;
      } // end loop p
      
      // setup tolerance for coarsening
#ifdef RINGDOMAIN
      double coarse_tol = (sigma - 1./M)*2.0*pow(rho,K)*mu*pow(2.0*pow(rho,K)*M*sigma,l-1)*0.1;
#else
      double coarse_tol = (sigma - 1./M)*2.0*pow(rho,K)*mu*pow(2.0*pow(rho,K)*M*sigma,l-1);
#endif
      cout << "tolerance for coarsening = " << coarse_tol << endl;
      cout << "norm of u_k = " << l2_norm(u_k) << endl;
      u_k.COARSE(coarse_tol, tmp2);
      u_k = tmp2;
      cout << "degrees of freedom after coarsening: " << u_k.size() << endl;

      // #####################################################################################
      // Approximate global EXACT residual and perform output.
      // #####################################################################################
      tend = clock();
      time += ((double) (tend-tstart))/((double) CLOCKS_PER_SEC);
      
      tmp = u_k;
      tmp.scale(&P,-1);

      //compute global residual
      tmp2.clear(); 
      for (int i = 0; i < m; i++) {
	P.RHS(1.0e-8, i, f);
	cout << "fsize exact res = " << f.size() << endl;
      	APPLY(P, i, u_k, 1.0e-8, w, jmax, CDD1);
      	tmp2 += f-w;
      }
      double residual_norm = l2_norm(tmp2);
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
	sprintf(name1, "%s%d%s%d%s", "./as_results22/as1D_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name2, "%s%d%s%d%s", "./as_results22/as1D_time_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name3, "%s%d%s%d%s", "./as_results22/as1D_weak_ell_tau_norms_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name4, "%s%d%s%d%s", "./as_results22/as1D_H1err_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name5, "%s%d%s%d%s", "./as_results22/as1D_H1err_time_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name6, "%s%d%s%d%s", "./as_results22/as1D_L2err_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name7, "%s%d%s%d%s", "./as_results22/as1D_L2err_time_P_jmax18_", d, "_dT", dT, ".m");
	break;
      }
      case 3: {
	weak_ell_tau_norms[k] = u_k.weak_norm(1./2.5);
	sprintf(name1, "%s%d%s%d%s", "./as_results33/as1D_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name2, "%s%d%s%d%s", "./as_results33/as1D_time_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name3, "%s%d%s%d%s", "./as_results33/as1D_weak_ell_tau_norms_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name4, "%s%d%s%d%s", "./as_results33/as1D_H1err_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name5, "%s%d%s%d%s", "./as_results33/as1D_H1err_time_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name6, "%s%d%s%d%s", "./as_results33/as1D_L2err_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name7, "%s%d%s%d%s", "./as_results33/as1D_L2err_time_P_jmax18_", d, "_dT", dT, ".m");
	break;
      }
      case 4: {
	weak_ell_tau_norms[k] = u_k.weak_norm(1./3.5);
	sprintf(name1, "%s%d%s%d%s", "./as_results46/as1D_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name2, "%s%d%s%d%s", "./as_results46/as1D_time_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name3, "%s%d%s%d%s", "./as_results46/as1D_weak_ell_tau_norms_P_jmax18_", d, "_dT", dT, ".m");  
	sprintf(name4, "%s%d%s%d%s", "./as_results46/as1D_H1err_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name5, "%s%d%s%d%s", "./as_results46/as1D_H1err_time_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name6, "%s%d%s%d%s", "./as_results46/as1D_L2err_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name7, "%s%d%s%d%s", "./as_results46/as1D_L2err_time_P_jmax18_", d, "_dT", dT, ".m");
	break;
      }
      };
      log_10_residual_error[log10(u_k.size())] = log10(residual_norm);
      std::ofstream os2c(name1);
      matlab_output(log_10_residual_error,os2c);
      os2c.close();

      log_10_residual_error_time[log10(time)] = log10(residual_norm);
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
	sprintf(name1, "%s%d%s%d%s", "./as_results2D_22/as2D_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name2, "%s%d%s%d%s", "./as_results2D_22/as2D_time_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name3, "%s%d%s%d%s", "./as_results2D_22/as2D_weak_ell_tau_norms_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name4, "%s%d%s%d%s", "./as_results2D_22/as2D_H1err_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name5, "%s%d%s%d%s", "./as_results2D_22/as2D_H1err_time_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name6, "%s%d%s%d%s", "./as_results2D_22/as2D_L2err_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name7, "%s%d%s%d%s", "./as_results2D_22/as2D_L2err_time_P_jmax18_", d, "_dT", dT, ".m");
	break;
      }
      case 3: {
	weak_ell_tau_norms[k] = u_k.weak_norm(1./1.5);
	sprintf(name1, "%s%d%s%d%s", "./as_results2D_33_3patch/as2D_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name2, "%s%d%s%d%s", "./as_results2D_33_3patch/as2D_time_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name3, "%s%d%s%d%s", "./as_results2D_33_3patch/as2D_weak_ell_tau_norms_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name4, "%s%d%s%d%s", "./as_results2D_33_3patch/as2D_H1err_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name5, "%s%d%s%d%s", "./as_results2D_33_3patch/as2D_H1err_time_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name6, "%s%d%s%d%s", "./as_results2D_33_3patch/as2D_L2err_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name7, "%s%d%s%d%s", "./as_results2D_33_3patch/as2D_L2err_time_P_jmax18_", d, "_dT", dT, ".m");
	break;
      }
      case 4: {
	weak_ell_tau_norms[k] = u_k.weak_norm(1./2.0);
	sprintf(name1, "%s%d%s%d%s", "./as_results2D_46/as2D_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name2, "%s%d%s%d%s", "./as_results2D_46/as2D_time_asymptotic_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name3, "%s%d%s%d%s", "./as_results2D_46/as2D_weak_ell_tau_norms_P_jmax18_", d, "_dT", dT, ".m");  
	sprintf(name4, "%s%d%s%d%s", "./as_results2D_46/as2D_H1err_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name5, "%s%d%s%d%s", "./as_results2D_46/as2D_H1err_time_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name6, "%s%d%s%d%s", "./as_results2D_46/as2D_L2err_P_jmax18_", d, "_dT", dT, ".m");
	sprintf(name7, "%s%d%s%d%s", "./as_results2D_46/as2D_L2err_time_P_jmax18_", d, "_dT", dT, ".m");
	break;
      }
      };
      log_10_residual_error[log10(u_k.size())] = log10(residual_norm);
      std::ofstream os2c(name1);
      matlab_output(log_10_residual_error,os2c);
      os2c.close();

      log_10_residual_error_time[log10(time)] = log10(residual_norm);
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
    
    // collect final approximation and its local parts
    approximations[P.basis().n_p()] = u_k;
    //approximations[P.basis().n_p()] = uks[0];
    
    for (int i = 0; i < P.basis().n_p(); i++) {
      approximations[i].clear();
      for (typename InfiniteVector<double, Index>::const_iterator it = u_k.begin(), itend = u_k.end();
	   it != itend; ++it)
	if (it.index().p() == i)
	  approximations[i].set_coefficient(it.index(),*it);
    }
    

  }
}
