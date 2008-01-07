// implementation for steepest_descent.h

#include <cmath>
#include <set>
#include <utils/plot_tools.h>
#include <adaptive/apply.h>
#include <numerics/corner_singularity.h>
#include <frame_evaluate.h>
#include <cdd1_local.h>
#include <error_H_scale.h>

using std::set;

namespace FrameTL
{

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

  /*!
  */
  template<class VALUE = double>
  class Singularity1D_2_prime
    : public Function<1, VALUE>
  {
  public:
    Singularity1D_2_prime() {};
    virtual ~Singularity1D_2_prime() {};
    VALUE value(const Point<1>& p,
		const unsigned int component = 0) const
    {

      if (0. <= p[0] && p[0] < 0.5)
	return -cos(3.*M_PI*p[0])*3*M_PI + 4.*p[0];

      if (0.5 <= p[0] && p[0] <= 1.0)
	return -cos(3.*M_PI*p[0])*3*M_PI - 4.*(1-p[0]);

      return 0.;

    }
  
    void vector_value(const Point<1> &p,
		      Vector<VALUE>& values) const { ; }
  
  };


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


// 	if (it.index().p() == 0)
// 	  u_very_sparse.set_coefficient(it.index(), *it);

	if ((it.index().p() == 1) && !((supp->a[0] > -OVERLAP) && (supp->b[1] < 0.)) ) {
	  u_sparse.set_coefficient(it.index(), *it);
	}


//	if ((it.index().p() == 1) && ( (supp->b[1] > 0.)) ) {
//	  u_sparse.set_coefficient(it.index(), *it);
//	}
      }
    }
    else if (i==1) {
      for (; it != u.end(); ++it) {
	const SuppType* supp = &(P.basis().all_patch_supports[it.index().number()]);
	if ((it.index().p() == 0) && (supp->a[0] < 0.) && (0. < supp->b[0])) {
	  u_very_sparse.set_coefficient(it.index(), *it);
	}
// 	if (it.index().p() == 1)
// 	  u_very_sparse.set_coefficient(it.index(), *it);

	
	if ( (it.index().p() == 0) && ( supp->b[0] > 0. ) ) {
	  //cout << supp->a[0] << " " << supp->b[0] << " " << supp->a[1] << " " << supp->b[1] << endl;
	  u_sparse.set_coefficient(it.index(), *it);
	}

      }
    }
#endif

  }

  template <class PROBLEM>
  void  MultSchw(const PROBLEM& P, const double epsilon,
		 Array1D<InfiniteVector<double, typename PROBLEM::Index> >& approximations)
  {
#ifdef ONE_D
    Singularity1D_2<double> exact1D;
    Singularity1D_2_prime<double> exact1D_prime;
#endif
    
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
    typedef typename PROBLEM::WaveletBasis Frame;

#ifdef TWO_D    
    const double mu = 1.6544; // = energy norm of the exact solution
    const double M = 1.; // spectral radius of stiffness matrix
    //const double rho = sqrt(1 - 1.0/4.0); // [Xu92] Theorem 4.4 with \omega_1=1, K_0 = K_1=1
    const double rho = 0.2; // tuned
    const double sigma = 4.;
#endif
#ifdef ONE_D
    const double mu = 7.44609; // = energy norm of the exact solution
    const double M = 1.; // spectral radius of stiffness matrix
    //const double rho = sqrt(1 - 1.0/4.0); // [Xu92] Theorem 4.4 with \omega_1=1, K_0 = K_1=1
    const double rho = 0.2; // tuned
    const double sigma = 4.;
#endif
 
    const int K = std::max(1,(int)ceil(log(1.0/(2.0 * M * sigma)) / log(rho)));

    const int L = (int)ceil(log(epsilon / mu ) / log(2.0*pow(rho,K)*M*sigma));

    cout << "mu = " << mu << endl
	 << "M = " << M << endl
	 << "rho = " <<  rho << endl
	 << "sigma = " << sigma << endl
	 << "K = " << K << endl
	 << "L = " << L << endl;

    InfiniteVector<double, Index> f, w, r, tmp, tmp_w;
    InfiniteVector<double, Index> u_k, u_k_sparse,u_k_very_sparse;
    InfiniteVector<double, Index> precond_r_i;

    Array1D<InfiniteVector<double, Index> > xks(P.basis().n_p()); // stores local approximations which are used as starting vectors
                                                                  // for the next cycle of local solves
    map<double,double> log_10_H1_error;    

    double time = 0.;
    clock_t tstart, tend;
    tstart = clock();
    const int m = P.basis().n_p();
    int k = 0;

    for (int l = 1; l < L; l++) {
      for (int p = 1; p <= K; p++) {
	for (int i = 0; i < P.basis().n_p(); i++) {
	  k = (l-1)*m*K+(p-1)*m+i+1;
	  cout << "################################" << endl;
	  cout << "number of iteration = " << k << endl;
	  cout << "################################" << endl;
	  double local_eps = mu*pow(2.0*pow(rho,K)*M*sigma,l-1)*pow(rho,p)/(m*K);
	  cout << "tolerance for solution of local problem = " << local_eps << endl;

	  // approximate right hand side for local problem
	  P.RHS(local_eps, i, f);
	  thin_out(P, i, u_k, u_k_sparse, u_k_very_sparse);
	  APPLY(P, i, u_k_very_sparse, local_eps, w, jmax, CDD1);
	  r = f-w;

	  precond_r_i.clear();

	  // ####### solution of local problem #######
	  set<Index> Lambda_i;
	  r.support(Lambda_i);
	  
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
	  // #########################################
	  
	  u_k = precond_r_i + u_k_sparse;
	  
	  cout << "degrees of freedom: " << u_k.size() << endl;
	  xks[i].clear();
	  xks[i] = precond_r_i;
	  
	  tend = clock();
	  time += ((double) (tend-tstart))/((double) CLOCKS_PER_SEC);
	  
	  if (k % 1 == 0) {
	    tmp = u_k;
	    tmp.scale(&P,-1);

#ifdef ONE_D	    
	    double H1err = error_H_scale_interval<Basis1D>(1, P.basis(), tmp, exact1D_prime);
	    cout << "H_1 error = " <<  H1err << endl;
	    
	    log_10_H1_error[log10(u_k.size())] = log10(H1err);
	    std::ofstream os2a("log_10_H1_error1D.m");
	    matlab_output(log_10_H1_error,os2a);
	    os2a.close();
#endif
	

#ifdef TWO_D
	   double H1err = error_H_scale_Lshaped<Basis1D>(1, P.basis(), tmp, singGrad);
	   cout << "H_1 error = " <<  H1err << endl;
	   
	   log_10_H1_error[log10(u_k.size())] = log10(H1err);
	   std::ofstream os2a("log_10_H1_error2D.m");
	   matlab_output(log_10_H1_error,os2a);
	   os2a.close();
#endif
	  }
	  tstart = clock();
	  

      	}// end loop over patches
	cout << "############## full cycle of local solves completed ##############"<< endl;
	
      }// end loop K
      double coarse_tol = (sigma - 1/M)*2.0*pow(rho,K)*M*sigma*mu*pow((2*pow(rho,K)*M*sigma),l-1);
      cout << "tolerance for coarsening = " << coarse_tol << endl;
      u_k.COARSE(coarse_tol, tmp_w);
      u_k = tmp_w;
      cout << "degrees of freedom after coarsening: " << u_k.size() << endl;
    }// end loop L


    // collect final approximation and its local parts
    approximations[P.basis().n_p()] = u_k;
    
    for (int i = 0; i < P.basis().n_p(); i++) {
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
