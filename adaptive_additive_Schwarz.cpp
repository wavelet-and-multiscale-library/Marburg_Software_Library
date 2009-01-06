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
      return (-sin(3.*M_PI*p[0])*9.*M_PI*M_PI - 4.);
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
  double compute_exact_residual_norm (const PROBLEM& P,
				      InfiniteVector<double, typename PROBLEM::Index>& u_k) {
    typedef typename PROBLEM::Index Index;
    typedef typename PROBLEM::WaveletBasis::IntervalBasis Basis1D;
    typedef typename PROBLEM::WaveletBasis Frame;
    const int jmax = JMAX;
    set<Index> Lambda1;
    for (Index lambda = FrameTL::first_generator<Basis1D,1,1,Frame>(&P.basis(), P.basis().j0());
	 lambda <= FrameTL::last_wavelet<Basis1D,1,1,Frame>(&P.basis(), jmax); ++lambda) {
      Lambda1.insert(lambda);
    }
    set<Index> Lambda2;
    u_k.support(Lambda2);
    SparseMatrix<double> A(Lambda1.size(), Lambda2.size());
    WaveletTL::setup_stiffness_matrix(P, Lambda1, Lambda2, A);

    // copy u_k in vector
    Vector<double> U_k(u_k.size());
    Vector<double> A_U_k(Lambda1.size());
    unsigned int id = 0;
    typename set<Index>::const_iterator it = Lambda2.begin();
    for (; it != Lambda2.end(); ++it, ++id) {
      U_k[id] = u_k.get_coefficient(*it);
    }

    Vector<double> F(Lambda1.size());    
    WaveletTL::setup_righthand_side(P, Lambda1, F);

    A.apply(U_k, A_U_k);
    return l2_norm(F-A_U_k);

  }


  // delete all entries corresponding to patch i
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
#ifdef ONE_D
    Singularity1D_2<double> exact1D;
    Singularity1D_2_prime<double> exact1D_prime;
#endif
    
    typedef typename PROBLEM::WaveletBasis::IntervalBasis Basis1D;
    const int jmax = JMAX;
    typedef typename PROBLEM::Index Index;
    typedef typename PROBLEM::WaveletBasis Frame;

    InfiniteVector<double, Index> f, w, r, tmp, tmp_w;
    InfiniteVector<double, Index> u_k, u_k_sparse,u_k_very_sparse;
    InfiniteVector<double, Index> precond_r_i;
    
    Array1D<InfiniteVector<double, Index> > xks(P.basis().n_p()); // stores local approximations which are used as starting vectors
                                                                  // for the next cycle of local solves

    Array1D<InfiniteVector<double, Index> > uks(P.basis().n_p());


    double eps = 1.0;
    const double alpha = 1.0;
    const double rho = 0.8;

    for (int k = 0; k < 25; k++) {
      
      // loop over patches
      for (int i = 0; i < P.basis().n_p(); i++) {
#ifdef SPARSE
	thin_out(P, i, u_k, u_k_sparse, u_k_very_sparse);
#endif
	  // CDD1
	cout << "entering CDD solver..." << endl;
#ifdef SPARSE
	tmp.clear();
	tmp = u_k;
	//remove_i<PROBLEM>(i, tmp);
	tmp=u_k-((1./alpha)*(tmp-u_k_sparse));
	//cout << "rhs " << tmp << endl;
	//cout << "i = " << i << " size = " << tmp.size() << endl; 
	CDD1_LOCAL_SOLVE(P, i, eps, xks[i], precond_r_i, tmp, jmax, CDD1);
	//CDD1_LOCAL_SOLVE(P, i, eps, xks[i], precond_r_i, u_k, jmax, CDD1);
#endif
	cout << "CDD 1 solve completed, size of output is " << precond_r_i.size() << endl;
	
	xks[i] = precond_r_i;

	//cout << precond_r_i << endl;


#ifdef SPARSE
 	tmp.clear();
 	tmp = u_k;
 	for (int j = 0; j < P.basis().n_p(); j++) {
 	  if (j!=i)
 	    remove_i<PROBLEM>(j, tmp);
 	}
	uks[i].clear();
	// new
	uks[i] = u_k_sparse + (alpha*precond_r_i);
	//uks[i] = (tmp + u_k_sparse) + (alpha*precond_r_i);
	//uks[i] = u_k + (alpha*precond_r_i);

	//cout << uks[i] << endl;
#endif

	cout << "degrees of freedom: " << u_k.size() << endl;
	xks[i].clear();
	xks[i] = precond_r_i;
      }

      tmp.clear();
      for (int j = 0; j < P.basis().n_p(); j++) {
	tmp = tmp + uks[j];
      }
      double m = ((double) P.basis().n_p())-1;
      //u_k = tmp - (m * u_k);
      u_k = (1./(m+1))*tmp;

      tmp.clear();
      u_k.COARSE(eps, tmp);
      u_k = tmp;
      
      
      //cout << u_k << endl;
      eps*=rho;


    }
    //u_k=uks[1];
    //cout << u_k << endl;
    //cout << "######################" << endl;

      
    // collect final approximation and its local parts
    approximations[P.basis().n_p()] = u_k;

   
    for (int i = 0; i < P.basis().n_p(); i++) {
      approximations[i].clear();
      for (typename InfiniteVector<double, Index>::const_iterator it = u_k.begin(), itend = u_k.end();
	   it != itend; ++it)
	if (it.index().p() == i) {
	  approximations[i].set_coefficient(it.index(),*it);
	}
      //cout << approximations[i] << endl;
      //cout << "######################" << endl;
    }



  }

}
