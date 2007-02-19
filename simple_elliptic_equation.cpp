// implementation for simple_elliptic_equation.h

#include <list>
#include <numerics/gauss_data.h>
#include <frame_index.h>
#include <frame_support.h>

using WaveletTL::CubeBasis;

namespace FrameTL
{

//   template <class IBASIS>
//   Index1D<IBASIS>::Index1D(const IntervalIndex<IBASIS>& ind,
// 			   const unsigned int p, const unsigned int der)
//     : ind_(ind), p_(p), der_(der)
//   {
//   }

//   template <class IBASIS>
//   bool
//   Index1D<IBASIS>::operator < (const Index1D<IBASIS>& lambda) const
//   {
//     return ind_ < lambda.index() ||
//       (
//        ind_ == lambda.index() &&
//        (
// 	p_ < lambda.p() ||
// 	(
// 	 p_ == lambda.p() && der_ < lambda.derivative()   
// 	 )
// 	)
//        );
//   };
  
//   template <class IBASIS>
//   bool Index1D<IBASIS>::operator == (const Index1D<IBASIS>& lambda) const
//   {
//     return (ind_ == lambda.index()) && (p_ == lambda.p()) && (der_ == lambda.derivative());
//   };

//   template <class IBASIS>
//   bool Index1D<IBASIS>::operator != (const Index1D<IBASIS>& lambda) const
//   {
//     return !(*this == lambda);
//   };
  
//   template <class IBASIS>
//   bool Index1D<IBASIS>::operator <= (const Index1D<IBASIS>& lambda) const
//   {
//     return (*this < lambda) || (*this == lambda);
//   };
  


  template <class IBASIS, unsigned int DIM>
  SimpleEllipticEquation<IBASIS,DIM>::SimpleEllipticEquation(const EllipticBVP<DIM>* ell_bvp,
							     const AggregatedFrame<IBASIS,DIM>* frame,
							     const int jmax)
    : ell_bvp_(ell_bvp), frame_(frame), jmax_(jmax)
  {
    
    compute_diagonal();
    compute_rhs();
  }

  template <class IBASIS, unsigned int DIM>
  double 
  SimpleEllipticEquation<IBASIS,DIM>::D(const typename AggregatedFrame<IBASIS,DIM>::Index& lambda) const
  {
    //return 1 << lambda.j();
    return stiff_diagonal.get_coefficient(lambda);
  }

  template <class IBASIS, unsigned int DIM>
  void
  SimpleEllipticEquation<IBASIS,DIM>::rescale(InfiniteVector<double,
					      typename AggregatedFrame<IBASIS,DIM>::Index>& coeffs,
					      const int n) const
  {
    typedef AggregatedFrame<IBASIS,DIM> Frame;
    for (typename InfiniteVector<double, typename Frame::Index>::const_iterator it(coeffs.begin());
	 it != coeffs.end(); ++it)
      {
	// TODO: implement an InfiniteVector::iterator to speed up this hack!
	coeffs.set_coefficient(it.index(), *it * pow(D(it.index()), n));
      }
  }


  template <class IBASIS, unsigned int DIM>
  void
  SimpleEllipticEquation<IBASIS,DIM>::compute_rhs()
  {
    cout << "SimpleEllipticEquation(): precompute right-hand side..." << endl;

    typedef AggregatedFrame<IBASIS,DIM> Frame;
    typedef typename Frame::Index Index;

    // precompute the right-hand side on a fine level
    InfiniteVector<double,Index> fhelp;
    const int j0   = frame_->j0();
    for (Index lambda(FrameTL::first_generator<IBASIS,DIM,DIM,Frame>(frame_,j0));; ++lambda)
      {
	const double coeff = f(lambda)/D(lambda);
	if (fabs(coeff)>1e-15) {
	  fhelp.set_coefficient(lambda, coeff);
	  //cout << lambda << endl;
	}
  	if (lambda == last_wavelet<IBASIS,DIM,DIM,Frame>(frame_,jmax_))
	  break;
      }
    fnorm_sqr = l2_norm_sqr(fhelp);

    cout << "... done, all integrals for right-hand side computed" << endl;

    // sort the coefficients into fcoeffs
    fcoeffs.resize(0); // clear eventual old values
    fcoeffs.resize(fhelp.size());
    unsigned int id(0);
    for (typename InfiniteVector<double,Index>::const_iterator it(fhelp.begin()), itend(fhelp.end());
	 it != itend; ++it, ++id)
      fcoeffs[id] = std::pair<Index,double>(it.index(), *it);
    sort(fcoeffs.begin(), fcoeffs.end(), typename InfiniteVector<double,Index>::decreasing_order());
  }

  template <class IBASIS, unsigned int DIM>
  void
  SimpleEllipticEquation<IBASIS,DIM>::compute_diagonal()
  {
    cout << "SimpleEllipticEquation(): precompute diagonal of stiffness matrix..." << endl;

    typedef AggregatedFrame<IBASIS,DIM> Frame;
    typedef typename Frame::Index Index;

    // precompute the right-hand side on a fine level
    const int j0   = frame_->j0();
    for (Index lambda(FrameTL::first_generator<IBASIS,DIM,DIM,Frame>(frame_,j0));; ++lambda)
      {
	stiff_diagonal.set_coefficient(lambda, sqrt(a(lambda,lambda)));
	if (lambda == last_wavelet<IBASIS,DIM,DIM,Frame>(frame_,jmax_))
	  break;
      }

    cout << "... done, digonal of stiffness matrix computed" << endl;
  }

  template <class IBASIS, unsigned int DIM>
  inline
  double
  SimpleEllipticEquation<IBASIS,DIM>:: integrate(const Index1D<IBASIS>& lambda,
						 const Index1D<IBASIS>& mu,
						 const int N_Gauss,
						 const int dir,
						 const typename CubeBasis<IBASIS,DIM>::Support* supp_lambda,
						 const typename CubeBasis<IBASIS,DIM>::Support* supp_mu) const
  {
    
    double res = 0;
    
    typename One_D_IntegralCache::iterator col_lb(one_d_integrals.lower_bound(lambda));
    typename One_D_IntegralCache::iterator col_it(col_lb);
    if (col_lb == one_d_integrals.end() ||
	one_d_integrals.key_comp()(lambda,col_lb->first))
      {
	// insert a new column
	typedef typename One_D_IntegralCache::value_type value_type;
	col_it = one_d_integrals.insert(col_lb, value_type(lambda, Column1D()));
      }
    
    Column1D& col(col_it->second);
    
    typename Column1D::iterator lb(col.lower_bound(mu));
    typename Column1D::iterator it(lb);
    if (lb == col.end() ||
	col.key_comp()(mu, lb->first))
      {

	// compute 1D irregular grid
	Array1D<double> irregular_grid;

	bool b = intersect_supports_1D<IBASIS,DIM>(*frame_, lambda, mu, supp_lambda, supp_mu, dir, irregular_grid);
	if (!b)
	  return 0.0;

	Array1D<double> gauss_points_la, gauss_points_mu, gauss_weights,
	  values_lambda, values_mu;
	
	//	cout << "N Gauss " << N_Gauss << endl;
 	gauss_points_la.resize(N_Gauss*(irregular_grid.size()-1));
 	gauss_points_mu.resize(N_Gauss*(irregular_grid.size()-1));
 	gauss_weights.resize(N_Gauss*(irregular_grid.size()-1));
  	for (unsigned int k = 0; k < irregular_grid.size()-1; k++)
 	  for (int n = 0; n < N_Gauss; n++) {
 	    gauss_points_la[ k*N_Gauss+n  ]
 	      = 0.5 * (irregular_grid[k+1]-irregular_grid[k]) * (GaussPoints[N_Gauss-1][n]+1)
 	      + irregular_grid[k];
	    
	    //	    	    cout << "gp = " << gauss_points_la[ k*N_Gauss+n  ] << endl;

 	    gauss_weights[ k*N_Gauss+n ]
 	      = (irregular_grid[k+1]-irregular_grid[k])*GaussWeights[N_Gauss-1][n];

	    //	    	    cout << "gw = " << gauss_weights[ k*N_Gauss+n  ] << endl;
 	  }
	
	IBASIS* basis1D_lambda = frame_->bases()[lambda.p()]->bases()[dir];// IMPORTANT: choose 1D basis of right direction !!!!!
	IBASIS* basis1D_mu     = frame_->bases()[mu.p()]->bases()[dir];
	
	const Chart<DIM>* chart_la = frame_->atlas()->charts()[lambda.p()];
	const Chart<DIM>* chart_mu = frame_->atlas()->charts()[mu.p()];
	

	//	cout << "deriv_lambda = " << lambda.derivative() << endl;
	WaveletTL::evaluate(*(basis1D_lambda), lambda.derivative(),
			    lambda.index(),
			    gauss_points_la, values_lambda);
	
	const Array1D<double> gouss_points_mu;
	Point<DIM> x; // = 0;
	Point<DIM> x_patch;
	Point<DIM> y;
	// setup mapped gauss points
	for (unsigned int i = 0; i < gauss_points_la.size(); i++) {
	  x[dir] = gauss_points_la[i];
	  chart_la->map_point(x,x_patch);
	  chart_mu->map_point_inv(x_patch,y);
	  gauss_points_mu[i] = y[dir];
	}

	//	cout << "deriv_mu = " << mu.derivative() << endl;
	WaveletTL::evaluate(*(basis1D_mu), mu.derivative(),
			    mu.index(),
			    gauss_points_mu, values_mu);
		
	//cout << "doing the job.." << endl;
	for (unsigned int i = 0; i < values_lambda.size(); i++)
	  res += gauss_weights[i] * values_lambda[i] * values_mu[i];

	typedef typename Column1D::value_type value_type;
	it = col.insert(lb, value_type(mu, res));
      }
    else {
      res = it->second;
    }
	
    return res;
  }


  template <class IBASIS, unsigned int DIM>
  double
  SimpleEllipticEquation<IBASIS,DIM>::a(const typename AggregatedFrame<IBASIS,DIM>::Index& la,
					const typename AggregatedFrame<IBASIS,DIM>::Index& nu) const
  {
    double r = 0.0;
    
    Index lambda = la;
    Index mu = nu;
    Index tmp_ind;

    typedef WaveletTL::CubeBasis<IBASIS,DIM> CUBEBASIS;

    typedef typename CUBEBASIS::Index CubeIndex;
 
    typename CUBEBASIS::Support supp_lambda_ = frame_->all_supports[lambda.number()];
    typename CUBEBASIS::Support supp_mu_ = frame_->all_supports[mu.number()];

    typename CUBEBASIS::Support tmp_supp;


    // swap indices and supports if necessary
    if (supp_mu_.j > supp_lambda_.j) {
      tmp_ind = lambda;
      lambda = mu;
      mu = tmp_ind;

      tmp_supp = supp_lambda_;
      supp_lambda_ = supp_mu_,
	supp_mu_ = tmp_supp;
    }

    const typename CUBEBASIS::Support* supp_lambda = &supp_lambda_;
    const typename CUBEBASIS::Support* supp_mu = &supp_mu_;

    const int N_Gauss = 3;

    bool b = intersect_supports<IBASIS,DIM,DIM>(*frame_, lambda, mu, supp_lambda, supp_mu);
    if ( !b )
      return 0.0;
              
    typedef typename IBASIS::Index Index_1D;
    
    const Chart<DIM>* chart_la = frame_->atlas()->charts()[lambda.p()];
    const Chart<DIM>* chart_mu = frame_->atlas()->charts()[mu.p()];

#if 1
    // loop over spatial direction
    for (int i = 0; i < (int) DIM; i++) {
      double t = 1.;
      
      for (int j = 0; j < (int) DIM; j++) {
	if (j == i)
	  continue;
	Index1D<IBASIS> i1(IntervalIndex<IBASIS> (
						  lambda.j(),lambda.e()[j],lambda.k()[j],
						  frame_->bases()[lambda.p()]->bases()[j]
						  ),
			   lambda.p(),j,0
			   );
	Index1D<IBASIS> i2(IntervalIndex<IBASIS> (mu.j(),mu.e()[j],mu.k()[j],
						  frame_->bases()[mu.p()]->bases()[j]
						  ),
			   mu.p(),j,0
			   );
	
	
	t *= integrate(i1, i2, N_Gauss, j, supp_lambda, supp_mu);
      }
      
      
      Index1D<IBASIS> i1(IntervalIndex<IBASIS> (
						lambda.j(),lambda.e()[i],lambda.k()[i],
						frame_->bases()[lambda.p()]->bases()[i]
						),
			 lambda.p(),i,1
			 );
      Index1D<IBASIS> i2(IntervalIndex<IBASIS> (mu.j(),mu.e()[i],mu.k()[i],
						frame_->bases()[mu.p()]->bases()[i]
						),
			 mu.p(),i,1
			 );
      
      t *= integrate(i1, i2, N_Gauss, i, supp_lambda, supp_mu);
      
      t *= 1./(chart_la->a_i(i) * chart_mu->a_i(i));
      
      r += t;
    }
    
    double tmp1 = 1., tmp2 = 1.;
    
    for (unsigned int i = 0; i < DIM; i++) {
      tmp1 *= chart_la->a_i(i);
      tmp2 *= chart_mu->a_i(i);
    }
    tmp1 = sqrt(fabs(tmp1)) / sqrt(fabs(tmp2));
    r *= tmp1;
    
#else
    // at this point the integration of the principal part of the operator has been
    // finished
    // now follows integration of q(x)*psi_lambda(x)*psi_mu(x)
    
    double s = 1.;
    
    // loop over spatial direction     
    for (int i = 0; i < (int) DIM; i++) {
      Index1D<IBASIS> i1(IntervalIndex<IBASIS> (
						lambda.j(),lambda.e()[i],lambda.k()[i],
						frame_->bases()[lambda.p()]->bases()[i]
						),
			 lambda.p(),i,0
			 );
      Index1D<IBASIS> i2(IntervalIndex<IBASIS> (mu.j(),mu.e()[i],mu.k()[i],
						frame_->bases()[mu.p()]->bases()[i]
						),
			 mu.p(),i,0
			 );
      s *= integrate(i1, i2, N_Gauss, i, supp_lambda, supp_mu);
    }
    
    double tmp1 = 1., tmp2 = 1.;
    
    for (unsigned int i = 0; i < DIM; i++) {
      tmp1 *= chart_la->a_i(i);
      tmp2 *= chart_mu->a_i(i);
    }
    tmp1 = sqrt(fabs(tmp1)) / sqrt(fabs(tmp2));
    
    s *= tmp1;
    
    r += s;
#endif      
    return r;
    
  }
  
  template <class IBASIS, unsigned int DIM>
  double
  SimpleEllipticEquation<IBASIS,DIM>::norm_A() const
  {
    if (normA == 0.0) {
      typedef typename AggregatedFrame<IBASIS,DIM>::Index Index;
      std::set<Index> Lambda;
      const int j0 = frame_->j0();
      const int jmax = j0+1;
      for (Index lambda = FrameTL::first_generator<IBASIS,DIM,DIM,Frame>(frame_,j0);; ++lambda) {
	Lambda.insert(lambda);
	if (lambda == FrameTL::last_wavelet<IBASIS,DIM,DIM,Frame>(frame_,jmax)) break;
      }
      SparseMatrix<double> A_Lambda;
      setup_stiffness_matrix(*this, Lambda, A_Lambda);
      
      Vector<double> xk(Lambda.size(), false);
      xk = 1;
      unsigned int iterations;
      normA = PowerIteration(A_Lambda, xk, 1e-6, 100, iterations);
    }

    return normA;
  }

  template <class IBASIS, unsigned int DIM>
  double
  SimpleEllipticEquation<IBASIS,DIM>::s_star() const
  {
    // notation from [St04a]
    const int t = operator_order();
    const int n = DIM;
    const int dT = frame_->bases()[0]->primal_vanishing_moments(); // we assume to have the same 'kind'
                                                                   // of wavelets on each patch, so use
                                                                   // patch 0 as reference case
    const double gamma = frame_->bases()[0]->primal_regularity();
    
    return (n == 1
	    ? t+dT // [St04a], Th. 2.3 for n=1
	    : std::min((t+dT)/(double)n, (gamma-t)/1./*(n-1.)*/)); // [St04a, Th. 2.3]
  }

  template <class IBASIS, unsigned int DIM>
  double
  SimpleEllipticEquation<IBASIS,DIM>::f(const typename AggregatedFrame<IBASIS,DIM>::Index& lambda) const
  {

    //\int_\Omega f(Kappa(x)) \psi^\Box (x) \sqrt(\sqrt(det ((D Kappa)^T(x) (D Kappa)(x)) ))
    //recall: \sqrt(\sqrt(...)) = Gram_factor

    double r = 0;

    const unsigned int p = lambda.p();
 
    typedef WaveletTL::CubeBasis<IBASIS,DIM> CUBEBASIS;
 
    // first compute supp(psi_lambda)
 
    typename CUBEBASIS::Support supp;
    
    typename CubeBasis<IBASIS,DIM>::Index lambda_c(lambda.j(),
						   lambda.e(),
						   lambda.k(),frame_->bases()[p]);
    
    WaveletTL::support<IBASIS,DIM>(*(frame_->bases()[p]), lambda_c, supp);

    const Chart<DIM>* chart = frame_->atlas()->charts()[p];

    // setup Gauss points and weights for a composite quadrature formula:
    const int N_Gauss = 6;
    //const double h = ldexp(1.0, -supp.j); // granularity for the quadrature
    const double h = 1.0 / (1 << supp.j); // granularity for the quadrature
    FixedArray1D<Array1D<double>,DIM> gauss_points, gauss_weights, v_values;
    for (unsigned int i = 0; i < DIM; i++) {
      gauss_points[i].resize(N_Gauss*(supp.b[i]-supp.a[i]));
      gauss_weights[i].resize(N_Gauss*(supp.b[i]-supp.a[i]));
      for (int patch = supp.a[i]; patch < supp.b[i]; patch++)
	for (int n = 0; n < N_Gauss; n++) {
	  gauss_points[i][(patch-supp.a[i])*N_Gauss+n]
	    = h*(2*patch+1+GaussPoints[N_Gauss-1][n])/2.;
	  gauss_weights[i][(patch-supp.a[i])*N_Gauss+n]
	    = h*GaussWeights[N_Gauss-1][n];
	}
    }

    // compute the point values of the integrand (where we use that it is a tensor product)
    for (unsigned int i = 0; i < DIM; i++) {
      WaveletTL::evaluate(*(frame_->bases()[p]->bases()[i]), 0,
			  typename IBASIS::Index(lambda.j(),
						 lambda.e()[i],
						 lambda.k()[i],
						 frame_->bases()[p]->bases()[i]),
			  gauss_points[i], v_values[i]);
    }

    // iterate over all points and sum up the integral shares
    int index[DIM]; // current multiindex for the point values
    for (unsigned int i = 0; i < DIM; i++)
      index[i] = 0;
    
    Point<DIM> x;
    Point<DIM> x_patch;
    
    while (true) {
      for (unsigned int i = 0; i < DIM; i++)
	x[i] = gauss_points[i][index[i]];

      chart->map_point(x,x_patch);

      double share = ell_bvp_->f(x_patch) * chart->Gram_factor(x);
      for (unsigned int i = 0; i < DIM; i++)
	share *= gauss_weights[i][index[i]] * v_values[i][index[i]];
      r += share;

      // "++index"
      bool exit = false;
      for (unsigned int i = 0; i < DIM; i++) {
	if (index[i] == N_Gauss*(supp.b[i]-supp.a[i])-1) {
	  index[i] = 0;
	  exit = (i == DIM-1);
	} else {
	  index[i]++;
	  break;
	}
      }
      if (exit) break;
    }

#if 1
    return r;
#endif
#if 0
    assert(DIM == 1);
    double tmp = 1;
    Point<DIM> p1;
    p1[0] = 0.5;
    Point<DIM> p2;
    chart->map_point_inv(p1,p2);
    tmp =  WaveletTL::evaluate(*(frame_->bases()[p]->bases()[0]), 0,
			       typename IBASIS::Index(lambda.j(),
						      lambda.e()[0],
						      lambda.k()[0],
						      frame_->bases()[p]->bases()[0]),
			       p2[0]);
    tmp /= chart->Gram_factor(p2);
  
  
    return 4.0*tmp + r;
#endif
  }

  template <class IBASIS, unsigned int DIM>
  void
  SimpleEllipticEquation<IBASIS,DIM>::RHS
  (const double eta,
   InfiniteVector<double, typename AggregatedFrame<IBASIS,DIM>::Index>& coeffs) const
  {
    coeffs.clear();
    double coarsenorm(0);
    double bound(fnorm_sqr - eta*eta);
    typedef typename AggregatedFrame<IBASIS,DIM>::Index Index;
    typename Array1D<std::pair<Index, double> >::const_iterator it(fcoeffs.begin());
    do {
      coarsenorm += it->second * it->second;
      coeffs.set_coefficient(it->first, it->second);
      ++it;
    } while (it != fcoeffs.end() && coarsenorm < bound);
  }

  template <class IBASIS, unsigned int DIM>
  void
  SimpleEllipticEquation<IBASIS,DIM>::set_bvp(const EllipticBVP<DIM>* bvp)
  {
    ell_bvp_ = bvp;
    compute_diagonal();
    compute_rhs();

  }

  template <class IBASIS, unsigned int DIM>
  void
  SimpleEllipticEquation<IBASIS,DIM>::add_level (const Index& lambda,
						 InfiniteVector<double, Index>& w, const int j,
						 const double factor,
						 const int J,
						 const CompressionStrategy strategy) const
  {
    typedef std::list<Index> IntersectingList;
    IntersectingList nus;
    if (strategy == CDD1) {
      // compute all wavelets on level j, such that supp(psi_lambda) and supp(psi_nu) intersect
      intersecting_wavelets(basis(), lambda,
			    std::max(j, basis().j0()),
			    j == (basis().j0()-1),
			    nus);
      
      // traverse the matrix block and update the result
      const double d1 = D(lambda);
      for (typename IntersectingList::const_iterator it(nus.begin()), itend(nus.end());
	   it != itend; ++it) {
	const double entry = a(*it, lambda) / (d1*D(*it));
	w.add_coefficient(*it, entry * factor);
      }
    }
    else if (strategy == St04a) {
      // compute all wavelets on level j, such that supp(psi_lambda) and supp(psi_nu) intersect
      intersecting_wavelets(basis(), lambda,
			    std::max(j, basis().j0()),
			    j == (basis().j0()-1),
			    nus);
      // traverse the matrix block and update the result
      const double d1 = D(lambda);
      for (typename IntersectingList::const_iterator it(nus.begin()), itend(nus.end());
	   it != itend; ++it)
	if (abs(lambda.j()-j) <= J/((double) space_dimension) ||
	    intersect_singular_support(basis(), lambda, *it)) {
	  const double entry = a(*it, lambda) / (d1*D(*it));
	  w.add_coefficient(*it, entry * factor);
	}
    }
  }
  
}
