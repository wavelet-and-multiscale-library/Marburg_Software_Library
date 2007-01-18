// implementation for biharmonic_equation.h

#include <list>
#include <numerics/gauss_data.h>
#include <frame_index.h>
#include <frame_support.h>
//#include <cuba.h>

//#include <cube/cube_support.h>

using WaveletTL::CubeBasis;

namespace FrameTL
{

  template <class IBASIS>
  Index1D<IBASIS>::Index1D(const IntervalIndex<IBASIS>& ind,
			   const unsigned int p, const unsigned int dir,
			   const unsigned int der)
    : ind_(ind), p_(p), dir_(dir), der_(der)
  {
  }

  template <class IBASIS>
  bool
  Index1D<IBASIS>::operator < (const Index1D<IBASIS>& lambda) const
  {
    return ind_ < lambda.index() ||
      (
       ind_ == lambda.index() &&
       (
	p_ < lambda.p() ||
	(
	 p_ == lambda.p() &&
	 (
	  dir_ < lambda.direction() ||
	  (
	   dir_ == lambda.direction() && der_ < lambda.derivative()   
	   )
	  )
	 )
	)
       );
  };
  
  template <class IBASIS>
  bool Index1D<IBASIS>::operator == (const Index1D<IBASIS>& lambda) const
  {
    return (ind_ == lambda.index()) && (p_ == lambda.p()) && (der_ == lambda.derivative());
  };

  template <class IBASIS>
  bool Index1D<IBASIS>::operator != (const Index1D<IBASIS>& lambda) const
  {
    return !(*this == lambda);
  };
  
  template <class IBASIS>
  bool Index1D<IBASIS>::operator <= (const Index1D<IBASIS>& lambda) const
  {
    return (*this < lambda) || (*this == lambda);
  };
  


  template <class IBASIS, unsigned int DIM>
  BiharmonicEquation<IBASIS,DIM>::BiharmonicEquation(const BiharmonicBVP<DIM>* bih_bvp,
						     const AggregatedFrame<IBASIS,DIM>* frame,
						     QuadratureStrategy qstrat)
    : bih_bvp_(bih_bvp), frame_(frame), qstrat_(qstrat)
  {
    
    //compute_diagonal();
    compute_rhs();
  }

  template <class IBASIS, unsigned int DIM>
  double 
  BiharmonicEquation<IBASIS,DIM>::D(const typename AggregatedFrame<IBASIS,DIM>::Index& lambda) const
  {
    return ldexp(1.0,lambda.j() * operator_order());
    //    return sqrt(a(lambda,lambda));
    //cout << stiff_diagonal.find(lambda).second << endl;
    //    cout << "1111111111" << endl;
    //return stiff_diagonal.get_coefficient(lambda);
  }

  template <class IBASIS, unsigned int DIM>
  void
  BiharmonicEquation<IBASIS,DIM>::rescale(InfiniteVector<double,
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
  BiharmonicEquation<IBASIS,DIM>::compute_rhs()
  {
    cout << "BiharmonicEquation(): precompute right-hand side..." << endl;

    typedef AggregatedFrame<IBASIS,DIM> Frame;
    typedef typename Frame::Index Index;

    // precompute the right-hand side on a fine level
    InfiniteVector<double,Index> fhelp;
    const int j0   = frame_->j0();
    const int jmax = 7; // for a first quick hack
    for (Index lambda(FrameTL::first_generator<IBASIS,DIM,DIM,Frame>(frame_,j0));; ++lambda)
      {
	const double coeff = f(lambda)/D(lambda);
	if (fabs(coeff)>1e-15) {
	  fhelp.set_coefficient(lambda, coeff);
	  //cout << lambda << endl;
	}
  	if (lambda == last_wavelet<IBASIS,DIM,DIM,Frame>(frame_,jmax))
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
  BiharmonicEquation<IBASIS,DIM>::compute_diagonal()
  {
    cout << "BiharmonicEquation(): precompute diagonal of stiffness matrix..." << endl;

    typedef AggregatedFrame<IBASIS,DIM> Frame;
    typedef typename Frame::Index Index;

    // precompute the right-hand side on a fine level
    const int j0   = frame_->j0();
    const int jmax = 7;  // for a first quick hack
    for (Index lambda(FrameTL::first_generator<IBASIS,DIM,DIM,Frame>(frame_,j0));; ++lambda)
      {
	stiff_diagonal.set_coefficient(lambda, sqrt(a(lambda,lambda)));
	//	cout << "lambda: " <<lambda<< endl;
	if (lambda == last_wavelet<IBASIS,DIM,DIM,Frame>(frame_,jmax))
	  break;
	//cout << " ok ";
      }

    cout << "... done, diagonal of stiffness matrix computed" << endl;
  }



  template <class IBASIS, unsigned int DIM>
  inline
  double
  BiharmonicEquation<IBASIS,DIM>:: integrate(const Index1D<IBASIS>& lambda,
					   const Index1D<IBASIS>& mu,
					   const FixedArray1D<Array1D<double>,DIM >& irregular_grid,
					   const int N_Gauss,
					   const int dir) const
  {
    
    double res = 0;
    
//     typename One_D_IntegralCache::iterator col_lb(one_d_integrals.lower_bound(lambda));
//     typename One_D_IntegralCache::iterator col_it(col_lb);
//     if (col_lb == one_d_integrals.end() ||
// 	one_d_integrals.key_comp()(lambda,col_lb->first))
//       {
// 	// insert a new column
// 	typedef typename One_D_IntegralCache::value_type value_type;
// 	col_it = one_d_integrals.insert(col_lb, value_type(lambda, Column1D()));
//       }
    
//     Column1D& col(col_it->second);
    
//     typename Column1D::iterator lb(col.lower_bound(mu));
//     typename Column1D::iterator it(lb);
//     if (lb == col.end() ||
// 	col.key_comp()(mu, lb->first))
//       {
//	const unsigned int dir = lambda.direction();
	Array1D<double> gauss_points_la, gauss_points_mu, gauss_weights,
	  values_lambda, values_mu;
	
	//	cout << "N Gauss " << N_Gauss << endl;
 	gauss_points_la.resize(N_Gauss*(irregular_grid[dir].size()-1));
 	gauss_points_mu.resize(N_Gauss*(irregular_grid[dir].size()-1));
 	gauss_weights.resize(N_Gauss*(irregular_grid[dir].size()-1));
  	for (unsigned int k = 0; k < irregular_grid[dir].size()-1; k++)
 	  for (int n = 0; n < N_Gauss; n++) {
 	    gauss_points_la[ k*N_Gauss+n  ]
 	      = 0.5 * (irregular_grid[dir][k+1]-irregular_grid[dir][k]) * (GaussPoints[N_Gauss-1][n]+1)
 	      + irregular_grid[dir][k];
	    
	    //	    	    cout << "gp = " << gauss_points_la[ k*N_Gauss+n  ] << endl;

 	    gauss_weights[ k*N_Gauss+n ]
 	      = (irregular_grid[dir][k+1]-irregular_grid[dir][k])*GaussWeights[N_Gauss-1][n];

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

// 	typedef typename Column1D::value_type value_type;
// 	it = col.insert(lb, value_type(mu, res));
//       }
//     else {
//       res = it->second;
//     }
	
    return res;
  }


  template <class IBASIS, unsigned int DIM>
  double
  BiharmonicEquation<IBASIS,DIM>::a_different_patches(const typename AggregatedFrame<IBASIS,DIM>::Index& la,
						      const typename AggregatedFrame<IBASIS,DIM>::Index& nu,
						      const unsigned int q, const unsigned int N) const
  {
    double r = 0.0;
    // cout<<"a begin";
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

    FixedArray1D<Array1D<double>,DIM > irregular_grid;
    
    const int N_Gauss = q;

    bool b = 0;

      b = intersect_supports<IBASIS,DIM,DIM>(*frame_, lambda, mu,
					     supp_lambda, supp_mu, irregular_grid);
      if ( !b )
	return 0.0;

    
#if 0
    for (unsigned int i = 0; i < DIM; i++) {
      cout << "intersect = " << b << endl;
      cout << "dim = " << i << " " << irregular_grid[i] << endl;
    }
#endif

    typedef typename IBASIS::Index Index_1D;
    
    const Chart<DIM>* chart_la = frame_->atlas()->charts()[lambda.p()];
    const Chart<DIM>* chart_mu = frame_->atlas()->charts()[mu.p()];

    FixedArray1D<IBASIS*,DIM> bases1D_lambda = frame_->bases()[lambda.p()]->bases();
    FixedArray1D<IBASIS*,DIM> bases1D_mu     = frame_->bases()[mu.p()]->bases();
    
    FixedArray1D<Array1D<double>,DIM> gauss_points, gauss_points_mu, gauss_weights,
      wav_values_lambda, wav_der_values_lambda, wav_values_mu, wav_der_values_mu;


#if 1
    int dim=(int)DIM;
    int terms=dim*dim;
    int d[2][dim]; 
    //loop over all terms
    for (int i = 0; i <terms; i++) {
	double t = 1.;
	// loop over spatial direction
	for(int l=0;l<dim;l++) {d[0][l]=0;d[1][l]=0;}
	d[0][i/dim]=2;
	d[1][i%dim]=2;
	for (int j = 0; j < (int)DIM; j++) {  
	  Index1D<IBASIS> i1(IntervalIndex<IBASIS> (
						    lambda.j(),lambda.e()[j],lambda.k()[j],
						    bases1D_lambda[j]
						    ),
			     lambda.p(),j,d[0][j]    
			     );
	  
	  Index1D<IBASIS> i2(IntervalIndex<IBASIS> (mu.j(),mu.e()[j],mu.k()[j],
						    bases1D_mu[j]
						    ),
			     mu.p(),j,d[1][j]
			     );

	  t *= integrate(i1, i2, irregular_grid, N_Gauss, j);
	  
 	}



	t *=(1./(chart_la->a_i(i/dim) * chart_mu->a_i(i%dim)))*(1./(chart_la->a_i(i/dim) * chart_mu->a_i(i%dim)));
	
	r += t;
  }

      double tmp1 = 1., tmp2 = 1.;

      for (unsigned int i = 0; i < DIM; i++) {
	tmp1 *= chart_la->a_i(i);
	tmp2 *= chart_mu->a_i(i);
      }
      tmp1 = sqrt(fabs(tmp1)) / sqrt(fabs(tmp2));
      r *= tmp1;
      //cout<<"a end";
#endif      
      return r;



    // dummy return
    return 0.0;
  }







  template <class IBASIS, unsigned int DIM>
  double
  BiharmonicEquation<IBASIS,DIM>::a(const typename AggregatedFrame<IBASIS,DIM>::Index& lambda,
				  const typename AggregatedFrame<IBASIS,DIM>::Index& nu) const
  {
    /* //the cases lambda.p() == nu.p() and lambda.p() != nu.p() have
    //to be treated seperately
    return lambda.p() == nu.p() ? a_same_patches(lambda, nu) :*/
    return a_different_patches(lambda, nu);
    
  }







  template <class IBASIS, unsigned int DIM>
  double
  BiharmonicEquation<IBASIS,DIM>::f(const typename AggregatedFrame<IBASIS,DIM>::Index& lambda) const
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

      double share = bih_bvp_->f(x_patch) * chart->Gram_factor(x);
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

  }

}
