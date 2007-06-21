// implementation for biharmonic_equation.h

#include <list>
#include <numerics/gauss_data.h>
#include <algebra/sparse_matrix.h>
#include <numerics/eigenvalues.h>
#include <frame_index.h>
#include <frame_support.h>
//#include <cuba.h>

//#include <cube/cube_support.h>

using WaveletTL::CubeBasis;

namespace FrameTL
{

  template <class IBASIS, unsigned int DIM>
  SimpleBiharmonicEquation<IBASIS,DIM>::SimpleBiharmonicEquation(const BiharmonicBVP<DIM>* bih_bvp,
						     const AggregatedFrame<IBASIS,DIM>* frame,
						     const int jmax,
						     QuadratureStrategy qstrat)
    : bih_bvp_(bih_bvp), frame_(frame),jmax_(jmax), qstrat_(qstrat)
  {
    compute_diagonal();
    // compute_rhs(false);
    compute_rhs(true);
  }

  template <class IBASIS, unsigned int DIM>
  double 
  SimpleBiharmonicEquation<IBASIS,DIM>::D(const typename AggregatedFrame<IBASIS,DIM>::Index& lambda) const
  {
    return stiff_diagonal[lambda.number()];
  }



  template <class IBASIS, unsigned int DIM>
  void
  SimpleBiharmonicEquation<IBASIS,DIM>::compute_rhs(bool compute)
  {
    cout << "BiharmonicEquation(): precompute right-hand side ..." << endl;
    
    typedef AggregatedFrame<IBASIS,DIM> Frame;
    typedef typename Frame::Index Index;

    // precompute the right-hand side on a fine level
    InfiniteVector<double,Index> fhelp;
    const int j0   = frame_->j0();
    if(compute)
      {    
	for (Index lambda(FrameTL::first_generator<IBASIS,DIM,DIM,Frame>(frame_,j0));; ++lambda)
	  {
	    const double coeff = f(lambda)/D(lambda);
	    //cout << "f(lambda):" << f(lambda) << endl;
	    //cout << lambda << " " << coeff << endl; 
	    if (fabs(coeff)>1e-15) {
	      fhelp.set_coefficient(lambda, coeff);
	    }
	    if (lambda == last_wavelet<IBASIS,DIM,DIM,Frame>(frame_,jmax_))
	      break;
	    //cout << lambda << endl;
	  }
      }
#if 0
    else 
      {
	 Vector<double> f;
	 Rhs_input(f);

	 for(unsigned int i=0; i<f.size(); i++)
	  {
	    Index lambda=FrameIndex<IBASIS, DIM>(i, frame_);
	    fhelp.set_coefficient(lambda, f[i]);
	  }
      }
#endif
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
  SimpleBiharmonicEquation<IBASIS,DIM>::compute_diagonal()
  {
    cout << "SimpleBiharmonicEquation(): precompute diagonal of stiffness matrix ..." << endl;

    stiff_diagonal.resize(frame_->degrees_of_freedom());
    for (int i = 0; i < frame_->degrees_of_freedom(); i++) {
      stiff_diagonal[i] = sqrt(a(*(frame_->get_wavelet(i)),*(frame_->get_wavelet(i))));
      cout << "stiff_diag["<< i << "]:" << stiff_diagonal[i] << endl;
    }

   //  typedef AggregatedFrame<IBASIS,DIM> Frame;
//     typedef typename Frame::Index Index;

//     // precompute the right-hand side on a fine level
//     const int j0   = frame_->j0();

//     for (Index lambda(FrameTL::first_generator<IBASIS,DIM,DIM,Frame>(frame_,j0));; ++lambda)
//       {
// 	stiff_diagonal.set_coefficient(lambda, sqrt(a(lambda,lambda)));
// 	//	cout << "lambda: " <<lambda<< endl;
// 	if (lambda == last_wavelet<IBASIS,DIM,DIM,Frame>(frame_,jmax_))
// 	  break;
	//cout << " ok ";
    //}

    cout << "... done, diagonal of stiffness matrix computed" << endl;
  }



  template <class IBASIS, unsigned int DIM>
  inline
  double
  SimpleBiharmonicEquation<IBASIS,DIM>:: integrate(const Index1D<IBASIS>& lambda,
					   const Index1D<IBASIS>& mu,
					   const FixedArray1D<Array1D<double>,DIM >& irregular_grid,
					   const int N_Gauss,
					   const int dir) const
  {    
    double res = 0;
    Array1D<double> gauss_points_la, gauss_points_mu, gauss_weights,
      values_lambda, values_mu;
    gauss_points_la.resize(N_Gauss*(irregular_grid[dir].size()-1));
    gauss_points_mu.resize(N_Gauss*(irregular_grid[dir].size()-1));
    gauss_weights.resize(N_Gauss*(irregular_grid[dir].size()-1));
    for (unsigned int k = 0; k < irregular_grid[dir].size()-1; k++)
      for (int n = 0; n < N_Gauss; n++) {
	gauss_points_la[ k * N_Gauss+n  ]
	  = 0.5 * (irregular_grid[dir][k+1]-irregular_grid[dir][k]) * (GaussPoints[N_Gauss-1][n]+1)
	  + irregular_grid[dir][k];
	gauss_weights[k*N_Gauss+n ]
	  = (irregular_grid[dir][k+1]-irregular_grid[dir][k])*GaussWeights[N_Gauss-1][n];
      }
    
    IBASIS* basis1D_lambda = frame_->bases()[lambda.p()]->bases()[dir];
    IBASIS* basis1D_mu     = frame_->bases()[mu.p()]->bases()[dir];
    
    const Chart<DIM>* chart_la = frame_->atlas()->charts()[lambda.p()];
    const Chart<DIM>* chart_mu = frame_->atlas()->charts()[mu.p()];
    
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
    WaveletTL::evaluate(*(basis1D_mu), mu.derivative(),
			mu.index(),
			gauss_points_mu, values_mu);
  
    for (unsigned int i = 0; i < values_lambda.size(); i++)
      res += gauss_weights[i] * values_lambda[i] * values_mu[i];

    return res;
  }


  template <class IBASIS, unsigned int DIM>
  double
  SimpleBiharmonicEquation<IBASIS,DIM>::a(const typename AggregatedFrame<IBASIS,DIM>::Index& la,
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
      supp_lambda_ = supp_mu_;
      supp_mu_ = tmp_supp;
    }

    const typename CUBEBASIS::Support* supp_lambda = &supp_lambda_;
    const typename CUBEBASIS::Support* supp_mu = &supp_mu_;

    FixedArray1D<Array1D<double>,DIM > irregular_grid;

    //int q_order;    
    const int N_Gauss = 3;

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
    

    int dim=(int)DIM;
    int terms=dim*dim;
    int d[2][dim]; 
    //loop over all terms
    for (int i = 0; i <terms; i++) {
	double t = 1.;
	// loop over spatial direction
	for(int l=0; l < dim; l++) {
	  d[0][l]=0;d[1][l]=0;
	}
 	d[0][i/dim]=2;
 	d[1][i%dim]=2;

	for (int j = 0; j < (int)DIM; j++) {
	  Index1D<IBASIS> i1(typename IBASIS::Index (
						    lambda.j(),lambda.e()[j],lambda.k()[j],
						    bases1D_lambda[j]
						    ),
			     lambda.p(),j,d[0][j]    
			     );
	  Index1D<IBASIS> i2(typename IBASIS::Index (mu.j(),mu.e()[j],mu.k()[j],
						    bases1D_mu[j]
						    ),
			     mu.p(),j,d[1][j]
			     );

	  t *= integrate(i1, i2, irregular_grid, N_Gauss, j);//cout << integrate(i1, i2, irregular_grid, N_Gauss, j) << endl;
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
    return r;
  }

  template <class IBASIS, unsigned int DIM>
  void
  SimpleBiharmonicEquation<IBASIS,DIM>::RHS
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
  double
  SimpleBiharmonicEquation<IBASIS,DIM>::norm_A() const
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
  SimpleBiharmonicEquation<IBASIS,DIM>::s_star() const
  {
    // cout << "Stern" << endl;
    // notation from [St04a]
    const double t = operator_order();
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
  SimpleBiharmonicEquation<IBASIS,DIM>::f(const typename AggregatedFrame<IBASIS,DIM>::Index& lambda) const
  {

    //\int_\Omega f(Kappa(x)) \psi^\Box (x) \sqrt(\sqrt(det ((D Kappa)^T(x) (D Kappa)(x)) ))
    //recall: \sqrt(\sqrt(...)) = Gram_factor

    int N=50;


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
    const int N_Gauss = 8;
    //const double h = ldexp(1.0, -supp.j); // granularity for the quadrature
    const double h = 1.0 / (1 << supp.j); // granularity for the quadrature
    FixedArray1D<Array1D<double>,DIM> gauss_points, gauss_weights, v_values;
    for (unsigned int i = 0; i < DIM; i++) {
      gauss_points[i].resize(N*N_Gauss*(supp.b[i]-supp.a[i]));
      gauss_weights[i].resize(N*N_Gauss*(supp.b[i]-supp.a[i]));
      for (int patch = supp.a[i]; patch < supp.b[i]; patch++)
	for(int m=0;m<N;m++)
	  for (int n = 0; n < N_Gauss; n++) {
	    gauss_points[i][N*(patch-supp.a[i])*N_Gauss+m*N_Gauss+n]
	      = h*(N*2*patch+1+GaussPoints[N_Gauss-1][n]+2.0*m)/(2.*N);
	    gauss_weights[i][N*(patch-supp.a[i])*N_Gauss+m*N_Gauss+n]
	      = h*GaussWeights[N_Gauss-1][n]/N;
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
	if (index[i] ==N* N_Gauss*(supp.b[i]-supp.a[i])-1) {
	  index[i] = 0;
	  exit = (i == DIM-1);
	} else {
	  index[i]++;
	  break;
	}
      }
      if (exit) break;
    }
#if 0
    return r;
#endif
#if 1
    assert(DIM == 1);
    double tmp = 1;
    //double tmp2 = 1;
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
    /*tmp2 =  WaveletTL::evaluate(*(frame_->bases()[p]->bases()[0]), 2,
			       typename IBASIS::Index(lambda.j(),
						      lambda.e()[0],
						      lambda.k()[0],
						      frame_->bases()[p]->bases()[0]),p2[0]);
    */
    //cout << "r:"<< r<<endl;
    if(tmp !=0.0){
      cout << "tmp:" << tmp << "lambda:" << lambda << endl;
      cout << "r:" << r << endl;
      cout << "D(lambda):" << D(lambda) << endl; 
    }
    tmp /= chart->Gram_factor(p2);
    
    
  
    return  48*tmp + r;
#endif
  }

}
