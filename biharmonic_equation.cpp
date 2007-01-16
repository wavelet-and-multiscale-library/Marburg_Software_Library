// implementation for elliptic_equation.h

#include <list>
#include <numerics/gauss_data.h>
#include <frame_index.h>
#include <frame_support.h>
//#include <cuba.h>

//#include <cube/cube_support.h>

using WaveletTL::CubeBasis;

namespace FrameTL
{

//   template <class IBASIS>
//   Index1D<IBASIS>::Index1D(const IntervalIndex<IBASIS>& ind,
// 			   const unsigned int p, const unsigned int dir,
// 			   const unsigned int der)
//     : ind_(ind), p_(p), dir_(dir), der_(der)
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
// 	 p_ == lambda.p() &&
// 	 (
// 	  dir_ < lambda.direction() ||
// 	  (
// 	   dir_ == lambda.direction() && der_ < lambda.derivative()   
// 	   )
// 	  )
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
  BiharmonicEquation<IBASIS,DIM>::BiharmonicEquation(const BiharmonicBVP<DIM>* bih_bvp,
						     const AggregatedFrame<IBASIS,DIM>* frame,
						     QuadratureStrategy qstrat)
    : bih_bvp_(bih_bvp), frame_(frame), qstrat_(qstrat)
  {
    
    compute_diagonal();
    compute_rhs();
  }

  template <class IBASIS, unsigned int DIM>
  double 
  BiharmonicEquation<IBASIS,DIM>::D(const typename AggregatedFrame<IBASIS,DIM>::Index& lambda) const
  {
    //return ldexp(1.0,lambda.j());
    //return 1 << lambda.j();
    
    //    return sqrt(a(lambda,lambda));
    //cout << stiff_diagonal.find(lambda).second << endl;
    //    cout << "1111111111" << endl;
    return stiff_diagonal.get_coefficient(lambda);
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
    const int jmax = 6; // for a first quick hack
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
    const int jmax = 6;  // for a first quick hack
    for (Index lambda(FrameTL::first_generator<IBASIS,DIM,DIM,Frame>(frame_,j0));; ++lambda)
      {
	stiff_diagonal.set_coefficient(lambda, sqrt(a(lambda,lambda)));
	if (lambda == last_wavelet<IBASIS,DIM,DIM,Frame>(frame_,jmax))
	  break;
      }

    cout << "... done, digonal of stiffness matrix computed" << endl;
  }

// //   template <class IBASIS, unsigned int DIM>
// //   double
// //   BiharmonicEquation<IBASIS,DIM>::a_same_patches(const typename AggregatedFrame<IBASIS,DIM>::Index& lambda,
// // 					       const typename AggregatedFrame<IBASIS,DIM>::Index& mu,
// // 					       const unsigned int q_order) const
// //   {
// //     double r = 0;
    
// //     //patchnumbers are assumed to be equal
// //     const unsigned int p = lambda.p();
     
// //     typedef WaveletTL::CubeBasis<IBASIS,DIM> CUBEBASIS;
// //     typedef typename CUBEBASIS::Index CubeIndex;

// //     typename CUBEBASIS::Support supp_intersect;
 
// //     bool b = WaveletTL::intersect_supports<IBASIS,DIM>
// //       (
// //        *(frame_->bases()[p]), 
// //        CubeIndex(lambda.j(), lambda.e(), lambda.k(), frame_->bases()[p]),
// //        CubeIndex(mu.j(), mu.e(), mu.k(), frame_->bases()[p]),
// //        supp_intersect
// //        );
    
// //     if (! b)
// //       return 0.0;

// //     CUBEBASIS* local_cube_basis = frame_->bases()[p];
// //     const Chart<DIM>* chart = frame_->atlas()->charts()[p];

// //     const int N_Gauss = q_order;
// //     //cout << "N_Gauss = " << N_Gauss << endl;
// //     //const double h = ldexp(1.0, -supp_intersect.j); // granularity for the quadrature
// //     const double h = 1.0 / (1 << supp_intersect.j);
// //     //cout << "h=" << h << endl;
// //     FixedArray1D<Array1D<double>,DIM> gauss_points, gauss_weights,
// //       wav_values_lambda, wav_der_values_lambda, wav_values_mu, wav_der_values_mu;

 

// //     for (unsigned int i = 0; i < DIM; i++) {
// //       gauss_points[i].resize(N_Gauss*(supp_intersect.b[i]-supp_intersect.a[i]));
// //       gauss_weights[i].resize(N_Gauss*(supp_intersect.b[i]-supp_intersect.a[i]));
// //       for (int patch = supp_intersect.a[i]; patch < supp_intersect.b[i]; patch++)
// // 	for (int n = 0; n < N_Gauss; n++) {
// // 	  gauss_points[i][(patch-supp_intersect.a[i])*N_Gauss+n]
// // 	    = h*(2*patch+1+GaussPoints[N_Gauss-1][n])/2.;
// // 	  gauss_weights[i][(patch-supp_intersect.a[i])*N_Gauss+n]
// // 	    = h*GaussWeights[N_Gauss-1][n];
// // 	}
// //     }

// //     // compute the point values of the wavelet part of the integrand
// //     for (unsigned int i = 0; i < DIM; i++) {
// //       WaveletTL::evaluate(*(local_cube_basis->bases()[i]), 0,
// // 			  typename IBASIS::Index(lambda.j(),
// // 						 lambda.e()[i],
// // 						 lambda.k()[i],
// // 						 local_cube_basis->bases()[i]),
// // 			  gauss_points[i], wav_values_lambda[i]);
// //       WaveletTL::evaluate(*(local_cube_basis->bases()[i]), 1,
// // 			  typename IBASIS::Index(lambda.j(),
// // 						 lambda.e()[i],
// // 						 lambda.k()[i],
// // 						 local_cube_basis->bases()[i]),
// // 			  gauss_points[i], wav_der_values_lambda[i]);

// //       WaveletTL::evaluate(*(local_cube_basis->bases()[i]), 0,
// // 			  typename IBASIS::Index(mu.j(),
// // 						 mu.e()[i],
// // 						 mu.k()[i],
// // 						 local_cube_basis->bases()[i]),
// // 			  gauss_points[i], wav_values_mu[i]);

// //       WaveletTL::evaluate(*(local_cube_basis->bases()[i]), 1,
// // 			  typename IBASIS::Index(mu.j(),
// // 						 mu.e()[i],
// // 						 mu.k()[i],
// // 						 local_cube_basis->bases()[i]),
// // 			  gauss_points[i], wav_der_values_mu[i]);
// //       //cout << wav_der_values_mu[i] << endl;
// //     }

// //     int index[DIM]; // current multiindex for the point values
// //     for (unsigned int i = 0; i < DIM; i++)
// //       index[i] = 0;

// //     Point<DIM> x;
// //     Point<DIM> x_patch;

// //     // loop over all quadrature knots
// //     while (true) {
// //       for (unsigned int i = 0; i < DIM; i++)
// // 	x[i] = gauss_points[i][index[i]];

// //       chart->map_point(x,x_patch);
// //       double sq_gram = chart->Gram_factor(x);

// //       Vector<double> values1(DIM);
// //       Vector<double> values2(DIM);
// //       double weights=1., psi_lambda=1., psi_mu=1.;

// //       for (unsigned int i = 0; i < DIM; i++) {
// // 	weights *= gauss_weights[i][index[i]];
// // 	psi_lambda *= wav_values_lambda[i][index[i]];
// // 	psi_mu *= wav_values_mu[i][index[i]];
// //       }
      
// //       if ( !(psi_lambda == 0. || psi_mu == 0.) ) {

// // 	for (unsigned int s = 0; s < DIM; s++) {
// // 	  double psi_der_lambda=1., psi_der_mu=1.;
// // 	  psi_der_lambda = (psi_lambda / wav_values_lambda[s][index[s]]) * wav_der_values_lambda[s][index[s]];
// // 	  psi_der_mu = (psi_mu / wav_values_mu[s][index[s]]) * wav_der_values_mu[s][index[s]];

// // 	  // for first part of the integral: \int_\Omega <a \Nabla u,\Nabla v> dx
// // 	  double tmp = chart->Gram_D_factor(s,x);
// // 	  values1[s] = ell_bvp_->a(x_patch) *
// // 	    (psi_der_lambda*sq_gram - (psi_lambda*tmp)) / (sq_gram*sq_gram);
// // 	  values2[s] =
// // 	    (psi_der_mu*sq_gram - (psi_mu*tmp)) / (sq_gram*sq_gram);

// // 	}//end for s
	
// // 	double t = 0.;
// // 	for (unsigned int i1 = 0; i1 < DIM; i1++) {
// // 	  double d1 = 0.;
// // 	  double d2 = 0.;
// // 	  for (unsigned int i2 = 0; i2 < DIM; i2++) {
// // 	    double tmp = chart->Dkappa_inv(i2, i1, x);
// // 	    d1 += values1[i2]*tmp;
// // 	    d2 += values2[i2]*tmp;
// // 	  }
// // 	  t += d1 * d2;

// // 	}
// // 	r += (t * (sq_gram*sq_gram) + ell_bvp_->q(x_patch) * psi_lambda * psi_mu)
// // 	  * weights;
	
// //       }
// //       // "++index"
// //       bool exit = false;
// //       for (unsigned int i = 0; i < DIM; i++) {
// // 	if (index[i] == N_Gauss*(supp_intersect.b[i]-supp_intersect.a[i])-1) {
// // 	  index[i] = 0;
// // 	  exit = (i == DIM-1);
// // 	} else {
// // 	  index[i]++;
// // 	  break;
// // 	}
// //       }
// //       if (exit) break;
	
// //     }//end while

// //     return r;
// //   }

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
  
    Index lambda = la;
    Index mu = nu;
    Index tmp_ind;

    typedef WaveletTL::CubeBasis<IBASIS,DIM> CUBEBASIS;

    typedef typename CUBEBASIS::Index CubeIndex;
 
    typename CUBEBASIS::Support supp_lambda_ = frame_->all_supports[lambda.number()];
    typename CUBEBASIS::Support supp_mu_ = frame_->all_supports[mu.number()];

    //     WaveletTL::support<IBASIS,DIM>(*frame_->bases()[lambda.p()], 
    // 				   CubeIndex(lambda.j(),
    // 					     lambda.e(),
    // 					     lambda.k(),
    // 					     frame_->bases()[lambda.p()]),
    // 				   supp_lambda);
    //     WaveletTL::support<IBASIS,DIM>(*frame_->bases()[mu.p()],
    // 				   CubeIndex(mu.j(),
    // 					     mu.e(),
    // 					     mu.k(),frame_->bases()[mu.p()]),
    // 				   supp_mu);
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


    
#if 0
    cout << "##########################" << endl;
    cout << "lambda = " << lambda << endl;
    for (unsigned int i = 0; i < 2; i++) {
      double dx1 = 1./(1 << supp_lambda->j);
      cout << "j = " << supp_lambda->j << " p = " << lambda.p() << endl;
      cout << "a = " << supp_lambda->a[i] <<  " b= " << supp_lambda->b[i] << " p = " << lambda.p() << endl;
      //cout << "a_lambda = " << supp_lambda->a[i]*dx1 <<  " b= " << supp_lambda->b[i]*dx1 << endl;
    }
    cout << "mu = " << mu << endl;
    for (unsigned int i = 0; i < 2; i++) {
      double dx2 = 1./(1 << supp_mu->j);
      cout << "j = " << supp_mu->j << " p = " << mu.p() << endl;
      cout << "a = " << supp_mu->a[i] <<  " b= " << supp_mu->b[i] << " p = " << mu.p() << endl;
      //cout << "a_mu = " << supp_mu->a[i]*dx2 <<  " b= " << supp_mu->b[i]*dx2 << endl;
    }

#endif

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
    int anzSum=dim*dim;
    int d[2][dim]; //Ordnung der Ableitung
    for (int i = 0; i < (int) anzSum; i++) {
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


#endif      
      return r;



    // dummy return
    return 0.0;
  }




















// #define NDIM 2
// #define NCOMP 1
// #define EPSREL 1e-8
// #define EPSABS 1e-15
// #define VERBOSE 3
// #define LAST 4
// #define MINEVAL 0
// #define MAXEVAL 400000

// #define NSTART 1000
// #define NINCREASE 500

// #define NNEW 1000
// #define FLATNESS 25.

// #define KEY1 47
// #define KEY2 1
// #define KEY3 1
// #define MAXPASS 5
// #define BORDER 0.
// #define MAXCHISQ 10.
// #define MINDEVIATION .25
// #define NGIVEN 0
// #define LDXGIVEN NDIM
// #define NEXTRA 0

// #define KEY 0

//   //typedef WaveletTL::DSBasis<2,2> I_BASIS;
//   typedef WaveletTL::PBasis<2,2> I_BASIS;

//   typedef WaveletTL::CubeBasis<I_BASIS,NDIM> CUBE_BASIS;

//   typedef I_BASIS::Index Index_1D;
  
//   typedef CUBE_BASIS::Index Cube_Index;
  
//   static CUBE_BASIS::Support* supp_lambda_;
//   static CUBE_BASIS::Support* supp_mu_;

//   static Chart<NDIM>* chart_la;
//   static Chart<NDIM>* chart_mu;

//   static FrameIndex<I_BASIS,NDIM>* lambda_;
//   static FrameIndex<I_BASIS,NDIM>* mu_;

//   static const AggregatedFrame<I_BASIS,NDIM>* frame_2;

//   static const FixedArray1D<I_BASIS*,NDIM>* bases1D_lambda;
//   static const FixedArray1D<I_BASIS*,NDIM>* bases1D_mu;

//   static const EllipticBVP<NDIM>* ell_bvp;

//   static double granularity;

//   template <class IBASIS, unsigned int DIM>
//   void
//   integration_setup(const AggregatedFrame<IBASIS,DIM>* frame,
// 		    const FrameIndex<IBASIS,DIM>* la, const FrameIndex<IBASIS,DIM>* nu,
// 		    const CUBE_BASIS::Support* supp_la, const CUBE_BASIS::Support* supp_nu,
// 		    const FixedArray1D<IBASIS*,DIM>* bases1D_la, const FixedArray1D<IBASIS*,DIM>* bases1D_nu)
//   {
    
//   }



 
//   template <class IBASIS, unsigned int DIM>
//   double
//   EllipticEquation<IBASIS,DIM>::a_different_patches_adaptive(const typename AggregatedFrame<IBASIS,DIM>::Index& la,
// 							     const typename AggregatedFrame<IBASIS,DIM>::Index& nu) const
//   {
//     int comp, nregions, neval, fail;
//     double integral[NCOMP], error[NCOMP], prob[NCOMP];
//     int verbose = 0;

//     double r = 0.0;
  
//     Index lambda = la;
//     Index mu = nu;
//     Index tmp_ind;

//     typedef WaveletTL::CubeBasis<IBASIS,DIM> CUBEBASIS;

//     typedef typename CUBEBASIS::Index CubeIndex;
 
//     typename CUBEBASIS::Support supp_lambda;
//     typename CUBEBASIS::Support supp_mu;

//     WaveletTL::support<IBASIS,DIM>(*frame_->bases()[lambda.p()], 
// 				   CubeIndex(lambda.j(),
// 					     lambda.e(),
// 					     lambda.k(),
// 					     frame_->bases()[lambda.p()]),
// 				   supp_lambda);
//     WaveletTL::support<IBASIS,DIM>(*frame_->bases()[mu.p()],
// 				   CubeIndex(mu.j(),
// 					     mu.e(),
// 					     mu.k(),frame_->bases()[mu.p()]),
// 				   supp_mu);
    

//     typename CUBEBASIS::Support tmp_supp;

//     // swap indices and supports if necessary
//     if (supp_mu.j > supp_lambda.j) {
//       tmp_ind = lambda;
//       lambda = mu;
//       mu = tmp_ind;

//       tmp_supp = supp_lambda;
//       supp_lambda = supp_mu;
//       supp_mu = tmp_supp;
//     }

//     bool b = 0;

//     b = intersect_supports<IBASIS,DIM,DIM>(*frame_, lambda, mu, supp_lambda, supp_mu);
    
//     if (! b)
//       return 0;

//     typedef typename IBASIS::Index Index_1D;
      
//     const double h = 1.0 / (1 << std::max(supp_lambda.j,supp_mu.j));// granularity for the quadrature

//     //    integration_setup<IBASIS,DIM>(frame_, lambda, mu, supp_lambda, supp_mu, bases1D_lambda, bases1D_mu);    

//     Index ind1(lambda);
//     Index ind2(mu);

//     frame_2 = &frame();
//     chart_la = frame_->atlas()->charts()[ind1.p()];
//     chart_mu = frame_->atlas()->charts()[ind2.p()];
//     lambda_ = &ind1;
//     mu_ = &ind2;
//     supp_lambda_ = &supp_lambda;
//     supp_mu_ = &supp_mu;
//     bases1D_lambda = &(frame_->bases()[ind1.p()]->bases());
//     bases1D_mu     = &(frame_->bases()[ind2.p()]->bases());
//     ell_bvp = &get_bvp();
//     granularity = h;
    
//     return r;

//   }

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


  
 }
