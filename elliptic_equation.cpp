// implementation for elliptic_equation.h

#include <numerics/gauss_data.h>
#include <frame_index.h>
#include <frame_support.h>
//#include <cube/cube_support.h>

using WaveletTL::CubeBasis;

namespace FrameTL
{

  template <class IBASIS, unsigned int DIM>
  EllipticEquation<IBASIS,DIM>::EllipticEquation(const EllipticBVP<DIM>* ell_bvp,
						 const AggregatedFrame<IBASIS,DIM>* frame,
						 QuadratureStrategy qstrat)
    : ell_bvp_(ell_bvp), frame_(frame), qstrat_(qstrat_)
  {
    compute_rhs();
  }

  template <class IBASIS, unsigned int DIM>
  double 
  EllipticEquation<IBASIS,DIM>::D(const typename AggregatedFrame<IBASIS,DIM>::Index& lambda) const
  {
    return ldexp(1.0,lambda.j());
  }

  template <class IBASIS, unsigned int DIM>
  void
  EllipticEquation<IBASIS,DIM>::rescale(InfiniteVector<double,
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
  EllipticEquation<IBASIS,DIM>::compute_rhs()
  {
    cout << "EllipticEquation(): precompute right-hand side..." << endl;

    typedef AggregatedFrame<IBASIS,DIM> Frame;
    typedef typename Frame::Index Index;

    // precompute the right-hand side on a fine level
    InfiniteVector<double,Index> fhelp;
    const int j0   = frame_->j0();
    const int jmax = 6; // for a first quick hack
    for (Index lambda(FrameTL::first_generator<IBASIS,DIM,DIM,Frame>(frame_,j0));; ++lambda)
      {
	const double coeff = f(lambda)/D(lambda);
	if (fabs(coeff)>1e-15)
	  fhelp.set_coefficient(lambda, coeff);
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
  double
  EllipticEquation<IBASIS,DIM>::a_same_patches(const typename AggregatedFrame<IBASIS,DIM>::Index& lambda,
					       const typename AggregatedFrame<IBASIS,DIM>::Index& mu,
					       const unsigned int q_order) const
  {
    double r = 0;
    
    //patchnumbers are assumed to be equal
    const unsigned int p = lambda.p();
     
    typedef WaveletTL::CubeBasis<IBASIS,DIM> CUBEBASIS;
 
    typename CUBEBASIS::Support supp_intersect;
    
    typename CubeBasis<IBASIS,DIM>::Index lambda_c(lambda.j(),
						   lambda.e(),
						   lambda.k(),frame_->bases()[p]);

    typename CubeBasis<IBASIS,DIM>::Index mu_c(mu.j(),
					       mu.e(),
					       mu.k(),frame_->bases()[p]);

    bool b = WaveletTL::intersect_supports<IBASIS,DIM>
      (
       *(frame_->bases()[p]), lambda_c, mu_c, supp_intersect
       );
    
    //    cout << "b =  " << b << endl;
//      for (unsigned int i = 0; i< DIM; i++)
//        cout << supp_intersect.j << " " << supp_intersect.a[i] << " " << supp_intersect.b[i] << endl;
    
    if (! b)
      return 0.0;
    
    const int N_Gauss = q_order;
    //cout << "N_Gauss = " << N_Gauss << endl;
    const double h = ldexp(1.0, -supp_intersect.j); // granularity for the quadrature
    //cout << "h=" << h << endl;
    FixedArray1D<Array1D<double>,DIM> gauss_points, gauss_weights,
      wav_values_lambda, wav_der_values_lambda, wav_values_mu, wav_der_values_mu;


    for (unsigned int i = 0; i < DIM; i++) {
      gauss_points[i].resize(N_Gauss*(supp_intersect.b[i]-supp_intersect.a[i]));
      gauss_weights[i].resize(N_Gauss*(supp_intersect.b[i]-supp_intersect.a[i]));
 	for (int patch = supp_intersect.a[i]; patch < supp_intersect.b[i]; patch++)
 	  for (int n = 0; n < N_Gauss; n++) {
 	    gauss_points[i][(patch-supp_intersect.a[i])*N_Gauss+n]
	      = h*(2*patch+1+GaussPoints[N_Gauss-1][n])/2.;
 	    gauss_weights[i][(patch-supp_intersect.a[i])*N_Gauss+n]
 	      = h*GaussWeights[N_Gauss-1][n];
 	  }
    }

    // compute the point values of the wavelet part of the integrand
    for (unsigned int i = 0; i < DIM; i++) {
      WaveletTL::evaluate(*(frame_->bases()[p]->bases()[i]), 0,
			  typename IBASIS::Index(lambda.j(),
						 lambda.e()[i],
						 lambda.k()[i],
						 frame_->bases()[p]->bases()[i]),
			  gauss_points[i], wav_values_lambda[i]);

      WaveletTL::evaluate(*(frame_->bases()[p]->bases()[i]), 1,
			  typename IBASIS::Index(lambda.j(),
						 lambda.e()[i],
						 lambda.k()[i],
						 frame_->bases()[p]->bases()[i]),
			  gauss_points[i], wav_der_values_lambda[i]);

      WaveletTL::evaluate(*(frame_->bases()[p]->bases()[i]), 0,
			  typename IBASIS::Index(mu.j(),
						 mu.e()[i],
						 mu.k()[i],
						 frame_->bases()[p]->bases()[i]),
			  gauss_points[i], wav_values_mu[i]);

      WaveletTL::evaluate(*(frame_->bases()[p]->bases()[i]), 1,
			  typename IBASIS::Index(mu.j(),
						 mu.e()[i],
						 mu.k()[i],
						 frame_->bases()[p]->bases()[i]),
			  gauss_points[i], wav_der_values_mu[i]);
      //cout << wav_der_values_mu[i] << endl;
    }

    int index[DIM]; // current multiindex for the point values
    for (unsigned int i = 0; i < DIM; i++)
      index[i] = 0;

    Point<DIM> x;
    Point<DIM> x_patch;

    // loop over all quadrature knots
    while (true) {
      for (unsigned int i = 0; i < DIM; i++)
	x[i] = gauss_points[i][index[i]];

      frame_->atlas()->charts()[p]->map_point(x,x_patch);
      double sq_gram = frame_->atlas()->charts()[p]->Gram_factor(x);

      Vector<double> values1(DIM);
      Vector<double> values2(DIM);
      double t1=1., t4=1., t5=1.;

      for (unsigned int i = 0; i < DIM; i++) {
	t1 *= gauss_weights[i][index[i]];
	t4 *= wav_values_lambda[i][index[i]];
	t5 *= wav_values_mu[i][index[i]];
      }

      for (unsigned int s = 0; s < DIM; s++) {
	double t2=1., t3=1.;
	for (unsigned int i = 0; i < DIM; i++) {
	  if (i != s) {
	    t2 *= wav_values_lambda[i][index[i]];
	    t3 *= wav_values_mu[i][index[i]];
	  }
	}
	t2 *= wav_der_values_lambda[s][index[s]];
	t3 *= wav_der_values_mu[s][index[s]];

	// for first part of the integral: \int_\Omega <a \Nabla u,\Nabla v> dx
	values1[s] = ell_bvp_->a(x_patch) *
	  (t2*sq_gram - (t4*frame_->atlas()->charts()[p]->Gram_D_factor(s,x))) / (sq_gram*sq_gram);
	values2[s] =
	  (t3*sq_gram - (t5*frame_->atlas()->charts()[p]->Gram_D_factor(s,x))) / (sq_gram*sq_gram);

      }//end for s

      Vector<double> tmp_values1(DIM);
      Vector<double> tmp_values2(DIM);

      for (unsigned int i1 = 0; i1 < DIM; i1++)
	for (unsigned int i2 = 0; i2 < DIM; i2++) {
	  //cout << frame_->atlas()->charts()[p]->Dkappa_inv(i2, i1, x) << endl;
	  tmp_values1[i1] += values1[i2]*frame_->atlas()->charts()[p]->Dkappa_inv(i2, i1, x);
	  tmp_values2[i1] += values2[i2]*frame_->atlas()->charts()[p]->Dkappa_inv(i2, i1, x);
	}

      // compute innerproduct
      double t = (tmp_values1 * tmp_values2) * (sq_gram*sq_gram) ;
  
      // for second part of the integral: \int_\Omega q u dx
      t += ell_bvp_->q(x_patch) * t4 * t5;
      
      t *= t1;
      
      r += t;

      // "++index"
      bool exit = false;
      for (unsigned int i = 0; i < DIM; i++) {
	if (index[i] == N_Gauss*(supp_intersect.b[i]-supp_intersect.a[i])-1) {
	  index[i] = 0;
	  exit = (i == DIM-1);
	} else {
	  index[i]++;
	  break;
	}
      }
      if (exit) break;
    }//end while

    return r;
  }

  template <class IBASIS, unsigned int DIM>
  double
  EllipticEquation<IBASIS,DIM>::a_different_patches(const typename AggregatedFrame<IBASIS,DIM>::Index& la,
						    const typename AggregatedFrame<IBASIS,DIM>::Index& nu,
						    const unsigned int q) const
  {

    double r = 0.0;

    typename AggregatedFrame<IBASIS,DIM>::Index lambda = la;
    typename AggregatedFrame<IBASIS,DIM>::Index mu = nu;
    typename AggregatedFrame<IBASIS,DIM>::Index tmp_ind;

    typedef WaveletTL::CubeBasis<IBASIS,DIM> CUBEBASIS;
 
    typename CUBEBASIS::Support supp_lambda;
    typename CUBEBASIS::Support supp_mu;

    const int N = 2;

    bool b = intersect_supports<IBASIS,DIM,DIM>(*frame_, lambda, mu, supp_lambda, supp_mu);
    if ( !b )
      return 0.0;

    const int N_Gauss = q;
    
    const double h = ldexp(1.0, -std::max(supp_lambda.j,supp_mu.j));// granularity for the quadrature

    typename CUBEBASIS::Support tmp_supp;

    // swap indices and supports if necessary
    if (supp_mu.j > supp_lambda.j) {
      tmp_ind = lambda;
      lambda = mu;
      mu = tmp_ind;

//       for (unsigned int i = 0; i < DIM; i++)
// 	cout << supp_lambda.a[i] << " " << supp_lambda.b[i]
// 	     << " " << supp_lambda.j << endl;
      tmp_supp = supp_lambda;
      supp_lambda = supp_mu;
      supp_mu = tmp_supp;
//       for (unsigned int i = 0; i < DIM; i++)
// 	cout << supp_lambda.a[i] << " " << supp_lambda.b[i]
// 	     << " " << supp_lambda.j << endl;

    }

    FixedArray1D<Array1D<double>,DIM> gauss_points, gauss_weights, gauss_points_other,
      wav_values_lambda, wav_der_values_lambda, wav_values_mu, wav_der_values_mu;


    for (unsigned int i = 0; i < DIM; i++) {
      gauss_points[i].resize(N * N_Gauss*(supp_lambda.b[i]-supp_lambda.a[i]));
      gauss_weights[i].resize(N * N_Gauss*(supp_lambda.b[i]-supp_lambda.a[i]));

      gauss_points_other[i].resize(N * N_Gauss*(supp_lambda.b[i]-supp_lambda.a[i]));

 	for (int patch = supp_lambda.a[i]; patch < supp_lambda.b[i]; patch++)
 	  for (int m = 0; m < N; m++)
	    for (int n = 0; n < N_Gauss; n++) {
	      //cout << "patch = " << patch << " N = " << N << " m = " << m << " n = " << n << endl; 
	      //cout << N*(patch-supp_lambda.a[i])*N_Gauss + m*N_Gauss+n << endl;
 	      gauss_points[i][ N*(patch-supp_lambda.a[i])*N_Gauss + m*N_Gauss+n ]
		= h*( 1.0/(2.*N)*(GaussPoints[N_Gauss-1][n]+1+2.0*m)+patch);
	      gauss_weights[i][ N*(patch-supp_lambda.a[i])*N_Gauss + m*N_Gauss+n ]
		= (h*GaussWeights[N_Gauss-1][n])/N;
 	  }
    }

    // compute the point values of the 'first' wavelet
    for (unsigned int i = 0; i < DIM; i++) {
      WaveletTL::evaluate(*(frame_->bases()[lambda.p()]->bases()[i]), 0,
			  typename IBASIS::Index(lambda.j(),
						 lambda.e()[i],
						 lambda.k()[i],
						 frame_->bases()[lambda.p()]->bases()[i]),
			  gauss_points[i], wav_values_lambda[i]);
      
      WaveletTL::evaluate(*(frame_->bases()[lambda.p()]->bases()[i]), 1,
			  typename IBASIS::Index(lambda.j(),
						 lambda.e()[i],
						 lambda.k()[i],
						 frame_->bases()[lambda.p()]->bases()[i]),
			  gauss_points[i], wav_der_values_lambda[i]);
    }
    //#####################################################
    // now get those gauss point that are relevant for
    // the 'coarser' wavelet
    //#####################################################
    int index[DIM]; // current multiindex for the point values
    for (unsigned int i = 0; i < DIM; i++)
      index[i] = 0;

    Point<DIM> x;
    Point<DIM> x_patch;
    Point<DIM> y;

    // loop over all quadrature knots
    while (true) {      
       for (unsigned int i = 0; i < DIM; i++)
	x[i] = gauss_points[i][index[i]];

       frame_->atlas()->charts()[lambda.p()]->map_point(x,x_patch);
       if ( frame_->atlas()->charts()[mu.p()]->in_patch(x_patch) ) {
	 frame_->atlas()->charts()[mu.p()]->map_point_inv(x_patch,y);
	 for (unsigned int i = 0; i < DIM; i++)
	   gauss_points_other[i][index[i]] = y[i];
       }
       else {
	 for (unsigned int i = 0; i < DIM; i++) 
	   gauss_points_other[i][index[i]] = -1.;
       }	 
      // "++index"
      bool exit = false;
      for (unsigned int i = 0; i < DIM; i++) {
	if (index[i] == N * N_Gauss * (supp_lambda.b[i]-supp_lambda.a[i]) - 1) {
	  index[i] = 0;
	  exit = (i == DIM-1);
	} else {
	  index[i]++;
	  break;
	}
      }
      if (exit) break;
    }
    //#####################################################

    // compute the point values of the 'second' wavelet
    for (unsigned int i = 0; i < DIM; i++) {
      WaveletTL::evaluate(*(frame_->bases()[mu.p()]->bases()[i]), 0,
			  typename IBASIS::Index(mu.j(),
						 mu.e()[i],
						 mu.k()[i],
						 frame_->bases()[mu.p()]->bases()[i]),
			  gauss_points_other[i], wav_values_mu[i]);
      
      WaveletTL::evaluate(*(frame_->bases()[mu.p()]->bases()[i]), 1,
			  typename IBASIS::Index(mu.j(),
						 mu.e()[i],
						 mu.k()[i],
						 frame_->bases()[mu.p()]->bases()[i]),
			  gauss_points_other[i], wav_der_values_mu[i]);
    }


    //------------------------------------------------------
    for (unsigned int i = 0; i < DIM; i++)
      index[i] = 0;

    // loop over all quadrature knots
    while (true) {
      for (unsigned int i = 0; i < DIM; i++) {
	x[i] = gauss_points[i][index[i]];
	y[i] = gauss_points_other[i][index[i]]; 
      }

      frame_->atlas()->charts()[lambda.p()]->map_point(x,x_patch);
      double sq_gram_la = frame_->atlas()->charts()[lambda.p()]->Gram_factor(x);
      double sq_gram_mu = frame_->atlas()->charts()[mu.p()]->Gram_factor(y);

      Vector<double> values1(DIM);
      Vector<double> values2(DIM);
      double t1=1., t4=1., t5=1.;

      for (unsigned int i = 0; i < DIM; i++) {
	t1 *= gauss_weights[i][index[i]];
	t4 *= wav_values_lambda[i][index[i]];
	t5 *= wav_values_mu[i][index[i]];
      }

      for (unsigned int s = 0; s < DIM; s++) {
	double t2=1., t3=1.;
	for (unsigned int i = 0; i < DIM; i++) {
	  if (i != s) {
	    t2 *= wav_values_lambda[i][index[i]];
	    t3 *= wav_values_mu[i][index[i]];
	  }
	}
	t2 *= wav_der_values_lambda[s][index[s]];
	t3 *= wav_der_values_mu[s][index[s]];

	// for first part of the integral: \int_\Omega <a \Nabla u,\Nabla v> dx
	values1[s] = ell_bvp_->a(x_patch) *
	  (t2*sq_gram_la - (t4*frame_->atlas()->charts()[lambda.p()]->Gram_D_factor(s,x))) / (sq_gram_la*sq_gram_la);
	values2[s] =
	  (t3*sq_gram_mu - (t5*frame_->atlas()->charts()[mu.p()]->Gram_D_factor(s,y))) / (sq_gram_mu*sq_gram_mu);

      }//end for s

      Vector<double> tmp_values1(DIM);
      Vector<double> tmp_values2(DIM);

      for (unsigned int i1 = 0; i1 < DIM; i1++)
	for (unsigned int i2 = 0; i2 < DIM; i2++) {
	  //cout << frame_->atlas()->charts()[p]->Dkappa_inv(i2, i1, x) << endl;
	  tmp_values1[i1] += values1[i2]*frame_->atlas()->charts()[lambda.p()]->Dkappa_inv(i2, i1, x);
	  tmp_values2[i1] += values2[i2]*frame_->atlas()->charts()[mu.p()]->Dkappa_inv(i2, i1, y);
	}

      // compute innerproduct
      double t = (tmp_values1 * tmp_values2) * (sq_gram_la*sq_gram_la) ;
  
      // for second part of the integral: \int_\Omega q u dx
      t += ell_bvp_->q(x_patch) * t4 * t5;
      
      t *= t1;
      
      r += t;

      // "++index"
      bool exit = false;
      for (unsigned int i = 0; i < DIM; i++) {
	if (index[i] == N * N_Gauss * (supp_lambda.b[i]-supp_lambda.a[i]) - 1) {
	  index[i] = 0;
	  exit = (i == DIM-1);
	} else {
	  index[i]++;
	  break;
	}
      }
      if (exit) break;
    }//end while

    return r;
  }


  template <class IBASIS, unsigned int DIM>
  double
  EllipticEquation<IBASIS,DIM>::a(const typename AggregatedFrame<IBASIS,DIM>::Index& lambda,
				  const typename AggregatedFrame<IBASIS,DIM>::Index& nu,
				  const unsigned int p) const
  {
    
    //the cases lambda.p() == nu.p() and lambda.p() != nu.p() have
    //to be treated seperately
    return lambda.p() == nu.p() ? a_same_patches(lambda, nu, p) : a_different_patches(lambda, nu, p);
    
  }

  template <class IBASIS, unsigned int DIM>
  double
  EllipticEquation<IBASIS,DIM>::norm_A() const
  {
    if (normA == 0.0) {
      typedef typename AggregatedFrame<IBASIS,DIM>::Index Index;
      std::set<Index> Lambda;
      const int j0 = frame->j0();
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
  EllipticEquation<IBASIS,DIM>::s_star() const
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
	    : std::min((t+dT)/(double)n, (gamma-t)/(n-1.))); // [St04a, Th. 2.3]
  }

  template <class IBASIS, unsigned int DIM>
  double
  EllipticEquation<IBASIS,DIM>::f(const typename AggregatedFrame<IBASIS,DIM>::Index& lambda) const
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

    // setup Gauss points and weights for a composite quadrature formula:
    const int N_Gauss = 1;
    const double h = ldexp(1.0, -supp.j); // granularity for the quadrature
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

      frame_->atlas()->charts()[p]->map_point(x,x_patch);

      double share = ell_bvp_->f(x_patch) * frame_->atlas()->charts()[p]->Gram_factor(x);

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
    
    return r;
  }

  template <class IBASIS, unsigned int DIM>
  void
  EllipticEquation<IBASIS,DIM>::RHS
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
  EllipticEquation<IBASIS,DIM>::set_bvp(const EllipticBVP<DIM>* bvp)
  {
    bvp_ = bvp;
    compute_rhs();
  }

  
}
