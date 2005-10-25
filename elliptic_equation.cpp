// implementation for elliptic_equation.h

#include <numerics/gauss_data.h>
#include <cube/cube_basis.h>
#include <cube/cube_support.h>

namespace FrameTL
{

  template <class IBASIS, unsigned int DIM>
  EllipticEquation<IBASIS,DIM>::EllipticEquation(const EllipticBVP<DIM>* ell_bvp,
						 const AggregatedFrame<IBASIS,DIM>* frame)
    : ell_bvp_(ell_bvp), frame_(frame)
  {
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

    bool b = WaveletTL::intersect_supports_<IBASIS,DIM,CUBEBASIS>
      (
       *(frame_->bases()[p]), lambda_c, mu_c, supp_intersect
       );
    
//     cout << "b =  " << b << endl;
//     for (unsigned int i = 0; i< DIM; i++)
//       cout << supp_intersect.j << " " << supp_intersect.a[i] << " " << supp_intersect.b[i] << endl;
    
    if (! b)
      return 0.0;
    
    const int N_Gauss = q_order;
    //cout << "N_Gauss = " << N_Gauss << endl;
    const double h = ldexp(1.0, -supp_intersect.j); // granularity for the quadrature
    cout << "h=" << h << endl;
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

    //compute the point values of the the wavelet part of the integrand
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
    }

    int index[DIM]; // current multiindex for the point values
    for (unsigned int i = 0; i < DIM; i++)
      index[i] = 0;

    Point<DIM> x;
    Point<DIM> x_patch;

    //loop over all quadrature knots
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

	//cout << "Dg= " << frame_->atlas()->charts()[p]->Gram_D_factor(s,x) << endl;
	//cout << "sq_g= " << sq_gram << endl;
	//cout << "a= " << ell_bvp_->a(x_patch) << endl;
	//cout << "q= " << ell_bvp_->q(x_patch) << endl;	

	//for first part of the integral: \int_\Omega <a \Nabla u,\Nabla v> dx
	values1[s] = ell_bvp_->a(x_patch) *
	  (t2*sq_gram - t4*frame_->atlas()->charts()[p]->Gram_D_factor(s,x)) / (sq_gram*sq_gram);
	values2[s] =
	  (t3*sq_gram - t5*frame_->atlas()->charts()[p]->Gram_D_factor(s,x)) / (sq_gram*sq_gram);

      }//end for s

      Vector<double> tmp_values1(DIM);
      Vector<double> tmp_values2(DIM);
      
      for (unsigned int i1 = 0; i1 < DIM; i1++)
	for (unsigned int i2 = 0; i2 < DIM; i2++) {
	  //cout << frame_->atlas()->charts()[p]->Dkappa_inv(i2, i1, x) << endl;
	  tmp_values1[i1] += values1[i2]*frame_->atlas()->charts()[p]->Dkappa_inv(i2, i1, x);
	  tmp_values2[i1] += values2[i2]*frame_->atlas()->charts()[p]->Dkappa_inv(i2, i1, x);
	}
      //cout << "#########" << endl;
      //compute innerproduct
      double t = (tmp_values1 * tmp_values2) * (sq_gram*sq_gram) ;
  
      //for second part of the integral: \int_\Omega q u dx
      t += ell_bvp_->q(x_patch) * t4 * t5;
      
      t *= t1;
      //cout << "t1 " << t1 << endl;
      
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
  EllipticEquation<IBASIS,DIM>::a_different_patches(const typename AggregatedFrame<IBASIS,DIM>::Index& lambda,
						    const typename AggregatedFrame<IBASIS,DIM>::Index& nu,
						    const unsigned int p) const
  {
    return 0.0;
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
    

    //    WaveletTL::support<IBASIS,DIM,MAPPEDCUBEBASIS>(*(frame_->bases()[p]), lambda_c, supp);
    WaveletTL::support<IBASIS,DIM,CUBEBASIS>(*(frame_->bases()[p]), lambda_c, supp);
    //cout << supp.j << endl;
//     for (int i = 0; i < DIM; i++)
//       cout << (supp.a)[i] << " " << (supp.b)[i] << endl;

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
//     cout << "f(), Gauss points and weights set up:" << endl;
//     for (unsigned int i = 0; i < DIM; i++)
//       cout << "i=" << i << " with points "
// 	   << gauss_points[i] << " and weights " << gauss_weights[i] << endl;

    // compute the point values of the integrand (where we use that it is a tensor product)
    for (unsigned int i = 0; i < DIM; i++) {
      WaveletTL::evaluate(*(frame_->bases()[p]->bases()[i]), 0,
	       typename IBASIS::Index(lambda.j(),
				      lambda.e()[i],
				      lambda.k()[i],
				      frame_->bases()[p]->bases()[i]),
	       gauss_points[i], v_values[i]);
//       cout << v_values[i] << endl;
//       cout << "#########################" << endl;
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
      //cout << x << endl;
      frame_->atlas()->charts()[p]->map_point(x,x_patch);
      //cout << x_patch << endl;
      double share = ell_bvp_->f(x_patch) * frame_->atlas()->charts()[p]->Gram_factor(x);
      //cout << "Gram = " << frame_->atlas()->charts()[p]->Gram_factor(x) << endl;
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
  
}
