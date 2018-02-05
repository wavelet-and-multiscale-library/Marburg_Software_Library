// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of FrameTL - the Frame Template Library          |
// |                                                                    |
// | Copyright (c) 2016-2019                                            |
// | Christoph Hartmann                                                 |
// +--------------------------------------------------------------------+

// implementation for pPoisson_equation.h

#include <list>
#include <numerics/gauss_data.h>
#include <frame_index.h>
#include <frame_support.h>
#include <index1D.h>

//#include <cube/cube_support.h>

using WaveletTL::CubeBasis;

namespace FrameTL
{
  template <class IBASIS, unsigned int DIM>
  pPoissonEquation<IBASIS,DIM>::pPoissonEquation(const PoissonBVP_Coeff<DIM>* ell_bvp,
						 const AggregatedFrame<IBASIS,DIM>* frame,
						 const double parameter_p,
                         const int max_lev_wav_rhs,
                         const int doe_quad_a,
                         const int min_res_quad_a,
                         const int doe_quad_f,
                         const int min_res_quad_f)
    : param_p(parameter_p), max_level_wavelet_basis(frame->jmax()), max_level_wavelets_rhs(max_lev_wav_rhs),
      degree_of_exactness_quadrature_a(doe_quad_a), min_res_quadrature_a(min_res_quad_a), degree_of_exactness_quadrature_f(doe_quad_f), 
      min_res_quadrature_f(min_res_quad_f), ell_bvp_(ell_bvp), frame_(frame)       
  {

    if (doe_quad_a == 0)
    {
      degree_of_exactness_quadrature_a = IBASIS::primal_polynomial_degree()*IBASIS::primal_polynomial_degree();
    }

    clock_t begin_compute_diag = clock();

    // precomputation of the right-hand side up to the maximal level
    compute_diagonal();

    clock_t end_compute_diag = clock();
    time_consumption_of_compute_diagonal = double(end_compute_diag - begin_compute_diag) / CLOCKS_PER_SEC;


    clock_t begin_compute_rhs = clock();

    // precomputation of the diagonal up to the maximal level
    compute_rhs();

    clock_t end_compute_rhs = clock();
    time_consumption_of_compute_rhs = double(end_compute_rhs - begin_compute_rhs) / CLOCKS_PER_SEC;
  }

  template <class IBASIS, unsigned int DIM>
  double
  pPoissonEquation<IBASIS,DIM>::D(const typename AggregatedFrame<IBASIS,DIM>::Index& lambda) const
  {
    //return stiff_diagonal.get_coefficient(lambda);
    return stiff_diagonal[lambda.number()];
  }

  template <class IBASIS, unsigned int DIM>
  void
  pPoissonEquation<IBASIS,DIM>::rescale(InfiniteVector<double,
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
  pPoissonEquation<IBASIS,DIM>::compute_rhs()
  {
    cout << "pPoissonEquation(): precompute right-hand side..." << endl;
    cout << std::setprecision(12);

    typedef AggregatedFrame<IBASIS,DIM> Frame;
    typedef typename Frame::Index Index;

    InfiniteVector<double,Index> fhelp;
    Array1D<InfiniteVector<double,Index> > fhelp_patch(frame_->n_p());

    // initialize array fnorms_sqr_patch
    fnorms_sqr_patch.resize(frame_->n_p());
    for (unsigned int i = 0; i < fnorms_sqr_patch.size(); i++)
    {
      fnorms_sqr_patch[i] = 0.;
    }


#if 1  // compute frame coordinates of f up to the level given by member 'max_level_wavelets_rhs'
    //! loop over all wavelets between minimal and maximal level
    int j = 0;
    for (Index lambda(frame_->first_generator(frame_->j0()));; ++lambda)
    {
	  double coeff = f(lambda)/D(lambda);
      if (fabs(coeff)>1e-15)
      {
	    fhelp.set_coefficient(lambda, coeff);
	    fhelp_patch[lambda.p()].set_coefficient(lambda, coeff);

	    fnorms_sqr_patch[lambda.p()] += coeff*coeff;
	    if (j % 100 == 0)
	    {
	      cout << lambda << " " << coeff << endl;
        }
	  }
      j++;
  	  if (lambda == frame_->last_wavelet(max_level_wavelets_rhs))
	    break;
    }
#endif

#if 0  // compute frame coordinates of f up to the maximal level, which is given by the max. level of 'frame_'
    //! loop over all wavelets between minimal and maximal level
    for (int i = 0; i < frame_->degrees_of_freedom(); i++)
    {
	// computation of one right-hand side coefficient
	  double coeff = f(*(frame_->get_wavelet(i)))/D(*(frame_->get_wavelet(i)));
	// put the coefficient into an InfiniteVector and successively
	// compute the squared \ell_2 norm
	  if (fabs(coeff)>1e-15)
      {
	    fhelp.set_coefficient(*(frame_->get_wavelet(i)), coeff);
	    fhelp_patch[frame_->get_wavelet(i)->p()].set_coefficient(*(frame_->get_wavelet(i)), coeff);

	    fnorms_sqr_patch[frame_->get_wavelet(i)->p()] += coeff*coeff;
	    if (i % 100 == 0)
	      cout << *(frame_->get_wavelet(i)) << " " << coeff << endl;
	  }
    }
#endif

    fnorm_sqr = l2_norm_sqr(fhelp);

    for (unsigned int i = 0; i < fnorms_sqr_patch.size(); i++)
      cout << fnorms_sqr_patch[i] << endl;

    cout << "norm rhs sqr = " << fnorm_sqr << endl;
    cout << "... done, all integrals for right-hand side computed" << endl;
    cout << std::setprecision(6);

    //! sort the coefficients into fcoeffs
    fcoeffs.resize(0); // clear eventual old values
    fcoeffs.resize(fhelp.size());
    unsigned int id(0);
    for (typename InfiniteVector<double,Index>::const_iterator it(fhelp.begin()), itend(fhelp.end());
	     it != itend; ++it, ++id)
    {
      fcoeffs[id] = std::pair<Index,double>(it.index(), *it);
    }
    sort(fcoeffs.begin(), fcoeffs.end(), typename InfiniteVector<double,Index>::decreasing_order());

    // the same patchwise
    fcoeffs_patch.resize(frame_->n_p());
    for (int i = 0; i < frame_->n_p(); i++)
    {
      //fcoeffs_patch[i].resize(0); // clear eventual old values
      fcoeffs_patch[i].resize(fhelp_patch[i].size());
      id = 0;
      for (typename InfiniteVector<double,Index>::const_iterator it(fhelp_patch[i].begin()), itend(fhelp_patch[i].end());
	       it != itend; ++it, ++id)
      {
	    (fcoeffs_patch[i])[id] = std::pair<Index,double>(it.index(), *it);
      }
      sort(fcoeffs_patch[i].begin(), fcoeffs_patch[i].end(), typename InfiniteVector<double,Index>::decreasing_order());
    }

  }

  template <class IBASIS, unsigned int DIM>
  void
  pPoissonEquation<IBASIS,DIM>::compute_diagonal()
  {
    cout << "pPoissonEquation(): precompute diagonal of stiffness matrix..." << endl;

    // precompute the right-hand side on a fine level

    stiff_diagonal.resize(frame_->degrees_of_freedom());
    for (int i = 0; i < frame_->degrees_of_freedom(); i++)
    {
      stiff_diagonal[i] = sqrt(a(*(frame_->get_wavelet(i)),*(frame_->get_wavelet(i))));
    }

    cout << "... done, diagonal of stiffness matrix computed" << endl;
  }

  template <class IBASIS, unsigned int DIM>
  double
  pPoissonEquation<IBASIS,DIM>::a_same_patches(const typename AggregatedFrame<IBASIS,DIM>::Index& lambda,
					       const typename AggregatedFrame<IBASIS,DIM>::Index& mu,
					       const unsigned int n_Gauss_knots) const
  {

    clock_t begin_a = clock();


    double r = 0;

    //patchnumbers are assumed to be equal
    const unsigned int p = lambda.p();

    typedef WaveletTL::CubeBasis<IBASIS,DIM> CUBEBASIS;
    typedef typename CUBEBASIS::Index CubeIndex;

    typename CUBEBASIS::Support supp_intersect;


    // check whether the supports of the reference wavelets intersect,
    // the intersection is returned in supp_intersect
    bool b = WaveletTL::intersect_supports<IBASIS,DIM>
      (
       *(frame_->bases()[p]),
       CubeIndex(lambda.j(), lambda.e(), lambda.k(), frame_->bases()[p]),
       CubeIndex(mu.j(), mu.e(), mu.k(), frame_->bases()[p]),
       supp_intersect
       );

    // no intersection ==> return zero
    if (! b)
    {
      return 0.0;
    }

    // pointer to loca basis
    CUBEBASIS* local_cube_basis = frame_->bases()[p];
    const Chart<DIM>* chart = frame_->atlas()->charts()[p];

    // number of Gauss quadrature knots
    const int N_Gauss = n_Gauss_knots;
    //cout << "N_Gauss = " << N_Gauss << endl;

    // granularity for the quadrature
    //const double h = ldexp(1.0, -supp_intersect.j);
    const double h = 1.0 / (1 << supp_intersect.j);
    //cout << "h=" << h << endl;

    // for each coordinate direction we shall need an array for
    // the quadrature knots and weights, and some others
    FixedArray1D<Array1D<double>,DIM> gauss_points, gauss_weights,
      wav_values_lambda, wav_der_values_lambda, wav_values_mu, wav_der_values_mu;


    for (unsigned int i = 0; i < DIM; i++)
    {
      gauss_points[i].resize(N_Gauss*(supp_intersect.b[i]-supp_intersect.a[i]));
      gauss_weights[i].resize(N_Gauss*(supp_intersect.b[i]-supp_intersect.a[i]));

      // set up Gauss knots and weights for coordinate direction i
      // on each interval where the the integrand is a polynomial
      for (int patch = supp_intersect.a[i]; patch < supp_intersect.b[i]; patch++)
      {
	    for (int n = 0; n < N_Gauss; n++)
	    {
	      gauss_points[i][(patch-supp_intersect.a[i])*N_Gauss+n]
	        = h*(2*patch+1+GaussPoints[N_Gauss-1][n])/2.;
	      gauss_weights[i][(patch-supp_intersect.a[i])*N_Gauss+n]
	        = h*GaussWeights[N_Gauss-1][n];
	    }
	  }

    }

    // compute the point values of the wavelet part of the integrand
    for (unsigned int i = 0; i < DIM; i++)
    {
      WaveletTL::evaluate(*(local_cube_basis->bases()[i]), 0,
			  typename IBASIS::Index(lambda.j(),
						 lambda.e()[i],
						 lambda.k()[i],
						 local_cube_basis->bases()[i]),
			  gauss_points[i], wav_values_lambda[i]);
      WaveletTL::evaluate(*(local_cube_basis->bases()[i]), 1,
			  typename IBASIS::Index(lambda.j(),
						 lambda.e()[i],
						 lambda.k()[i],
						 local_cube_basis->bases()[i]),
			  gauss_points[i], wav_der_values_lambda[i]);

      WaveletTL::evaluate(*(local_cube_basis->bases()[i]), 0,
			  typename IBASIS::Index(mu.j(),
						 mu.e()[i],
						 mu.k()[i],
						 local_cube_basis->bases()[i]),
			  gauss_points[i], wav_values_mu[i]);

      WaveletTL::evaluate(*(local_cube_basis->bases()[i]), 1,
			  typename IBASIS::Index(mu.j(),
						 mu.e()[i],
						 mu.k()[i],
						 local_cube_basis->bases()[i]),
			  gauss_points[i], wav_der_values_mu[i]);
      //cout << wav_der_values_mu[i] << endl;
    }

    int index[DIM]; // current multiindex for the point values
    for (unsigned int i = 0; i < DIM; i++)
    {
      index[i] = 0;
    }

    Point<DIM> x;
    Point<DIM> x_patch;

    // now we perform the quadrature
    // loop over all quadrature knots
    while (true)
    {
      for (unsigned int i = 0; i < DIM; i++)
      {
	    x[i] = gauss_points[i][index[i]];
      }
      chart->map_point(x,x_patch);
      double sq_gram = chart->Gram_factor(x);

      Vector<double> values1(DIM);
      Vector<double> values2(DIM);
      double weights=1., psi_lambda=1., psi_mu=1.;

      for (unsigned int i = 0; i < DIM; i++)
      {
	    weights *= gauss_weights[i][index[i]];
	    psi_lambda *= wav_values_lambda[i][index[i]];
	    psi_mu *= wav_values_mu[i][index[i]];
      }

      if ( !(psi_lambda == 0. || psi_mu == 0.) )
      {

        double value_coeff_a = ell_bvp_->a(x_patch);

	    for (unsigned int s = 0; s < DIM; s++)
	    {
	      double psi_der_lambda=1., psi_der_mu=1.;
	      psi_der_lambda = (psi_lambda / wav_values_lambda[s][index[s]]) * wav_der_values_lambda[s][index[s]];
	      psi_der_mu = (psi_mu / wav_values_mu[s][index[s]]) * wav_der_values_mu[s][index[s]];

	      // for first part of the integral: \int_\Omega <a \Nabla u,\Nabla v> dx
	      double tmp = chart->Gram_D_factor(s,x);

	      values1[s] = value_coeff_a *
	        (psi_der_lambda*sq_gram - (psi_lambda*tmp)) / (sq_gram*sq_gram);
	      values2[s] =
	        (psi_der_mu*sq_gram - (psi_mu*tmp)) / (sq_gram*sq_gram);

	    }//end for s


	    double t = 0.;
	    for (unsigned int i1 = 0; i1 < DIM; i1++)
	    {
	      double d1 = 0.;
	      double d2 = 0.;
	      for (unsigned int i2 = 0; i2 < DIM; i2++)
	      {
	        double tmp = chart->Dkappa_inv(i2, i1, x);
	        d1 += values1[i2]*tmp;
	        d2 += values2[i2]*tmp;
	      }
	      t += d1 * d2;

	    }

	    r += (t * (sq_gram*sq_gram) ) * weights;

      }


      // "++index"
      bool exit = false;
      for (unsigned int i = 0; i < DIM; i++)
      {
	    if (index[i] == N_Gauss*(supp_intersect.b[i]-supp_intersect.a[i])-1)
	    {
	      index[i] = 0;
	      exit = (i == DIM-1);
	    }
	    else
	    {
	      index[i]++;
	      break;
	    }
      }

      if (exit) {
        break;
      }

    }//end while


    clock_t end_a = clock();
    double elapsed_secs = double(end_a - begin_a) / CLOCKS_PER_SEC;
    time_consumption_of_a_same_patches += elapsed_secs;

    return r;
  }



  template <class IBASIS, unsigned int DIM>
  double
  pPoissonEquation<IBASIS,DIM>::a_different_patches(const typename AggregatedFrame<IBASIS,DIM>::Index& la,
						    const typename AggregatedFrame<IBASIS,DIM>::Index& nu,
						    const unsigned int n_Gauss_knots, const unsigned int N) const
  {

    clock_t begin_a = clock();

    double r = 0.0;

    Index lambda = la;
    Index mu = nu;
    Index tmp_ind;

    typedef WaveletTL::CubeBasis<IBASIS,DIM> CUBEBASIS;

//    typedef typename CUBEBASIS::Index CubeIndex;

    typename CUBEBASIS::Support supp_lambda_ = frame_->all_supports[lambda.number()];
    typename CUBEBASIS::Support supp_mu_ = frame_->all_supports[mu.number()];

    typename CUBEBASIS::Support tmp_supp;

    // we want to assume that the index lambda has the
    // larger level j of the two indices involved,
    // swap indices and supports if necessary
    if (supp_mu_.j > supp_lambda_.j)
    {
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

    // number of Gauss quadrature knots
    const int N_Gauss = n_Gauss_knots;

    // check whether the supports of the reference wavelets intersect,
    // the intersection is returned in supp_intersect
    bool b = 0;
    b = intersect_supports<IBASIS,DIM,DIM>(*frame_, lambda, mu, supp_lambda, supp_mu);

    // no intersection ==> return zero
    if ( !b )
      return 0.0;

    typedef typename IBASIS::Index Index_1D;

    // get pointers to the parametric mappings
    const Chart<DIM>* chart_la = frame_->atlas()->charts()[lambda.p()];
    const Chart<DIM>* chart_mu = frame_->atlas()->charts()[mu.p()];

    // get array of the one-dimensional bases from each spatial direction
    FixedArray1D<IBASIS*,DIM> bases1D_lambda = frame_->bases()[lambda.p()]->bases();
    FixedArray1D<IBASIS*,DIM> bases1D_mu     = frame_->bases()[mu.p()]->bases();

    // for each coordinate direction we shall need an array for
    // the quadrature knots and weights, and some others
    FixedArray1D<Array1D<double>,DIM> gauss_points, gauss_points_mu, gauss_weights,
      wav_values_lambda, wav_der_values_lambda, wav_values_mu, wav_der_values_mu;

    //const double h = ldexp(1.0, -std::max(supp_lambda->j,supp_mu->j)); // granularity for the quadrature
    const double h = 1.0 / (1 << std::max(supp_lambda->j,supp_mu->j)); // granularity for the quadrature
    //const int N = 4;

    for (unsigned int i = 0; i < DIM; i++)
    {
      gauss_points[i].resize(N * N_Gauss*(supp_lambda->b[i]-supp_lambda->a[i]));
      gauss_weights[i].resize(N * N_Gauss*(supp_lambda->b[i]-supp_lambda->a[i]));

      // set up Gauss knots and weights for coordinate direction i
      // We have a dyadic partition of the support of the higher level wavelet.
      // With respect to this partition the wavelet is a piecewise polynomial.
      // On each of the polynomial parts, we apply a composite quadrature rule.
      for (int patch = supp_lambda->a[i]; patch < supp_lambda->b[i]; patch++)
      {
	    for (unsigned int m = 0; m < N; m++)
	    {
	      for (int n = 0; n < N_Gauss; n++)
	      {
	        gauss_points[i][ N*(patch-supp_lambda->a[i])*N_Gauss + m*N_Gauss+n ]
	          = h*( 1.0/(2.*N)*(GaussPoints[N_Gauss-1][n]+1+2.0*m)+patch);

	        gauss_weights[i][ N*(patch-supp_lambda->a[i])*N_Gauss + m*N_Gauss+n ]
	          = (h*GaussWeights[N_Gauss-1][n])/N;
	      }
        }
      }

    }

    for (unsigned int i = 0; i < DIM; i++)
    {
      // compute function values at gauss_points[i] of the 1D wavelets given by the index lambda
      WaveletTL::evaluate(*(bases1D_lambda[i]), 0,
			  Index_1D(lambda.j(), lambda.e()[i], lambda.k()[i], bases1D_lambda[i]),
			  gauss_points[i], wav_values_lambda[i]);

      // compute function values at gauss_points[i] of the derivatives of the 1D wavelets given by the index lambda
      WaveletTL::evaluate(*(bases1D_lambda[i]), 1,
			  Index_1D(lambda.j(), lambda.e()[i], lambda.k()[i], bases1D_lambda[i]),
			  gauss_points[i], wav_der_values_lambda[i]);
    }

    int index[DIM]; // current multiindex for the point values
    for (unsigned int i = 0; i < DIM; i++)
      index[i] = 0;

    Point<DIM> x;
    Point<DIM> x_patch;
    Point<DIM> y;

    // now we perform the quadrature,
    // loop over all quadrature knots
    while (true)
    {
      for (unsigned int i = 0; i < DIM; i++)
      {
	    x[i] = gauss_points[i][index[i]];
      }
      //cout << x << endl;
      chart_la->map_point(x,x_patch);
      if ( in_support(*frame_,mu, supp_mu, x_patch) )
	  {
	    chart_mu->map_point_inv(x_patch,y);


	    double sq_gram_la = chart_la->Gram_factor(x);
	    double sq_gram_mu = chart_mu->Gram_factor(y);

	    double weight=1., psi_lambda=1., psi_mu=1.;
	    for (unsigned int i = 0; i < DIM; i++)
	    {
	      weight *= gauss_weights[i][index[i]];
	      psi_lambda *= wav_values_lambda[i][index[i]];
	      psi_mu *=  WaveletTL::evaluate(*(bases1D_mu[i]), 0,
		  			   Index_1D(mu.j(), mu.e()[i], mu.k()[i], bases1D_mu[i]), y[i]);

	    }
	    if ( !(psi_lambda == 0.|| psi_mu == 0.) )
	    {

	      Vector<double> values1(DIM);
	      Vector<double> values2(DIM);

	      for (unsigned int s = 0; s < DIM; s++)
	      {
	        double psi_der_lambda=1., psi_der_mu=1.;
	        psi_der_lambda = (psi_lambda / wav_values_lambda[s][index[s]]) * wav_der_values_lambda[s][index[s]];

	        psi_der_mu = (psi_mu / WaveletTL::evaluate(*(bases1D_mu[s]), 0,
							 Index_1D(mu.j(), mu.e()[s], mu.k()[s], bases1D_mu[s]), y[s])
			    ) * WaveletTL::evaluate(*(bases1D_mu[s]), 1,
						    Index_1D(mu.j(), mu.e()[s], mu.k()[s], bases1D_mu[s]), y[s]);

	        values1[s] = ell_bvp_->a(x_patch) *
		      (psi_der_lambda*sq_gram_la - (psi_lambda*chart_la->Gram_D_factor(s,x)))
		      / (sq_gram_la*sq_gram_la);
	        values2[s] =
              (psi_der_mu*sq_gram_mu - (psi_mu*chart_mu->Gram_D_factor(s,y)))
		      / (sq_gram_mu*sq_gram_mu);

	      } // end loop s

	      double tmp = 0.;
	      for (unsigned int i1 = 0; i1 < DIM; i1++)
	      {
	        double d1 = 0.;
	        double d2 = 0.;
	        for (unsigned int i2 = 0; i2 < DIM; i2++)
	        {
		      d1 += values1[i2]*chart_la->Dkappa_inv(i2, i1, x);
		      d2 += values2[i2]*chart_mu->Dkappa_inv(i2, i1, y);
	        }
	        tmp += d1 * d2;
	      }

	      r += (tmp * (sq_gram_la*sq_gram_la) +
		       ell_bvp_->q(x_patch) * psi_lambda * (psi_mu/sq_gram_mu) * sq_gram_la) * weight;
	    }
	  }
      // "++index"
      bool exit = false;
      for (unsigned int i = 0; i < DIM; i++)
      {
	    if (index[i] == (int) N * N_Gauss * (supp_lambda->b[i]-supp_lambda->a[i]) - 1)
	    {
	      index[i] = 0;
	      exit = (i == DIM-1);
	    }
	    else
	    {
	      index[i]++;
	      break;
	    }
      }
      if (exit) break;
    }


    clock_t end_a = clock();
    double elapsed_secs = double(end_a - begin_a) / CLOCKS_PER_SEC;
    time_consumption_of_a_different_patches += elapsed_secs;


    return r;
  }



//! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! Variant of a_same_patches(,) which caches the point evaluations of the
//! coefficient function bvp_->a    (by use of member 'cached_coeff_a')
//! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  template <class IBASIS, unsigned int DIM>
  double
  pPoissonEquation<IBASIS,DIM>::a_same_patches2(const typename AggregatedFrame<IBASIS,DIM>::Index& lambda,
					       const typename AggregatedFrame<IBASIS,DIM>::Index& mu,
					       const unsigned int doe_quadrature, const unsigned int min_resolution) const
  {

    clock_t begin_a = clock();


    double r = 0;

    //patchnumbers are assumed to be equal
    const unsigned int p = lambda.p();

    typedef WaveletTL::CubeBasis<IBASIS,DIM> CUBEBASIS;
    typedef typename CUBEBASIS::Index CubeIndex;

    typename CUBEBASIS::Support supp_intersect;


    // check whether the supports of the reference wavelets intersect,
    // the intersection is returned in supp_intersect
    bool b = WaveletTL::intersect_supports<IBASIS,DIM>
      (
       *(frame_->bases()[p]),
       CubeIndex(lambda.j(), lambda.e(), lambda.k(), frame_->bases()[p]),
       CubeIndex(mu.j(), mu.e(), mu.k(), frame_->bases()[p]),
       supp_intersect
       );

    // no intersection ==> return zero
    if (! b)
    {
      return 0.0;
    }

    // pointer to loca basis
    CUBEBASIS* local_cube_basis = frame_->bases()[p];
    const Chart<DIM>* chart = frame_->atlas()->charts()[p];

    // number of Gauss quadrature knots
    const int N_Gauss = (doe_quadrature+1)/2;
    //cout << "N_Gauss = " << N_Gauss << endl;

    // refine, if necessary, the resolution of the support
    if (supp_intersect.j < (int) min_resolution)
    {
      unsigned int diff = (1 << (min_resolution - supp_intersect.j));
      for (unsigned int i = 0; i < DIM; i++)
      {
        supp_intersect.a[i] *= diff;
        supp_intersect.b[i] *= diff;
      }
      supp_intersect.j = min_resolution;
    }

    // granularity for the quadrature
    //const double h = ldexp(1.0, -supp_intersect.j);
    const double h = 1.0 / (1 << supp_intersect.j);
    //cout << "h=" << h << endl;

    // for each coordinate direction we shall need an array for
    // the quadrature knots and weights, and some others
    FixedArray1D<Array1D<double>,DIM> gauss_points, gauss_weights,
      wav_values_lambda, wav_der_values_lambda, wav_values_mu, wav_der_values_mu;


    for (unsigned int i = 0; i < DIM; i++)
    {
      gauss_points[i].resize(N_Gauss*(supp_intersect.b[i]-supp_intersect.a[i]));
      gauss_weights[i].resize(N_Gauss*(supp_intersect.b[i]-supp_intersect.a[i]));

      // set up Gauss knots and weights for coordinate direction i
      // on each interval where the the integrand is a polynomial
      for (int patch = supp_intersect.a[i]; patch < supp_intersect.b[i]; patch++)
      {
	    for (int n = 0; n < N_Gauss; n++)
	    {
	      gauss_points[i][(patch-supp_intersect.a[i])*N_Gauss+n]
	        = h*(2*patch+1+GaussPoints[N_Gauss-1][n])/2.;
	      gauss_weights[i][(patch-supp_intersect.a[i])*N_Gauss+n]
	        = h*GaussWeights[N_Gauss-1][n];
	    }
	  }

    }

    // compute the point values of the wavelet part of the integrand
    for (unsigned int i = 0; i < DIM; i++)
    {
      WaveletTL::evaluate(*(local_cube_basis->bases()[i]), 0,
			  typename IBASIS::Index(lambda.j(),
						 lambda.e()[i],
						 lambda.k()[i],
						 local_cube_basis->bases()[i]),
			  gauss_points[i], wav_values_lambda[i]);
      WaveletTL::evaluate(*(local_cube_basis->bases()[i]), 1,
			  typename IBASIS::Index(lambda.j(),
						 lambda.e()[i],
						 lambda.k()[i],
						 local_cube_basis->bases()[i]),
			  gauss_points[i], wav_der_values_lambda[i]);

      WaveletTL::evaluate(*(local_cube_basis->bases()[i]), 0,
			  typename IBASIS::Index(mu.j(),
						 mu.e()[i],
						 mu.k()[i],
						 local_cube_basis->bases()[i]),
			  gauss_points[i], wav_values_mu[i]);

      WaveletTL::evaluate(*(local_cube_basis->bases()[i]), 1,
			  typename IBASIS::Index(mu.j(),
						 mu.e()[i],
						 mu.k()[i],
						 local_cube_basis->bases()[i]),
			  gauss_points[i], wav_der_values_mu[i]);
      //cout << wav_der_values_mu[i] << endl;
    }

    int index[DIM]; // current multiindex for the point values
    for (unsigned int i = 0; i < DIM; i++)
    {
      index[i] = 0;
    }

    Point<DIM> x;
    Point<DIM> x_patch;

    // now we perform the quadrature
    // loop over all quadrature knots
    while (true)
    {
      for (unsigned int i = 0; i < DIM; i++)
      {
	    x[i] = gauss_points[i][index[i]];
      }
      chart->map_point(x,x_patch);

      Vector<double> values1(DIM);
      Vector<double> values2(DIM);
      double weights=1.;

      for (unsigned int i = 0; i < DIM; i++)
      {
	    weights *= gauss_weights[i][index[i]];
      }


      // current patch:
      int current_patch = 0;
      int offset = 1;
      for (unsigned int i = 0; i < DIM; i++)
      {
        current_patch += offset*((index[i] / N_Gauss) + supp_intersect.a[i]);
        offset *= (1 << supp_intersect.j);
      }

      // gauss point in current patch
      int gauss_point_in_patch = 0;
      int offset2 = 1;
      for (unsigned int i = 0; i < DIM; i++)
      {
        gauss_point_in_patch += offset2*(index[i] % N_Gauss);
        offset2 *= N_Gauss;
      }


      // evaluate coeff. function ell_bvp_->a in x_patch (cached!)
      double value_coeff_a = cached_coeff_a(p, current_patch, gauss_point_in_patch, supp_intersect.j, N_Gauss);



#if 0   //! sanity check for cached_coeff_a
        double a_value_original = ell_bvp_->a(x_patch);

        if (abs(value_coeff_a - a_value_original) > 1e-16)
        {
          cout << std::setprecision(20);
          cout << "ERROR in a: cached_coeff_a = " << value_coeff_a << ", bvp_->a(x) = " << a_value_original << endl;
          cout << "global patch: " << p <<" current patch: " << current_patch << " gauss point in patch: " << gauss_point_in_patch << " j: " << supp_intersect.j << " N_Gauss: " << N_Gauss << endl;
          cout << std::setprecision(6);
          exit(1);
        }
#endif


      values1[0] = value_coeff_a * wav_values_lambda[1][index[1]] * wav_der_values_lambda[0][index[0]];
      values1[1] = value_coeff_a * wav_values_lambda[0][index[0]] * wav_der_values_lambda[1][index[1]];
      values2[0] = wav_values_mu[1][index[1]] * wav_der_values_mu[0][index[0]];
      values2[1] = wav_values_mu[0][index[0]] * wav_der_values_mu[1][index[1]];

	  double t = 0.;
	  for (unsigned int i = 0; i < DIM; i++)
	  {
	    t += ( values1[i] * values2[i] * chart->Dkappa_inv(i, i, x) * chart->Dkappa_inv(i, i, x) );
	  }

	  r += (t * weights);

      // "++index"
      bool exit = false;
      for (unsigned int i = 0; i < DIM; i++)
      {
	    if (index[i] == N_Gauss*(supp_intersect.b[i]-supp_intersect.a[i])-1)
	    {
	      index[i] = 0;
	      exit = (i == DIM-1);
	    }
	    else
	    {
	      index[i]++;
	      break;
	    }
      }

      if (exit) {
        break;
      }

    }//end while


    clock_t end_a = clock();
    double elapsed_secs = double(end_a - begin_a) / CLOCKS_PER_SEC;
    time_consumption_of_a_same_patches += elapsed_secs;


    return r;
  }


//! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! Variant of a_different_patches(,) which caches the point evaluations of the
//! coefficient function bvp_->a    (by use of member 'cached_coeff_a')
//! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  template <class IBASIS, unsigned int DIM>
  double
  pPoissonEquation<IBASIS,DIM>::a_different_patches2(const typename AggregatedFrame<IBASIS,DIM>::Index& la,
						    const typename AggregatedFrame<IBASIS,DIM>::Index& nu,
						    const unsigned int doe_quadrature, const unsigned int min_resolution,
						    const unsigned int N) const
  {

    clock_t begin_a = clock();

    double r = 0.0;

    Index lambda = la;
    Index mu = nu;
    Index tmp_ind;

    typedef WaveletTL::CubeBasis<IBASIS,DIM> CUBEBASIS;

//    typedef typename CUBEBASIS::Index CubeIndex;

    typename CUBEBASIS::Support supp_lambda_ = frame_->all_supports[lambda.number()];
    typename CUBEBASIS::Support supp_mu_ = frame_->all_supports[mu.number()];

    typename CUBEBASIS::Support tmp_supp;


    // check whether the supports of the reference wavelets intersect,
    bool b = 0;
    b = intersect_supports<IBASIS,DIM,DIM>(*frame_, lambda, mu, &supp_lambda_, &supp_mu_);

    // no intersection ==> return zero
    if ( !b )
      return 0.0;


    // we want to assume that the index lambda has the
    // larger level j of the two indices involved,
    // swap indices and supports if necessary
    if (supp_mu_.j > supp_lambda_.j)
    {
      tmp_ind = lambda;
      lambda = mu;
      mu = tmp_ind;

      tmp_supp = supp_lambda_;
      supp_lambda_ = supp_mu_,
	  supp_mu_ = tmp_supp;
    }

    const unsigned int p = lambda.p();


    // refine, if necessary, the resolution of the support
    if (supp_lambda_.j < (int) min_resolution)
    {
      unsigned int diff = (1 << (min_resolution - supp_lambda_.j));
      for (unsigned int i = 0; i < DIM; i++)
      {
        supp_lambda_.a[i] *= diff;
        supp_lambda_.b[i] *= diff;
      }
      supp_lambda_.j = min_resolution;
    }


    const typename CUBEBASIS::Support* supp_lambda = &supp_lambda_;
    const typename CUBEBASIS::Support* supp_mu = &supp_mu_;


    FixedArray1D<Array1D<double>,DIM > irregular_grid;

    // number of Gauss quadrature knots
    const int N_Gauss = (doe_quadrature+1)/2;


    typedef typename IBASIS::Index Index_1D;

    // get pointers to the parametric mappings
    const Chart<DIM>* chart_la = frame_->atlas()->charts()[lambda.p()];
    const Chart<DIM>* chart_mu = frame_->atlas()->charts()[mu.p()];

    // get array of the one-dimensional bases from each spatial direction
    FixedArray1D<IBASIS*,DIM> bases1D_lambda = frame_->bases()[lambda.p()]->bases();
    FixedArray1D<IBASIS*,DIM> bases1D_mu     = frame_->bases()[mu.p()]->bases();

    // for each coordinate direction we shall need an array for
    // the quadrature knots and weights, and some others
    FixedArray1D<Array1D<double>,DIM> gauss_points, gauss_points_mu, gauss_weights,
      wav_values_lambda, wav_der_values_lambda, wav_values_mu, wav_der_values_mu;

    //const double h = ldexp(1.0, -std::max(supp_lambda->j,supp_mu->j)); // granularity for the quadrature
    const double h = 1.0 / (1 << std::max(supp_lambda->j,supp_mu->j)); // granularity for the quadrature
    //const int N = 4;

    const int res_level = (supp_lambda->j +  (N > 1 ? log2(N) : 0 ) );

    for (unsigned int i = 0; i < DIM; i++)
    {
      gauss_points[i].resize(N * N_Gauss*(supp_lambda->b[i]-supp_lambda->a[i]));
      gauss_weights[i].resize(N * N_Gauss*(supp_lambda->b[i]-supp_lambda->a[i]));

      // set up Gauss knots and weights for coordinate direction i
      // We have a dyadic partition of the support of the higher level wavelet.
      // With respect to this partition the wavelet is a piecewise polynomial.
      // On each of the polynomial parts, we apply a composite quadrature rule.
      for (int patch = supp_lambda->a[i]; patch < supp_lambda->b[i]; patch++)
      {
	    for (unsigned int m = 0; m < N; m++)
	    {
	      for (int n = 0; n < N_Gauss; n++)
	      {
	        gauss_points[i][ N*(patch-supp_lambda->a[i])*N_Gauss + m*N_Gauss+n ]
	          = h*( 1.0/(2.*N)*(GaussPoints[N_Gauss-1][n]+1+2.0*m)+patch);

	        gauss_weights[i][ N*(patch-supp_lambda->a[i])*N_Gauss + m*N_Gauss+n ]
	          = (h*GaussWeights[N_Gauss-1][n])/N;
	      }
        }
      }

    }

    for (unsigned int i = 0; i < DIM; i++)
    {
      // compute function values at gauss_points[i] of the 1D wavelets given by the index lambda
      WaveletTL::evaluate(*(bases1D_lambda[i]), 0,
			  Index_1D(lambda.j(), lambda.e()[i], lambda.k()[i], bases1D_lambda[i]),
			  gauss_points[i], wav_values_lambda[i]);

      // compute function values at gauss_points[i] of the derivatives of the 1D wavelets given by the index lambda
      WaveletTL::evaluate(*(bases1D_lambda[i]), 1,
			  Index_1D(lambda.j(), lambda.e()[i], lambda.k()[i], bases1D_lambda[i]),
			  gauss_points[i], wav_der_values_lambda[i]);
    }

    int index[DIM]; // current multiindex for the point values
    for (unsigned int i = 0; i < DIM; i++)
      index[i] = 0;

    Point<DIM> x;
    Point<DIM> x_patch;
    Point<DIM> y;

    // now we perform the quadrature,
    // loop over all quadrature knots
    while (true)
    {
      for (unsigned int i = 0; i < DIM; i++)
      {
	    x[i] = gauss_points[i][index[i]];
      }
      //cout << x << endl;
      chart_la->map_point(x,x_patch);
      if ( in_support(*frame_,mu, supp_mu, x_patch) )
	  {
	    chart_mu->map_point_inv(x_patch,y);


	    double sq_gram_la = chart_la->Gram_factor(x);
	    double sq_gram_mu = chart_mu->Gram_factor(y);

	    double weight=1., psi_lambda=1., psi_mu=1.;
	    for (unsigned int i = 0; i < DIM; i++)
	    {
	      weight *= gauss_weights[i][index[i]];
	      psi_lambda *= wav_values_lambda[i][index[i]];
	      psi_mu *=  WaveletTL::evaluate(*(bases1D_mu[i]), 0,
		  			   Index_1D(mu.j(), mu.e()[i], mu.k()[i], bases1D_mu[i]), y[i]);

	    }
	    if ( !(psi_lambda == 0.|| psi_mu == 0.) )
	    {

	      Vector<double> values1(DIM);
	      Vector<double> values2(DIM);


          // current patch:
          int current_patch = 0;
          int offset = 1;
          for (unsigned int i = 0; i < DIM; i++)
          {
            current_patch += offset*((index[i] / N_Gauss) + supp_lambda->a[i] * N);
            offset *= (1 << res_level);
          }

          // gauss point in current patch
          int gauss_point_in_patch = 0;
          int offset2 = 1;
          for (unsigned int i = 0; i < DIM; i++)
          {
            gauss_point_in_patch += offset2*(index[i] % N_Gauss);
            offset2 *= N_Gauss;
          }


          // evaluate coeff. function ell_bvp_->a in x_patch (cached!)
          double value_coeff_a = cached_coeff_a(p, current_patch, gauss_point_in_patch, res_level, N_Gauss);


#if 0   //! sanity check for cached_coeff_a
        double a_value_original = ell_bvp_->a(x_patch);

        if (abs(value_coeff_a - a_value_original) > 1e-16)
        {
          cout << std::setprecision(20);
          cout << "ERROR in a: cached_coeff_a = " << value_coeff_a << ", bvp_->a(x) = " << a_value_original << endl;
          cout << "global patch: " << p <<" current patch: " << current_patch << " gauss point in patch: " << gauss_point_in_patch << " j: " << supp_lambda->j << " N_Gauss: " << N_Gauss << endl;
          cout << std::setprecision(6);
          exit(1);
        }
#endif


	      for (unsigned int s = 0; s < DIM; s++)
	      {
	        double psi_der_lambda=1., psi_der_mu=1.;
	        psi_der_lambda = (psi_lambda / wav_values_lambda[s][index[s]]) * wav_der_values_lambda[s][index[s]];

	        psi_der_mu = (psi_mu / WaveletTL::evaluate(*(bases1D_mu[s]), 0,
							 Index_1D(mu.j(), mu.e()[s], mu.k()[s], bases1D_mu[s]), y[s])
			    ) * WaveletTL::evaluate(*(bases1D_mu[s]), 1,
						    Index_1D(mu.j(), mu.e()[s], mu.k()[s], bases1D_mu[s]), y[s]);

	        values1[s] = value_coeff_a *
		      (psi_der_lambda*sq_gram_la - (psi_lambda*chart_la->Gram_D_factor(s,x)))
		      / (sq_gram_la*sq_gram_la);
	        values2[s] =
              (psi_der_mu*sq_gram_mu - (psi_mu*chart_mu->Gram_D_factor(s,y)))
		      / (sq_gram_mu*sq_gram_mu);

	      } // end loop s

	      double tmp = 0.;
	      for (unsigned int i1 = 0; i1 < DIM; i1++)
	      {
	        double d1 = 0.;
	        double d2 = 0.;
	        for (unsigned int i2 = 0; i2 < DIM; i2++)
	        {
		      d1 += values1[i2]*chart_la->Dkappa_inv(i2, i1, x);
		      d2 += values2[i2]*chart_mu->Dkappa_inv(i2, i1, y);
	        }
	        tmp += d1 * d2;
	      }

	      r += (tmp * (sq_gram_la*sq_gram_la) +
		       ell_bvp_->q(x_patch) * psi_lambda * (psi_mu/sq_gram_mu) * sq_gram_la) * weight;
	    }
	  }
      // "++index"
      bool exit = false;
      for (unsigned int i = 0; i < DIM; i++)
      {
	    if (index[i] == (int) N * N_Gauss * (supp_lambda->b[i]-supp_lambda->a[i]) - 1)
	    {
	      index[i] = 0;
	      exit = (i == DIM-1);
	    }
	    else
	    {
	      index[i]++;
	      break;
	    }
      }
      if (exit) break;
    }


    clock_t end_a = clock();
    double elapsed_secs = double(end_a - begin_a) / CLOCKS_PER_SEC;
    time_consumption_of_a_different_patches += elapsed_secs;


    return r;
  }



//! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! (cached!) point evaluation of the coefficient function bvp_->a
//! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  template <class IBASIS, unsigned int DIM>
  double
  pPoissonEquation<IBASIS,DIM>::cached_coeff_a(const int global_patch,
                                               const int current_patch,
                                               const int gauss_point_in_patch,
                                               const int j,
                                               const int N_Gauss) const
  {
      double r = 0.0;

      // search for Global_Patch 'global_patch'
      typename Global_Patch::iterator gp_lb(coeff_cache.lower_bound(global_patch));
      typename Global_Patch::iterator gp_it(gp_lb);

      if (gp_lb == coeff_cache.end() ||
	  coeff_cache.key_comp()(global_patch, gp_lb->first))
      {
	    // insert a new resolution
	    typedef typename Global_Patch::value_type value_type;
	    gp_it = coeff_cache.insert(gp_lb, value_type(global_patch, Resolution()));
	  }

      Resolution& res(gp_it->second);


      // search for partition with resolution j
      typename Resolution::iterator res_lb(res.lower_bound(j));
      typename Resolution::iterator res_it(res_lb);

      if (res_lb == res.end() ||
	  res.key_comp()(j, res_lb->first))
      {
	    // insert a new partition
	    typedef typename Resolution::value_type value_type;
	    res_it = res.insert(res_lb, value_type(j, Partition()));
	  }

      Partition& part(res_it->second);


      // search for the patch 'current_patch'
      typename Partition::iterator part_lb(part.lower_bound(current_patch));
      typename Partition::iterator part_it(part_lb);

      // no entries have ever been computed for this patch and this partition
      if (part_lb == part.end() ||
	  part.key_comp()(current_patch, part_lb->first))
	  {
        // insert a new patch
        typedef typename Partition::value_type value_type;
	    part_it = part.insert(part_lb, value_type(current_patch, Patch()));

      }

      Patch& patch(part_it->second);


      // search for Func_values 'N_Gauss'
      typename Patch::iterator patch_lb(patch.lower_bound(N_Gauss));
      typename Patch::iterator patch_it(patch_lb);

      // no entries have ever been computed for this patch and this partition and this Func_values
      if (patch_lb == patch.end() ||
	  patch.key_comp()(N_Gauss, patch_lb->first))
	  {
        // insert a Func_values
        typedef typename Patch::value_type value_type;
	    patch_it = patch.insert(patch_lb, value_type(N_Gauss, Func_values()));

	    Func_values& func_values(patch_it->second);


	    // compute the coefficient function values for ALL Gauss points in this patch

        Point<DIM> x;
        Point<DIM> x_patch;
        const double h = ldexp(1.0, -j);

        int offset[DIM];
        int dummy = current_patch;

        for (unsigned int i = 0; i < DIM; i++)
        {
          offset[i] = (dummy % (1 << j));
          dummy -= offset[i];
          dummy /= (1 << j);
        }


       // int gauss_points_per_patch = 1;
        int index[DIM];
        for (unsigned int i = 0; i < DIM; i++)
        {
         // gauss_points_per_patch *= N_Gauss;
          index[i] = 0;
        }

        int current_index = 0;

        while (true)
        {

          // setup gauss point x with index i
          for (unsigned int i = 0; i < DIM; i++)
          {
            x[i] = h*(2*offset[i]+1+GaussPoints[N_Gauss-1][index[i]])/2.;
          }

          // map x to patch
          frame_->atlas()->charts()[global_patch]->map_point(x,x_patch);

          // compute entry at gauss point x
          double entry = ell_bvp_->a(x_patch);


	      typedef typename Func_values::value_type value_type_Func_values;

	      func_values.insert(func_values.end(), value_type_Func_values(current_index, entry));
	      if (current_index == gauss_point_in_patch)
	      {
		    r = entry;
	      }

          // "++index"
	      bool exit = false;
	      for (unsigned int i = 0; i < DIM; i++)
	      {
	        if (index[i] == (N_Gauss-1))
	        {
		      index[i] = 0;
		      exit = (i == DIM-1);
	        }
	        else
	        {
		      index[i]++;
		      break;
	        }
	      }
	      if (exit) break;
          current_index++;

	    }

      }

      // Func_values = N_Gauss exists --> extract value corresponding to 'gauss_point_in_patch'
      else
      {
        Func_values& func_values(patch_it->second);

        r = func_values[gauss_point_in_patch];
      }

      return r;
  }



  template <class IBASIS, unsigned int DIM>
  double
  pPoissonEquation<IBASIS,DIM>::a(const typename AggregatedFrame<IBASIS,DIM>::Index& lambda,
				  const typename AggregatedFrame<IBASIS,DIM>::Index& nu) const
  {
    // we seperately treat the entries from diagonal blocks and non-diagonal blocks
    return lambda.p() == nu.p() ? a_same_patches2(lambda, nu, degree_of_exactness_quadrature_a, min_res_quadrature_a) : a_different_patches2(lambda, nu, degree_of_exactness_quadrature_a, min_res_quadrature_a);

  }

  template <class IBASIS, unsigned int DIM>
  double
  pPoissonEquation<IBASIS,DIM>::norm_A() const
  {
    if (normA == 0.0)
    {
      typedef typename AggregatedFrame<IBASIS,DIM>::Index Index;
      std::set<Index> Lambda;
      const int j0 = frame_->j0();
      const int jmax = j0+1;
      for (Index lambda = FrameTL::first_generator<IBASIS,DIM,DIM,Frame>(frame_,j0);; ++lambda)
      {
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
  pPoissonEquation<IBASIS,DIM>::s_star() const
  {
    const double t = operator_order();
    const int n = DIM;
    const int dT = frame_->bases()[0]->primal_vanishing_moments(); // we assume to have the same 'kind'
                                                                   // of wavelets on each patch, so use
                                                                   // patch 0 as reference case
    const double gamma = frame_->bases()[0]->primal_regularity(); // = spline order - 0.5

    // cf. Manuel's thesis Theorem 5.1 and Remark 5.2
    return (n == 1
	    ? t+dT
	    : std::min((t+dT)/(double)n, (gamma-t)/double(n-1)));
  }



  template <class IBASIS, unsigned int DIM>
  double
  pPoissonEquation<IBASIS,DIM>::f(const typename AggregatedFrame<IBASIS,DIM>::Index& lambda) const
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
    const int N_Gauss = (degree_of_exactness_quadrature_f+1)/2;

    if (min_res_quadrature_f > supp.j)  // adjust resolution of composite quadrature if necessary
    {
      unsigned int diff = (1 << (min_res_quadrature_f - supp.j));
      for (unsigned int i = 0; i < DIM; i++)
      {
        supp.a[i] *= diff;
        supp.b[i] *= diff;
      }
      supp.j = min_res_quadrature_f;
    }

    //const double h = ldexp(1.0, -supp.j); // granularity for the quadrature
    const double h = 1.0 / (1 << supp.j); // granularity for the quadrature
    FixedArray1D<Array1D<double>,DIM> gauss_points, gauss_weights;
    for (unsigned int i = 0; i < DIM; i++)
    {
      gauss_points[i].resize(N_Gauss*(supp.b[i]-supp.a[i]));
      gauss_weights[i].resize(N_Gauss*(supp.b[i]-supp.a[i]));
      for (int patch = supp.a[i]; patch < supp.b[i]; patch++)
      {
	    for (int n = 0; n < N_Gauss; n++)
	    {
	    gauss_points[i][(patch-supp.a[i])*N_Gauss+n]
	      = h*(2*patch+1+GaussPoints[N_Gauss-1][n])/2.;
	    gauss_weights[i][(patch-supp.a[i])*N_Gauss+n]
	      = h*GaussWeights[N_Gauss-1][n];
	    }
	  }
    }


#if !PPE_USE_ALTERNATIVE_F
    //! compute f(\lambda) = \int_{\Omega} f(x) * \Psi_{\lambda}(x) dx
    FixedArray1D<Array1D<double>,DIM> v_values;
    // compute the point values of the integrand (where we use that it is a tensor product)
    for (unsigned int i = 0; i < DIM; i++)
    {
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

    while (true)
    {
      for (unsigned int i = 0; i < DIM; i++)
      {
	    x[i] = gauss_points[i][index[i]];
      }
      chart->map_point(x,x_patch);

      double share = ell_bvp_->f(x_patch) * chart->Gram_factor(x);
      for (unsigned int i = 0; i < DIM; i++)
      {
	    share *= gauss_weights[i][index[i]] * v_values[i][index[i]];
      }
      r += share;

      // "++index"
      bool exit = false;
      for (unsigned int i = 0; i < DIM; i++)
      {
	    if (index[i] == N_Gauss*(supp.b[i]-supp.a[i])-1)
	    {
	      index[i] = 0;
	      exit = (i == DIM-1);
	    }
	    else
	    {
	      index[i]++;
	      break;
	    }
      }
      if (exit) break;
    }


#else
 //! compute f(\lambda) = \int_{\Omega} f_ar(x) * \nabla \Psi_{\lambda}(x) dx

// compute point values of the integrand (where we use that it is a tensor product)
	FixedArray1D<Array1D<double>,DIM>
	  psi_lambda_values,     // values of the components of psi_lambda at gauss_points[i]
	  psi_lambda_der_values; // values of the 1st deriv. of the components of psi_lambda at gauss_points[i]

	for (unsigned int i = 0; i < DIM; i++)
	{
	   WaveletTL::evaluate(*(frame_->bases()[p]->bases()[i]), 0,
			  typename IBASIS::Index(lambda.j(),
						 lambda.e()[i],
						 lambda.k()[i],
						 frame_->bases()[p]->bases()[i]),
			  gauss_points[i], psi_lambda_values[i]);


	  WaveletTL::evaluate(*(frame_->bases()[p]->bases()[i]), 1,
		   typename IBASIS::Index(lambda.j(),
					  lambda.e()[i],
					  lambda.k()[i],
					  frame_->bases()[p]->bases()[i]),
		   gauss_points[i], psi_lambda_der_values[i]);
	}


	double grad_psi_lambda[DIM], weights;
    Vector<double> f_values;
    f_values.resize(2);

    // iterate over all points and sum up the integral shares
    int index[DIM]; // current multiindex for the point values
    for (unsigned int i = 0; i < DIM; i++)
    {
      index[i] = 0;
    }

    Point<DIM> x;
    Point<DIM> x_patch;

    while (true) {
      // current Gauss point x
      for (unsigned int i = 0; i < DIM; i++)
      {
	    x[i] = gauss_points[i][index[i]];
      }

      chart->map_point(x,x_patch);

      // product of current Gauss weights
      weights = 1.0;
      for (unsigned int i = 0; i < DIM; i++)
      {
        weights *= gauss_weights[i][index[i]];
      }

      // compute grad Psi_lambda
      for (unsigned int i = 0; i < DIM; i++)
	    {
	      grad_psi_lambda[i] = 1.0;
	      for (unsigned int s = 0; s < DIM; s++)
	      {
		    if (i == s)
		    {
		      grad_psi_lambda[i] *= psi_lambda_der_values[i][index[i]];
		    }
		    else
		    {
		      grad_psi_lambda[i] *= psi_lambda_values[s][index[s]];
		    }
	      }
	    }


      double share = 0.0;
      ell_bvp_->f_ar(x_patch, f_values);

      for (unsigned int i = 0; i < DIM; i++)
      {
	    share += (grad_psi_lambda[i] * f_values[i] * chart->Dkappa_inv(i, i, x));
      }
      share *= weights;
      share *= chart->Gram_factor(x);
      r += share;


      // "++index"
      bool exit = false;
      for (unsigned int i = 0; i < DIM; i++)
      {
	    if (index[i] == N_Gauss*(supp.b[i]-supp.a[i])-1)
	    {
	      index[i] = 0;
	      exit = (i == DIM-1);
	    }
	    else
	    {
	      index[i]++;
	      break;
	    }
      }
      if (exit) break;
    }
#endif

#ifdef TWO_D
    return r;
#endif

// #####################################################################################
// Attention! This is a bad hack! It is assumed that in case macro ONE_D is defined,
// the right hand side is the special functional defines on page 105 of Manuels thesis.
// We should use the class Functional instead!
// #####################################################################################
#ifdef ONE_D
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
  pPoissonEquation<IBASIS,DIM>::RHS
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
  pPoissonEquation<IBASIS,DIM>::RHS(const double eta, const int p,
        InfiniteVector<double, typename AggregatedFrame<IBASIS,DIM>::Index>& coeffs) const
  {
    cout.precision(12);
    coeffs.clear();

    if (fnorms_sqr_patch[p] < 1.0e-15)
      return;

    double coarsenorm(0);
    double bound(fnorms_sqr_patch[p] - eta*eta);
    typedef typename AggregatedFrame<IBASIS,DIM>::Index Index;
    typename Array1D<std::pair<Index, double> >::const_iterator it(fcoeffs_patch[p].begin());
    do {
      coarsenorm += it->second * it->second;
      //cout << it->first << " " << it->second << " " << coarsenorm  << " " << bound << endl;
      coeffs.set_coefficient(it->first, it->second);
      ++it;
    } while (it != fcoeffs_patch[p].end() && coarsenorm < bound);
  }


  template <class IBASIS, unsigned int DIM>
  void
  pPoissonEquation<IBASIS,DIM>::set_bvp(const PoissonBVP_Coeff<DIM>* bvp)
  {
    ell_bvp_ = bvp;


    clock_t begin_compute_diag = clock();

    // precomputation of the right-hand side up to the maximal level
    compute_diagonal();

    clock_t end_compute_diag = clock();
    time_consumption_of_compute_diagonal = double(end_compute_diag - begin_compute_diag) / CLOCKS_PER_SEC;


    clock_t begin_compute_rhs = clock();

    // precomputation of the diagonal up to the maximal level
    compute_rhs();

    clock_t end_compute_rhs = clock();
    time_consumption_of_compute_rhs = double(end_compute_rhs - begin_compute_rhs) / CLOCKS_PER_SEC;


  }

  template <class IBASIS, unsigned int DIM>
  void
  pPoissonEquation<IBASIS,DIM>::add_level (const Index& lambda,
					   InfiniteVector<double, Index>& w, const int j,
					   const double factor,
					   const int J,
					   const CompressionStrategy strategy) const
  {
    typedef std::list<Index> IntersectingList;
    IntersectingList nus;
    if (strategy == CDD1)
    {
      // compute all wavelets on level j, such that supp(psi_lambda) and supp(psi_nu) intersect
      intersecting_wavelets(basis(), lambda,
			    std::max(j, basis().j0()),
			    j == (basis().j0()-1),
			    nus);

      // traverse the matrix block and update the result
      const double d1 = D(lambda);
      for (typename IntersectingList::const_iterator it(nus.begin()), itend(nus.end());
	       it != itend; ++it)
      {
	    const double entry = a(*it, lambda) / (d1*D(*it));
	    w.add_coefficient(*it, entry * factor);
      }
    }
    else if (strategy == St04a)
    {
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
