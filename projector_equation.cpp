// implementation for projector_equation.h

#include <list>
#include <numerics/gauss_data.h>
#include <frame_index.h>
#include <frame_support.h>

using WaveletTL::CubeBasis;

#define _FRAMETL_ADAPTIVE_COMPUTATION 0
#define  GAUSS_POINTS_WEIGHTED 6

namespace FrameTL
{
  template <class IBASIS, unsigned int DIM>
  ProjectorEquation<IBASIS,DIM>::ProjectorEquation(const EllipticBVP<DIM>* ell_bvp,
							     const AggregatedFrame<IBASIS,DIM>* frame,
							     const int jmax, const bool weighted, const int patch = 0)
    : ell_bvp_(ell_bvp), frame_(frame), jmax_(jmax), weighted_(weighted), patch_(patch)
  {
    // precomputation of the right-hand side up to the maximal level
    compute_diagonal();

    // precomputation of the diagonal up to the maximal level
    compute_rhs();
  }

  template <class IBASIS, unsigned int DIM>
  double 
  ProjectorEquation<IBASIS,DIM>::D(const typename AggregatedFrame<IBASIS,DIM>::Index& lambda) const
  {

    return stiff_diagonal[lambda.number()];
  }

  template <class IBASIS, unsigned int DIM>
  void
  ProjectorEquation<IBASIS,DIM>::rescale(InfiniteVector<double,
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


  // we do not compute the RHS here
  template <class IBASIS, unsigned int DIM>
  void
  ProjectorEquation<IBASIS,DIM>::compute_rhs()
  {
  }

  // For the computation of the diagonal side, we include the option
  // to read in the data from file, because the computation of all coefficients
  // up to a high level may take a while.
  // The precomputation will only be performed in case the preprocessor macro
  // PRECOMP_DIAG is defined. Otherwise the coefficients are computed.
  // First we read the data into a sparse matrix
  // and copy it into stiff_diagonal afterwards.
  // We should write file io routines directly for the Infinite vector classes to get
  // rid of this hack.
  template <class IBASIS, unsigned int DIM>
  void
  ProjectorEquation<IBASIS,DIM>::compute_diagonal()
  {
    cout << "ProjectorEquation(): precompute diagonal of stiffness matrix..." << endl;

    SparseMatrix<double> diag(1,frame_->degrees_of_freedom());
    char filename[50];
    char matrixname[50];

    int d = IBASIS::primal_polynomial_degree();
    int dT = IBASIS::primal_vanishing_moments();
    
    // prepare filenames for 1D and 2D case
#ifdef ONE_D
    sprintf(filename, "%s%d%s%d", "stiff_diagonal_poisson_interval_lap07_d", d, "_dT", dT);
    sprintf(matrixname, "%s%d%s%d", "stiff_diagonal_poisson_1D_lap07_d", d, "_dT", dT);
#endif
#ifdef TWO_D
    sprintf(filename, "%s%d%s%d", "stiff_diagonal_poisson_lshaped_lap1_d", d, "_dT", dT);
    sprintf(matrixname, "%s%d%s%d", "stiff_diagonal_poisson_2D_lap1_d", d, "_dT", dT);
#endif

#ifndef PRECOMP_DIAG
    std::list<Vector<double>::size_type> indices;
    std::list<double> entries;
#endif
    
#ifdef PRECOMP_DIAG
    cout << "reading in diagonal of unpreconditioned stiffness matrix from file "
	 << filename << "..." << endl;
    diag.matlab_input(filename);
    cout << "...ready" << endl;
#endif


    stiff_diagonal.resize(frame_->degrees_of_freedom());
    for (int i = 0; i < frame_->degrees_of_freedom(); i++) {
#ifdef PRECOMP_DIAG
      stiff_diagonal[i] = diag.get_entry(0,i);
#endif
#ifndef PRECOMP_DIAG
      stiff_diagonal[i] = sqrt(a(*(frame_->get_wavelet(i)),*(frame_->get_wavelet(i))));


      // avoid zeros on main diagonal!!!
      // if (abs(stiff_diagonal[i]) < 10e-24) stiff_diagonal[i]=10e-15;


      indices.push_back(i);
      entries.push_back(stiff_diagonal[i]);
#endif
      //cout << stiff_diagonal[i] << " " << *(frame_->get_wavelet(i)) << endl;
    }

#ifndef PRECOMP_DIAG
    diag.set_row(0,indices, entries);
    diag.matlab_output(filename, matrixname, 1);
#endif

    cout << "... done, diagonal of stiffness matrix computed" << endl;
  }


  // integration in the one-d case
  template <class IBASIS, unsigned int DIM>
  inline
  double
  ProjectorEquation<IBASIS,DIM>:: integrate(const Index1D<IBASIS>& lambda,
						 const Index1D<IBASIS>& mu,
						 const int N_Gauss,
						 const int dir,
						 const typename CubeBasis<IBASIS,DIM>::Support* supp_lambda,
						 const typename CubeBasis<IBASIS,DIM>::Support* supp_mu) const
  {
    
     double res = 0;
     
	// compute 1D irregular grid
	Array1D<double> irregular_grid;


	// Check whether the supports of the two functions intersect and compute the intersection
	// of the two singular supports. We obtain a non-uniform grid with respect to which the
	// present integrand is a piecewise polynomial; see �6.3, page 57-62 in Manuel diploma thesis.
	bool b = intersect_supports_1D<IBASIS,DIM>(*frame_, lambda, mu, supp_lambda, supp_mu, dir, irregular_grid);
	if (!b)
	  return 0.0;


	Array1D<double> gauss_points_la, gauss_points_mu, gauss_weights,
	  values_lambda, values_mu;
	
	// Setup gauss knots and weights for each of the little pieces where the integrand is smooth.
	// The gauss knots and weights for the interval [-1,1] are given in <numerics/gauss_data.h>.
 	gauss_points_la.resize(N_Gauss*(irregular_grid.size()-1));
 	gauss_points_mu.resize(N_Gauss*(irregular_grid.size()-1));
 	gauss_weights.resize(N_Gauss*(irregular_grid.size()-1));
  	for (unsigned int k = 0; k < irregular_grid.size()-1; k++)
 	  for (int n = 0; n < N_Gauss; n++) {
 	    gauss_points_la[ k*N_Gauss+n  ]
 	      = 0.5 * (irregular_grid[k+1]-irregular_grid[k]) * (GaussPoints[N_Gauss-1][n]+1)
 	      + irregular_grid[k];
	    
	    gauss_weights[ k*N_Gauss+n ]
 	      = (irregular_grid[k+1]-irregular_grid[k])*GaussWeights[N_Gauss-1][n];
 	  }
	
	IBASIS* basis1D_lambda = frame_->bases()[lambda.p()]->bases()[dir];
	IBASIS* basis1D_mu     = frame_->bases()[mu.p()]->bases()[dir];
	
	const Chart<DIM>* chart_la = frame_->atlas()->charts()[lambda.p()];
	const Chart<DIM>* chart_mu = frame_->atlas()->charts()[mu.p()];
	

	//	cout << "deriv_lambda = " << lambda.derivative() << endl;
	WaveletTL::evaluate(*(basis1D_lambda), lambda.derivative(),
			    lambda.index(),
			    gauss_points_la, values_lambda);
	
	const Array1D<double> gouss_points_mu;

	// setup mapped gauss points
	for (unsigned int i = 0; i < gauss_points_la.size(); i++) {
	  gauss_points_mu[i] = chart_mu->map_point_inv(chart_la->map_point(gauss_points_la[i], dir), dir);
	}

	//	cout << "deriv_mu = " << mu.derivative() << endl;
	WaveletTL::evaluate(*(basis1D_mu), mu.derivative(),
			    mu.index(),
			    gauss_points_mu, values_mu);
		
	//cout << "doing the job.." << endl;

        double factor;
        double x_i;
        Point<2> p_i;
	for (unsigned int i = 0; i < values_lambda.size(); i++)
        {
          if (!weighted_) factor = 1.0;
          else
          {
            // calculate mapped point

            x_i = chart_la->map_point(gauss_points_la[i],dir);
            factor = FrameTL::C2(x_i, 0.3, 0.7);

            if (patch_ == 1) factor = 1 - factor;

          }
                  
	  res += gauss_weights[i] * values_lambda[i] * values_mu[i] * factor;
        }

    return res;
  }


  template <class IBASIS, unsigned int DIM>
  double
  ProjectorEquation<IBASIS,DIM>::a(const typename AggregatedFrame<IBASIS,DIM>::Index& la,
					const typename AggregatedFrame<IBASIS,DIM>::Index& nu) const
  {
    double r = 0.0;

    typedef WaveletTL::CubeBasis<IBASIS,DIM> CUBEBASIS;

    const int jla = (frame_->all_supports[la.number()]).j;
    const int jnu = (frame_->all_supports[nu.number()]).j;

    switched_indices = (jnu > jla);
    
    const Index* lambda = switched_indices ? &nu : &la;
    const Index* mu     = switched_indices ? &la : &nu;

    const typename CUBEBASIS::Support* supp_lambda =
      switched_indices ? &(frame_->all_supports[nu.number()]) : &(frame_->all_supports[la.number()]);

    const typename CUBEBASIS::Support* supp_mu =
      switched_indices ? &(frame_->all_supports[la.number()]) : &(frame_->all_supports[nu.number()]);

    //const int N_Gauss = 3;
    //const int N_Gauss = IBASIS::primal_polynomial_degree();


    // use more Gauss-points in case the integrand is weighted
    const int N_Gauss = (weighted_) ? GAUSS_POINTS_WEIGHTED : IBASIS::primal_polynomial_degree();


    typedef typename IBASIS::Index Index_1D;
    
    const Chart<DIM>* chart_la = frame_->atlas()->charts()[lambda->p()];
    const Chart<DIM>* chart_mu = frame_->atlas()->charts()[mu->p()];

    // now follows integration of q(x)*psi_lambda(x)*psi_mu(x)

#ifdef ONE_D

    double s = 1.;
    // loop over spatial direction     
    for (int i = 0; i < (int) DIM; i++) {
      Index1D<IBASIS> i1(IntervalIndex<IBASIS> (
						lambda->j(),lambda->e()[i],lambda->k()[i],
						frame_->bases()[lambda->p()]->bases()[i]
						),
			 lambda->p(),i,0
			 );
      Index1D<IBASIS> i2(IntervalIndex<IBASIS> (mu->j(),mu->e()[i],mu->k()[i],
						frame_->bases()[mu->p()]->bases()[i]
						),
			 mu->p(),i,0
			 );
      s *= integrate(i1, i2, N_Gauss, i, supp_lambda, supp_mu);
    }
    
#else
    double s;
    if (!weighted_)
    {
        s = 1.;
        // loop over spatial direction
        for (int i = 0; i < (int) DIM; i++) {
          Index1D<IBASIS> i1(IntervalIndex<IBASIS> (
                                                    lambda->j(),lambda->e()[i],lambda->k()[i],
                                                    frame_->bases()[lambda->p()]->bases()[i]
                                                    ),
                             lambda->p(),i,0
                             );
          Index1D<IBASIS> i2(IntervalIndex<IBASIS> (mu->j(),mu->e()[i],mu->k()[i],
                                                    frame_->bases()[mu->p()]->bases()[i]
                                                    ),
                             mu->p(),i,0
                             );
          s *= integrate(i1, i2, N_Gauss, i, supp_lambda, supp_mu);
        }
    }
    else
    {

       s=0;

       Array1D<double> x_values, y_values, x_points, y_points;
       for (int dir=0; dir <= 1; dir++) {

           Index1D<IBASIS> lambda_int(IntervalIndex<IBASIS> (
                                                        lambda->j(),lambda->e()[dir],lambda->k()[dir],
                                                        frame_->bases()[lambda->p()]->bases()[dir]
                                                        ),
                                 lambda->p(),dir,0
                                 );
           Index1D<IBASIS> mu_int(IntervalIndex<IBASIS> (mu->j(),mu->e()[dir],mu->k()[dir],
                                                        frame_->bases()[mu->p()]->bases()[dir]
                                                        ),
                                 mu->p(),dir,0
                                 );


            // compute 1D irregular grid
            Array1D<double> irregular_grid;

            // Check whether the supports of the two functions intersect and compute the intersection
            // of the two singular supports. We obtain a non-uniform grid with respect to which the
            // present integrand is a piecewise polynomial; see �6.3, page 57-62 in Manuel diploma thesis.
            bool b = intersect_supports_1D<IBASIS,DIM>(*frame_, lambda_int, mu_int, supp_lambda, supp_mu, dir, irregular_grid);
            if (!b)
              return 0.0;


            Array1D<double> gauss_points_la, gauss_points_mu, gauss_weights,
              values_lambda, values_mu;



            // Setup gauss knots and weights for each of the little pieces where the integrand is smooth.
            // The gauss knots and weights for the interval [-1,1] are given in <numerics/gauss_data.h>.
            gauss_points_la.resize(N_Gauss*(irregular_grid.size()-1));
            gauss_points_mu.resize(N_Gauss*(irregular_grid.size()-1));
            gauss_weights.resize(N_Gauss*(irregular_grid.size()-1));
            for (unsigned int k = 0; k < irregular_grid.size()-1; k++)
              for (int n = 0; n < N_Gauss; n++) {
                gauss_points_la[ k*N_Gauss+n ]
                  = 0.5 * (irregular_grid[k+1]-irregular_grid[k]) * (GaussPoints[N_Gauss-1][n]+1)
                  + irregular_grid[k];

                gauss_weights[ k*N_Gauss+n]
                  = (irregular_grid[k+1]-irregular_grid[k])*GaussWeights[N_Gauss-1][n];
              }

            IBASIS* basis1D_lambda = frame_->bases()[lambda_int.p()]->bases()[dir];
            IBASIS* basis1D_mu     = frame_->bases()[mu_int.p()]->bases()[dir];

            const Chart<DIM>* chart_la = frame_->atlas()->charts()[lambda_int.p()];
            const Chart<DIM>* chart_mu = frame_->atlas()->charts()[mu_int.p()];


            //	cout << "deriv_lambda = " << lambda.derivative() << endl;
            WaveletTL::evaluate(*(basis1D_lambda), lambda_int.derivative(),
                                lambda_int.index(),
                                gauss_points_la, values_lambda);

            // setup mapped gauss points
            for (unsigned int i = 0; i < gauss_points_la.size(); i++) {
              gauss_points_mu[i] = chart_mu->map_point_inv(chart_la->map_point(gauss_points_la[i], dir), dir);
            }

            WaveletTL::evaluate(*(basis1D_mu), mu_int.derivative(),
                                mu_int.index(),
                                gauss_points_mu, values_mu);

            if(dir == 0)
            {
                x_values.resize(N_Gauss*(irregular_grid.size()-1));
                x_points.resize(N_Gauss*(irregular_grid.size()-1));

    //            cout << "x_values: " << x_values.size() << endl;
    //            cout << "gauss_weights: " << gauss_weights.size() << endl;
    //            cout << "values_lambda: " << values_lambda.size() << endl;
    //            cout << "gauss_points_la: " << gauss_points_la.size() << endl;
    //            cout << "values_mu: " << values_mu.size() << endl;

                for (unsigned int i = 0; i<x_values.size(); i++){

                    x_values[i] = gauss_weights[i]*values_lambda[i]*values_mu[i];
                    x_points[i] = chart_la->map_point(gauss_points_la[i],dir);
                }

            }
            else
            {
                y_values.resize(N_Gauss*(irregular_grid.size()-1));
                y_points.resize(N_Gauss*(irregular_grid.size()-1));
                for (unsigned int  i = 0; i<y_values.size(); i++){
                    y_values[i] = gauss_weights[i]*values_lambda[i]*values_mu[i];
                    y_points[i] = chart_la->map_point(gauss_points_la[i],dir);
                }
            }
            //cout << "doing the job.." << endl;
       }

       // sum up weighted values

//       cout << *lambda << " ; " << *mu << endl;
//       for(unsigned int  i=0; i < x_values.size(); i++)
//       {
//           cout << i << ": " << x_points[i] << endl;
//       }
//       for(unsigned int  i=0; i < y_values.size(); i++)
//       {
//           cout << i << ": " << y_points[i] << endl;
//       }

       for(unsigned int  i=0; i<x_values.size(); i++)
            for(unsigned int j=0; j<y_values.size(); j++)
            {
                Point<2> p;
                p[0] = x_points[i];
                p[1] = y_points[j];
                double factor = FrameTL::sigmaL(p,'2');
                if (patch_ == 1) factor = 1.0 -factor;
                factor = (weighted_) ? factor : 1.0;
                s += x_values[i]*y_values[j]*factor;
            }


    }

#endif

    double tmp1 = 1., tmp2 = 1.;
    
    for (unsigned int i = 0; i < DIM; i++) {
      tmp1 *= chart_la->a_i(i);
      tmp2 *= chart_mu->a_i(i);
    }
    tmp1 = sqrt(fabs(tmp1)) / sqrt(fabs(tmp2));

  
    s *= tmp1;
    
    r += s;

    assert(r == r);
    
    if (abs(r) > 10e3) cout << "a(" << la << ", " << nu << ")= " << r << endl;

    return r;

  }
  
  template <class IBASIS, unsigned int DIM>
  double
  ProjectorEquation<IBASIS,DIM>::norm_A() const
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
  ProjectorEquation<IBASIS,DIM>::s_star() const
  {

    const double t = operator_order();
    const int n = DIM;
    const int dT = frame_->bases()[0]->primal_vanishing_moments(); // we assume to have the same 'kind'
                                                                   // of wavelets on each patch, so use
                                                                   // patch 0 as reference case
    const double gamma = frame_->bases()[0]->primal_regularity();
    
    // cf. Manuel's thesis Theorem 5.1 and Remark 5.2
    return (n == 1
	    ? t+dT 
	    : std::min((t+dT)/(double)n, (gamma-t)/double(n-1)));
  }

  template <class IBASIS, unsigned int DIM>
  double
  ProjectorEquation<IBASIS,DIM>::f(const typename AggregatedFrame<IBASIS,DIM>::Index& lambda) const
  {
      return 0;
  }

  template <class IBASIS, unsigned int DIM>
  void
  ProjectorEquation<IBASIS,DIM>::RHS
  (const double eta,
   InfiniteVector<double, typename AggregatedFrame<IBASIS,DIM>::Index>& coeffs) const
  {
    
  }

  template <class IBASIS, unsigned int DIM>
  void
  ProjectorEquation<IBASIS,DIM>::RHS
  (const double eta,
   const int p,
   InfiniteVector<double, typename AggregatedFrame<IBASIS,DIM>::Index>& coeffs) const
  {
  }



  template <class IBASIS, unsigned int DIM>
  void
  ProjectorEquation<IBASIS,DIM>::set_bvp(const EllipticBVP<DIM>* bvp)
  {
    ell_bvp_ = bvp;
    compute_diagonal();
    compute_rhs();

  }

  template <class IBASIS, unsigned int DIM>
  void
  ProjectorEquation<IBASIS,DIM>::add_level (const Index& lambda,
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
#ifdef BUGFIX
          assert(d1*D(*it) != 0);
#endif
	const double entry = a(*it, lambda) / (d1*D(*it));
#ifdef BUGFIX
        assert(entry == entry);
#endif
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
            
#ifdef BUGFIX
          assert(d1*D(*it) != 0);
#endif
	  const double entry = a(*it, lambda) / (d1*D(*it));

	  w.add_coefficient(*it, entry * factor);
	}
    }
  }
  
}


