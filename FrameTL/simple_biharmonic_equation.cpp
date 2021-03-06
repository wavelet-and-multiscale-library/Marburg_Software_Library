// implementation for biharmonic_equation.h

#include <list>
#include <utils/array1d.h>
#include <numerics/gauss_data.h>
#include <algebra/sparse_matrix.h>
#include <numerics/eigenvalues.h>
#include <frame_index.h>
#include <frame_support.h>
//#include <cuba.h>

//#include <cube/cube_support.h>

using WaveletTL::CubeBasis;
using MathTL::Array1D;

namespace FrameTL
{
  template <class IBASIS, unsigned int DIM>
  SimpleBiharmonicEquation<IBASIS,DIM>::SimpleBiharmonicEquation(const Functional<IBASIS,DIM>* rhs,
								 const AggregatedFrame<IBASIS,DIM>* frame,
								 const int jmax)
    : rhs_(rhs), frame_(frame),jmax_(jmax)
  {
    // precomputation of the right-hand side up to the maximal level
    compute_diagonal();
    // precomputation of the diagonal up to the maximal level
    compute_rhs();
  }

  template <class IBASIS, unsigned int DIM>
  double 
  SimpleBiharmonicEquation<IBASIS,DIM>::D(const typename AggregatedFrame<IBASIS,DIM>::Index& lambda) const
  {
    return stiff_diagonal[lambda.number()];
  }


  // For the computation of the right-hand side, we include the option
  // to read in the data from file, because the computation of all coefficients
  // up to a high level may take a while.
  // The precomputation will only be performed in case the preprocessor macro
  // PRECOMP_RHS is defined. Otherwise the coefficients are computed.
  // First we read the data into a sparse matrix
  // and copy it into fcoeffs and fcoeffs_patch afterwards.
  // We should write file io routines directly for the Infinite vector classes to get
  // rid of this hack.
  template <class IBASIS, unsigned int DIM>
  void
  SimpleBiharmonicEquation<IBASIS,DIM>::compute_rhs()
  {
    cout << "BiharmonicEquation(): precompute right-hand side ..." << endl;
    typedef AggregatedFrame<IBASIS,DIM> Frame;
    typedef typename Frame::Index Index;

    // the sparse matrix in which we shall put the right-hand side coefficients
    SparseMatrix<double> rhs(1,frame_->degrees_of_freedom());
    char filename[50];
    char matrixname[50];

    
    
    // prepare filenames for 2D case
#ifdef TWO_D
    int d = IBASIS::primal_polynomial_degree();
    int dT = IBASIS::primal_vanishing_moments();
    sprintf(filename, "%s%d%s%d", "rhs_biharm_lshaped_lap1_d", d, "_dT", dT);
    sprintf(matrixname, "%s%d%s%d", "rhs_biharm_2D_lap1_d", d, "_dT", dT);
#endif

    // initialize array fnorms_sqr_patch
    fnorms_sqr_patch.resize(frame_->n_p());
    for (unsigned int i = 0; i < fnorms_sqr_patch.size(); i++)
      fnorms_sqr_patch[i] = 0.;
     
    
    // precompute the right-hand side on a fine level
    InfiniteVector<double,Index> fhelp;
    Array1D<InfiniteVector<double,Index> > fhelp_patch(frame_->n_p());

#ifdef PRECOMP_RHS
    // we read in the right hand side from file
    // we assume that ot had been precomputed on a sufficiently high level
    cout << "reading in right hand side from file " << filename << "..." << endl;
    rhs.matlab_input(filename);
    cout << "...ready" << endl;
#endif
#ifndef PRECOMP_RHS
    std::list<Vector<double>::size_type> indices;
    std::list<double> entries;
#endif

    // loop over all wavelets between minimal and maximal level
    for (int i = 0; i < frame_->degrees_of_freedom(); i++)
      {
#ifdef PRECOMP_RHS
	double coeff = rhs.get_entry(0,i);
#else	
	// computation of one right-hand side coefficient
	double coeff = f(*(frame_->get_wavelet(i)))/D(*(frame_->get_wavelet(i)));
	if (fabs(coeff)>1e-15) {
	  indices.push_back(i);
	  entries.push_back(coeff);
	  //rhs.set_entry(0, i, coeff);
	}
#endif
	// put the coefficient into an InfiniteVector and successively
	// compute the squared \ell_2 norm
	if (fabs(coeff)>1e-15) {
	  fhelp.set_coefficient(*(frame_->get_wavelet(i)), coeff);
	  fhelp_patch[frame_->get_wavelet(i)->p()].set_coefficient(*(frame_->get_wavelet(i)), coeff);

	  fnorms_sqr_patch[frame_->get_wavelet(i)->p()] += coeff*coeff;
	  if (i % 100 == 0)
	    cout << *(frame_->get_wavelet(i)) << " " << coeff << endl;
	}
      }

    // write the right-hand side into file in case ot has just been computed
#ifndef PRECOMP_RHS
#ifndef MSL_GUI
    rhs.set_row(0, indices, entries);
    // write right hand side to file
    cout << "writing right hand side into file..." << filename << "..." << endl;
    rhs.matlab_output(filename, matrixname, 1);
    cout << "...ready" << endl;
#endif
#endif

    fnorm_sqr = l2_norm_sqr(fhelp);
    for (unsigned int i = 0; i < fnorms_sqr_patch.size(); i++)
      cout << fnorms_sqr_patch[i] << endl;
    

    cout << "norm rhs sqr = " << fnorm_sqr << endl;
    cout << "... done, all integrals for right-hand side computed" << endl;

    // sort the coefficients into fcoeffs
    fcoeffs.resize(0); // clear eventual old values
    fcoeffs.resize(fhelp.size());
    unsigned int id(0);
    for (typename InfiniteVector<double,Index>::const_iterator it(fhelp.begin()), itend(fhelp.end());
	 it != itend; ++it, ++id)
      fcoeffs[id] = std::pair<Index,double>(it.index(), *it);
    sort(fcoeffs.begin(), fcoeffs.end(), typename InfiniteVector<double,Index>::decreasing_order());

    //     for (int i = 0; i < fcoeffs.size(); i++) {
    //       cout << fcoeffs[i].first << " " << fcoeffs[i].second << endl;
    //     }


    fcoeffs_patch.resize(frame_->n_p());
    for (int i = 0; i < frame_->n_p(); i++) {
      //fcoeffs_patch[i].resize(0); // clear eventual old values
      fcoeffs_patch[i].resize(fhelp_patch[i].size());
      id = 0;
      for (typename InfiniteVector<double,Index>::const_iterator it(fhelp_patch[i].begin()), itend(fhelp_patch[i].end());
	   it != itend; ++it, ++id) {
	(fcoeffs_patch[i])[id] = std::pair<Index,double>(it.index(), *it);
      }
      sort(fcoeffs_patch[i].begin(), fcoeffs_patch[i].end(), typename InfiniteVector<double,Index>::decreasing_order());
    } 
  }


  //   template <class IBASIS, unsigned int DIM>
  //   void
  //   SimpleBiharmonicEquation<IBASIS,DIM>::compute_rhs(bool compute)
  //   {
  //     cout << "BiharmonicEquation(): precompute right-hand side ..." << endl;
    
  //     typedef AggregatedFrame<IBASIS,DIM> Frame;
  //     typedef typename Frame::Index Index;

  //     // precompute the right-hand side on a fine level
  //     InfiniteVector<double,Index> fhelp;
  //     const int j0   = frame_->j0();
  //     if(compute)
  //       {    
  // 	for (Index lambda(FrameTL::first_generator<IBASIS,DIM,DIM,Frame>(frame_,j0));; ++lambda)
  // 	  {
  // 	    const double coeff = f(lambda)/D(lambda);
  // 	    //cout << "f(lambda):" << f(lambda) << endl;
  // 	    //cout << lambda << " " << coeff << endl; 
  // 	    if (fabs(coeff)>1e-15) {
  // 	      fhelp.set_coefficient(lambda, coeff);
  // 	    }
  // 	    if (lambda == last_wavelet<IBASIS,DIM,DIM,Frame>(frame_,jmax_))
  // 	      break;
  // 	    //cout << lambda << endl;
  // 	  }
  //       }
  // #if 0
  //     else 
  //       {
  // 	 Vector<double> f;
  // 	 Rhs_input(f);

  // 	 for(unsigned int i=0; i<f.size(); i++)
  // 	  {
  // 	    Index lambda=FrameIndex<IBASIS, DIM>(i, frame_);
  // 	    fhelp.set_coefficient(lambda, f[i]);
  // 	  }
  //       }
  // #endif
  //     fnorm_sqr = l2_norm_sqr(fhelp);
      
  //     cout << "... done, all integrals for right-hand side computed" << endl;

  //     // sort the coefficients into fcoeffs
  //     fcoeffs.resize(0); // clear eventual old values
  //     fcoeffs.resize(fhelp.size());
  //     unsigned int id(0);
  //     for (typename InfiniteVector<double,Index>::const_iterator it(fhelp.begin()), itend(fhelp.end());
  // 	 it != itend; ++it, ++id)
  //       fcoeffs[id] = std::pair<Index,double>(it.index(), *it);
  //     sort(fcoeffs.begin(), fcoeffs.end(), typename InfiniteVector<double,Index>::decreasing_order());
  //   }

  template <class IBASIS, unsigned int DIM>
  void
  SimpleBiharmonicEquation<IBASIS,DIM>::compute_diagonal()
  {
    cout << "SimpleBiharmonicEquation(): precompute diagonal of stiffness matrix ..." << endl;

    stiff_diagonal.resize(frame_->degrees_of_freedom());
    for (int i = 0; i < frame_->degrees_of_freedom(); i++) {
      stiff_diagonal[i] = sqrt(a(*(frame_->get_wavelet(i)),*(frame_->get_wavelet(i))));
      if (i%1000 == 0)
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

     // If the dimension is larger that just 1, it makes sense to store the one dimensional
     // integrals arising when we make use of the tensor product structure. This costs quite
     // some memory, but really speeds up the algorithm!
#ifdef TWO_D
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
#endif

   
    Array1D<double> gauss_points_la, gauss_points_mu, gauss_weights,
      values_lambda, values_mu;

    // Setup gauss knots and weights for each of the little pieces where the integrand is smooth.
    // The gauss knots and weights for the interval [-1,1] are given in <numerics/gauss_data.h>.
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

// in the 2D case store the calculated value
#ifdef TWO_D
	typedef typename Column1D::value_type value_type;
	it = col.insert(lb, value_type(mu, res));
      }
    else {
      res = it->second;
    }
#endif

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

    //typedef typename CUBEBASIS::Index CubeIndex;
 
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
    const int N_Gauss = IBASIS::primal_polynomial_degree();

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

    //typedef typename IBASIS::Index Index_1D;
    
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
  void
  SimpleBiharmonicEquation<IBASIS,DIM>::RHS
  (const double eta,
   const int p,
   InfiniteVector<double, typename AggregatedFrame<IBASIS,DIM>::Index>& coeffs) const
  {
    cout.precision(12);
    coeffs.clear();
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
  SimpleBiharmonicEquation<IBASIS,DIM>::f(const typename AggregatedFrame<IBASIS,DIM>::Index& lambda) const
  {
    // evaluate righthand side functional
    return rhs_->evaluate(lambda);
  }

}
