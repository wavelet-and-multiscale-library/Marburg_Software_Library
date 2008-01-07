// implementation for cube_equation.h

namespace WaveletTL
{
  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  CubeEquation<IBASIS,DIM,CUBEBASIS>::CubeEquation(const EllipticBVP<DIM>* bvp,
						   const FixedArray1D<bool,2*DIM>& bc)
    : bvp_(bvp), basis_(bc), normA(0.0), normAinv(0.0)
  {
    compute_rhs();
    const int jmax = 5; // for a first quick hack
    basis_.set_jmax(jmax);
  }

  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  CubeEquation<IBASIS,DIM,CUBEBASIS>::CubeEquation(const EllipticBVP<DIM>* bvp,
						   const FixedArray1D<int,2*DIM>& bc)
    : bvp_(bvp), basis_(bc), normA(0.0), normAinv(0.0)
  {
    compute_rhs();
    const int jmax = 5; // for a first quick hack
    basis_.set_jmax(jmax);
 
  }

  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  CubeEquation<IBASIS,DIM,CUBEBASIS>::CubeEquation(const CubeEquation& eq)
    : bvp_(eq.bvp_), basis_(eq.basis_),
      fcoeffs(eq.fcoeffs), fnorm_sqr(eq.fnorm_sqr),
      normA(eq.normA), normAinv(eq.normAinv)
  {
    const int jmax = 5; // for a first quick hack
    basis_.set_jmax(jmax);
  }

  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  void
  CubeEquation<IBASIS,DIM,CUBEBASIS>::compute_rhs()
  {
    cout << "CubeEquation(): precompute right-hand side..." << endl;

    typedef typename WaveletBasis::Index Index;

    // precompute the right-hand side on a fine level
    InfiniteVector<double,Index> fhelp;
    const int j0   = basis().j0();
    const int jmax = 5; // for a first quick hack
    for (Index lambda(basis_.first_generator(j0));; ++lambda)
      {
	const double coeff = f(lambda)/D(lambda);
	if (fabs(coeff)>1e-15)
	  fhelp.set_coefficient(lambda, coeff);
  	if (lambda == basis_.last_wavelet(jmax))
	  break;
      }
    fnorm_sqr = l2_norm_sqr(fhelp);

    cout << "... done, sort the entries in modulus..." << endl;

    // sort the coefficients into fcoeffs
    fcoeffs.resize(0); // clear eventual old values
    fcoeffs.resize(fhelp.size());
    unsigned int id(0);
    for (typename InfiniteVector<double,Index>::const_iterator it(fhelp.begin()), itend(fhelp.end());
	 it != itend; ++it, ++id)
      fcoeffs[id] = std::pair<Index,double>(it.index(), *it);
    sort(fcoeffs.begin(), fcoeffs.end(), typename InfiniteVector<double,Index>::decreasing_order());
    cout << "... done, all integrals for right-hand side computed!" << endl;
  }
  
  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  inline
  double
  CubeEquation<IBASIS,DIM,CUBEBASIS>::D(const typename WaveletBasis::Index& lambda) const
  {
//     return ldexp(1.0, lambda.j());
    return sqrt(a(lambda, lambda));
  }

  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  inline
  double
  CubeEquation<IBASIS,DIM,CUBEBASIS>::a(const typename WaveletBasis::Index& lambda,
					const typename WaveletBasis::Index& mu) const
  {
    return a(lambda, mu, IBASIS::primal_polynomial_degree()*IBASIS::primal_polynomial_degree());
  }
  
  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  double
  CubeEquation<IBASIS,DIM,CUBEBASIS>::a(const typename WaveletBasis::Index& lambda,
					const typename WaveletBasis::Index& mu,
					const unsigned int p) const
  {
    // a(u,v) = \int_Omega [a(x)grad u(x)grad v(x)+q(x)u(x)v(x)] dx

    double r = 0;

    // first decide whether the supports of psi_lambda and psi_mu intersect
    typedef typename CUBEBASIS::Support Support;
    Support supp;
    
    if (intersect_supports(basis_, lambda, mu, supp))
      {
	// setup Gauss points and weights for a composite quadrature formula:
	const int N_Gauss = (p+1)/2;
	const double h = ldexp(1.0, -supp.j); // granularity for the quadrature
	FixedArray1D<Array1D<double>,DIM> gauss_points, gauss_weights;
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
	
	// compute point values of the integrand (where we use that it is a tensor product)
	FixedArray1D<Array1D<double>,DIM>
	  psi_lambda_values,     // values of the components of psi_lambda at gauss_points[i]
	  psi_mu_values,         // -"-, for psi_mu
	  psi_lambda_der_values, // values of the 1st deriv. of the components of psi_lambda at gauss_points[i]
	  psi_mu_der_values;     // -"-, for psi_mu
	for (unsigned int i = 0; i < DIM; i++) {
	  evaluate(*basis_.bases()[i], 0,
		   typename IBASIS::Index(lambda.j(),
					  lambda.e()[i],
					  lambda.k()[i],
					  basis_.bases()[i]),
		   gauss_points[i], psi_lambda_values[i]);
	  evaluate(*basis_.bases()[i], 1,
		   typename IBASIS::Index(lambda.j(),
					  lambda.e()[i],
					  lambda.k()[i],
					  basis_.bases()[i]),
		   gauss_points[i], psi_lambda_der_values[i]);
	  evaluate(*basis_.bases()[i], 0,
		   typename IBASIS::Index(mu.j(),
					  mu.e()[i],
					  mu.k()[i],
					  basis_.bases()[i]),
		   gauss_points[i], psi_mu_values[i]);
	  evaluate(*basis_.bases()[i], 1,
		   typename IBASIS::Index(mu.j(),
					  mu.e()[i],
					  mu.k()[i],
					  basis_.bases()[i]),
		   gauss_points[i], psi_mu_der_values[i]);
	}
	
	// iterate over all points and sum up the integral shares
	int index[DIM]; // current multiindex for the point values
	for (unsigned int i = 0; i < DIM; i++)
	  index[i] = 0;
	
	Point<DIM> x;
	const double ax = bvp_->constant_coefficients() ? bvp_->a(x) : 0.0;
	const double qx = bvp_->constant_coefficients() ? bvp_->q(x) : 0.0;
	double grad_psi_lambda[DIM], grad_psi_mu[DIM], weights;
	if (bvp_->constant_coefficients()) {
	  while (true) {
	    for (unsigned int i = 0; i < DIM; i++)
	      x[i] = gauss_points[i][index[i]];
	    
	    // product of current Gauss weights
	    weights = 1.0;
	    for (unsigned int i = 0; i < DIM; i++)
	      weights *= gauss_weights[i][index[i]];
	    
	    // compute the share a(x)(grad psi_lambda)(x)(grad psi_mu)(x)
	    for (unsigned int i = 0; i < DIM; i++) {
	      grad_psi_lambda[i] = 1.0;
	      grad_psi_mu[i] = 1.0;
	      for (unsigned int s = 0; s < DIM; s++) {
		if (i == s) {
		  grad_psi_lambda[i] *= psi_lambda_der_values[i][index[i]];
		  grad_psi_mu[i]     *= psi_mu_der_values[i][index[i]];
		} else {
		  grad_psi_lambda[i] *= psi_lambda_values[s][index[s]];
		  grad_psi_mu[i] *= psi_mu_values[s][index[s]];
		}
	      }
	    }
	    double share = 0;
	    for (unsigned int i = 0; i < DIM; i++)
	      share += grad_psi_lambda[i]*grad_psi_mu[i];
 	    r += ax * weights * share;
	    
	    // compute the share q(x)psi_lambda(x)psi_mu(x)
 	    share = qx * weights;
	    for (unsigned int i = 0; i < DIM; i++)
	      share *= psi_lambda_values[i][index[i]] * psi_mu_values[i][index[i]];
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
	} else {
	  while (true) {
	    for (unsigned int i = 0; i < DIM; i++)
	      x[i] = gauss_points[i][index[i]];
	    
	    // product of current Gauss weights
	    weights = 1.0;
	    for (unsigned int i = 0; i < DIM; i++)
	      weights *= gauss_weights[i][index[i]];
	    
	    // compute the share a(x)(grad psi_lambda)(x)(grad psi_mu)(x)
	    for (unsigned int i = 0; i < DIM; i++) {
	      grad_psi_lambda[i] = 1.0;
	      grad_psi_mu[i] = 1.0;
	      for (unsigned int s = 0; s < DIM; s++) {
		if (i == s) {
		  grad_psi_lambda[i] *= psi_lambda_der_values[i][index[i]];
		  grad_psi_mu[i]     *= psi_mu_der_values[i][index[i]];
		} else {
		  grad_psi_lambda[i] *= psi_lambda_values[s][index[s]];
		  grad_psi_mu[i] *= psi_mu_values[s][index[s]];
		}
	      }
	    }
	    double share = 0;
	    for (unsigned int i = 0; i < DIM; i++)
	      share += grad_psi_lambda[i]*grad_psi_mu[i];
	    r += bvp_->a(x) * weights * share;
	    
	    // compute the share q(x)psi_lambda(x)psi_mu(x)
	    share = bvp_->q(x) * weights;
	    for (unsigned int i = 0; i < DIM; i++)
	      share *= psi_lambda_values[i][index[i]] * psi_mu_values[i][index[i]];
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
	}
	
      }

    return r;
  }
  
  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  double
  CubeEquation<IBASIS,DIM,CUBEBASIS>::f(const typename WaveletBasis::Index& lambda) const
  {
    // f(v) = \int_0^1 g(t)v(t) dt


    double r = 1.;
    for (unsigned int i = 0; i < DIM; i++)
      r *= evaluate(*basis_.bases()[i], 0,
		    typename IBASIS::Index(lambda.j(),
					   lambda.e()[i],
					   lambda.k()[i],
					   basis_.bases()[i]),
		    0.5);

    return r;


#if 0
    double r = 0;

    // first compute supp(psi_lambda)
    typename CUBEBASIS::Support supp;
    support(basis_, lambda, supp);

    // setup Gauss points and weights for a composite quadrature formula:
    const int N_Gauss = 5;
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
    for (unsigned int i = 0; i < DIM; i++)
      evaluate(*basis_.bases()[i], 0,
	       typename IBASIS::Index(lambda.j(),
				      lambda.e()[i],
				      lambda.k()[i],
				      basis_.bases()[i]),
	       gauss_points[i], v_values[i]);

    // iterate over all points and sum up the integral shares
    int index[DIM]; // current multiindex for the point values
    for (unsigned int i = 0; i < DIM; i++)
      index[i] = 0;
    
    Point<DIM> x;
    while (true) {
      for (unsigned int i = 0; i < DIM; i++)
	x[i] = gauss_points[i][index[i]];
      double share = bvp_->f(x);
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
#endif
  }

  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  void
  CubeEquation<IBASIS,DIM,CUBEBASIS>::RHS
  (const double eta,
   InfiniteVector<double, typename WaveletBasis::Index>& coeffs) const
  {
    coeffs.clear();
    double coarsenorm(0);
    double bound(fnorm_sqr - eta*eta);
    typedef typename WaveletBasis::Index Index;
    typename Array1D<std::pair<Index, double> >::const_iterator it(fcoeffs.begin());
    do {
      coarsenorm += it->second * it->second;
      coeffs.set_coefficient(it->first, it->second);
      ++it;
    } while (it != fcoeffs.end() && coarsenorm < bound);
  }

  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  void
  CubeEquation<IBASIS,DIM,CUBEBASIS>::set_bvp(const EllipticBVP<DIM>* bvp)
  {
    bvp_ = bvp;
    compute_rhs();
  }

  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  double
  CubeEquation<IBASIS,DIM,CUBEBASIS>::s_star() const
  {
    // notation from [St04a]
    const double t = operator_order();
    const int n = DIM;
    const int dT = WaveletBasis::primal_vanishing_moments();
    const double gamma = WaveletBasis::primal_regularity();
    
    return (n == 1
	    ? t+dT // [St04a], Th. 2.3 for n=1
	    : std::min((t+dT)/(double)n, (gamma-t)/(n-1.))); // [St04a, Th. 2.3]
  }

  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  double
  CubeEquation<IBASIS,DIM,CUBEBASIS>::norm_A() const
  {
    if (normA == 0.0) {
      typedef typename WaveletBasis::Index Index;
      std::set<Index> Lambda;
      const int j0 = basis().j0();
      const int jmax = j0+3;
      for (Index lambda = basis().first_generator(j0);; ++lambda) {
	Lambda.insert(lambda);
	if (lambda == basis().last_wavelet(jmax)) break;
      }
      SparseMatrix<double> A_Lambda;
      setup_stiffness_matrix(*this, Lambda, A_Lambda);
      
#if 1
      double help;
      unsigned int iterations;
      LanczosIteration(A_Lambda, 1e-6, help, normA, 200, iterations);
      normAinv = 1./help;
#else
      Vector<double> xk(Lambda.size(), false);
      xk = 1;
      unsigned int iterations;
      normA = PowerIteration(A_Lambda, xk, 1e-6, 100, iterations);
#endif
    }

    return normA;
  }
   
  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  double
  CubeEquation<IBASIS,DIM,CUBEBASIS>::norm_Ainv() const
  {
    if (normAinv == 0.0) {
      typedef typename WaveletBasis::Index Index;
      std::set<Index> Lambda;
      const int j0 = basis().j0();
      const int jmax = j0+3;
      for (Index lambda = basis().first_generator(j0);; ++lambda) {
	Lambda.insert(lambda);
	if (lambda == basis().last_wavelet(jmax)) break;
      }
      SparseMatrix<double> A_Lambda;
      setup_stiffness_matrix(*this, Lambda, A_Lambda);
      
#if 1
      double help;
      unsigned int iterations;
      LanczosIteration(A_Lambda, 1e-6, help, normA, 200, iterations);
      normAinv = 1./help;
#else
      Vector<double> xk(Lambda.size(), false);
      xk = 1;
      unsigned int iterations;
      normAinv = InversePowerIteration(A_Lambda, xk, 1e-6, 200, iterations);
#endif
    }

    return normAinv;
  }

}
