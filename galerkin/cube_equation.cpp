// implementation for cube_equation.h

namespace WaveletTL
{
  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  CubeEquation<IBASIS,DIM,CUBEBASIS>::CubeEquation(const EllipticBVP<DIM>* bvp,
						   const FixedArray1D<bool,2*DIM>& bc)
    : bvp_(bvp), basis_(bc)
  {
    typedef typename WaveletBasis::Index Index;
    
    cout << "CubeEquation() setup ..." << endl;

    // precompute the right-hand side on a fine level
    InfiniteVector<double,Index> fhelp;
    const int j0   = basis().j0();
    const int jmax = 6; // for a first quick hack
    for (Index lambda(first_generator<IBASIS,DIM,CUBEBASIS>(&basis_, j0));; ++lambda)
      {
// 	cout << "CubeEquation() setup, integration of psi_lambda, lambda=" << lambda << endl;
	const double coeff = f(lambda)/D(lambda);
	if (fabs(coeff)>1e-15)
	  fhelp.set_coefficient(lambda, coeff);
  	if (lambda == last_wavelet<IBASIS,DIM,CUBEBASIS>(&basis_, jmax))
	  break;
      }
    fnorm_sqr = l2_norm_sqr(fhelp);

    cout << "CubeEquation() setup, all integrals for right-hand side computed" << endl;

    // sort the coefficients into fcoeffs
    fcoeffs.resize(fhelp.size());
    unsigned int id(0);
    for (typename InfiniteVector<double,Index>::const_iterator it(fhelp.begin()), itend(fhelp.end());
	 it != itend; ++it, ++id)
      fcoeffs[id] = std::pair<Index,double>(it.index(), *it);
    sort(fcoeffs.begin(), fcoeffs.end(), typename InfiniteVector<double,Index>::decreasing_order());
  }
  
  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  inline
  double
  CubeEquation<IBASIS,DIM,CUBEBASIS>::D(const typename WaveletBasis::Index& lambda) const
  {
    return ldexp(1.0, lambda.j());
  }

  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  inline
  void
  CubeEquation<IBASIS,DIM,CUBEBASIS>::rescale
  (InfiniteVector<double,typename WaveletBasis::Index>& coeffs,
   const int n) const
  {
    for (typename InfiniteVector<double,typename WaveletBasis::Index>::const_iterator it(coeffs.begin());
	 it != coeffs.end(); ++it)
      {
	// TODO: implement an InfiniteVector::iterator to speed up this hack!
	coeffs.set_coefficient(it.index(), *it * pow(D(it.index()), n));
      }
  }

  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  double
  CubeEquation<IBASIS,DIM,CUBEBASIS>::a(const typename WaveletBasis::Index& lambda,
					const typename WaveletBasis::Index& nu,
					const unsigned int p) const
  {
    // a(u,v) = \int_Omega [a(x)grad u(x)grad v(x)+q(x)u(x)v(x)] dx

    double r = 0;

    // first decide whether the supports of psi_lambda and psi_nu intersect
    typedef typename CUBEBASIS::Support Support;
    Support supp;
    
    if (intersect_supports<IBASIS,DIM,CUBEBASIS>(basis_, lambda, nu, supp))
      {
      }
						

//     // First we compute the support intersection of \psi_\lambda and \psi_\nu:
//     typedef typename WBASIS::Support Support;

//     Support supp;

//     if (intersect_supports(basis_, lambda, nu, supp))
//       {
// 	// Set up Gauss points and weights for a composite quadrature formula:
// 	// (TODO: maybe use an instance of MathTL::QuadratureRule instead of computing
// 	// the Gauss points and weights)
// 	const unsigned int N_Gauss = (p+1)/2;
// 	const double h = ldexp(1.0, -supp.j);
// 	Array1D<double> gauss_points (N_Gauss*(supp.k2-supp.k1)), func1values, func2values, der1values, der2values;
// 	for (int patch = supp.k1, id = 0; patch < supp.k2; patch++) // refers to 2^{-j}[patch,patch+1]
// 	  for (unsigned int n = 0; n < N_Gauss; n++, id++)
// 	    gauss_points[id] = h*(2*patch+1+GaussPoints[N_Gauss-1][n])/2.;

// 	// - compute point values of the integrands
// 	evaluate(basis_, lambda, gauss_points, func1values, der1values);
// 	evaluate(basis_, nu, gauss_points, func2values, der2values);

// 	// - add all integral shares
// 	for (int patch = supp.k1, id = 0; patch < supp.k2; patch++)
// 	  for (unsigned int n = 0; n < N_Gauss; n++, id++) {
// 	    const double t = gauss_points[id];
// 	    const double gauss_weight = GaussWeights[N_Gauss-1][n] * h;
	    
// 	    const double pt = bvp_.p(t);
//   	    if (pt != 0)
// 	      r += pt * der1values[id] * der2values[id] * gauss_weight;
	    
// 	    const double qt = bvp_.q(t);
//   	    if (qt != 0)
// 	      r += qt * func1values[id] * func2values[id] * gauss_weight;
// 	  }
//       }

    return r;
  }
  
  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  double
  CubeEquation<IBASIS,DIM,CUBEBASIS>::f(const typename WaveletBasis::Index& lambda) const
  {
    // f(v) = \int_0^1 g(t)v(t) dt

    double r = 0;

    // first compute supp(psi_lambda)
    typename CUBEBASIS::Support supp;
    support<IBASIS,DIM,CUBEBASIS>(basis_, lambda, supp);

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
}
