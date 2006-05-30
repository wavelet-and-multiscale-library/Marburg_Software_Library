// implementation for ldomain_equation.h

namespace WaveletTL
{
  template <class IBASIS>
  LDomainEquation<IBASIS>::LDomainEquation(const EllipticBVP<2>* bvp)
    : bvp_(bvp), basis_(), normA(0.0), normAinv(0.0)
  {
    compute_rhs();
  }

  template <class IBASIS>
  LDomainEquation<IBASIS>::LDomainEquation(const LDomainEquation& eq)
    : bvp_(eq.bvp_), basis_(eq.basis_),
      fcoeffs(eq.fcoeffs), fnorm_sqr(eq.fnorm_sqr),
      normA(eq.normA), normAinv(eq.normAinv)
  {
  }

  template <class IBASIS>
  inline
  double
  LDomainEquation<IBASIS>::D(const typename WaveletBasis::Index& lambda) const
  {
    return ldexp(1.0, lambda.j());
//     return sqrt(a(lambda, lambda));
//     return lambda.e() == 0 ? 1.0 : ldexp(1.0, lambda.j()); // do not scale the generators
//     return lambda.e() == 0 ? 1.0 : sqrt(a(lambda, lambda)); // do not scale the generators
  }

  template <class IBASIS>
  double
  LDomainEquation<IBASIS>::f(const typename WaveletBasis::Index& lambda) const
  {
    // f(v) = \int_0^1 g(t)v(t) dt

    double r = 0;

    // first compute supp(psi_lambda)
    typename WaveletBasis::Support supp;
    support(basis_, lambda, supp);

    // setup Gauss points and weights for a composite quadrature formula:
    const int N_Gauss = 5;
    const double h = ldexp(1.0, -supp.j); // granularity for the quadrature
//     FixedArray1D<Array1D<double>,2> gauss_points, gauss_weights, v_values;
//     for (unsigned int i = 0; i <=1; i++) {
//       gauss_points[i].resize(N_Gauss*(supp.b[i]-supp.a[i]));
//       gauss_weights[i].resize(N_Gauss*(supp.b[i]-supp.a[i]));
//       for (int patch = supp.a[i]; patch < supp.b[i]; patch++)
// 	for (int n = 0; n < N_Gauss; n++) {
// 	  gauss_points[i][(patch-supp.a[i])*N_Gauss+n]
// 	    = h*(2*patch+1+GaussPoints[N_Gauss-1][n])/2.;
// 	  gauss_weights[i][(patch-supp.a[i])*N_Gauss+n]
// 	    = h*GaussWeights[N_Gauss-1][n];
// 	}
//     }

//     // compute the point values of the integrand (where we use that it is a tensor product)
//     for (unsigned int i = 0; i < DIM; i++)
//       evaluate(*basis_.bases()[i], 0,
// 	       typename IBASIS::Index(lambda.j(),
// 				      lambda.e()[i],
// 				      lambda.k()[i],
// 				      basis_.bases()[i]),
// 	       gauss_points[i], v_values[i]);

//     // iterate over all points and sum up the integral shares
//     int index[DIM]; // current multiindex for the point values
//     for (unsigned int i = 0; i < DIM; i++)
//       index[i] = 0;
    
//     Point<DIM> x;
//     while (true) {
//       for (unsigned int i = 0; i < DIM; i++)
// 	x[i] = gauss_points[i][index[i]];
//       double share = bvp_->f(x);
//       for (unsigned int i = 0; i < DIM; i++)
// 	share *= gauss_weights[i][index[i]] * v_values[i][index[i]];
//       r += share;

//       // "++index"
//       bool exit = false;
//       for (unsigned int i = 0; i < DIM; i++) {
// 	if (index[i] == N_Gauss*(supp.b[i]-supp.a[i])-1) {
// 	  index[i] = 0;
// 	  exit = (i == DIM-1);
// 	} else {
// 	  index[i]++;
// 	  break;
// 	}
//       }
//       if (exit) break;
//     }
    
    return r;
  }

  template <class IBASIS>
  void
  LDomainEquation<IBASIS>::compute_rhs()
  {
    cout << "LDomainEquation(): precompute right-hand side..." << endl;

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

  template <class IBASIS>
  void
  LDomainEquation<IBASIS>::RHS
  (const double eta,
   InfiniteVector<double, typename WaveletBasis::Index>& coeffs) const
  {
//     coeffs.clear();
//     double coarsenorm(0);
//     double bound(fnorm_sqr - eta*eta);
//     typedef typename WaveletBasis::Index Index;
//     typename Array1D<std::pair<Index, double> >::const_iterator it(fcoeffs.begin());
//     do {
//       coarsenorm += it->second * it->second;
//       coeffs.set_coefficient(it->first, it->second);
//       ++it;
//     } while (it != fcoeffs.end() && coarsenorm < bound);
  }


}
