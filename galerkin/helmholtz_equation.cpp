// implementation for helmholtz_equation.h

#include <numerics/quadrature.h>
#include <numerics/schoenberg_splines.h>
#include <interval/interval_bspline.h>
#include <numerics/eigenvalues.h>

using namespace MathTL;

namespace WaveletTL
{
  template <int d, int dT>
  HelmholtzEquation1D<d,dT>::HelmholtzEquation1D(const Function<1>* f,
						 const double alpha,
						 const bool precompute_f)
    : f_(f), alpha_(alpha),
      basis_("P","",1,1,0,0), // PBasis, complementary b.c.'s
      A_(basis_, alpha, no_precond),
      normA(0.0), normAinv(0.0)
  {
    if (precompute_f && f != 0) precompute_rhs();
  }

  template <int d, int dT>
  void
  HelmholtzEquation1D<d,dT>::set_rhs(const InfiniteVector<double,Index>& rhs) const
  {
    fcoeffs_unsorted = rhs;
    fcoeffs_unsorted.scale(this, -1);
    fnorm_sqr = l2_norm_sqr(fcoeffs_unsorted);
    
    // sort the coefficients into fcoeffs
    fcoeffs.resize(fcoeffs_unsorted.size());
    size_type i = 0;
    for (typename InfiniteVector<double,Index>::const_iterator it(fcoeffs_unsorted.begin()),
	   itend(fcoeffs_unsorted.end()); it != itend; ++it, ++i)
      fcoeffs[i] = std::pair<Index,double>(it.index(), *it);
    sort(fcoeffs.begin(), fcoeffs.end(), typename InfiniteVector<double,Index>::decreasing_order());
    
    rhs_precomputed = true;
  }

  template <int d, int dT>
  void
  HelmholtzEquation1D<d,dT>::precompute_rhs() const
  {
    typedef typename WaveletBasis::Index Index;
    // precompute the right-hand side on a fine level
    const int jmax = 15;
    A_.set_level(jmax);
    fcoeffs_unsorted.clear();

    // setup rhs in the phi_{j,k} basis,
    Vector<double> rhs_phijk(A_.row_dimension(), false);
    // perform quadrature with a composite rule on [0,1]
    SimpsonRule simpson;
    CompositeRule<1> composite(simpson, 72);
    SchoenbergIntervalBSpline_td<d> sbs(jmax,0);
    for (int k = basis_.DeltaLmin(); k <= basis_.DeltaRmax(jmax); k++) {
      sbs.set_k(k);
      ProductFunction<1> integrand(f_, &sbs);
      rhs_phijk[k-basis_.DeltaLmin()]
	= composite.integrate(integrand,
			      Point<1>(std::max(0.0, (k+ell1<d>())*ldexp(1.0, -jmax))),
			      Point<1>(std::min(1.0, (k+ell2<d>())*ldexp(1.0, -jmax))));
    }

    // transform rhs into that of psi_{j,k} basis:
    // 1. apply T_{j-1}^T
    Vector<double> rhs(A_.row_dimension(), false);
    assert(jmax > basis_.j0());
    basis_.apply_Tj_transposed(jmax-1, rhs_phijk, rhs);
    
    InfiniteVector<double,Index> rhs_helper;
    size_type i(0);
    for (Index lambda(basis_.first_generator(basis_.j0())); i < rhs.size(); ++lambda, i++)
      {
 	const double coeff = rhs[i];
 	if (fabs(coeff)>1e-15)
 	  rhs_helper.set_coefficient(lambda, coeff);
      }

    set_rhs(rhs_helper);
  }

  template <int d, int dT>
  void
  HelmholtzEquation1D<d,dT>::set_alpha(const double alpha) const
  {
    assert(alpha >= 0);
    alpha_ = alpha;
    A_.set_alpha(alpha);
  }

  template <int d, int dT>
  inline
  double
  HelmholtzEquation1D<d,dT>::D(const typename WaveletBasis::Index& lambda) const
  {
#if 1
    // determine number of index lambda
    size_type number = 0;
    if (lambda.e() == 0) {
      number = lambda.k()-basis_.DeltaLmin();
    } else {
      number = basis_.Deltasize(lambda.j())+lambda.k()-basis_.Nablamin();
    }
    
    return sqrt(A_.diagonal(number));
#else
    return sqrt(a(lambda, lambda));
#endif
  }
  
  template <int d, int dT>
  inline
  double
  HelmholtzEquation1D<d,dT>::a(const typename WaveletBasis::Index& lambda,
			       const typename WaveletBasis::Index& nu) const
  {
    return a(lambda, nu, WaveletBasis::primal_polynomial_degree()*WaveletBasis::primal_polynomial_degree());
  }

  template <int d, int dT>
  double
  HelmholtzEquation1D<d,dT>::a(const typename WaveletBasis::Index& lambda,
			       const typename WaveletBasis::Index& nu,
			       const unsigned int p) const
  {
    double r = 0;

    // first compute the support intersection of \psi_\lambda and \psi_\nu:
    typedef typename WaveletBasis::Support Support;
    Support supp;
    if (intersect_supports(basis_, lambda, nu, supp))
      {
#if 0
	// Set up Gauss points and weights for a composite quadrature formula:
	// (TODO: maybe use an instance of MathTL::QuadratureRule instead of computing
	// the Gauss points and weights)
	const unsigned int N_Gauss = (p+1)/2;
	const double h = ldexp(1.0, -supp.j);
	Array1D<double> gauss_points (N_Gauss*(supp.k2-supp.k1)), func1values, func2values, der1values, der2values;
	for (int patch = supp.k1, id = 0; patch < supp.k2; patch++) // refers to 2^{-j}[patch,patch+1]
	  for (unsigned int n = 0; n < N_Gauss; n++, id++)
	    gauss_points[id] = h*(2*patch+1+GaussPoints[N_Gauss-1][n])/2.;
	
	// - compute point values of the integrands
	evaluate(basis_, lambda, gauss_points, func1values, der1values);
	evaluate(basis_, nu, gauss_points, func2values, der2values);

 	// - add all integral shares
 	for (int patch = supp.k1, id = 0; patch < supp.k2; patch++)
 	  for (unsigned int n = 0; n < N_Gauss; n++, id++) {
//  	    const double t = gauss_points[id];
 	    const double gauss_weight = GaussWeights[N_Gauss-1][n] * h;
	    r += (der1values[id] * der2values[id]
		  + alpha_ * func1values[id] * func2values[id]) * gauss_weight;
 	  }
#else
 	// determine numbers of indices
	size_type number_lambda = 0, number_nu = 0;
	if (lambda.e() == 0) {
	  number_lambda = lambda.k()-basis_.DeltaLmin();
	} else {
	  number_lambda = basis_.Deltasize(lambda.j())+lambda.k()-basis_.Nablamin();
	}
	if (nu.e() == 0) {
	  number_nu = nu.k()-basis_.DeltaLmin();
	} else {
	  number_nu = basis_.Deltasize(nu.j())+nu.k()-basis_.Nablamin();
	}
	
	A_.set_level(std::max(lambda.j()+lambda.e(),nu.j()+nu.e()));
	return A_.get_entry(number_nu, number_lambda);
#endif
      }
    
    return r;
    
  }
  
  template <int d, int dT>
  inline
  void
  HelmholtzEquation1D<d,dT>::RHS(const double eta,
				 InfiniteVector<double, typename WaveletBasis::Index>& coeffs) const
  {
    if (!rhs_precomputed) precompute_rhs();
    
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

  template <int d, int dT>
  double
  HelmholtzEquation1D<d,dT>::f(const typename WaveletBasis::Index& lambda) const
  {
    if (!rhs_precomputed) precompute_rhs();
    return fcoeffs_unsorted.get_coefficient(lambda) * D(lambda);
  }

  template <int d, int dT>
  double
  HelmholtzEquation1D<d,dT>::norm_A() const
  {
    if (normA == 0.0) {
      FullHelmholtz<d,dT> A(basis_, alpha_, energy);
      A.set_level(basis().j0()+4);
      double help;
      unsigned int iterations;
      LanczosIteration(A, 1e-6, help, normA, 200, iterations);
      normAinv = 1./help;
    }

    return normA;
  }
   
  template <int d, int dT>
  double
  HelmholtzEquation1D<d,dT>::norm_Ainv() const
  {
    if (normAinv == 0.0) {
      FullHelmholtz<d,dT> A(basis_, alpha_, energy);
      A.set_level(basis().j0()+4);
      double help;
      unsigned int iterations;
      LanczosIteration(A, 1e-6, help, normA, 200, iterations);
      normAinv = 1./help;
    }

    return normAinv;
  }

  template <int d, int dT>
  void
  HelmholtzEquation1D<d,dT>::add_level (const Index& lambda,
					InfiniteVector<double, Index>& w, const int j,
					const double factor,
					const int J,
					const CompressionStrategy strategy) const
  {
#if 1
    // quick and dirty:
    // compute a full column of the stiffness matrix
    FullHelmholtz<d,dT> A(basis_, alpha_, energy);
    const int jmax = std::max(j+1, lambda.j()+lambda.e());
    A.set_level(jmax);
    std::map<size_type,double> e_lambda, col_lambda;
    size_type number_lambda = 0;
    if (lambda.e() == 0) {
      number_lambda = lambda.k()-basis_.DeltaLmin();
    } else {
      number_lambda = basis_.Deltasize(lambda.j())+lambda.k()-basis_.Nablamin();
    }
//     cout << "add_level(): lambda=" << lambda << ", nr.=" << number_lambda << endl;
    e_lambda[number_lambda] = 1.0;
    A.apply(e_lambda, col_lambda, true);
    
    // extract the entries from level j
    if (j == basis_.j0()-1) {
      // "generator block"
      size_type startrow = 0;
      size_type endrow   = basis_.Deltasize(basis_.j0())-1;
      std::map<size_type,double>::const_iterator it(col_lambda.lower_bound(startrow));
//       assert(startrow <= it->first);
//       assert(it->first <= endrow);
      for (; it != col_lambda.end() && it->first <= endrow; ++it) {
	w.add_coefficient(Index(basis_.j0(), 0, basis_.DeltaLmin()+it->first, &basis_),
			  it->second * factor);
      }
    } else {
      // j>=j0, a "wavelet block"
      size_type startrow = basis_.Deltasize(j);
      size_type endrow   = basis_.Deltasize(j+1)-1;
//       cout << "add_level(): startrow=" << startrow << ", endrow=" << endrow << endl;
      std::map<size_type,double>::const_iterator it(col_lambda.lower_bound(startrow));
//       cout << "add_level(): it->first=" << it->first << endl;
//       assert(startrow <= it->first);
//       assert(it->first <= endrow);
      for (; it != col_lambda.end() && it->first <= endrow; ++it) {
	w.add_coefficient(Index(j, 1, basis_.Nablamin()+it->first-startrow, &basis_),
			  it->second * factor);
      }
    }
#else
    typedef std::list<Index> IntersectingList;
    
    IntersectingList nus;
    
    intersecting_wavelets(basis(), lambda,
			  std::max(j, basis().j0()),
			  j == (basis().j0()-1),
			  nus);
    
    // do the rest of the job
    const double d1 = D(lambda);
    if (strategy == St04a) {
      for (typename IntersectingList::iterator it(nus.begin()), itend(nus.end());
	   it != itend; ++it) {
	if (abs(lambda.j()-j) <= J/((double) space_dimension) ||
	    intersect_singular_support(basis(), lambda, *it)) {
	  const double entry = a(*it, lambda);
	  w.add_coefficient(*it, (entry / (d1*D(*it))) * factor);
	}
      }
    }
    else if (strategy == CDD1) {
      for (typename IntersectingList::const_iterator it(nus.begin()), itend(nus.end());
	   it != itend; ++it) {
	const double entry = a(*it, lambda);
	w.add_coefficient(*it, (entry / (d1 * D(*it))) * factor);
      }
    }   
#endif
  }

  
}
