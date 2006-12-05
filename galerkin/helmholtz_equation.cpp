// implementation for helmholtz_equation.h

#include <numerics/quadrature.h>
#include <numerics/schoenberg_splines.h>
#include <interval/interval_bspline.h>

using namespace MathTL;

namespace WaveletTL
{
  template <int d, int dT>
  HelmholtzEquation1D<d,dT>::HelmholtzEquation1D(const Function<1>* f,
						 const double alpha,
						 const bool precompute_f)
    : f_(f), alpha_(alpha),
      basis_("P","",1,1,0,0), // PBasis, complementary b.c.'s
      A_(basis_, alpha, no_precond)
  {
    if (precompute_f) precompute_rhs();
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

    const int j0 = basis_.j0();
    unsigned int i = 0;
    for (Index lambda(basis_.first_generator(j0));; ++lambda, i++)
      {
 	const double coeff = rhs[i]/sqrt(A_.diagonal(i));
 	if (fabs(coeff)>1e-15)
 	  fcoeffs_unsorted.set_coefficient(lambda, coeff);
 	if (lambda == basis_.last_wavelet(jmax-1))
 	  break;
      }
    fnorm_sqr = l2_norm_sqr(fcoeffs_unsorted);
    
    // sort the coefficients into fcoeffs
    fcoeffs.resize(fcoeffs_unsorted.size());
    i = 0;
    for (typename InfiniteVector<double,Index>::const_iterator it(fcoeffs_unsorted.begin()),
	   itend(fcoeffs_unsorted.end()); it != itend; ++it, ++i)
      fcoeffs[i] = std::pair<Index,double>(it.index(), *it);
    sort(fcoeffs.begin(), fcoeffs.end(), typename InfiniteVector<double,Index>::decreasing_order());
    
    rhs_precomputed = true;
  }

  template <int d, int dT>
  inline
  double
  HelmholtzEquation1D<d,dT>::D(const typename WaveletBasis::Index& lambda) const
  {
    // determine number of index lambda
    size_type number = 0;
    if (lambda.e() == 0) {
      number = lambda.k()-basis_.DeltaLmin();
    } else {
      number = basis_.Deltasize(lambda.j())+lambda.k()-basis_.Nablamin();
    }
    
    return sqrt(A_.diagonal(number));
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

}
