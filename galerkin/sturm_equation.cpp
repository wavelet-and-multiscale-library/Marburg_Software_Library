// implementation for sturm_equation.h

#include <cmath>
#include <utils/array1d.h>
#include <numerics/gauss_data.h>

namespace WaveletTL
{
  template <class WBASIS>
  SturmEquation<WBASIS>::SturmEquation(const simpleSturmBVP& bvp)
    : bvp_(bvp), wbasis_(bvp.bc_left(), bvp.bc_right())
  {
  }

  template <class WBASIS>
  double
  SturmEquation<WBASIS>::a(const typename WBASIS::Index& lambda,
			   const typename WBASIS::Index& nu,
			   const unsigned int p) const
  {
    // a(u,v) = \int_0^1 [p(t)u'(t)v'(t)+q(t)u(t)v(t)] dt

    double r = 0;

    // Remark: There are of course many possibilities to evaluate
    // a(u,v) numerically.
    // In this implementation, we rely on the fact that the primal functions in
    // WBASIS are splines with respect to a dyadic subgrid.
    // We can then apply an appropriate composite quadrature rule.
    // In the scope of WBASIS, the routines intersect_supports() and evaluate()
    // must exist, which is the case for DSBasis<d,dT>.

    // First we compute the support intersection of \psi_\lambda and \psi_\nu:
    int j, k1, k2;
    bool inter = intersect_supports(wbasis_, lambda, nu, j, k1, k2);

    if (inter)
      {
	// Set up Gauss points and weights for a composite quadrature formula:
	// (TODO: maybe use an instance of MathTL::QuadratureRule instead of computing
	// the Gauss points and weights)
	const unsigned int N_Gauss = (p+1)/2;
	const double h = ldexp(1.0, -j);
	Array1D<double> gauss_points (N_Gauss*(k2-k1)), func1values, func2values, der1values, der2values;
	for (int patch = k1, id = 0; patch < k2; patch++) // refers to 2^{-j}[patch,patch+1]
	  for (unsigned int n = 0; n < N_Gauss; n++, id++)
	    gauss_points[id] = h*(2*patch+1+GaussPoints[N_Gauss-1][n])/2.;

	// - compute point values of the integrands
	evaluate(wbasis_, lambda, gauss_points, func1values, der1values);
	evaluate(wbasis_, nu, gauss_points, func2values, der2values);

	// - add all integral shares
	for (int patch = k1, id = 0; patch < k2; patch++)
	  for (unsigned int n = 0; n < N_Gauss; n++, id++) {
	    const double t = gauss_points[id];
	    const double gauss_weight = GaussWeights[N_Gauss-1][n] * h;
	    
	    const double pt = bvp_.p(t);
	    if (pt != 0)
	      r += pt * der1values[id] * der2values[id] * gauss_weight;
	    
	    const double qt = bvp_.q(t);
	    if (qt != 0)
	      r += qt * func1values[id] * func2values[id] * gauss_weight;
	  }
      }

    return r;
  }

  template <class WBASIS>
  double
  SturmEquation<WBASIS>::f(const typename WBASIS::Index& lambda) const
  {
    // f(v) = \int_0^1 g(t)v(t) dt

    double r = 0;

    const int j = lambda.j()+lambda.e();
    int k1, k2;
    support(wbasis_, lambda, k1, k2);

    // Set up Gauss points and weights for a composite quadrature formula:
    const unsigned int N_Gauss = 5;
    const double h = ldexp(1.0, -j);
    Array1D<double> gauss_points (N_Gauss*(k2-k1));
    for (int patch = k1; patch < k2; patch++) // refers to 2^{-j}[patch,patch+1]
      for (unsigned int n = 0; n < N_Gauss; n++)
	gauss_points[(patch-k1)*N_Gauss+n] = h*(2*patch+1+GaussPoints[N_Gauss-1][n])/2;
    
    // - add all integral shares
    for (unsigned int n = 0; n < N_Gauss; n++)
      {
	const double gauss_weight = GaussWeights[N_Gauss-1][n] * h;
	for (int patch = k1; patch < k2; patch++)
	  {
	    const double t = gauss_points[(patch-k1)*N_Gauss+n];
	    
	    const double gt = bvp_.g(t);
	    if (gt != 0)
	      r += gt
		* evaluate(wbasis_, 0, lambda, t)
		* gauss_weight;
	  }
      }
    
    return r;
  }

  template <class WBASIS>
  void
  SturmEquation<WBASIS>::RHS(InfiniteVector<double, typename WBASIS::Index>& coeffs,
			     const double eta) const
  {
    coeffs.clear();

    // remark: for a quick hack, we use a projection of f onto a space V_{jmax}
    // of the given multiresolution analysis

    const int j0 = wbasis_.j0();
    const int jmax = j0+1;
    for (typename WBASIS::Index lambda(first_generator(&wbasis_, j0));; ++lambda)
      {
	coeffs.set_coefficient(lambda, f(lambda)/D(lambda));
  	if (lambda == last_wavelet(&wbasis_, jmax))
	  break;
      }
  }

  template <class WBASIS>
  inline
  double
  SturmEquation<WBASIS>::D(const typename WBASIS::Index& lambda) const
  {
    return ldexp(1.0, lambda.j());
  }

  template <class WBASIS>
  inline
  void
  SturmEquation<WBASIS>::rescale(InfiniteVector<double, typename WBASIS::Index>& coeffs,
				 const int n) const
  {
    for (typename InfiniteVector<double, typename WBASIS::Index>::const_iterator it(coeffs.begin());
	 it != coeffs.end(); ++it)
      {
	// TODO: implement an InfiniteVector::iterator to speed up this hack!
	coeffs.set_coefficient(it.index(), *it * pow(D(it.index()), n));
      }
  }
}
