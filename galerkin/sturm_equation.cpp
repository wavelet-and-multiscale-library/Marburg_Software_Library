// implementation for sturm_equation.h

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
			   const typename WBASIS::Index& nu) const
  {
    // a(u,v) = \int_0^1 [p(t)u'(t)v'(t)+q(t)u(t)v(t)] dt

    double r = 0;

    // Remark: There are of course many possibilities to evaluate
    // a(u,v) numerically.
    // In this implementation, we rely on the fact that the primal functions in
    // WBASIS are splines with respect to a dyadic subgrid.
    // We can then apply an appropriate composite quadrature rule.
    // In the scope of WBASIS, the routines intersect_supports() and evaluate()
    // must exist, which is the case for DKUBasis<d,dT>.

    // First we compute the support intersection of \psi_\lambda and \psi_\nu:
    int j, k1, k2;
    bool inter(intersect_supports(wbasis_, lambda, nu, j, k1, k2));

    if (inter)
      {
	// Set up Gauss points and weights for a composite quadrature formula:
	// (TODO: maybe use an instance of MathTL::QuadratureRule instead of computing
	// the Gauss points and weights)
	const unsigned int N_Gauss = 4;
	const double h = ldexp(1.0, -j);
	Array1D<double> gauss_points (N_Gauss*(k2-k1));
	for (int patch = k1; patch < k2; patch++) // refers to 2^{-j}[patch,patch+1]
	  for (unsigned int n = 0; n < N_Gauss; n++)
	    gauss_points[(patch-k1)*N_Gauss+n] = h*(2*patch+1+GaussPoints[N_Gauss-1][n])/2;

	// - add all integral shares
	for (unsigned int n = 0; n < N_Gauss; n++)
	  for (int patch = k1; patch < k2; patch++)
	    {
	      const double gauss_weight = GaussWeights[N_Gauss-1][n] * h;
	      const double t = gauss_points[(patch-k1)*N_Gauss+n];

	      const double pt = bvp_.p(t);
	      if (pt != 0)
		r += pt
		  * evaluate(wbasis_, 1, lambda, t)
		  * evaluate(wbasis_, 1, nu, t)
		  * gauss_weight;

	      const double qt = bvp_.q(t);
	      if (qt != 0)
		r += qt
		  * evaluate(wbasis_, 0, lambda, t)
		  * evaluate(wbasis_, 0, nu, t)
		  * gauss_weight;
	    }
      }

    return r;
  }

  template <class WBASIS>
  inline
  double
  SturmEquation<WBASIS>::D(const typename WBASIS::Index& lambda) const
  {
    return ldexp(1.0, lambda.j());
  }
}
