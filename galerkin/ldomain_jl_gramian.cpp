// implementation for ldomain_jl_gramian.h

#include <Ldomain/ldomain_jl_support.h>
#include <numerics/quadrature.h>
#include <numerics/gauss_quadrature.h>
#include <numerics/bezier.h>

namespace WaveletTL
{
  LDomainJLGramian::LDomainJLGramian(const WaveletBasis& basis,
				     const InfiniteVector<double,Index>& y)
    : basis_(basis), y_(y), normA(0.0), normAinv(0.0)
  {
  }
  
  double
  LDomainJLGramian::a(const Index& lambda,
		      const Index& mu) const
  {
    double r = 0;
    
    // First we compute the support intersection of \psi_\lambda and \psi_\mu:
    typedef WaveletBasis::Support Support;
    Support supp;
    
    if (intersect_supports(basis_, lambda, mu, supp))
      {
#if 0
	// first variant:
	// compute the 2D generator expansions of \psi_\lambda and \psi_\mu
	InfiniteVector<double,Index> psilambdacoeffs, psimucoeffs;
	const int psilambdalevel = lambda.j()+ (lambda.e()[0]==1 || lambda.e()[1]==1 ? 1 : 0);
	basis_.reconstruct_1(lambda, psilambdalevel, psilambdacoeffs);
	const int psimulevel = mu.j()+ (mu.e()[0]==1 || mu.e()[1]==1 ? 1 : 0);
	basis_.reconstruct_1(lambda, psimulevel, psimucoeffs);
	
	const double h = ldexp(1.0, -supp.j); // granularity for the quadrature
	MathTL::GaussLegendreRule gauss1d(5);
	MathTL::QuadratureRule<2> gauss2d(gauss1d, gauss1d); // tensor product quadrature rule

	// iterate through the involved generators and collect the integral shares
	for (InfiniteVector<double,Index>::const_iterator itlambda = psilambdacoeffs.begin();
	     itlambda != psilambdacoeffs.end(); ++itlambda)
	  for (InfiniteVector<double,Index>::const_iterator itmu = psimucoeffs.begin();
	       itmu != psimucoeffs.end(); ++itmu) {
	    // integrate two generators over supp
	    MathTL::CubicHermiteInterpolant2D_td glambda(psilambdalevel,
							 itlambda.index().c()[0],
							 itlambda.index().c()[1],
							 itlambda.index().k()[0],
							 itlambda.index().k()[1]);
	    MathTL::CubicHermiteInterpolant2D_td gmu(psimulevel,
						     itmu.index().c()[0],
						     itmu.index().c()[1],
						     itmu.index().k()[0],
						     itmu.index().k()[1]);
	    MathTL::ProductFunction<2> integrand(&glambda, &gmu);

	    // iterate through the squares of supp
	    for (int k0 = supp.xmin; k0 < supp.xmax; k0++)
	      for (int k1 = supp.ymin; k1 < supp.ymax; k1++) {
		// check whether 2^{-supp.j}[k0,k0+1]x[k1,k1+1] is contained in Omega
		if ((k0 >= -(1<<supp.j) && k0 < (1<<supp.j) && k1 >= -(1<<supp.j) && k1 < 0)
		    || (k0 >= -(1<<supp.j) && k0 < 0 && k1 >= 0 && k1 < (1<<supp.j))) {
		  r += *itlambda * *itmu *
		    gauss2d.integrate(integrand,
				      MathTL::Point<2>(    k0*h,    k1*h),
				      MathTL::Point<2>((k0+1)*h,(k1+1)*h));
		}
	      }
	    
	  }
#else
	// use that both \psi_\lambda and \psi_\mu are a tensor product of 1D bases;
	// first perform 2 1D reconstruct_1D() calls:
	


	r = 42;
#endif
      }

    return r;
  }

}
