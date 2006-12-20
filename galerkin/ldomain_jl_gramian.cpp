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
		  r += *itlambda * *itmu
		    * gauss2d.integrate(integrand,
					MathTL::Point<2>(    k0*h,    k1*h),
					MathTL::Point<2>((k0+1)*h,(k1+1)*h));
		}
	      }
	    
	  }
#else
	// use that both \psi_\lambda and \psi_\mu are a tensor product of 1D bases;
	// the entry of the Gramian on each square subpatch is a product of 2 1D integrals

	// first compute the 1D generator coefficients;
	FixedArray1D<int,2> psilambdalevel, psimulevel;
	FixedArray1D<InfiniteVector<double,RMWIndex>,2> coeffs_lambda, coeffs_mu;
	for (int n = 0; n <= 1; n++) {
	  psilambdalevel[n] = lambda.j()+(lambda.e()[n]==1 ? 1 : 0);
	  basis_.reconstruct_1_1d(RMWIndex(lambda.j(),lambda.e()[n],lambda.c()[n],lambda.k()[n]),
				  psilambdalevel[n], coeffs_lambda[n]);
	  psimulevel[n] = mu.j()+(mu.e()[n]==1 ? 1 : 0);
	  basis_.reconstruct_1_1d(RMWIndex(mu.j(),mu.e()[n],mu.c()[n],mu.k()[n]),
				  psimulevel[n], coeffs_mu[n]);
	}

	// then iterate through the subsquares of supp and compute the integral shares
	const double h = ldexp(1.0, -supp.j); // granularity for the quadrature
	MathTL::GaussLegendreRule gauss1d(8);
	FixedArray1D<int,2> k;
	for (k[0] = supp.xmin; k[0] < supp.xmax; k[0]++)
	  for (k[1] = supp.ymin; k[1] < supp.ymax; k[1]++) {
	    // check whether 2^{-supp.j}[k0,k0+1]x[k1,k1+1] is contained in Omega
	    if ((k[0] >= -(1<<supp.j) && k[0] < (1<<supp.j) && k[1] >= -(1<<supp.j) && k[1] < 0)
		|| (k[0] >= -(1<<supp.j) && k[0] < 0 && k[1] >= 0 && k[1] < (1<<supp.j))) {
	      // perform the integration
	      FixedArray1D<double,2> factor;
	      for (int n=0; n <= 1; n++) {
		factor[n] = 0;
		for (InfiniteVector<double,RMWIndex>::const_iterator itlambda(coeffs_lambda[n].begin());
		     itlambda != coeffs_lambda[n].end(); ++itlambda) {
		  MathTL::CubicHermiteInterpolant_td glambda(psilambdalevel[n],
							     itlambda.index().c(),
							     itlambda.index().k());
		  for (InfiniteVector<double,RMWIndex>::const_iterator itmu(coeffs_mu[n].begin());
		       itmu != coeffs_mu[n].end(); ++itmu) {
		    MathTL::CubicHermiteInterpolant_td gmu(psimulevel[n],
							   itmu.index().c(),
							   itmu.index().k());
		    ProductFunction<1> integrand(&glambda, &gmu);
		    factor[n] += *itlambda * *itmu
		      * gauss1d.integrate(integrand,
					  MathTL::Point<1>(    k[n]*h),
					  MathTL::Point<1>((k[n]+1)*h));
		  }
		}
	      }
	      r += factor[0] * factor[1];
	    }
	  }
#endif
      }

    return r;
  }

}
