// implementation for ldomain_jl_expansion.h

#include <numerics/quadrature.h>
#include <numerics/gauss_quadrature.h>
#include <numerics/bezier.h>

namespace WaveletTL
{
  double
  integrate(const Function<2>* f,
	    const LDomainJLBasis& basis,
	    const Index& lambda)
  {
    double r = 0;
    
    // compute the generator expansion of psi_lambda
    InfiniteVector<double, Index> gcoeffs;
    const int level = lambda.j()+ (lambda.e()[0]==1 || lambda.e()[1]==1 ? 1 : 0);
    basis.reconstruct_1(lambda, level, gcoeffs);
    
    // iterate through the involved generators and add the shares
    const double h = ldexp(1.0, -level); // granularity for the quadrature
    MathTL::GaussLegendreRule gauss1d(5);
    MathTL::QuadratureRule<2> gauss2d(gauss1d, gauss1d); // tensor product quadrature rule
    MathTL::CubicHermiteInterpolant2D_td generator(level, 0, 0, 0, 0);

    for (InfiniteVector<double,Index>::const_iterator it(gcoeffs.begin()),
  	   itend(gcoeffs.end()); it != itend; ++it)
      {
	// We know that the support of the involved generator is the intersection
	// of a square 2^{-j}[k0-1,k0+1]x[k1-1,k1+1] with the L-shaped domain.
	// On each of the four subpatches, the function is smooth, so we can simply employ
	// a quadrature rule four (or less) times.
	generator.set_c0(it.index().c()[0]);
	generator.set_c1(it.index().c()[1]);
	const int k0 = it.index().k()[0];
	const int k1 = it.index().k()[1];
	generator.set_k0(k0);
	generator.set_k1(k1);
	MathTL::ProductFunction<2> integrand(f, &generator);
	
 	if (k0 >= 1-(1<<level) && k1 >= 1-(1<<level)) { // lower left subsquare
 	  r += gauss2d.integrate(integrand,
				 MathTL::Point<2>((k0-1)*h,(k1-1)*h),
				 MathTL::Point<2>(    k0*h,    k1*h));
	}

	if (!(k1 == -(1<<level) || k0 == (1<<level) || (k0 == 0 && k1 >= 1))) { // lower right subsquare
 	  r += gauss2d.integrate(integrand,
				 MathTL::Point<2>(    k0*h,(k1-1)*h),
				 MathTL::Point<2>((k0+1)*h,    k1*h));
	}

	if (!(k0 == -(1<<level) || k1 == (1<<level) || (k1 == 0 && k0 >= 1))) { // upper left subsquare
 	  r += gauss2d.integrate(integrand,
				 MathTL::Point<2>((k0-1)*h,    k1*h),
				 MathTL::Point<2>(    k0*h,(k1+1)*h));
	}

	if (!(k1 == (1<<level) || (k0 == 0 && k1 >= 0) || (k1 == 0 && k0 >= 0) || k0 == (1<<level))) { // upper right subsquare
 	  r += gauss2d.integrate(integrand,
				 MathTL::Point<2>(    k0*h,    k1*h),
				 MathTL::Point<2>((k0+1)*h,(k1+1)*h));
	}
      }
    
    return r;
  }
  
  void
  expand(const Function<2>* f,
	 const LDomainJLBasis& basis,
	 const bool primal,
	 const int jmax,
	 InfiniteVector<double,Index>& coeffs)
  {
    const int j0 = basis.j0();
    for (Index lambda = first_generator(j0);;++lambda)
      {
	const double coeff = integrate(f, basis, lambda);
	if (fabs(coeff) > 1e-15)
	  coeffs.set_coefficient(lambda, coeff);
 	if (lambda == last_wavelet(jmax))
 	  break;
      } 

    if (!primal) {
      // TODO!!!
    }
  }

}
