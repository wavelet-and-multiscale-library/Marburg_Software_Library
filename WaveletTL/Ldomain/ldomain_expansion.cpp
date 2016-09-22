// implementation for ldomain_expansion.h

#include <cmath>
#include <time.h>
#include <utils/fixed_array1d.h>
#include <numerics/gauss_data.h>

namespace WaveletTL
{
  template <class IBASIS>
  double integrate(const Function<2>* f,
		   const LDomainBasis<IBASIS>& basis,
		   const typename LDomainBasis<IBASIS>::Index& lambda)
  {
    double r = 0;

#if 1
    // new variant, integrate over supp(psi_lambda)
    typename LDomainBasis<IBASIS>::Support supp;
    basis.support(lambda, supp);
    
    const int N_Gauss = 5; // number of Gauss points in x- and y-direction (per sing. supp. subpatch)
    const double h = ldexp(1.0, -supp.j); // granularity for the quadrature
    
    Array1D<double> gauss_points0, gauss_points1, gauss_weights0, gauss_weights1,
      psi_lambda_values; // point values of psi_lambda
    
    // per patch, collect all point values
    for (int p = 0; p <= 2; p++) {
      if (supp.xmin[p] != -1) { // psi_lambda is nontrivial on patch p
	// prepare Gauss points and weights in x direction
	gauss_points0.resize(N_Gauss*(supp.xmax[p]-supp.xmin[p]));
	gauss_weights0.resize(N_Gauss*(supp.xmax[p]-supp.xmin[p]));
	for (int subpatch = supp.xmin[p]; subpatch < supp.xmax[p]; subpatch++)
	  for (int n = 0; n < N_Gauss; n++) {
	    gauss_points0[(subpatch-supp.xmin[p])*N_Gauss+n]
	      = h*(2*subpatch+1+GaussPoints[N_Gauss-1][n])/2.;
	    gauss_weights0[(subpatch-supp.xmin[p])*N_Gauss+n]
	      = h*GaussWeights[N_Gauss-1][n];
	  }
	
	// prepare Gauss points and weights in y direction
	gauss_points1.resize(N_Gauss*(supp.ymax[p]-supp.ymin[p]));
	gauss_weights1.resize(N_Gauss*(supp.ymax[p]-supp.ymin[p]));
	for (int subpatch = supp.ymin[p]; subpatch < supp.ymax[p]; subpatch++)
	  for (int n = 0; n < N_Gauss; n++) {
	    gauss_points1[(subpatch-supp.ymin[p])*N_Gauss+n]
	      = h*(2*subpatch+1+GaussPoints[N_Gauss-1][n])/2.;
	    gauss_weights1[(subpatch-supp.ymin[p])*N_Gauss+n]
	      = h*GaussWeights[N_Gauss-1][n];
	  }
	
	// evaluate pullback of psi_lambda
	basis.evaluate(lambda, p, gauss_points0, gauss_points1, psi_lambda_values);
	
	// determine chart action
	const double xoffset = (p == 2 ? 0 : -1);
	const double yoffset = (p == 0 ? 0 : -1);

	// aggregate the integral shares
	Point<2> x;
	for (unsigned int i0 = 0, id = 0; i0 < gauss_points0.size(); i0++) {
	  x[0] = gauss_points0[i0] + xoffset; // apply chart
	  for (unsigned int i1 = 0; i1 < gauss_points1.size(); i1++) {
	    x[1] = gauss_points1[i1] + yoffset; // apply chart

	    r += gauss_weights0[i0] * gauss_weights1[i1]
	      * psi_lambda_values[id++] * f->value(x);
	  }
	}
      }
    }
#else
    // compute the generator expansion of psi_lambda
    InfiniteVector<double, Index> gcoeffs;
    typename InfiniteVector<double,Index>::const_iterator it(gcoeffs.begin()), gcoeffs_end(gcoeffs.begin());

    const int ecode(lambda.e()[0]+2*lambda.e()[1]);
    const int level = lambda.j() + (ecode == 0 ? 0 : 1);
    if (ecode == 0) {
      gcoeffs.set_coefficient(lambda, 1.0);
      it = gcoeffs.begin();
      gcoeffs_end = gcoeffs.end();
    }
    else
      basis.reconstruct_1(lambda, level, it, gcoeffs_end);

    // iterate through the involved generators and collect the point evaluations
    const int N_Gauss = 5; // number of Gauss points in x- and y-direction (per sing. supp. subpatch)
    const double h = ldexp(1.0, -level); // granularity for the quadrature

    typename LDomainBasis<IBASIS>::Support supp;
    Array1D<double> gauss_points0, gauss_points1, gauss_weights0, gauss_weights1,
      psi_mu_0_values, psi_mu_1_values;
    
    for (; it != gcoeffs_end; ++it)
      {
 	// compute the support of the generator corresponding to mu=it.index()
  	support(basis, it.index(), supp);
	
 	// per patch, collect all point values
 	for (int p = 0; p <= 2; p++) {
 	  if (supp.xmin[p] != -1) { // psi_mu is nontrivial on patch p
 	    // prepare Gauss points and weights in x direction
 	    gauss_points0.resize(N_Gauss*(supp.xmax[p]-supp.xmin[p]));
  	    gauss_weights0.resize(N_Gauss*(supp.xmax[p]-supp.xmin[p]));
 	    for (int subpatch = supp.xmin[p]; subpatch < supp.xmax[p]; subpatch++)
 	      for (int n = 0; n < N_Gauss; n++) {
		gauss_points0[(subpatch-supp.xmin[p])*N_Gauss+n]
		  = h*(2*subpatch+1+GaussPoints[N_Gauss-1][n])/2.;
		gauss_weights0[(subpatch-supp.xmin[p])*N_Gauss+n]
		  = h*GaussWeights[N_Gauss-1][n];
	      }

	    // prepare Gauss points and weights in y direction
	    gauss_points1.resize(N_Gauss*(supp.ymax[p]-supp.ymin[p]));
	    gauss_weights1.resize(N_Gauss*(supp.ymax[p]-supp.ymin[p]));
	    for (int subpatch = supp.ymin[p]; subpatch < supp.ymax[p]; subpatch++)
	      for (int n = 0; n < N_Gauss; n++) {
		gauss_points1[(subpatch-supp.ymin[p])*N_Gauss+n]
		  = h*(2*subpatch+1+GaussPoints[N_Gauss-1][n])/2.;
		gauss_weights1[(subpatch-supp.ymin[p])*N_Gauss+n]
		  = h*GaussWeights[N_Gauss-1][n];
	      }

	    // evaluate the two Kronecker components of psi_mu
	    switch (it.index().p()) {
	    case 0:
	    case 1:
	    case 2:
	      // psi_mu completely lives on patch 0/1/2,
	      // this can only happen if p == it.index().p()
	      assert(p == it.index().p());
	      
	      evaluate(basis.basis1d(), 0,
		       typename IBASIS::Index(level,
					      0,
					      it.index().k()[0],
					      &basis.basis1d()),
		       gauss_points0, psi_mu_0_values);
	      
	      evaluate(basis.basis1d(), 0,
		       typename IBASIS::Index(level,
					      0,
					      it.index().k()[1],
					      &basis.basis1d()),
		       gauss_points1, psi_mu_1_values);
	      break;
	    case 3:
	      // psi_mu lives on patches 0 and 1
	      if (p == 0) {
		evaluate(basis.basis1d(), 0,
			 typename IBASIS::Index(level,
						0,
						it.index().k()[0],
						&basis.basis1d()),
			 gauss_points0, psi_mu_0_values);
		
		evaluate(basis.basis1d(), 0,
			 typename IBASIS::Index(level,
						0,
						basis.basis1d().DeltaLmin(),
						&basis.basis1d()),
			 gauss_points1, psi_mu_1_values);
	      } else {
		assert(p == 1);
		
		evaluate(basis.basis1d(), 0,
			 typename IBASIS::Index(level,
						0,
						it.index().k()[0],
						&basis.basis1d()),
			 gauss_points0, psi_mu_0_values);
		
		evaluate(basis.basis1d(), 0,
			 typename IBASIS::Index(level,
						0,
						basis.basis1d().DeltaRmax(level),
						&basis.basis1d()),
			 gauss_points1, psi_mu_1_values);
	      }
	      break;
	    case 4:
	      // psi_mu lives on patches 1 and 2
	      if (p == 1) {
		evaluate(basis.basis1d(), 0,
			 typename IBASIS::Index(level,
						0,
						basis.basis1d().DeltaRmax(level),
						&basis.basis1d()),
			 gauss_points0, psi_mu_0_values);
		
		evaluate(basis.basis1d(), 0,
			 typename IBASIS::Index(level,
						0,
						it.index().k()[1],
						&basis.basis1d()),
			 gauss_points1, psi_mu_1_values);
	      } else {
		assert(p == 2);
		
		evaluate(basis.basis1d(), 0,
			 typename IBASIS::Index(level,
						0,
						basis.basis1d().DeltaLmin(),
						&basis.basis1d()),
			 gauss_points0, psi_mu_0_values);
		
		evaluate(basis.basis1d(), 0,
			 typename IBASIS::Index(level,
						0,
						it.index().k()[1],
						&basis.basis1d()),
			 gauss_points1, psi_mu_1_values);
	      }
	      break;
	    default:
	      break;
	    }
	    
	    // determine weight factor (interface generators have to be scaled)
	    const double psi_mu_factor = it.index().p() <= 2 ? 1.0 : M_SQRT1_2;
	    
	    // determine chart action
	    const double xoffset = (p == 2 ? 0 : -1);
	    const double yoffset = (p == 0 ? 0 : -1);
	    
	    // Now we know that on patch p, the generator psi_mu is nontrivial.
	    // We have to collect all the integral shares:
	    
	    Point<2> x;
	    for (unsigned int i0 = 0; i0 < gauss_points0.size(); i0++) {
	      x[0] = gauss_points0[i0] + xoffset; // apply chart

	      const double temp = *it * gauss_weights0[i0] * psi_mu_0_values[i0] * psi_mu_factor;
	      
	      if (temp != 0) {
		for (unsigned int i1 = 0; i1 < gauss_points1.size(); i1++) {
		  x[1] = gauss_points1[i1] + yoffset; // apply chart
		  
		  r += gauss_weights1[i1]
		    * psi_mu_1_values[i1]
		    * temp
		    * f->value(x);
		}
	      }
	    }
	  }
 	}
      }
#endif    

    return r;
  }

  template <class IBASIS>
  void
  expand(const Function<2>* f,
	 const LDomainBasis<IBASIS>& basis,
	 const bool primal,
	 const int jmax,
	 InfiniteVector<double, typename LDomainBasis<IBASIS>::Index>& coeffs)
  {
    typedef typename LDomainBasis<IBASIS>::Index Index;
    const int j0 = basis.j0();
    
    for (Index lambda = first_generator(&basis, j0);;++lambda)
      {
 	coeffs.set_coefficient(lambda, integrate(f, basis, lambda));
 	if (lambda == last_wavelet(&basis, jmax))
 	  break;
      }
  }
  
}
