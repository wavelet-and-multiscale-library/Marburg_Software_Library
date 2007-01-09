// implementation for ldomain_gramian.h

namespace WaveletTL
{
  template <int d, int dT>
  LDomainGramian<SplineBasis<d,dT,DS_construction> >::LDomainGramian
  (const WaveletBasis& basis,
   const InfiniteVector<double,Index>& y)
    : basis_(basis), y_(y), normA(0.0), normAinv(0.0)
  {
  }

  template <int d, int dT>
  double
  LDomainGramian<SplineBasis<d,dT,DS_construction> >::a
  (const Index& lambda1,
   const Index& lambda2) const
  {
    double r = 0;

    // first compute the support intersection of psi_lambda1 and psi_lambda2
    typename LDomainBasis<IntervalBasis>::Support supp;
    if (intersect_supports(basis_, lambda1, lambda2, supp)) {

#if 0
      // new variant, integrate over supp
      const int N_Gauss = d+1; // number of Gauss points in x- and y-direction (per sing. supp. subpatch)
      const double h = ldexp(1.0, -supp.j); // granularity for the quadrature
      
      Array1D<double> gauss_points0, gauss_points1, gauss_weights0, gauss_weights1,
	psi_lambda1_values, psi_lambda2_values; // point values of the integrands

      // per patch, collect all point values
      for (int p = 0; p <= 2; p++) {
	if (supp.xmin[p] != -1) { // psi_lambda1 and psi_lambda2 are nontrivial on patch p

	  // TODO!!!

	}
      }
#else
      // compute the generator expansions of psi_lambda1 and psi_lambda2
      InfiniteVector<double, Index> gcoeffs1, gcoeffs2;
      const int ecode1 = lambda1.e()[0]+2*lambda1.e()[1];
      const int level1 = lambda1.j() + (ecode1 == 0 ? 0 : 1);
      if (ecode1 == 0)
	gcoeffs1.set_coefficient(lambda1, 1.0);
      else
	basis_.reconstruct_1(lambda1, level1, gcoeffs1);

      const int ecode2 = lambda2.e()[0]+2*lambda2.e()[1];
      const int level2 = lambda2.j() + (ecode2 == 0 ? 0 : 1);
      if (ecode2 == 0)
	gcoeffs2.set_coefficient(lambda2, 1.0);
      else
	basis_.reconstruct_1(lambda2, level2, gcoeffs2);

      // iterate through the involved generators and collect the point evaluations
      const int N_Gauss = d+1; // number of Gauss points in x- and y-direction (per sing. supp. subpatch)
      const double h = ldexp(1.0, -supp.j); // granularity for the quadrature
      
      Array1D<double> gauss_points0, gauss_points1, gauss_weights0, gauss_weights1,
	psi_mu1_0_values, psi_mu1_1_values,     // point values of the Kron. factors of psi_mu1 (mui are set later)
	psi_mu2_0_values, psi_mu2_1_values;     // point values of the two Kronecker factors of psi_mu2

      // first variant: compute support intersection of psi_mu1 and psi_mu2, integration over that domain 
      typename LDomainBasis<IntervalBasis>::Support supp1, supp1cap2;
      for (typename InfiniteVector<double,Index>::const_iterator it1(gcoeffs1.begin()),
	     itend1(gcoeffs1.end()); it1 != itend1; ++it1) {
 	// compute the support of the generator corresponding to mu1=it1.index()
  	support(basis_, it1.index(), supp1);
	for (typename InfiniteVector<double,Index>::const_iterator it2(gcoeffs2.begin()),
	       itend2(gcoeffs2.end()); it2 != itend2; ++it2) {
	  // compute the intersection of supp(psi_mu2) with supp1, mu2=it2.index()
	  if (intersect_supports(basis_, it2.index(), supp1, supp1cap2)) {
	    // per patch, collect all point values
	    for (int p = 0; p <= 2; p++) {
	      if (supp1cap2.xmin[p] != -1) { // psi_mu1 and psi_mu2 are nontrivial on patch p
		// prepare Gauss points and weights in x direction
		gauss_points0.resize(N_Gauss*(supp1cap2.xmax[p]-supp1cap2.xmin[p]));
		gauss_weights0.resize(N_Gauss*(supp1cap2.xmax[p]-supp1cap2.xmin[p]));
		for (int subpatch = supp1cap2.xmin[p]; subpatch < supp1cap2.xmax[p]; subpatch++)
		  for (int n = 0; n < N_Gauss; n++) {
		    gauss_points0[(subpatch-supp1cap2.xmin[p])*N_Gauss+n]
		      = h*(2*subpatch+1+GaussPoints[N_Gauss-1][n])/2.;
		    gauss_weights0[(subpatch-supp1cap2.xmin[p])*N_Gauss+n]
		      = h*GaussWeights[N_Gauss-1][n];
		  }
		
		// prepare Gauss points and weights in y direction
		gauss_points1.resize(N_Gauss*(supp1cap2.ymax[p]-supp1cap2.ymin[p]));
		gauss_weights1.resize(N_Gauss*(supp1cap2.ymax[p]-supp1cap2.ymin[p]));
		for (int subpatch = supp1cap2.ymin[p]; subpatch < supp1cap2.ymax[p]; subpatch++)
		  for (int n = 0; n < N_Gauss; n++) {
		    gauss_points1[(subpatch-supp1cap2.ymin[p])*N_Gauss+n]
		      = h*(2*subpatch+1+GaussPoints[N_Gauss-1][n])/2.;
		    gauss_weights1[(subpatch-supp1cap2.ymin[p])*N_Gauss+n]
		      = h*GaussWeights[N_Gauss-1][n];
		  }
		
		// evaluate the two Kronecker components of psi_mu1
		switch (it1.index().p()) {
		case 0:
		case 1:
		case 2:
		  // psi_mu1 completely lives on patch 0/1/2,
		  // this can only happen if p == it.index().p()
		  assert(p == it1.index().p());
		  
		  basis_.basis1d().evaluate(0,
					    level1, 0, it1.index().k()[0],
					    gauss_points0,
					    psi_mu1_0_values);
		  
		  basis_.basis1d().evaluate(0,
					    level1, 0, it1.index().k()[1],
					    gauss_points1,
					    psi_mu1_1_values);
		  break;
		case 3:
		  // psi_mu1 lives on patches 0 and 1
		  if (p == 0) {
		    basis_.basis1d().evaluate(0,
					      level1, 0, it1.index().k()[0],
					      gauss_points0,
					      psi_mu1_0_values);
		    
		    basis_.basis1d().evaluate(0,
					      level1, 0, basis_.basis1d().DeltaLmin(),
					      gauss_points1,
					      psi_mu1_1_values);
		  } else {
		    assert(p == 1);
		    
		    basis_.basis1d().evaluate(0,
					      level1, 0, it1.index().k()[0],
					      gauss_points0,
					      psi_mu1_0_values);
		    
		    basis_.basis1d().evaluate(0,
					      level1, 0, basis_.basis1d().DeltaRmax(level1),
					      gauss_points1,
					      psi_mu1_1_values);
		  }
		  break;
		case 4:
		  // psi_mu1 lives on patches 1 and 2
		  if (p == 1) {
		    basis_.basis1d().evaluate(0,
					      level1, 0, basis_.basis1d().DeltaRmax(level1),
					      gauss_points0,
					      psi_mu1_0_values);
		    
		    basis_.basis1d().evaluate(0,
					      level1, 0, it1.index().k()[1],
					      gauss_points1,
					      psi_mu1_1_values);
		  } else {
		    assert(p == 2);
		    
		    basis_.basis1d().evaluate(0,
					      level1, 0, basis_.basis1d().DeltaLmin(),
					      gauss_points0,
					      psi_mu1_0_values);
		    
		    basis_.basis1d().evaluate(0,
					      level1, 0, it1.index().k()[1],
					      gauss_points1,
					      psi_mu1_1_values);
		  }
		  break;
		default:
		  break;
		}
		
		// evaluate the two Kronecker components of psi_mu2
		switch (it2.index().p()) {
		case 0:
		case 1:
		case 2:
		  // psi_mu2 completely lives on patch 0/1/2,
		  // this can only happen if p == it.index().p()
		  assert(p == it2.index().p());
		  
		  basis_.basis1d().evaluate(0,
					    level2, 0, it2.index().k()[0],
					    gauss_points0,
					    psi_mu2_0_values);
		  
		  basis_.basis1d().evaluate(0,
					    level2, 0, it2.index().k()[1],
					    gauss_points1,
					    psi_mu2_1_values);
		  break;
		case 3:
		  // psi_mu1 lives on patches 0 and 1
		  if (p == 0) {
		    basis_.basis1d().evaluate(0,
					      level2, 0,it2.index().k()[0],
					      gauss_points0,
					      psi_mu2_0_values);
		    
		    basis_.basis1d().evaluate(0,
					      level2, 0, basis_.basis1d().DeltaLmin(),
					      gauss_points1,
					      psi_mu2_1_values);
		  } else {
		    assert(p == 1);
		    
		    basis_.basis1d().evaluate(0,
					      level2, 0, it2.index().k()[0],
					      gauss_points0,
					      psi_mu2_0_values);
		    
		    basis_.basis1d().evaluate(0,
					      level2, 0, basis_.basis1d().DeltaRmax(level2),
					      gauss_points1,
					      psi_mu2_1_values);
		  }
		  break;
		case 4:
		  // psi_mu2 lives on patches 1 and 2
		  if (p == 1) {
		    basis_.basis1d().evaluate(0,
					      level2, 0, basis_.basis1d().DeltaRmax(level2),
					      gauss_points0,
					      psi_mu2_0_values);
		    
		    basis_.basis1d().evaluate(0,
					      level2, 0, it2.index().k()[1],
					      gauss_points1,
					      psi_mu2_1_values);
		  } else {
		    assert(p == 2);
		    
		    basis_.basis1d().evaluate(0,
					      level2, 0, basis_.basis1d().DeltaLmin(),
					      gauss_points0,
					      psi_mu2_0_values);
		    
		    basis_.basis1d().evaluate(0,
					      level2, 0, it2.index().k()[1],
					      gauss_points1,
					      psi_mu2_1_values);
		  }
		  break;
		default:
		  break;
		}
		
		// determine weight factors (interface generators have to be scaled)
		const double psi_mu1_factor = it1.index().p() <= 2 ? 1.0 : M_SQRT1_2;
		const double psi_mu2_factor = it2.index().p() <= 2 ? 1.0 : M_SQRT1_2;

		// Now we know that on patch p, the generator psi_mu is nontrivial.
		// We have to collect all the integral shares:
		for (unsigned int i0 = 0; i0 < gauss_points0.size(); i0++) {
		  const double temp = gauss_weights0[i0] * psi_mu1_factor * psi_mu2_factor * *it1 * *it2;
		  const double temp1= psi_mu1_0_values[i0]   * psi_mu2_0_values[i0];

		  for (unsigned int i1 = 0; i1 < gauss_points1.size(); i1++) {
		    // compute the share psi_mu1(x)psi_mu2(x)
		    r += temp1
		      * gauss_weights1[i1] * temp
		      * psi_mu1_1_values[i1] * psi_mu2_1_values[i1];
 		  }
		}
	      }
	    }
	  }
	}
      }
#endif
    } // if intersect_supports(...

    return r;
  }
}
