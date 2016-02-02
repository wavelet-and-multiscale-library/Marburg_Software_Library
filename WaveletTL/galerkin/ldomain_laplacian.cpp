// implementation for ldomain_laplacian.h

namespace WaveletTL
{
//   template <int d, int dT>
//   LDomainLaplacian<SplineBasis<d,dT,DS_construction> >::LDomainLaplacian
//   (const WaveletBasis& basis,
//    const InfiniteVector<double,Index>& y)
//     : basis_(basis), y_(y), normA(0.0), normAinv(0.0)
//   {
//   }
//   
//   template <int d, int dT>
//   double
//   LDomainLaplacian<SplineBasis<d,dT,DS_construction> >::a
//   (const Index& lambda1,
//    const Index& lambda2) const
//   {
//     double r = 0;
// 
//     static typename LDomainBasis<IntervalBasis>::Support supp;
// 
//     // first exlude the case lambda1==lambda2
//     if (lambda1 == lambda2) {
//       basis_.support(lambda1, supp);
//  
//       // integrate over supp
//       const int N_Gauss = d; // number of Gauss points in x- and y-direction (per sing. supp. subpatch)
//       const double h = ldexp(1.0, -supp.j); // granularity for the quadrature
//       
//       Array1D<double> gauss_points0, gauss_points1, gauss_weights0, gauss_weights1,
// 	psi_lambda1_derx_values, // partial x-derivatives
// 	psi_lambda1_dery_values; // partial y-derivatives
//       
//       // per patch, collect all point values
//       for (int p = 0; p <= 2; p++) {
// 	if (supp.xmin[p] != -1) { // psi_lambda1 is nontrivial on patch p
// 	  // prepare Gauss points and weights in x direction
// 	  gauss_points0.resize(N_Gauss*(supp.xmax[p]-supp.xmin[p]));
// 	  gauss_weights0.resize(N_Gauss*(supp.xmax[p]-supp.xmin[p]));
// 	  for (int subpatch = supp.xmin[p]; subpatch < supp.xmax[p]; subpatch++)
// 	    for (int n = 0; n < N_Gauss; n++) {
// 	      gauss_points0[(subpatch-supp.xmin[p])*N_Gauss+n]
// 		= h*(2*subpatch+1+GaussPoints[N_Gauss-1][n])/2.;
// 	      gauss_weights0[(subpatch-supp.xmin[p])*N_Gauss+n]
// 		= h*GaussWeights[N_Gauss-1][n];
// 	    }
// 	  
// 	  // prepare Gauss points and weights in y direction
// 	  gauss_points1.resize(N_Gauss*(supp.ymax[p]-supp.ymin[p]));
// 	  gauss_weights1.resize(N_Gauss*(supp.ymax[p]-supp.ymin[p]));
// 	  for (int subpatch = supp.ymin[p]; subpatch < supp.ymax[p]; subpatch++)
// 	    for (int n = 0; n < N_Gauss; n++) {
// 	      gauss_points1[(subpatch-supp.ymin[p])*N_Gauss+n]
// 		= h*(2*subpatch+1+GaussPoints[N_Gauss-1][n])/2.;
// 	      gauss_weights1[(subpatch-supp.ymin[p])*N_Gauss+n]
// 		= h*GaussWeights[N_Gauss-1][n];
// 	    }
// 	  
//  	  // evaluate pullbacks of the integrand
//  	  basis_.evaluate(lambda1, p, gauss_points0, gauss_points1,
// 			  psi_lambda1_derx_values, psi_lambda1_dery_values);
// 	  
//  	  // aggregate the integral shares
//  	  for (unsigned int i0 = 0, id = 0; i0 < gauss_points0.size(); i0++) {
// 	    double help = 0;
//  	    for (unsigned int i1 = 0; i1 < gauss_points1.size(); i1++, id++) {
// 	      help += gauss_weights1[i1]
// 		* (psi_lambda1_derx_values[id] * psi_lambda1_derx_values[id]
// 		   + psi_lambda1_dery_values[id] * psi_lambda1_dery_values[id]);
//  	    }
// 	    r += gauss_weights0[i0] * help;
//  	  }
// 	}
//       }
//     } else {
//       // compute the support intersection of psi_lambda1 and psi_lambda2
//       if (intersect_supports(basis_, lambda1, lambda2, supp)) {
// 	// integrate over supp
// 	const int N_Gauss = d; // number of Gauss points in x- and y-direction (per sing. supp. subpatch)
// 	const double h = ldexp(1.0, -supp.j); // granularity for the quadrature
// 	
// 	Array1D<double> gauss_points0, gauss_points1, gauss_weights0, gauss_weights1,
// 	  psi_lambda1_derx_values, psi_lambda2_derx_values, // partial x-derivatives
// 	  psi_lambda1_dery_values, psi_lambda2_dery_values; // partial y-derivatives
// 	
// 	// per patch, collect all point values
// 	for (int p = 0; p <= 2; p++) {
// 	  if (supp.xmin[p] != -1) { // psi_lambda1 and psi_lambda2 are nontrivial on patch p
// 	    // prepare Gauss points and weights in x direction
// 	    gauss_points0.resize(N_Gauss*(supp.xmax[p]-supp.xmin[p]));
// 	    gauss_weights0.resize(N_Gauss*(supp.xmax[p]-supp.xmin[p]));
// 	    for (int subpatch = supp.xmin[p]; subpatch < supp.xmax[p]; subpatch++)
// 	      for (int n = 0; n < N_Gauss; n++) {
// 		gauss_points0[(subpatch-supp.xmin[p])*N_Gauss+n]
// 		  = h*(2*subpatch+1+GaussPoints[N_Gauss-1][n])/2.;
// 		gauss_weights0[(subpatch-supp.xmin[p])*N_Gauss+n]
// 		  = h*GaussWeights[N_Gauss-1][n];
// 	      }
// 	    
// 	    // prepare Gauss points and weights in y direction
// 	    gauss_points1.resize(N_Gauss*(supp.ymax[p]-supp.ymin[p]));
// 	    gauss_weights1.resize(N_Gauss*(supp.ymax[p]-supp.ymin[p]));
// 	    for (int subpatch = supp.ymin[p]; subpatch < supp.ymax[p]; subpatch++)
// 	      for (int n = 0; n < N_Gauss; n++) {
// 		gauss_points1[(subpatch-supp.ymin[p])*N_Gauss+n]
// 		  = h*(2*subpatch+1+GaussPoints[N_Gauss-1][n])/2.;
// 		gauss_weights1[(subpatch-supp.ymin[p])*N_Gauss+n]
// 		  = h*GaussWeights[N_Gauss-1][n];
// 	      }
// 	    
// 	    // evaluate pullbacks of the integrands
// 	    basis_.evaluate(lambda1, p, gauss_points0, gauss_points1,
// 			    psi_lambda1_derx_values, psi_lambda1_dery_values);
// 	    basis_.evaluate(lambda2, p, gauss_points0, gauss_points1,
// 			    psi_lambda2_derx_values, psi_lambda2_dery_values);	  
// 	    
// 	    // aggregate the integral shares
// 	    for (unsigned int i0 = 0, id = 0; i0 < gauss_points0.size(); i0++) {
// 	      double help = 0;
// 	      for (unsigned int i1 = 0; i1 < gauss_points1.size(); i1++, id++) {
// 		help += gauss_weights1[i1]
// 		  * (psi_lambda1_derx_values[id] * psi_lambda2_derx_values[id]
// 		     + psi_lambda1_dery_values[id] * psi_lambda2_dery_values[id]);
// 	      }
// 	      r += gauss_weights0[i0] * help;
// 	    }
// 	  }
// 	}
//       }
//     }
// 
//     return r;
//   }
}
