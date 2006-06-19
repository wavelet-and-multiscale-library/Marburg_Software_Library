// implementation for ldomain_equation.h

#include <cmath>
#include <time.h>
#include <utils/fixed_array1d.h>
#include <numerics/gauss_data.h>

namespace WaveletTL
{
  template <class IBASIS>
  LDomainEquation<IBASIS>::LDomainEquation(const EllipticBVP<2>* bvp)
    : bvp_(bvp), basis_(), normA(0.0), normAinv(0.0)
  {
    compute_rhs();
  }

  template <class IBASIS>
  LDomainEquation<IBASIS>::LDomainEquation(const LDomainEquation& eq)
    : bvp_(eq.bvp_), basis_(eq.basis_),
      fcoeffs(eq.fcoeffs), fnorm_sqr(eq.fnorm_sqr),
      normA(eq.normA), normAinv(eq.normAinv)
  {
  }

  template <class IBASIS>
  inline
  double
  LDomainEquation<IBASIS>::D(const Index& lambda) const
  {
    return ldexp(1.0, lambda.j());
//     return sqrt(a(lambda, lambda));
  }

  template <class IBASIS>
  inline
  double
  LDomainEquation<IBASIS>::a(const Index& lambda,
			     const Index& nu) const
  {
    return a(lambda, nu, 4);
  }

  template <class IBASIS>
  double
  LDomainEquation<IBASIS>::a(const Index& lambda1,
			     const Index& lambda2,
			     const unsigned int pno) const
  {
    // a(u,v) = \int_Omega [a(x)grad u(x)grad v(x)+q(x)u(x)v(x)] dx

    double r = 0;

    // first compute the support intersection of psi_lambda1 and psi_lambda2
    typename LDomainBasis<IBASIS>::Support supp, supp1, supp1cap2;
    if (intersect_supports(basis_, lambda1, lambda2, supp)) {
      // compute the generator expansions of psi_lambda1 and psi_lambda2
      InfiniteVector<double, Index> gcoeffs1, gcoeffs2;
      const int level1 = lambda1.j()+ (lambda1.e()[0]==1 || lambda1.e()[1]==1 ? 1 : 0);
      const int level2 = lambda2.j()+ (lambda2.e()[0]==1 || lambda2.e()[1]==1 ? 1 : 0);
      basis_.reconstruct_1(lambda1, level1, gcoeffs1);
      basis_.reconstruct_1(lambda2, level2, gcoeffs2);

      // iterate through the involved generators and collect the point evaluations
      const int N_Gauss = 5; // number of Gauss points in x- and y-direction (per sing. supp. subpatch)
      const double h = ldexp(1.0, -supp.j); // granularity for the quadrature
      
      Array1D<double> gauss_points0, gauss_points1, gauss_weights0, gauss_weights1,
	psi_mu1_0_values, psi_mu1_1_values,     // point values of the Kron. factors of psi_mu1 (mui are set later)
	psi_mu1_0_x_values, psi_mu1_1_x_values, // point values of the derivatives
	psi_mu2_0_values, psi_mu2_1_values,     // point values of the two Kronecker factors of psi_mu2
	psi_mu2_0_x_values, psi_mu2_1_x_values; // point values of the derivatives

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
		  
		  evaluate(basis_.basis1d(),
			   typename IBASIS::Index(level1,
						  0,
						  it1.index().k()[0],
						  &basis_.basis1d()),
			   gauss_points0,
			   psi_mu1_0_values, psi_mu1_0_x_values);
		  
		  evaluate(basis_.basis1d(),
			   typename IBASIS::Index(level1,
						  0,
						  it1.index().k()[1],
						  &basis_.basis1d()),
			   gauss_points1,
			   psi_mu1_1_values, psi_mu1_1_x_values);
		  break;
		case 3:
		  // psi_mu1 lives on patches 0 and 1
		  if (p == 0) {
		    evaluate(basis_.basis1d(),
			     typename IBASIS::Index(level1,
						    0,
						    it1.index().k()[0],
						    &basis_.basis1d()),
			     gauss_points0,
			     psi_mu1_0_values, psi_mu1_0_x_values);
		    
		    evaluate(basis_.basis1d(),
			     typename IBASIS::Index(level1,
						    0,
						    basis_.basis1d().DeltaLmin(),
						    &basis_.basis1d()),
			     gauss_points1,
			     psi_mu1_1_values, psi_mu1_1_x_values);
		  } else {
		    assert(p == 1);
		    
		    evaluate(basis_.basis1d(),
			     typename IBASIS::Index(level1,
						    0,
						    it1.index().k()[0],
						    &basis_.basis1d()),
			     gauss_points0,
			     psi_mu1_0_values, psi_mu1_0_x_values);
		    
		    evaluate(basis_.basis1d(),
			     typename IBASIS::Index(level1,
						    0,
						    basis_.basis1d().DeltaRmax(level1),
						    &basis_.basis1d()),
			     gauss_points1,
			     psi_mu1_1_values, psi_mu1_1_x_values);
		  }
		  break;
		case 4:
		  // psi_mu1 lives on patches 1 and 2
		  if (p == 1) {
		    evaluate(basis_.basis1d(),
			     typename IBASIS::Index(level1,
						    0,
						    basis_.basis1d().DeltaRmax(level1),
						    &basis_.basis1d()),
			     gauss_points0,
			     psi_mu1_0_values, psi_mu1_0_x_values);
		    
		    evaluate(basis_.basis1d(),
			     typename IBASIS::Index(level1,
						    0,
						    it1.index().k()[1],
						    &basis_.basis1d()),
			     gauss_points1,
			     psi_mu1_1_values, psi_mu1_1_x_values);
		  } else {
		    assert(p == 2);
		    
		    evaluate(basis_.basis1d(),
			     typename IBASIS::Index(level1,
						    0,
						    basis_.basis1d().DeltaLmin(),
						    &basis_.basis1d()),
			     gauss_points0,
			     psi_mu1_0_values, psi_mu1_0_x_values);
		    
		    evaluate(basis_.basis1d(),
			     typename IBASIS::Index(level1,
						    0,
						    it1.index().k()[1],
						    &basis_.basis1d()),
			     gauss_points1,
			     psi_mu1_1_values, psi_mu1_1_x_values);
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
		  
		  evaluate(basis_.basis1d(),
			   typename IBASIS::Index(level2,
						  0,
						  it2.index().k()[0],
						  &basis_.basis1d()),
			   gauss_points0,
			   psi_mu2_0_values, psi_mu2_0_x_values);
		  
		  evaluate(basis_.basis1d(),
			   typename IBASIS::Index(level2,
						  0,
						  it2.index().k()[1],
						  &basis_.basis1d()),
			   gauss_points1,
			   psi_mu2_1_values, psi_mu2_1_x_values);
		  break;
		case 3:
		  // psi_mu1 lives on patches 0 and 1
		  if (p == 0) {
		    evaluate(basis_.basis1d(),
			     typename IBASIS::Index(level2,
						    0,
						    it2.index().k()[0],
						    &basis_.basis1d()),
			     gauss_points0,
			     psi_mu2_0_values, psi_mu2_0_x_values);
		    
		    evaluate(basis_.basis1d(),
			     typename IBASIS::Index(level2,
						    0,
						    basis_.basis1d().DeltaLmin(),
						    &basis_.basis1d()),
			     gauss_points1,
			     psi_mu2_1_values, psi_mu2_1_x_values);
		  } else {
		    assert(p == 1);
		    
		    evaluate(basis_.basis1d(),
			     typename IBASIS::Index(level2,
						    0,
						    it2.index().k()[0],
						    &basis_.basis1d()),
			     gauss_points0,
			     psi_mu2_0_values, psi_mu2_0_x_values);
		    
		    evaluate(basis_.basis1d(),
			     typename IBASIS::Index(level2,
						    0,
						    basis_.basis1d().DeltaRmax(level2),
						    &basis_.basis1d()),
			     gauss_points1,
			     psi_mu2_1_values, psi_mu2_1_x_values);
		  }
		  break;
		case 4:
		  // psi_mu2 lives on patches 1 and 2
		  if (p == 1) {
		    evaluate(basis_.basis1d(),
			     typename IBASIS::Index(level2,
						    0,
						    basis_.basis1d().DeltaRmax(level2),
						    &basis_.basis1d()),
			     gauss_points0,
			     psi_mu2_0_values, psi_mu2_0_x_values);
		    
		    evaluate(basis_.basis1d(),
			     typename IBASIS::Index(level2,
						    0,
						    it2.index().k()[1],
						    &basis_.basis1d()),
			     gauss_points1,
			     psi_mu2_1_values, psi_mu2_1_x_values);
		  } else {
		    assert(p == 2);
		    
		    evaluate(basis_.basis1d(),
			     typename IBASIS::Index(level2,
						    0,
						    basis_.basis1d().DeltaLmin(),
						    &basis_.basis1d()),
			     gauss_points0,
			     psi_mu2_0_values, psi_mu2_0_x_values);
		    
		    evaluate(basis_.basis1d(),
			     typename IBASIS::Index(level2,
						    0,
						    it2.index().k()[1],
						    &basis_.basis1d()),
			     gauss_points1,
			     psi_mu2_1_values, psi_mu2_1_x_values);
		  }
		  break;
		default:
		  break;
		}
		
		// determine weight factor (interface generators have to be scaled)
		const double psi_mu1_factor = it1.index().p() <= 2 ? 1.0 : M_SQRT1_2;
		const double psi_mu2_factor = it2.index().p() <= 2 ? 1.0 : M_SQRT1_2;
		
		// determine chart action
		const double xoffset = (p == 2 ? 0 : -1);
		const double yoffset = (p == 0 ? 0 : -1);


		// Now we know that on patch p, the generator psi_mu is nontrivial.
		// We have to collect all the integral shares:
		
		Point<2> x;
		for (unsigned int i0 = 0; i0 < gauss_points0.size(); i0++) {
		  x[0] = gauss_points0[i0] + xoffset; // apply chart

		  for (unsigned int i1 = 0; i1 < gauss_points1.size(); i1++) {
		    x[1] = gauss_points1[i1] + yoffset; // apply chart
		    
		    // compute the share a(x)(grad psi_mu1)(x)(grad psi_mu2)(x)
		    r += bvp_->a(x)
		      * gauss_weights0[i0]
		      * gauss_weights1[i1]
		      * psi_mu1_factor
		      * psi_mu2_factor
		      * ((psi_mu1_0_x_values[i0] // d/dx1 psi_mu1 * d/dx1 psi_mu2
			  * psi_mu1_1_values[i1]
			  * psi_mu2_0_x_values[i0]
			  * psi_mu2_1_values[i1])
			 + (psi_mu1_0_values[i0] // d/dx2 psi_mu1 * d/dx2 psi_mu2
			    * psi_mu1_1_x_values[i1]
			    * psi_mu2_0_values[i0]
			    * psi_mu2_1_x_values[i1]))
		      * *it1
		      * *it2;
		    
		    // compute the share q(x)psi_mu1(x)psi_mu2(x)
		    r += bvp_->q(x)
		      * gauss_weights0[i0]
		      * gauss_weights1[i1]
		      * psi_mu1_factor
		      * psi_mu2_factor
		      * psi_mu1_0_values[i0]
		      * psi_mu1_1_values[i1]
		      * psi_mu2_0_values[i0]
		      * psi_mu2_1_values[i1]
		      * *it1
		      * *it2;
 		  }
		}
	      }
	    }
	  }
	}
      }
    }

    return r;
  }

  template <class IBASIS>
  double
  LDomainEquation<IBASIS>::f(const Index& lambda) const
  {
    // f(v) = \int_0^1 g(t)v(t) dt

//     cout << "LDomainEquation::f() called with lambda=" << lambda << endl;

    double r = 0;

    // compute the generator expansion of lambda
    InfiniteVector<double, Index> gcoeffs;
    const int level = lambda.j()+ (lambda.e()[0]==1 || lambda.e()[1]==1 ? 1 : 0);
    basis_.reconstruct_1(lambda, level, gcoeffs);

    // iterate through the involved generators and collect the point evaluations
    const int N_Gauss = 5; // number of Gauss points in x- and y-direction (per sing. supp. subpatch)
    const double h = ldexp(1.0, -level); // granularity for the quadrature

    typename LDomainBasis<IBASIS>::Support supp;
    Array1D<double> gauss_points0, gauss_points1, gauss_weights0, gauss_weights1,
      psi_mu_0_values, psi_mu_1_values;
    
    for (typename InfiniteVector<double,Index>::const_iterator it(gcoeffs.begin()),
  	   itend(gcoeffs.end()); it != itend; ++it)
      {
// 	cout << "LDomainEquation::f(), consider involved generator " << it.index() << endl;

 	// compute the support of the generator corresponding to mu=it.index()
  	support(basis_, it.index(), supp);
	
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
	      
	      evaluate(basis_.basis1d(), 0,
		       typename IBASIS::Index(level,
					      0,
					      it.index().k()[0],
					      &basis_.basis1d()),
		       gauss_points0, psi_mu_0_values);
	      
	      evaluate(basis_.basis1d(), 0,
		       typename IBASIS::Index(level,
					      0,
					      it.index().k()[1],
					      &basis_.basis1d()),
		       gauss_points1, psi_mu_1_values);
	      break;
	    case 3:
	      // psi_mu lives on patches 0 and 1
	      if (p == 0) {
		evaluate(basis_.basis1d(), 0,
			 typename IBASIS::Index(level,
						0,
						it.index().k()[0],
						&basis_.basis1d()),
			 gauss_points0, psi_mu_0_values);
		
		evaluate(basis_.basis1d(), 0,
			 typename IBASIS::Index(level,
						0,
						basis_.basis1d().DeltaLmin(),
						&basis_.basis1d()),
			 gauss_points1, psi_mu_1_values);
	      } else {
		assert(p == 1);
		
		evaluate(basis_.basis1d(), 0,
			 typename IBASIS::Index(level,
						0,
						it.index().k()[0],
						&basis_.basis1d()),
			 gauss_points0, psi_mu_0_values);
		
		evaluate(basis_.basis1d(), 0,
			 typename IBASIS::Index(level,
						0,
						basis_.basis1d().DeltaRmax(level),
						&basis_.basis1d()),
			 gauss_points1, psi_mu_1_values);
	      }
	      break;
	    case 4:
	      // psi_mu lives on patches 1 and 2
	      if (p == 1) {
		evaluate(basis_.basis1d(), 0,
			 typename IBASIS::Index(level,
						0,
						basis_.basis1d().DeltaRmax(level),
						&basis_.basis1d()),
			 gauss_points0, psi_mu_0_values);
		
		evaluate(basis_.basis1d(), 0,
			 typename IBASIS::Index(level,
						0,
						it.index().k()[1],
						&basis_.basis1d()),
			 gauss_points1, psi_mu_1_values);
	      } else {
		assert(p == 2);
		
		evaluate(basis_.basis1d(), 0,
			 typename IBASIS::Index(level,
						0,
						basis_.basis1d().DeltaLmin(),
						&basis_.basis1d()),
			 gauss_points0, psi_mu_0_values);
		
		evaluate(basis_.basis1d(), 0,
			 typename IBASIS::Index(level,
						0,
						it.index().k()[1],
						&basis_.basis1d()),
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
		    * bvp_->f(x);
		}
	      }
	    }
	  }
 	}
      }
    
//     cout << "LDomainEquation::f() result is r=" << r << endl;
    
    return r;
  }

  template <class IBASIS>
  void
  LDomainEquation<IBASIS>::compute_rhs()
  {
    cout << "LDomainEquation(): precompute right-hand side..." << endl;

    // precompute the right-hand side on a fine level
    InfiniteVector<double,Index> fhelp;
    const int j0   = basis().j0();
//     const int jmax = 5; // for a first quick hack
    const int jmax = j0; // for a first quick hack
    for (Index lambda(basis_.first_generator(j0));; ++lambda)
      {
	const double coeff = f(lambda)/D(lambda);
// 	if (fabs(coeff)>1e-15)
	  fhelp.set_coefficient(lambda, coeff);
  	if (lambda == basis_.last_wavelet(jmax))
	  break;
      }
    fnorm_sqr = l2_norm_sqr(fhelp);

    cout << "... done, sort the entries in modulus..." << endl;

    // sort the coefficients into fcoeffs
    fcoeffs.resize(0); // clear eventual old values
    fcoeffs.resize(fhelp.size());
    unsigned int id(0);
    for (typename InfiniteVector<double,Index>::const_iterator it(fhelp.begin()), itend(fhelp.end());
	 it != itend; ++it, ++id)
      fcoeffs[id] = std::pair<Index,double>(it.index(), *it);
    sort(fcoeffs.begin(), fcoeffs.end(), typename InfiniteVector<double,Index>::decreasing_order());
    cout << "... done, all integrals for right-hand side computed!" << endl;
  }

  template <class IBASIS>
  void
  LDomainEquation<IBASIS>::RHS
  (const double eta,
   InfiniteVector<double,Index>& coeffs) const
  {
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


}
