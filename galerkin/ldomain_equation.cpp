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
  double
  LDomainEquation<IBASIS>::f(const Index& lambda) const
  {
    // f(v) = \int_0^1 g(t)v(t) dt

    cout << "LDomainEquation::f() called with lambda=" << lambda << endl;

    double r = 0;

    // compute the generator expansion of lambda
    InfiniteVector<double, Index> gcoeffs;
    const int level = lambda.j()+ (lambda.e()[0]==1 || lambda.e()[1]==1 ? 1 : 0);
    basis_.reconstruct_1(lambda, level, gcoeffs);

//     // iterate through the involved generators and collect the point evaluations
//     const int N_Gauss = 5; // number of Gauss points in x- and y-direction (per sing. supp. subpatch)
//     const double h = ldexp(1.0, -level); // granularity for the quadrature

//     typename LDomainBasis<IBASIS>::Support supp;
//     for (typename InfiniteVector<double,Index>::const_iterator it(gcoeffs.begin()),
//  	   itend(gcoeffs.end()); it != itend; ++it)
//       {
// 	FixedArray1D<Array1D<double>,2> gauss_points, gauss_weights, v_values;

// 	// compute the support of the generator corresponding to it.index()
//  	support(basis_, it.index(), supp);

// 	// per patch, collect all point values
// 	for (int p = 0; p <= 3; p++) {
// 	  if (supp.xmin[p] != -1) {
// 	    // prepare Gauss points in x direction
// 	    gauss_points [0].resize(N_Gauss*(supp.xmax[p]-supp.xmin[p]));
// 	    gauss_weights[0].resize(N_Gauss*(supp.xmax[p]-supp.xmin[p]));
// 	    for (int subpatch = supp.xmin[p]; subpatch < supp.xmax[p]; subpatch++)
// 	      for (int n = 0; n < N_Gauss; n++) {
// 		gauss_points[0][(subpatch-supp.xmin[p])*N_Gauss+n]
// 		  = h*(2*subpatch+1+GaussPoints[N_Gauss-1][n])/2.;
// 		gauss_weights[0][(subpatch-supp.xmin[p])*N_Gauss+n]
// 		  = h*GaussWeights[N_Gauss-1][n];
// 	      }

// 	    // prepare Gauss points in y direction
// 	    gauss_points [1].resize(N_Gauss*(supp.ymax[p]-supp.ymin[p]));
// 	    gauss_weights[1].resize(N_Gauss*(supp.ymax[p]-supp.ymin[p]));
// 	    for (int subpatch = supp.ymin[p]; subpatch < supp.ymax[p]; subpatch++)
// 	      for (int n = 0; n < N_Gauss; n++) {
// 		gauss_points[1][(subpatch-supp.ymin[p])*N_Gauss+n]
// 		  = h*(2*subpatch+1+GaussPoints[N_Gauss-1][n])/2.;
// 		gauss_weights[1][(subpatch-supp.ymin[p])*N_Gauss+n]
// 		  = h*GaussWeights[N_Gauss-1][n];
// 	      }


// 	  }
// 	}

// 	r += 0;
//       }

    return r;

	
// 	// point values on patch 1 [-1,0]x[-1,0]: relevant generators live on patches 1, 3 or 4
// 	switch (it.index().p()) {
// 	case 1:
// 	  // patch generator

// 	  // compute support bounds
// 	  support(basis_.basis1d(),
// 		  typename IBASIS::Index(it.index().j(),
// 					 0,
// 					 it.index().k()[0],
// 					 &basis_.basis1d()),
// 		  xmin,
// 		  xmax);
	  
// 	  support(basis_.basis1d(),
// 		  typename IBASIS::Index(it.index().j(),
// 					 0,
// 					 it.index().k()[1],
// 					 &basis_.basis1d()),
// 		  ymin,
// 		  ymax);
	  
// 	  break;
// 	case 3:
// 	  // interface generator

// 	  // compute support bounds
// 	  support(basis_.basis1d(),
// 		  typename IBASIS::Index(it.index().j(),
// 					 0,
// 					 it.index().k()[0],
// 					 &basis_.basis1d()),
// 		  xmin,
// 		  xmax);
	  
// 	  support(basis_.basis1d(),
// 		  typename IBASIS::Index(it.index().j(),
// 					 0,
// 					 basis_.basis1d().DeltaRmax(it.index().j()),
// 					 &basis_.basis1d()),
// 		  ymin,
// 		  ymax);
	  
// 	  break;
// 	case 4:
// 	  // interface generator

// 	  // compute support bounds
// 	  support(basis_.basis1d(),
// 		  typename IBASIS::Index(it.index().j(),
// 					 0,
// 					 basis_.basis1d().DeltaRmax(it.index().j()),
// 					 &basis_.basis1d()),
// 		  xmin,
// 		  xmax);
	  
// 	  support(basis_.basis1d(),
// 		  typename IBASIS::Index(it.index().j(),
// 					 0,
// 					 it.index().k()[1],
// 					 &basis_.basis1d()),
// 		  ymin,
// 		  ymax);
	  
// 	  break;
// 	}

// 	// point values on patch 2 [ 0,1]x[-1,0]: relevant generators live on patches 2 or 4
// 	switch (it.index().p()) {
// 	case 2:
// 	  // patch generator

// 	  // compute support bounds
// 	  support(basis_.basis1d(),
// 		  typename IBASIS::Index(it.index().j(),
// 					 0,
// 					 it.index().k()[0],
// 					 &basis_.basis1d()),
// 		  xmin,
// 		  xmax);
	  
// 	  support(basis_.basis1d(),
// 		  typename IBASIS::Index(it.index().j(),
// 					 0,
// 					 it.index().k()[1],
// 					 &basis_.basis1d()),
// 		  ymin,
// 		  ymax);
	  
// 	  break;
// 	case 4:
// 	  // interface generator

// 	  // compute support bounds
// 	  support(basis_.basis1d(),
// 		  typename IBASIS::Index(it.index().j(),
// 					 0,
// 					 basis_.basis1d().DeltaLmin(),
// 					 &basis_.basis1d()),
// 		  xmin,
// 		  xmax);
	  
// 	  support(basis_.basis1d(),
// 		  typename IBASIS::Index(it.index().j(),
// 					 0,
// 					 it.index().k()[1],
// 					 &basis_.basis1d()),
// 		  ymin,
// 		  ymax);
// 	  break;
// 	}
//      }
  }

//     // compute the point values of the integrand (where we use that it is a tensor product)
//     for (unsigned int i = 0; i < DIM; i++)
//       evaluate(*basis_.bases()[i], 0,
// 	       typename IBASIS::Index(lambda.j(),
// 				      lambda.e()[i],
// 				      lambda.k()[i],
// 				      basis_.bases()[i]),
// 	       gauss_points[i], v_values[i]);

//     // iterate over all points and sum up the integral shares
//     int index[DIM]; // current multiindex for the point values
//     for (unsigned int i = 0; i < DIM; i++)
//       index[i] = 0;
    
//     Point<DIM> x;
//     while (true) {
//       for (unsigned int i = 0; i < DIM; i++)
// 	x[i] = gauss_points[i][index[i]];
//       double share = bvp_->f(x);
//       for (unsigned int i = 0; i < DIM; i++)
// 	share *= gauss_weights[i][index[i]] * v_values[i][index[i]];
//       r += share;

//       // "++index"
//       bool exit = false;
//       for (unsigned int i = 0; i < DIM; i++) {
// 	if (index[i] == N_Gauss*(supp.b[i]-supp.a[i])-1) {
// 	  index[i] = 0;
// 	  exit = (i == DIM-1);
// 	} else {
// 	  index[i]++;
// 	  break;
// 	}
//       }
//       if (exit) break;
//     }













  template <class IBASIS>
  void
  LDomainEquation<IBASIS>::compute_rhs()
  {
    cout << "LDomainEquation(): precompute right-hand side..." << endl;

    // precompute the right-hand side on a fine level
    InfiniteVector<double,Index> fhelp;
    const int j0   = basis().j0();
    const int jmax = 5; // for a first quick hack
    for (Index lambda(basis_.first_generator(j0));; ++lambda)
      {
	const double coeff = f(lambda)/D(lambda);
	if (fabs(coeff)>1e-15)
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
//     coeffs.clear();
//     double coarsenorm(0);
//     double bound(fnorm_sqr - eta*eta);
//     typedef typename WaveletBasis::Index Index;
//     typename Array1D<std::pair<Index, double> >::const_iterator it(fcoeffs.begin());
//     do {
//       coarsenorm += it->second * it->second;
//       coeffs.set_coefficient(it->first, it->second);
//       ++it;
//     } while (it != fcoeffs.end() && coarsenorm < bound);
  }


}
