// implementation for refinable.h

#include <algorithm>
#include <cmath>
#include <iostream>
#include <set>

#include <utils/array1d.h>
#include <utils/tiny_tools.h>
#include <algebra/matrix.h>
#include <algebra/triangular_matrix.h>
#include <algebra/multi_laurent_polynomial.h>
#include <numerics/matrix_decomp.h>
#include <numerics/eigenvalues.h>
#include <numerics/differences.h>

namespace WaveletTL
{
  template <class MASK, unsigned int DIMENSION>
  inline
  InfiniteVector<double, MultiIndex<int, DIMENSION> >
  MultivariateRefinableFunction<MASK, DIMENSION>::evaluate
  (const int resolution) const
  {
    return evaluate(MultiIndex<unsigned int, DIMENSION>(),
		    resolution);
  }

  template <class MASK, unsigned int DIMENSION>
  InfiniteVector<double, MultiIndex<int, DIMENSION> >
  MultivariateRefinableFunction<MASK, DIMENSION>::evaluate
  (const MultiIndex<unsigned int, DIMENSION>& mu,
   const int resolution) const
  {
    InfiniteVector<double, MultiIndex<int, DIMENSION> > r;
    
    // First we calculate the values on \mathbb Z^d and write the result to r.

    // compute a support cube
    int suppleft(MultivariateLaurentPolynomial<double, DIMENSION>::begin().index()[0]);
    int suppright(suppleft);
    
    for (typename MASK::const_iterator it(MultivariateLaurentPolynomial<double, DIMENSION>::begin());
	 it != MultivariateLaurentPolynomial<double, DIMENSION>::end(); ++it) {
      for (unsigned int i(0); i < DIMENSION; i++) {
	suppleft = std::min(suppleft, it.index()[i]);
	suppright = std::max(suppright, it.index()[i]);
      }
    }
    
    // exclude the special case of \chi_{[0,1)^d}
    if (suppleft == suppright-1) {
      r.set_coefficient(MultiIndex<int, DIMENSION>(), 1);
    } else {
      // for convenience, collect all integer points from the interior of the support cube
      MultiIndex<int, DIMENSION> alpha, beta;
      for (unsigned int i(0); i < DIMENSION; i++) {
	alpha[i] = suppleft+1;
	beta[i] = suppright-1;
      }
      std::set<MultiIndex<int, DIMENSION> > indices
	(cuboid_indices<int, DIMENSION>(alpha, beta));
	
      // In the following, we set up the eigenvalue problem for the values
      //   V_\alpha := D^\mu\phi(\alpha), \alpha\in\mathbb Z^d
      // The eigenvector is determined uniquely by the following equations,
      // see [DM] for details (so we don't need an iterative scheme):
      //
      // (3.22) eigenvalue condition
      //   2^{-|\mu|}V_\alpha = \sum_\beta a_{2\alpha-\beta}V_\beta, \alpha\in\mathbb Z^d
      //
      // (3.23) orthogonality condition
      //   \sum_\alpha (-\alpha)^\nu V_\alpha = \mu!\delta_{\mu,\nu}, |\nu|\le|\mu|
	
      unsigned int degmu(multi_degree(mu)), facmu(multi_faculty(mu));
	
      // we also collect all \nu\in\mathbb N^d, such that |\nu|\le\mu|
      std::set<MultiIndex<unsigned int, DIMENSION> > nus;
      for (unsigned int degnu(0); degnu <= degmu; degnu++) {
	std::set<MultiIndex<unsigned int, DIMENSION> > sofar(nus);
	std::set<MultiIndex<unsigned int, DIMENSION> > plus(degree_indices<DIMENSION>(degnu));
	std::set_union(sofar.begin(), sofar.end(),
		       plus.begin(), plus.end(),
		       inserter(nus, nus.begin()));
      }
	
      Matrix<double> A(indices.size() + nus.size(), indices.size());
      Vector<double> b(indices.size() + nus.size());
	
      unsigned int m, n;
      typename std::set<MultiIndex<int, DIMENSION> >::const_iterator rowit1, colit;
      for (rowit1 = indices.begin(), m = 0; rowit1 != indices.end(); ++rowit1, m++) {
	// (3.22)
	for (colit = indices.begin(), n = 0; colit != indices.end(); ++colit, n++) {
	  MultiIndex<int, DIMENSION> index;
	  for (unsigned int i(0); i < DIMENSION; i++)
	    index[i] = 2*(*rowit1)[i] - (*colit)[i];
	    
	  A(m, n) = MultivariateLaurentPolynomial<double, DIMENSION>::get_coefficient(index);
	}
	A(m, m) -= ldexp(1.0, -degmu);
	// b[m] = 0;
      }
	
      typename std::set<MultiIndex<unsigned int, DIMENSION> >::const_iterator rowit2;
      for (rowit2 = nus.begin(), m = indices.size(); rowit2 != nus.end(); ++rowit2, m++) {
	// (3.23)
	unsigned int degnu(multi_degree(*rowit2));
	for (colit = indices.begin(), n = 0; colit != indices.end(); ++colit, n++)
	  {
	    A(m, n) = minus1power(degnu) * multi_power(*colit, *rowit2);
	  }
	b[m] = (*rowit2 == mu ? facmu : 0);
      }
	
      // the system matrix is rectangular, but Ax=b is solvable via a QR decomposition
      QRDecomposition<double> qr(A);
      assert(qr.hasFullRank());
      Vector<double> x;
      qr.solve(b, x);
      x.compress(1e-15);
	
      // reinterpret the entries of x
      for (colit = indices.begin(), n = 0; colit != indices.end(); ++colit, n++)
	r.set_coefficient(*colit, x[n]);
	
    }   

    // For the remaining points we use the refinement relation of phi
    if (resolution > 0) {
      for (int newres(1); newres <= resolution; newres++) {
	// copy the coarse values \phi(2^{-j}m) = \phi(2^{-(j+1)}2m)
	InfiniteVector<double, MultiIndex<int, DIMENSION> > coarse(r);
	r.clear();
	for (typename InfiniteVector<double, MultiIndex<int, DIMENSION> >::const_iterator it(coarse.begin());
	     it != coarse.end(); ++it) {
	  MultiIndex<int, DIMENSION> index;
	  for (unsigned int i(0); i < DIMENSION; i++)
	    index[i] = 2 * it.index()[i];
	  r.set_coefficient(index, *it);
	}
	
	// \phi(2^{-(j+1)}m) = \sum_k a_k\phi(2^{-j}(m-2^jk))
	
	// compute the set of all entries 2^{-(j+1)}m, such that m has at least one odd entry
	MultiIndex<int, DIMENSION> a, b;
	for (unsigned int i(0); i < DIMENSION; i++) {
	  a[i] = (1<<(newres-1)) * suppleft;
	  b[i] = (1<<(newres-1)) * suppright - 1;
	}
	std::set<MultiIndex<int, DIMENSION> > help(cuboid_indices<int, DIMENSION>(a, b));
	for (typename std::set<MultiIndex<int, DIMENSION> >::const_iterator it(help.begin());
	     it != help.end(); ++it)
	  {
	    MultiIndex<int, DIMENSION> m, l;
	    for (unsigned int i(0); i < DIMENSION; i++)
	      m[i] = 2*(*it)[i]+1;

	    // compute phi(2^{-newres}*m)
	    for (typename MASK::const_iterator maskit(MultivariateLaurentPolynomial<double, DIMENSION>::begin());
		 maskit != MultivariateLaurentPolynomial<double, DIMENSION>::end(); ++maskit) {
	      for (unsigned int i(0); i < DIMENSION; i++)
		l[i] = m[i] - (1<<(newres-1)) * maskit.index()[i];
	      r[m] += *maskit * coarse.get_coefficient(l);
	    }
	  }
      }
    }

    return r;
  }

  template <class MASK, unsigned int DIMENSION>
  SampledMapping<DIMENSION>
  MultivariateRefinableFunction<MASK, DIMENSION>::evaluate
  (const MultiIndex<unsigned int, DIMENSION>& mu,
   const int j,
   const MultiIndex<int, DIMENSION>& k,
   const MultiIndex<int, DIMENSION>& a,
   const MultiIndex<int, DIMENSION>& b,
   const int resolution) const
  {
    InfiniteVector<double, MultiIndex<int, DIMENSION> > help(evaluate(mu, resolution-j)), v;
    const double factor(sqrt(ldexp(1.0, j*DIMENSION)) * ldexp(1.0, j*multi_degree(mu)));
    for (typename InfiniteVector<double, MultiIndex<int, DIMENSION> >::const_iterator it(help.begin());
	 it != help.end(); ++it) {
      MultiIndex<int, DIMENSION> index;
      for (unsigned int i(0); i < DIMENSION; i++)
	index[i] = it.index()[i] + (1<<(resolution-j))*k[i];
      v.set_coefficient(index, *it * factor);
    }

    return SampledMapping<DIMENSION>(a, b, v, resolution);
  }

  template <class MASK, unsigned int DIMENSION>
  double
  MultivariateRefinableFunction<MASK, DIMENSION>::moment
  (const MultiIndex<unsigned int, DIMENSION>& alpha) const
  {
    double r(1.0);
    
    int degree(multi_degree(alpha));
    if (degree > 0)
      {
	r = 0.0;

	// collect all multiindices \beta, such that (0,...,0)\le\beta\le\alpha
	std::set<MultiIndex<unsigned int, DIMENSION> > indices
	  (cuboid_indices<unsigned int, DIMENSION>(MultiIndex<unsigned int, DIMENSION>(), alpha));
	indices.erase(alpha);

	for (typename std::set<MultiIndex<unsigned int, DIMENSION> >::const_iterator it(indices.begin());
	     it != indices.end(); ++it)
	  {
	    MultiIndex<unsigned int, DIMENSION> alphaminusbeta;
	    for (unsigned int i(0); i < DIMENSION; i++)
	      alphaminusbeta[i] = alpha[i] - (*it)[i];
	    
	    double calphabeta(0.0);
 	    for (typename MASK::const_iterator maskit(MultivariateLaurentPolynomial<double, DIMENSION>::begin());
 		 maskit != MultivariateLaurentPolynomial<double, DIMENSION>::end(); ++maskit) {
 	      calphabeta += *maskit * multi_power(maskit.index(), alphaminusbeta);
 	    }
	    
	    r += ldexp(1.0, -(int)(degree+DIMENSION)) * calphabeta * multi_binomial(alpha, *it)
	      * moment(*it);
	  }
	r /= 1.0 - ldexp(1.0, -degree);
      }

    return r;
  }
}
