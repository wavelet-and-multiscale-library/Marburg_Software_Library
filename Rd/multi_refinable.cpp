// implementation for multi_refinable.h

#include <algorithm>
#include <cmath>
#include <iostream>
#include <set>
#include <utils/fixed_array1d.h>
#include <utils/tiny_tools.h>
#include <algebra/matrix.h>
#include <algebra/triangular_matrix.h>
#include <numerics/eigenvalues.h>
#include <numerics/differences.h>
#include <numerics/matrix_decomp.h>

using std::cout;
using std::endl;
using namespace MathTL;

namespace WaveletTL
{
  template <class MASK, unsigned int DIMENSION>
  inline
  InfiniteVector<double, MultiIndex<int, DIMENSION> >
  MultivariateRefinableFunction<MASK, DIMENSION>::evaluate() const
  {
    return evaluate(MultiIndex<unsigned int, DIMENSION>());
  }

  template <class MASK, unsigned int DIMENSION>
  InfiniteVector<double, MultiIndex<int, DIMENSION> >
  MultivariateRefinableFunction<MASK, DIMENSION>::evaluate
  (const MultiIndex<unsigned int, DIMENSION>& mu) const
  {
    InfiniteVector<double, MultiIndex<int, DIMENSION> > r;

    // compute a support cube
    int suppleft(MultivariateLaurentPolynomial<double, DIMENSION>::begin().index()[0]);
    int suppright(suppleft);
    
    for (typename MASK::const_iterator it(MultivariateLaurentPolynomial<double, DIMENSION>::begin());
	 it != MultivariateLaurentPolynomial<double, DIMENSION>::end(); ++it)
      {
 	for (unsigned int i(0); i < DIMENSION; i++)
 	  {
 	    suppleft = std::min(suppleft, it.index()[i]);
 	    suppright = std::max(suppright, it.index()[i]);
 	  }
      }

    // for convenience, collect all integer points from the interior of the support cube
    MultiIndex<int, DIMENSION> alpha, beta;
    for (unsigned int i(0); i < DIMENSION; i++)
      {
	alpha[i] = suppleft+1;
	beta[i] = suppright-1;
      }
    std::set<MultiIndex<int, DIMENSION> > indices
      (cuboid_indices<int, DIMENSION>(alpha, beta));

    // In the following, we set up the eigenvalue problem for the values
    //   V_\alpha := D^\mu\phi(\alpha), \alpha\in\mathbb Z^d
    // The eigenvector is determined uniquely by the following equations,
    // see [DM] for details:
    //
    // (3.22) eigenvalue condition
    //   2^{-|\mu|}V_\alpha = \sum_\beta a_{2\alpha-\beta}V_\beta, \alpha\in\mathbb Z^d
    //
    // (3.23) orthogonality condition
    //   \sum_\alpha (-\alpha)^\nu V_\alpha = \mu!\delta_{\mu,\nu}, |\nu|\le|\mu|

    unsigned int degmu(multi_degree(mu)), facmu(multi_faculty(mu));

    // we also collect all \nu\in\mathbb N^d, such that |\nu|\le\mu|
    std::set<MultiIndex<unsigned int, DIMENSION> > nus;
    for (unsigned int degnu(0); degnu <= degmu; degnu++)
      {
	std::set<MultiIndex<unsigned int, DIMENSION> > sofar(nus);
	std::set<MultiIndex<unsigned int, DIMENSION> > plus(degree_indices<DIMENSION>(degnu));
	std::set_union(sofar.begin(), sofar.end(),
		       plus.begin(), plus.end(),
		       inserter(nus, nus.begin()));
      }

    // for safety, re-compute the number of multiindices \nu\in\mathbb N^d, such that |\nu|\le|\mu|
    unsigned int extra_rows(0);
    for (unsigned int k(0); k <= degmu; k++)
      extra_rows += binomial(DIMENSION+k-1, k); // \#\{\nu\in\mathbb N^d: |\nu|=k\}

    Matrix<double> A(indices.size() + nus.size(), indices.size());
    Vector<double> b(indices.size() + nus.size());

    unsigned int m, n;
    typename std::set<MultiIndex<int, DIMENSION> >::const_iterator rowit1, colit;
    for (rowit1 = indices.begin(), m = 0; rowit1 != indices.end(); ++rowit1, m++)
      {
	// (3.22)
	for (colit = indices.begin(), n = 0; colit != indices.end(); ++colit, n++)
	  {
	    MultiIndex<int, DIMENSION> index;
	    for (unsigned int i(0); i < DIMENSION; i++)
	      index[i] = 2*(*rowit1)[i] - (*colit)[i];
	    
	    A(m, n) = MultivariateLaurentPolynomial<double, DIMENSION>::get_coefficient(index);
	  }
	A(m, m) -= ldexp(1.0, -degmu);
	// b[m] = 0;
      }
    
    typename std::set<MultiIndex<unsigned int, DIMENSION> >::const_iterator rowit2;
    for (rowit2 = nus.begin(), m = indices.size(); rowit2 != nus.end(); ++rowit2, m++)
      {
	// (3.23)
	unsigned int degnu(multi_degree(*rowit2));
	for (colit = indices.begin(), n = 0; colit != indices.end(); ++colit, n++)
	  {
	    A(m, n) = minus1power(degnu) * multi_power(alpha, *rowit2);
	  }
	b[m] = (*rowit2 == mu ? facmu : 0);
      }

    QRDecomposition<double> qr(A);
    assert(qr.hasFullRank());
    Vector<double> x;
    qr.solve(b, x);
    x.compress(1e-15);

    for (colit = indices.begin(), n = 0; colit != indices.end(); ++colit, n++)
      r.set_coefficient(*colit, x[n]);

    return r;
  }
}
