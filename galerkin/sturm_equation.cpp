// implementation for sturm_equation.h

#include <cmath>
#include <algorithm>
#include <utils/array1d.h>
#include <algebra/vector.h>
#include <algebra/sparse_matrix.h>
#include <numerics/eigenvalues.h>
#include <numerics/gauss_data.h>

namespace WaveletTL
{
  template <class WBASIS>
  SturmEquation<WBASIS>::SturmEquation(const simpleSturmBVP& bvp)
    : bvp_(bvp), basis_(bvp.bc_left(), bvp.bc_right())
  {
    // estimate ||A||
    typedef typename WBASIS::Index Index;
    std::set<Index> Lambda;
    const int j0 = basis_.j0();
    const int jmax = 8;
    for (Index lambda = first_generator(&basis_, j0);; ++lambda) {
      Lambda.insert(lambda);
      if (lambda == last_wavelet(&basis_, jmax)) break;
    }
    SparseMatrix<double> A_Lambda;
    setup_stiffness_matrix(Lambda, A_Lambda);

    Vector<double> xk(Lambda.size(), false);
    xk = 1;
    unsigned int iterations;
    normA = PowerIteration(A_Lambda, xk, 1e-6, 100, iterations);

    // estimate ||A^{-1}||
    xk = 1;
    normAinv = InversePowerIteration(A_Lambda, xk, 1e-6, 200, iterations);
  }

  template <class WBASIS>
  inline
  double
  SturmEquation<WBASIS>::D(const typename WBASIS::Index& lambda) const
  {
    return ldexp(1.0, lambda.j());
  }

  template <class WBASIS>
  inline
  void
  SturmEquation<WBASIS>::rescale(InfiniteVector<double, typename WBASIS::Index>& coeffs,
				 const int n) const
  {
    for (typename InfiniteVector<double, typename WBASIS::Index>::const_iterator it(coeffs.begin());
	 it != coeffs.end(); ++it)
      {
	// TODO: implement an InfiniteVector::iterator to speed up this hack!
	coeffs.set_coefficient(it.index(), *it * pow(D(it.index()), n));
      }
  }

  template <class WBASIS>
  double
  SturmEquation<WBASIS>::a(const typename WBASIS::Index& lambda,
			   const typename WBASIS::Index& nu,
			   const unsigned int p) const
  {
    // a(u,v) = \int_0^1 [p(t)u'(t)v'(t)+q(t)u(t)v(t)] dt

    double r = 0;

    // Remark: There are of course many possibilities to evaluate
    // a(u,v) numerically.
    // In this implementation, we rely on the fact that the primal functions in
    // WBASIS are splines with respect to a dyadic subgrid.
    // We can then apply an appropriate composite quadrature rule.
    // In the scope of WBASIS, the routines intersect_supports() and evaluate()
    // must exist, which is the case for DSBasis<d,dT>.

    // First we compute the support intersection of \psi_\lambda and \psi_\nu:
    int j, k1, k2;
    bool inter = intersect_supports(basis_, lambda, nu, j, k1, k2);

    if (inter)
      {
	// Set up Gauss points and weights for a composite quadrature formula:
	// (TODO: maybe use an instance of MathTL::QuadratureRule instead of computing
	// the Gauss points and weights)
	const unsigned int N_Gauss = (p+1)/2;
	const double h = ldexp(1.0, -j);
	Array1D<double> gauss_points (N_Gauss*(k2-k1)), func1values, func2values, der1values, der2values;
	for (int patch = k1, id = 0; patch < k2; patch++) // refers to 2^{-j}[patch,patch+1]
	  for (unsigned int n = 0; n < N_Gauss; n++, id++)
	    gauss_points[id] = h*(2*patch+1+GaussPoints[N_Gauss-1][n])/2.;

	// - compute point values of the integrands
	evaluate(basis_, lambda, gauss_points, func1values, der1values);
	evaluate(basis_, nu, gauss_points, func2values, der2values);

	// - add all integral shares
	for (int patch = k1, id = 0; patch < k2; patch++)
	  for (unsigned int n = 0; n < N_Gauss; n++, id++) {
	    const double t = gauss_points[id];
	    const double gauss_weight = GaussWeights[N_Gauss-1][n] * h;
	    
	    const double pt = bvp_.p(t);
	    if (pt != 0)
	      r += pt * der1values[id] * der2values[id] * gauss_weight;
	    
	    const double qt = bvp_.q(t);
	    if (qt != 0)
	      r += qt * func1values[id] * func2values[id] * gauss_weight;
	  }
      }

    return r;
  }

  template <class WBASIS>
  void
  SturmEquation<WBASIS>::setup_stiffness_matrix(const std::set<typename WBASIS::Index>& Lambda,
						SparseMatrix<double>& A_Lambda) const {
    A_Lambda.resize(Lambda.size(), Lambda.size());
    unsigned int i = 0;
    typedef typename WBASIS::Index Index;
    for (typename std::set<Index>::const_iterator it1 = Lambda.begin(), itend = Lambda.end();
	 it1 != itend; ++it1, ++i)
      {
	const double d1 = D(*it1);
	typename std::set<Index>::const_iterator it2 = Lambda.begin();
	for (; *it2 != *it1; ++it2);
	unsigned int j = i;
	for (; it2 != itend; ++it2, ++j)
	  {
	    const double d1d2 = d1 * D(*it2);
	    double entry = a(*it2, *it1);
	    if (entry != 0)
	      {
		entry /= d1d2;
		A_Lambda.set_entry(i, j, entry);
		A_Lambda.set_entry(j, i, entry); // symmetry
	      }
	  }
      }
  }

  template <class WBASIS>
  void
  SturmEquation<WBASIS>::add_column(const double factor,
				    const typename WBASIS::Index& lambda,
				    const int J,
				    InfiniteVector<double, typename WBASIS::Index>& w) const
  {
    // check whether the corresponding column already exists, create it if necessary
    typename MatrixColumnCache::iterator col_it_lb(cache_.lower_bound(lambda));
    typename MatrixColumnCache::iterator col_it(col_it_lb);
    if (col_it_lb == cache_.end() || cache_.key_comp()(lambda, col_it_lb->first))
      col_it = cache_.insert(col_it_lb, typename MatrixColumnCache::value_type(lambda, MatrixBlockCache()));

    MatrixBlockCache& col = col_it->second;

    // traverse all necessary level blocks for the FMVM, compute them if necessary

    // generator block for the coarsest level (always present)
    int level = basis_.j0()-1;
    typename MatrixBlockCache::iterator col_block_it_lb(col.lower_bound(level));
    typename MatrixBlockCache::iterator col_block_it(col_block_it_lb);
    if (col_block_it_lb == col.end() || col.key_comp()(level, col_block_it_lb->first))
      compute_matrix_block(lambda, level, col, col_block_it);
    MatrixBlock& col_block = col_block_it->second;
    for (unsigned int id = 0; id < col_block.indices.size(); ++id) {
      w.set_coefficient(col_block.indices[id],
			w.get_coefficient(col_block.indices[id])
			+ col_block.entries[id] * factor);
    }

    // wavelet blocks
    for (level = std::max(basis_.j0(), lambda.j()-J);
// 	 level <= lambda.j()+J; level++) {
	 level <= 8; level++) {
      col_block_it_lb = col.lower_bound(level);
      col_block_it    = col_block_it_lb;
      if (col_block_it_lb == col.end() || col.key_comp()(level, col_block_it_lb->first))
	compute_matrix_block(lambda, level, col, col_block_it);
      col_block = col_block_it->second;
      for (unsigned int id = 0; id < col_block.indices.size(); ++id) {
	w.set_coefficient(col_block.indices[id],
			  w.get_coefficient(col_block.indices[id])
			  + col_block.entries[id] * factor);
      }
    }
  }
  
  template <class WBASIS>
  void
  SturmEquation<WBASIS>::compute_matrix_block(const typename WBASIS::Index& lambda,
					      const int level,
					      MatrixBlockCache& col,
					      typename MatrixBlockCache::iterator& hint_and_result) const
  {
    typedef std::list<std::pair<typename WBASIS::Index, typename WBASIS::Support> > SupportList;
    SupportList nus;
    intersecting_wavelets(basis_, lambda, std::max(level, basis_.j0()), level == basis_.j0()-1, nus);

    MatrixBlock block;
    const unsigned int N = nus.size();
    block.indices.resize(N);
    block.entries.resize(N);
    
    unsigned int id = 0;
    const double d1 = D(lambda);
    for (typename SupportList::const_iterator it(nus.begin()); id < N; ++it, ++id) {
      block.indices[id] = it->first;
      block.entries[id] = a(it->first, lambda) / (d1*D(it->first));
    }

    hint_and_result = col.insert(hint_and_result, typename MatrixBlockCache::value_type(level, block));
  }

  template <class WBASIS>
  double
  SturmEquation<WBASIS>::f(const typename WBASIS::Index& lambda) const
  {
    // f(v) = \int_0^1 g(t)v(t) dt

    double r = 0;

    const int j = lambda.j()+lambda.e();
    int k1, k2;
    support(basis_, lambda, k1, k2);

    // Set up Gauss points and weights for a composite quadrature formula:
    const unsigned int N_Gauss = 5;
    const double h = ldexp(1.0, -j);
    Array1D<double> gauss_points (N_Gauss*(k2-k1)), vvalues;
    for (int patch = k1; patch < k2; patch++) // refers to 2^{-j}[patch,patch+1]
      for (unsigned int n = 0; n < N_Gauss; n++)
	gauss_points[(patch-k1)*N_Gauss+n] = h*(2*patch+1+GaussPoints[N_Gauss-1][n])/2;

    // - compute point values of the integrand
    evaluate(basis_, 0, lambda, gauss_points, vvalues);
    
    // - add all integral shares
    for (int patch = k1, id = 0; patch < k2; patch++)
      for (unsigned int n = 0; n < N_Gauss; n++, id++) {
	const double t = gauss_points[id];
	const double gauss_weight = GaussWeights[N_Gauss-1][n] * h;
	    
	const double gt = bvp_.g(t);
	if (gt != 0)
	  r += gt
	    * vvalues[id]
	    * gauss_weight;
      }
    
    return r;
  }
  
  template <class WBASIS>
  void
  SturmEquation<WBASIS>::setup_righthand_side(const std::set<typename WBASIS::Index>& Lambda,
					      Vector<double>& F_Lambda) const {
    F_Lambda.resize(Lambda.size());
    unsigned int i = 0;
    typedef typename WBASIS::Index Index;
    for (typename std::set<Index>::const_iterator it = Lambda.begin(), itend = Lambda.end();
	 it != itend; ++it, ++i)
      F_Lambda[i] = f(*it)/D(*it);
  }

  template <class WBASIS>
  void
  SturmEquation<WBASIS>::RHS(const double eta,
			     InfiniteVector<double, typename WBASIS::Index>& coeffs) const
  {
    coeffs.clear();

    // remark: for a quick hack, we use a projection of f onto a fine space V_{jmax}
    // of the given multiresolution analysis

    const int j0 = basis_.j0();
    const int jmax = 8;
    for (typename WBASIS::Index lambda(first_generator(&basis_, j0));; ++lambda)
      {
	coeffs.set_coefficient(lambda, f(lambda)/D(lambda));
  	if (lambda == last_wavelet(&basis_, jmax))
	  break;
      }
  }
}
