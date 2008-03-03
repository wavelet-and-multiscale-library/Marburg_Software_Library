// implementation of PeriodicBasis methods

#include <cmath>
#include <algebra/sparse_matrix.h>
#include <geometry/point.h>
#include <geometry/grid.h>
#include <utils/tiny_tools.h>
#include <numerics/gauss_data.h>
#include <numerics/iteratsolv.h>

using MathTL::Grid;
using MathTL::SampledMapping;
using MathTL::Point;
using MathTL::SparseMatrix;

namespace WaveletTL
{
  template <class RBASIS>
  PeriodicBasis<RBASIS>::PeriodicBasis()
    : r_basis(),
      Mj0_(j0(), r_basis.offset_a, r_basis.band_a, M_SQRT1_2),
      Mj1_(j0(), r_basis.offset_b, r_basis.band_b, M_SQRT1_2),
      Mj0T_(j0(), r_basis.offset_aT, r_basis.band_aT, M_SQRT1_2),
      Mj1T_(j0(), r_basis.offset_bT, r_basis.band_bT, M_SQRT1_2)
  {
//     cout << "Mj0=" << endl << Mj0_ << endl;
//     cout << "Mj0T=" << endl << Mj0T_ << endl;
//     cout << "Mj1=" << endl << Mj1_ << endl;
//     cout << "Mj1T=" << endl << Mj1T_ << endl;
  }
  
  template <class RBASIS>
  inline
  typename PeriodicBasis<RBASIS>::Index
  PeriodicBasis<RBASIS>::first_generator(const int j)
  {
    assert(j >= j0());
    return Index(j, 0, PeriodicBasis<RBASIS>::DeltaLmin());
  }
  
  template <class RBASIS>
  inline
  typename PeriodicBasis<RBASIS>::Index
  PeriodicBasis<RBASIS>::last_generator(const int j)
  {
    assert(j >= j0());
    return Index(j, 0, PeriodicBasis<RBASIS>::DeltaRmax(j));
  }

  template <class RBASIS>
  inline
  typename PeriodicBasis<RBASIS>::Index
  PeriodicBasis<RBASIS>::first_wavelet(const int j)
  {
    assert(j >= j0());
    return Index(j, 1, PeriodicBasis<RBASIS>::Nablamin());
  }
  
  template <class RBASIS>
  inline
  typename PeriodicBasis<RBASIS>::Index
  PeriodicBasis<RBASIS>::last_wavelet(const int j)
  {
    assert(j >= j0());
    return Index(j, 1, PeriodicBasis<RBASIS>::Nablamax(j));
  }

  template <class RBASIS>
  inline
  typename PeriodicBasis<RBASIS>::Index
  PeriodicBasis<RBASIS>::first_index(const int j, const int e)
  {
    return (e == 0 ? first_generator(j) : first_wavelet(j));
  }
  
  template <class RBASIS>
  inline
  typename PeriodicBasis<RBASIS>::Index
  PeriodicBasis<RBASIS>::last_index(const int j, const int e)
  {
    return (e == 0 ? last_generator(j) : last_wavelet(j));
  }

  template <class RBASIS>
  template <class V>
  inline
  void
  PeriodicBasis<RBASIS>::apply_Mj0(const int j, const V& x, V& y,
				   const size_type x_offset, const size_type y_offset,
				   const bool add_to) const
  {
    Mj0_.set_level(j);
    Mj0_.apply(x, y, x_offset, y_offset, add_to);
  }

  template <class RBASIS>
  template <class V>
  inline
  void
  PeriodicBasis<RBASIS>::apply_Mj1(const int j, const V& x, V& y,
				   const size_type x_offset, const size_type y_offset,
				   const bool add_to) const
  {
    Mj1_.set_level(j);
    Mj1_.apply(x, y, x_offset, y_offset, add_to);
  }

  template <class RBASIS>
  template <class V>
  void
  PeriodicBasis<RBASIS>::apply_Mj(const int j, const V& x, V& y) const
  {
    Mj0_.set_level(j);
    Mj1_.set_level(j);
    
    // decompose x=(x1 x2) appropriately
    Mj0_.apply(x, y, 0, 0, false);                      // apply Mj0 to first block x1
    Mj1_.apply(x, y, Mj0_.column_dimension(), 0, true); // apply Mj1 to second block x2 and add result
  }

  template <class RBASIS>
  template <class V>
  void
  PeriodicBasis<RBASIS>::apply_Mj_transposed(const int j, const V& x, V& y) const
  {
    Mj0_.set_level(j);
    Mj1_.set_level(j);

    // y=(y1 y2) is a block vector
    Mj0_.apply_transposed(x, y, 0, 0, false);                       // write into first block y1
    Mj1_.apply_transposed(x, y, 0, Mj0_.column_dimension(), false); // write into second block y2
  }

  template <class RBASIS>
  template <class V>
  void
  PeriodicBasis<RBASIS>::apply_Gj(const int j, const V& x, V& y) const
  {
    Mj0T_.set_level(j);
    Mj1T_.set_level(j);

    // y=(y1 y2) is a block vector
    Mj0T_.apply_transposed(x, y, 0, 0, false);                        // write into first block y1
    Mj1T_.apply_transposed(x, y, 0, Mj0T_.column_dimension(), false); // write into second block y2
  }

  template <class RBASIS>
  template <class V>
  void
  PeriodicBasis<RBASIS>::apply_Gj_transposed(const int j, const V& x, V& y) const
  {
    Mj0T_.set_level(j);
    Mj1T_.set_level(j);
    
    // decompose x=(x1 x2) appropriately
    Mj0T_.apply(x, y, 0, 0, false);                       // apply Mj0T to first block x1
    Mj1T_.apply(x, y, Mj0T_.column_dimension(), 0, true); // apply Mj1T to second block x2 and add result
  }
  
  template <class RBASIS>
  template <class V>
  void
  PeriodicBasis<RBASIS>::apply_Tj(const int j, const V& x, V& y) const
  { 
    y = x;
    V z(x);
    apply_Mj(j0(), z, y);
    for (int k = j0()+1; k <= j; k++) {
      apply_Mj(k, y, z);
      y.swap(z);
    }
  }

  template <class RBASIS>
  template <class V>
  void
  PeriodicBasis<RBASIS>::apply_Tjinv(const int j, const V& x, V& y) const
  { 
    // T_j^{-1}=diag(G_{j0},I)*...*diag(G_{j-1},I)*G_j
    V z(x);
    apply_Gj(j, x, y);
    for (int k = j-1; k >= j0(); k--) {
      z.swap(y);
      apply_Gj(k, z, y);
      for (int i = Deltasize(k+1); i < Deltasize(j+1); i++)
	y[i] = z[i];
    }
  }
  
  template <class RBASIS>
  void
  PeriodicBasis<RBASIS>::support(const Index& lambda, int& k1, int& k2)
  {
    if (lambda.e() == 0) // generator
      {
	// For the generators on the real line, we have
        //   \supp\phi_{j,k} = 2^{-j}[ell1+k,ell2+k]
	k1 = dyadic_modulo(RBASIS::primal_mask::begin()+lambda.k(), lambda.j());
	k2 = dyadic_modulo(RBASIS::primal_mask::end()  +lambda.k(), lambda.j());
      }
    else // wavelet
      {
	// For the wavelets on the real line, we have
        //   \supp\phi_{j,k} = 2^{-(j+1)}[ell1+1-ell2T+2*k,ell2+1-ell1T+2*k]
	k1 = dyadic_modulo(RBASIS::primal_mask::begin()+1-RBASIS::dual_mask::end()+2*lambda.k(), lambda.j()+1);
	k2 = dyadic_modulo(RBASIS::primal_mask::end()+1-RBASIS::dual_mask::begin()+2*lambda.k(), lambda.j()+1);
      }
  }

  template <class RBASIS>
  bool
  PeriodicBasis<RBASIS>::intersect_supports(const Index& lambda, const Index& mu)
  {
    // A note on the strategy:
    // Both supp(psi_lambda) and supp(psi_mu) are subintervals of the circle.
    // In order to decide whether they intersect or not, we compute the distance
    // between the respective midpoints.

    // compute dyadic coordinate of the midpoint of supp(psi_lambda)
    const int j_lambda = lambda.j()+1; // independent from lambda.e()!

    int k1_lambda, k2_lambda;
    if (lambda.e() == 0) {
      k1_lambda = 2*dyadic_modulo(RBASIS::primal_mask::begin()+lambda.k(), lambda.j());
      k2_lambda = 2*dyadic_modulo(RBASIS::primal_mask::end()  +lambda.k(), lambda.j());
    } else {
      k1_lambda = dyadic_modulo(RBASIS::primal_mask::begin()+1-RBASIS::dual_mask::end()+2*lambda.k(), j_lambda);
      k2_lambda = dyadic_modulo(RBASIS::primal_mask::end()+1-RBASIS::dual_mask::begin()+2*lambda.k(), j_lambda);
    }
    int length_lambda = k2_lambda-k1_lambda+(k2_lambda>k1_lambda ? 0 : 1<<j_lambda);
    int mid_lambda = dyadic_modulo((k2_lambda+k1_lambda+(k2_lambda>k1_lambda ? 0 : 1<<j_lambda))/2, j_lambda);

    // do the same for supp(psi_mu)
    const int j_mu = mu.j()+1;
    int k1_mu, k2_mu;
    if (mu.e() == 0) {
      k1_mu = 2*dyadic_modulo(RBASIS::primal_mask::begin()+mu.k(), mu.j());
      k2_mu = 2*dyadic_modulo(RBASIS::primal_mask::end()  +mu.k(), mu.j());
    } else {
      k1_mu = dyadic_modulo(RBASIS::primal_mask::begin()+1-RBASIS::dual_mask::end()+2*mu.k(), mu.j()+1);
      k2_mu = dyadic_modulo(RBASIS::primal_mask::end()+1-RBASIS::dual_mask::begin()+2*mu.k(), mu.j()+1);
    }
    int length_mu = k2_mu-k1_mu+(k2_mu>k1_mu ? 0 : 1<<j_mu);
    int mid_mu = dyadic_modulo((k2_mu+k1_mu+(k2_mu>k1_mu ? 0 : 1<<j_mu))/2, j_mu);
    
    const int j = std::max(j_lambda, j_mu);
    mid_lambda    <<= (j-j_lambda);
    length_lambda <<= (j-j_lambda);
    mid_mu        <<= (j-j_mu);
    length_mu     <<= (j-j_mu);

    const int dist_midpoints =
      mid_lambda > mid_mu
      ? std::min(mid_lambda-mid_mu, mid_mu+(1<<j)-mid_lambda)
      : std::min(mid_mu-mid_lambda, mid_lambda+(1<<j)-mid_mu);
    
    return (2*dist_midpoints < length_lambda+length_mu);
  }

  //
  //
  // point evaluation subroutines

  template <class RBASIS>
  SampledMapping<1>
  PeriodicBasis<RBASIS>::evaluate
  (const typename PeriodicBasis<RBASIS>::Index& lambda,
   const int resolution) const
  {
    Grid<1> grid(0, 1, 1<<resolution);
    Array1D<double> values((1<<resolution)+1);
    for (unsigned int i(0); i < values.size(); i++) {
      const double x = i*ldexp(1.0, -resolution);
      values[i] = evaluate(0, lambda, x);
    }     
    return SampledMapping<1>(grid, values);
  }

//   template <int d, int dT, SplineBasisFlavor flavor, int s0, int s1, int sT0, int sT1>
//   SampledMapping<1>
//   SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::evaluate
//   (const InfiniteVector<double, typename SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::Index>& coeffs,
//    const int resolution) const
//   {
//     Grid<1> grid(0, 1, 1<<resolution);
//     SampledMapping<1> result(grid); // zero
//     if (coeffs.size() > 0) {
//       // determine maximal level
//       int jmax(0);
//       typedef typename SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::Index Index;
//       for (typename InfiniteVector<double,Index>::const_iterator it(coeffs.begin()),
//   	     itend(coeffs.end()); it != itend; ++it)
//   	jmax = std::max(it.index().j()+it.index().e(), jmax);
      
//       // insert coefficients into a dense vector
//       Vector<double> wcoeffs(Deltasize(jmax));
//       for (typename InfiniteVector<double,Index>::const_iterator it(coeffs.begin()),
//  	     itend(coeffs.end()); it != itend; ++it) {
//  	// determine number of the wavelet
//  	typedef typename Vector<double>::size_type size_type;
//  	size_type number = 0;
//  	if (it.index().e() == 0) {
//  	  number = it.index().k()-DeltaLmin();
//  	} else {
//  	  number = Deltasize(it.index().j())+it.index().k()-Nablamin();
//  	}
//  	wcoeffs[number] = *it;
//       }
      
//       // switch to generator representation
//       Vector<double> gcoeffs(wcoeffs.size(), false);
//       if (jmax == j0())
//  	gcoeffs = wcoeffs;
//       else
//  	apply_Tj(jmax-1, wcoeffs, gcoeffs);
      
//       Array1D<double> values((1<<resolution)+1);
//       for (unsigned int i(0); i < values.size(); i++) {
//  	values[i] = 0;
//  	const double x = i*ldexp(1.0, -resolution);
//  	for (unsigned int k = 0; k < gcoeffs.size(); k++) {
// 	  values[i] += gcoeffs[k] * evaluate(0, Index(jmax, 0, DeltaLmin()+k), x);
//  	}
//       }
      
//       return SampledMapping<1>(grid, values);
//     }
    
//     return result;
//   }


  
  template <class RBASIS>
  inline
  double
  PeriodicBasis<RBASIS>::evaluate
  (const unsigned int derivative, const Index& lambda, const double x) const
  {
    return r_basis.evaluate(derivative, lambda,
			    x-floor(x-ldexp(1.0,-lambda.j())
				    *((lambda.e() == 0
				       ? RBASIS::primal_mask::begin()
				       : (RBASIS::primal_mask::begin()+1-RBASIS::dual_mask::end())/2)
				      +lambda.k())));
  }
  
  template <class RBASIS>
  inline
  void
  PeriodicBasis<RBASIS>::evaluate(const unsigned int derivative,
				  const Index& lambda,
				  const Array1D<double>& points,
				  Array1D<double>& values) const
  {
    values.resize(points.size());
    for (unsigned int i = 0; i < points.size(); i++)
      values[i] = evaluate(derivative, lambda, points[i]);
  }
    
  template <class RBASIS>
  inline
  void
  PeriodicBasis<RBASIS>::evaluate(const Index& lambda,
				  const Array1D<double>& points,
				  Array1D<double>& funcvalues,
				  Array1D<double>& dervalues) const
  {
    funcvalues.resize(points.size());
    dervalues.resize(points.size());
    evaluate(0, lambda, points, funcvalues);
    evaluate(1, lambda, points, dervalues );
  }

  template <class RBASIS>
  void
  PeriodicBasis<RBASIS>::expand(const Function<1>* f,
				const bool primal,
				const int jmax,
				InfiniteVector<double, Index>& coeffs) const
  {
    for (Index lambda = first_generator(j0());;++lambda)
      {
	coeffs.set_coefficient(lambda, integrate(f, lambda));
	if (lambda == last_wavelet(jmax))
	  break;
      }

    if (!primal) {
      // setup active index set
      std::set<Index> Lambda;
      for (Index lambda = first_generator(j0());; ++lambda) {
 	Lambda.insert(lambda);
	if (lambda == last_wavelet(jmax)) break;
      }
      
      // setup Gramian A_Lambda
      SparseMatrix<double> A_Lambda(Lambda.size());
      typedef typename SparseMatrix<double>::size_type size_type;     
      size_type row = 0;
      for (typename std::set<Index>::const_iterator it1(Lambda.begin()), itend(Lambda.end());
	   it1 != itend; ++it1, ++row)
	{
	  std::list<size_type> indices;
	  std::list<double> entries;
	  
	  size_type column = 0;
	  for (typename std::set<Index>::const_iterator it2(Lambda.begin());
	       it2 != itend; ++it2, ++column)
	    {
	      double entry = integrate(*it2, *it1);
	      
	      if (entry != 0) {
		indices.push_back(column);
		entries.push_back(entry);
	      }
	    }
	  A_Lambda.set_row(row, indices, entries);
	} 
      
      // solve A_Lambda*x = b
      Vector<double> b(Lambda.size());
      row = 0;
      for (typename std::set<Index>::const_iterator it(Lambda.begin()), itend(Lambda.end());
	   it != itend; ++it, ++row)
	b[row] = coeffs.get_coefficient(*it);
      
      Vector<double> x(b);
      unsigned int iterations;
      CG(A_Lambda, b, x, 1e-15, 500, iterations);
      
      coeffs.clear();
      row = 0;
      for (typename std::set<Index>::const_iterator it(Lambda.begin()), itend(Lambda.end());
	   it != itend; ++it, ++row)
	coeffs.set_coefficient(*it, x[row]);
    }
  }

  template <class RBASIS>
  void
  PeriodicBasis<RBASIS>::expand(const Function<1>* f,
				const bool primal,
				const int jmax,
				Vector<double>& coeffs) const
  {
    assert(primal);

    coeffs.resize(Deltasize(jmax+1));

    int id = 0;
    for (Index lambda = first_generator(j0());;++lambda, ++id)
      {
	coeffs[id] = integrate(f, lambda);
	if (lambda == last_wavelet(jmax))
	  break;
      }   
  }
  
  template <class RBASIS>
  double
  PeriodicBasis<RBASIS>::integrate(const Function<1>* f,
				   const Index& lambda) const
  {
    double r = 0;

    // first we compute the support of psi_lambda
    const int j = lambda.j()+lambda.e();
    int k1, k2;
    support(lambda, k1, k2); // note: k2 may be less or equal to k1 in the case of overlap
    const int length = (k2 > k1 ? k2-k1 : k2+(1<<j)-k1); // number of subintervals
    
    // setup Gauss points and weights for a composite quadrature formula:
    const unsigned int N_Gauss = 5;
    const double h = 1.0/(1<<j);

    Array1D<double> gauss_points (N_Gauss*length);
    int k = k1;
    for (int patch = 0; patch < length; patch++, k = dyadic_modulo(++k,j)) // work on 2^{-j}[k,k+1]
      for (unsigned int n = 0; n < N_Gauss; n++)
 	gauss_points[patch*N_Gauss+n] = h*(2*k+1+GaussPoints[N_Gauss-1][n])/2;
    
    // add all integral shares
    for (unsigned int n = 0; n < N_Gauss; n++)
      {
 	const double gauss_weight = GaussWeights[N_Gauss-1][n] * h;
 	for (int patch = k1; patch < k1+length; patch++)
 	  {
 	    const double t = gauss_points[(patch-k1)*N_Gauss+n];
	    
 	    const double ft = f->value(MathTL::Point<1>(t));
 	    if (ft != 0)
 	      r += ft
 		* evaluate(0, lambda, t)
 		* gauss_weight;
 	  }
      }
    
    return r;
  }
  
  template <class RBASIS>
  double
  PeriodicBasis<RBASIS>::integrate(const Index& lambda,
				   const Index& mu) const
  {
    double r = 0;
    
    if (intersect_supports(lambda, mu))
      {
	// first we determine the support over which to integrate
	int j, k1, k2, length;
	if (lambda.j()+lambda.e() >= mu.j()+mu.e()) {
	  j = lambda.j()+lambda.e();
	  support(lambda, k1, k2); // note: k2 may be less or equal to k1 in the case of overlap
	} else {
	  j = mu.j()+mu.e();
	  support(mu, k1, k2); // note: k2 may be less or equal to k1 in the case of overlap
	}
	length = (k2 > k1 ? k2-k1 : k2+(1<<j)-k1); // number of subintervals
	
	// setup Gauss points and weights for a composite quadrature formula:
 	const unsigned int N_Gauss = primal_polynomial_degree();
	const double h = 1.0/(1<<j);
	
	Array1D<double> gauss_points (N_Gauss*(length)), func1values, func2values;
	int k = k1;
	for (int patch = 0; patch < length; patch++, k = dyadic_modulo(++k,j)) // work on 2^{-j}[k,k+1]
	  for (unsigned int n = 0; n < N_Gauss; n++)
	    gauss_points[patch*N_Gauss+n] = h*(2*k+1+GaussPoints[N_Gauss-1][n])/2;

  	// - compute point values of the integrands
   	evaluate(0, lambda, gauss_points, func1values);
  	evaluate(0, mu, gauss_points, func2values);
	
  	// - add all integral shares
	for (int patch = k1, id = 0; patch < k1+length; patch++)
	  for (unsigned int n = 0; n < N_Gauss; n++, id++) {
	    r += func1values[id] * func2values[id] * GaussWeights[N_Gauss-1][n] * h;
	  }

	return r;
      }

    return 0;
  }

  template <class RBASIS>
  void
  PeriodicBasis<RBASIS>::decompose(const InfiniteVector<double, Index>& c,
				   const int j0,
				   InfiniteVector<double, Index>& d) const {
    InfiniteVector<double, Index> help;
    for (typename InfiniteVector<double, Index>::const_iterator it(c.begin()), itend(c.end());
	 it != itend; ++it) {
      decompose_1(it.index(), j0, help); // calls help.clear() first
      d.add(*it, help);
    }
  }
  
  template <class RBASIS>
  void
  PeriodicBasis<RBASIS>::decompose_t(const InfiniteVector<double, Index>& c,
				     const int j0,
				     InfiniteVector<double, Index>& d) const {
    InfiniteVector<double, Index> help;
    for (typename InfiniteVector<double, Index>::const_iterator it(c.begin()), itend(c.end());
	 it != itend; ++it) {
      decompose_t_1(it.index(), j0, help); // calls help.clear() first
      d.add(*it, help);
    }
  }

  template <class RBASIS>
  void
  PeriodicBasis<RBASIS>::reconstruct(const InfiniteVector<double, Index>& c,
				     const int j,
				     InfiniteVector<double, Index>& d) const {
    for (typename InfiniteVector<double, Index>::const_iterator it(c.begin()), itend(c.end());
	 it != itend; ++it) {
      InfiniteVector<double, Index> help;
      reconstruct_1(it.index(), j, help);
      d.add(*it, help);
    }
  }
  
  template <class RBASIS>
  void
  PeriodicBasis<RBASIS>::reconstruct_t(const InfiniteVector<double, Index>& c,
				       const int j,
				       InfiniteVector<double, Index>& d) const {
    for (typename InfiniteVector<double, Index>::const_iterator it(c.begin()), itend(c.end());
	 it != itend; ++it) {
      InfiniteVector<double, Index> help;
      reconstruct_t_1(it.index(), j, help);
      d.add(*it, help);
    }
  }
  
  template <class RBASIS>
  void
  PeriodicBasis<RBASIS>::decompose_1(const Index& lambda,
				     const int j0,
				     InfiniteVector<double, Index>& c) const {
    assert(lambda.j() >= j0);
    c.clear();
    if (lambda.e() == 1) {
      // the true wavelet coefficients don't have to be modified
      c.set_coefficient(lambda, 1.0);
    } else {
      // a generator on a (possibly) fine level
      if (lambda.j() == j0) {
	// generators from the coarsest level can be copied
	c.set_coefficient(lambda, 1.0);
      }	else {
	// j>j0, perform multiscale decomposition
	
	const int aTbegin = RBASIS::dual_mask::begin();
	const int aTend   = RBASIS::dual_mask::end();
	const int bTbegin = 1-RBASIS::primal_mask::end();
	const int bTend   = 1-RBASIS::primal_mask::begin();
	
	// compute d_{j-1}
	for (int n = 0; n < 1<<(lambda.j()-1); n++) {
	  double cn = 0;
	  for (int m = (int) ceil(ldexp(1.0, 1-lambda.j()) * ((lambda.k() - bTend) / 2.0 - n));
	       m <= (int) floor(ldexp(1.0, 1-lambda.j()) * ((lambda.k() - bTbegin) / 2.0 - n));
	       m++)
	    cn += r_basis.bT(lambda.k()-2*((1<<(lambda.j()-1))*m+n));
	  if (cn != 0)
	    c[Index(lambda.j()-1, 1, n)] = M_SQRT1_2 * cn;
	}
	
	// compute c_{j_0} via recursion
	for (int n = 0; n < 1<<(lambda.j()-1); n++) {
	  double cn = 0;
	  for (int m = (int) ceil(ldexp(1.0, 1-lambda.j()) * ((lambda.k() - aTend) / 2.0 - n));
	       m <= (int) floor(ldexp(1.0, 1-lambda.j()) * ((lambda.k() - aTbegin) / 2.0 - n));
	       m++)
	    cn += r_basis.aT(lambda.k()-2*((1<<(lambda.j()-1))*m+n));
	  InfiniteVector<double,Index> d;
	  decompose_1(Index(lambda.j()-1, 0, n), j0, d);
	  c.add(M_SQRT1_2 * cn, d);
	}
      }
    }
  }
  
  template <class RBASIS>
  void
  PeriodicBasis<RBASIS>::decompose_t_1(const Index& lambda,
				       const int j0,
				       InfiniteVector<double, Index>& c) const {
    assert(lambda.j() >= j0);
    c.clear();
    if (lambda.e() == 1) {
      // the true wavelet coefficients don't have to be modified
      c.set_coefficient(lambda, 1.0);
    } else {
      // a generator on a (possibly) fine level
      if (lambda.j() == j0) {
	// generators from the coarsest level can be copied
	c.set_coefficient(lambda, 1.0);
      }	else {
	// j>j0, perform multiscale decomposition
	
	const int abegin = RBASIS::primal_mask::begin();
	const int aend   = RBASIS::primal_mask::end();
	const int bbegin = 1-RBASIS::dual_mask::end();
	const int bend   = 1-RBASIS::dual_mask::begin();
	
	// compute d_{j-1}
	for (int n = 0; n < 1<<(lambda.j()-1); n++) {
	  double cn = 0;
	  for (int m = (int) ceil(ldexp(1.0, 1-lambda.j()) * ((lambda.k() - bend) / 2.0 - n));
	       m <= (int) floor(ldexp(1.0, 1-lambda.j()) * ((lambda.k() - bbegin) / 2.0 - n));
	       m++)
	    cn += r_basis.b(lambda.k()-2*((1<<(lambda.j()-1))*m+n));
	  if (cn != 0)
	    c[Index(lambda.j()-1, 1, n)] = M_SQRT1_2 * cn;
	}
	
	// compute c_{j_0} via recursion
	for (int n = 0; n < 1<<(lambda.j()-1); n++) {
	  double cn = 0;
	  for (int m = (int) ceil(ldexp(1.0, 1-lambda.j()) * ((lambda.k() - aend) / 2.0 - n));
	       m <= (int) floor(ldexp(1.0, 1-lambda.j()) * ((lambda.k() - abegin) / 2.0 - n));
	       m++)
	    cn += r_basis.a(lambda.k()-2*((1<<(lambda.j()-1))*m+n));
	  InfiniteVector<double,Index> d;
	  decompose_t_1(Index(lambda.j()-1, 0, n), j0, d);
	  c.add(M_SQRT1_2 * cn, d);
	}
      }
    }
  }

  template <class RBASIS>
  void
  PeriodicBasis<RBASIS>::reconstruct_1(const Index& lambda,
				       const int j,
				       InfiniteVector<double, Index>& c) const {
    if (lambda.j() >= j) {
      // then we can just copy \psi_\lambda
      c.set_coefficient(lambda, 1.0);
    } else {
      // reconstruct by recursion
      
      const int abegin = RBASIS::primal_mask::begin();
      const int aend   = RBASIS::primal_mask::end();
      const int bbegin = 1-RBASIS::dual_mask::end();
      const int bend   = 1-RBASIS::dual_mask::begin();
      
      if (lambda.e() == 0) {
	for (int n = 0; n < 1<<(lambda.j()+1); n++) {
	  double cn = 0;
	  for (int m = (int) ceil(ldexp(1.0, -lambda.j()-1) * (abegin+2*lambda.k()-n));
	       m <= (int) floor(ldexp(1.0, -lambda.j()-1) * (aend+2*lambda.k()-n));
	       m++)
	    cn += r_basis.a((1<<(lambda.j()+1))*m+n-2*lambda.k());
	  if (cn != 0) {
	    InfiniteVector<double,Index> d;
	    reconstruct_1(Index(lambda.j()+1, 0, n), j, d);
	    c.add(M_SQRT1_2 * cn, d);
	  }
	}
      } else {
	for (int n = 0; n < 1<<(lambda.j()+1); n++) {
	  double cn = 0;
	  for (int m = (int) ceil(ldexp(1.0, -lambda.j()-1) * (bbegin+2*lambda.k()-n));
	       m <= (int) floor(ldexp(1.0, -lambda.j()-1) * (bend+2*lambda.k()-n));
	       m++)
	    cn += r_basis.b((1<<(lambda.j()+1))*m+n-2*lambda.k());
	  if (cn != 0) {
	    InfiniteVector<double,Index> d;
	    reconstruct_1(Index(lambda.j()+1, 0, n), j, d);
	    c.add(M_SQRT1_2 * cn, d);
	  }
	}
      }
    }
  }
  
  template <class RBASIS>
  void
  PeriodicBasis<RBASIS>::reconstruct_t_1(const Index& lambda,
					 const int j,
					 InfiniteVector<double, Index>& c) const {
    if (lambda.j() >= j) {
      // then we can just copy \psi_\lambda
      c.set_coefficient(lambda, 1.0);
    } else {
      // reconstruct by recursion
   
      const int aTbegin = RBASIS::dual_mask::begin();
      const int aTend   = RBASIS::dual_mask::end();
      const int bTbegin = 1-RBASIS::primal_mask::end();
      const int bTend   = 1-RBASIS::primal_mask::begin();
   
      if (lambda.e() == 0) {
	for (int n = 0; n < 1<<(lambda.j()+1); n++) {
	  double cn = 0;
	  for (int m = (int) ceil(ldexp(1.0, -lambda.j()-1) * (aTbegin+2*lambda.k()-n));
	       m <= (int) floor(ldexp(1.0, -lambda.j()-1) * (aTend+2*lambda.k()-n));
	       m++)
	    cn += r_basis.aT((1<<(lambda.j()+1))*m+n-2*lambda.k());
	  if (cn != 0) {
	    InfiniteVector<double,Index> d;
	    reconstruct_t_1(Index(lambda.j()+1, 0, n), j, d);
	    c.add(M_SQRT1_2 * cn, d);
	  }
	}
      } else {
	for (int n = 0; n < 1<<(lambda.j()+1); n++) {
	  double cn = 0;
	  for (int m = (int) ceil(ldexp(1.0, -lambda.j()-1) * (bTbegin+2*lambda.k()-n));
	       m <= (int) floor(ldexp(1.0, -lambda.j()-1) * (bTend+2*lambda.k()-n));
	       m++)
	    cn += r_basis.bT((1<<(lambda.j()+1))*m+n-2*lambda.k());
	  if (cn != 0) {
	    InfiniteVector<double,Index> d;
	    reconstruct_t_1(Index(lambda.j()+1, 0, n), j, d);
	    c.add(M_SQRT1_2 * cn, d);
	  }
	}
      }
    }
  }

}
