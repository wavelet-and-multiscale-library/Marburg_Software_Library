// implementation of PeriodicBasis methods

#include <cmath>
#include <utils/tiny_tools.h>
#include <numerics/gauss_data.h>

namespace WaveletTL
{
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
  inline
  double
  PeriodicBasis<RBASIS>::evaluate
  (const unsigned int derivative, const Index& lambda, const double x) const
  {
    return r_basis.evaluate(derivative, lambda,
			    x-(int)floor(x-ldexp(1.0,-lambda.j())
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
    assert(!primal);
    
    for (Index lambda = first_generator(j0());;++lambda)
      {
	coeffs.set_coefficient(lambda, integrate(f, lambda));
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
    const double h = ldexp(1.0, -j);

    Array1D<double> gauss_points (N_Gauss*(length));
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
	    
 	    const double ft = f->value(Point<1>(t));
 	    if (ft != 0)
 	      r += ft
 		* evaluate(0, lambda, t)
 		* gauss_weight;
 	  }
      }
    
    return r;
  }
  

  template <class RBASIS>
  void
  PeriodicBasis<RBASIS>::decompose(const InfiniteVector<double, Index>& c,
				   const int j0,
				   InfiniteVector<double, Index>& d) {
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
				     InfiniteVector<double, Index>& d) {
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
				     InfiniteVector<double, Index>& d) {
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
				       InfiniteVector<double, Index>& d) {
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
				     InfiniteVector<double, Index>& c) {
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
	
	static const int aTbegin = RBASIS::dual_mask::begin();
	static const int aTend   = RBASIS::dual_mask::end();
	static const int bTbegin = 1-RBASIS::primal_mask::end();
	static const int bTend   = 1-RBASIS::primal_mask::begin();
	
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
				       InfiniteVector<double, Index>& c) {
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
				       InfiniteVector<double, Index>& c) {
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
					 InfiniteVector<double, Index>& c) {
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
