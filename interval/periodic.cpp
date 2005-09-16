// implementation of PeriodicBasis methods

#include <cmath>

namespace WaveletTL
{
  template <class RBASIS>
  PeriodicBasis<RBASIS>::PeriodicBasis() {
    const int alength  = rbasis.a().rbegin().index()[0] - rbasis.a().begin().index()[0];
    const int aTlength = rbasis.aT().rbegin().index()[0] - rbasis.aT().begin().index()[0];
    j0_ = (int) ceil(log((double) std::max(alength, aTlength))/M_LN2);
  }

  template <class RBASIS>
  void
  PeriodicBasis<RBASIS>::decompose(const InfiniteVector<double, Index>& c,
				   const int j0,
				   InfiniteVector<double, Index>& d) const {
    for (typename InfiniteVector<double, Index>::const_iterator it(c.begin()), itend(c.end());
	 it != itend; ++it) {
      InfiniteVector<double, Index> help;
      decompose_1(it.index(), j0, help);
      d += *it * help;
    }
  }
  
  template <class RBASIS>
  void
  PeriodicBasis<RBASIS>::decompose_t(const InfiniteVector<double, Index>& c,
				     const int j0,
				     InfiniteVector<double, Index>& d) const {
    for (typename InfiniteVector<double, Index>::const_iterator it(c.begin()), itend(c.end());
	 it != itend; ++it) {
      InfiniteVector<double, Index> help;
      decompose_t_1(it.index(), j0, help);
      d += *it * help;
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
      d += *it * help;
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
      d += *it * help;
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
	
	const int aTbegin = rbasis.aT().begin().index()[0];
	const int aTend   = rbasis.aT().rbegin().index()[0];
	const int bTbegin = rbasis.bT().begin().index()[0];
	const int bTend   = rbasis.bT().rbegin().index()[0];
	
	// compute d_{j-1}
	for (int n = 0; n < 1<<(lambda.j()-1); n++) {
	  double cn = 0;
	  for (int m = (int) ceil(ldexp(1.0, 1-lambda.j()) * ((lambda.k() - bTend) / 2.0 - n));
	       m <= (int) floor(ldexp(1.0, 1-lambda.j()) * ((lambda.k() - bTbegin) / 2.0 - n));
	       m++)
	    cn += rbasis.bT().get_coefficient(MultiIndex<int, 1>(lambda.k()-2*((1<<(lambda.j()-1))*m+n)));
	  if (cn != 0)
	    c[Index(lambda.j()-1, 1, n, this)] = M_SQRT1_2 * cn;
	}
	
	// compute c_{j_0} via recursion
	for (int n = 0; n < 1<<(lambda.j()-1); n++) {
	  double cn = 0;
	  for (int m = (int) ceil(ldexp(1.0, 1-lambda.j()) * ((lambda.k() - aTend) / 2.0 - n));
	       m <= (int) floor(ldexp(1.0, 1-lambda.j()) * ((lambda.k() - aTbegin) / 2.0 - n));
	       m++)
	    cn += rbasis.aT().get_coefficient(MultiIndex<int, 1>(lambda.k()-2*((1<<(lambda.j()-1))*m+n)));
	  InfiniteVector<double,Index> d;
	  decompose_1(Index(lambda.j()-1, 0, n, this), j0, d);
	  c += M_SQRT1_2 * cn * d;
	}
      }
    }
  }
  
  template <class RBASIS>
  void
  PeriodicBasis<RBASIS>::decompose_t_1(const Index& lambda,
				       const int j0,
				       InfiniteVector<double, Index>& c) const {
//     assert(lambda.j() >= j0);
//     c.clear();
//     if (lambda.e() == 1)
//       {
// 	// the true wavelet coefficients don't have to be modified
// 	c[lambda] = 1.0;
//       }
//     else
//       {
// 	// a generator on a (possibly) fine level
// 	if (lambda.j() == j0)
// 	  {
// 	    // generators from the coarsest level can be copied
// 	    c[lambda] = 1.0;
// 	  }
// 	else
// 	  {
// 	    // j>j0, perform multiscale decomposition
	    
// 	    const int abegin = a().begin().index()[0];
// 	    const int aend   = a().rbegin().index()[0];
// 	    const int bbegin = b().begin().index()[0];
// 	    const int bend   = b().rbegin().index()[0];
	    
// 	    // compute d_{j-1}
//  	    for (int l((int)ceil((lambda.k()-bend)/2.0));
// 		 l <= (int)floor((lambda.k()-bbegin)/2.0); l++)
//  	      c[Index(lambda.j()-1, 1, l)] =
// 		M_SQRT1_2 * b().get_coefficient(MultiIndex<int, 1>(lambda.k()-2*l));
	    
// 	    // compute c_{j_0} via recursion
// 	    for (int l((int)ceil((lambda.k()-aend)/2.0));
// 		 l <= (int)floor((lambda.k()-abegin)/2.0); l++)
// 	      {
// 		InfiniteVector<double, Index> d;
// 		decompose_t_1(Index(lambda.j()-1, 0, l), j0, d);
// 		c += M_SQRT1_2 * a().get_coefficient(MultiIndex<int, 1>(lambda.k()-2*l)) * d;
// 	      }
// 	  }
//       }
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
      
      const int abegin = rbasis.a().begin().index()[0];
      const int aend   = rbasis.a().rbegin().index()[0];
      const int bbegin = rbasis.b().begin().index()[0];
      const int bend   = rbasis.b().rbegin().index()[0];
      
      if (lambda.e() == 0) {
	for (int n = 0; n < 1<<(lambda.j()+1); n++) {
	  double cn = 0;
	  for (int m = (int) ceil(ldexp(1.0, -lambda.j()-1) * (abegin+2*lambda.k()-n));
	       m <= (int) floor(ldexp(1.0, -lambda.j()-1) * (aend+2*lambda.k()-n));
	       m++)
	    cn += rbasis.a().get_coefficient(MultiIndex<int, 1>((1<<(lambda.j()+1))*m+n-2*lambda.k()));
	  if (cn != 0) {
	    InfiniteVector<double,Index> d;
	    reconstruct_1(Index(lambda.j()+1, 0, n, this), j, d);
	    c += M_SQRT1_2 * cn * d;
	  }
	}
      } else {
	for (int n = 0; n < 1<<(lambda.j()+1); n++) {
	  double cn = 0;
	  for (int m = (int) ceil(ldexp(1.0, -lambda.j()-1) * (bbegin+2*lambda.k()-n));
	       m <= (int) floor(ldexp(1.0, -lambda.j()-1) * (bend+2*lambda.k()-n));
	       m++)
	    cn += rbasis.b().get_coefficient(MultiIndex<int, 1>((1<<(lambda.j()+1))*m+n-2*lambda.k()));
	  if (cn != 0) {
	    InfiniteVector<double,Index> d;
	    reconstruct_1(Index(lambda.j()+1, 0, n, this), j, d);
	    c += M_SQRT1_2 * cn * d;
	  }
	}
      }
    }
  }
  
  template <class RBASIS>
  void
  PeriodicBasis<RBASIS>::reconstruct_t_1(const Index& lambda,
					 const int j,
					 InfiniteVector<double, Index>& c) const
  {
//     if (lambda.j() >= j)
//       {
// 	// then we can just copy \psi_\lambda
// 	c[lambda] += 1.0;
//       }
//     else
//       {
// 	// reconstruct by recursion

// 	const int aTbegin = aT().begin().index()[0];
// 	const int aTend   = aT().rbegin().index()[0];
// 	const int bTbegin = bT().begin().index()[0];
// 	const int bTend   = bT().rbegin().index()[0];
	    
// 	if (lambda.e() == 0)
// 	  {
//  	    for (int l(2*lambda.k()+aTbegin); l <= 2*lambda.k()+aTend; l++)
//  	      {
//  		InfiniteVector<double, Index> d;
//  		reconstruct_t_1(Index(lambda.j()+1, 0, l), j, d);
//  		c += M_SQRT1_2 * aT().get_coefficient(MultiIndex<int, 1>(l-2*lambda.k())) * d;
//  	      }
//  	  }
//  	else
//  	  {
//  	    for (int l(2*lambda.k()+bTbegin); l <= 2*lambda.k()+bTend; l++)
//  	      {
//  		InfiniteVector<double, Index> d;
//  		reconstruct_t_1(Index(lambda.j()+1, 0, l), j, d);
//  		c += M_SQRT1_2 * bT().get_coefficient(MultiIndex<int, 1>(l-2*lambda.k())) * d;
//  	      }
// 	  }
//       }
  }

}
