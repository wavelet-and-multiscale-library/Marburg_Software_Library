// implementation of PeriodicBasis methods

#include <cmath>

namespace WaveletTL
{
  template <class RBASIS>
  void
  PeriodicBasis<RBASIS>::support(const Index& lambda, int& k1, int& k2) const
  {
    if (lambda.e() == 0) // generator
      {
	// For the generators on the real line, we have
        //   \supp\phi_{j,k} = 2^{-j}[ell1+k,ell2+k]
	// To obtain the support of the periodized generators, we assume that 0 <= ell_1+2^j
	const int help = 1<<lambda.j();
	k1 = (RBASIS::primal_mask::begin()+lambda.k()+help) % help;
	k2 = (RBASIS::primal_mask::end()  +lambda.k()+help) % help;
      }
    else // wavelet
      {
	// For the wavelets on the real line, we have
        //   \supp\phi_{j,k} = 2^{-(j+1)}[ell1+1-ell2T+2*k,ell2+1-ell1T+2*k] =: 2^{-(j+1)}[c1,c2]
	// To obtain the support of the periodized generators, we assume that 0 <= c1+2^(j+1)
	const int help = 1<<(lambda.j()+1);
	k1 = (RBASIS::primal_mask::begin()+1-RBASIS::dual_mask::end()+2*lambda.k()+help) % help;
	k2 = (RBASIS::primal_mask::end()+1-RBASIS::dual_mask::begin()+2*lambda.k()+help) % help;
      }
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
