// implementation of r_basis.h

#include <cmath>
#include <cassert>

namespace WaveletTL
{
  template <class PRIMALMASK, class DUALMASK>
  RBasis<PRIMALMASK, DUALMASK>::RBasis()
  {
  }

  template <class PRIMALMASK, class DUALMASK>
  void
  RBasis<PRIMALMASK, DUALMASK>::decompose(const InfiniteVector<double, Index>& c,
					  const int j0,
					  InfiniteVector<double, Index>& d) const
  {
    for (InfiniteVector<double, Index>::const_iterator it(c.begin()), itend(c.end());
	 it != itend; ++it)
      {
 	InfiniteVector<double, Index> help;
 	decompose_1(it.index(), j0, help);
 	d.add(*it, help);
      }
  }
  
  template <class PRIMALMASK, class DUALMASK>
  void
  RBasis<PRIMALMASK, DUALMASK>::decompose_t(const InfiniteVector<double, Index>& c,
					    const int j0,
					    InfiniteVector<double, Index>& d) const
  {
    for (InfiniteVector<double, Index>::const_iterator it(c.begin()), itend(c.end());
	 it != itend; ++it)
      {
	InfiniteVector<double, Index> help;
	decompose_t_1(it.index(), j0, help);
	d.add(*it, help);
      }
  }

  template <class PRIMALMASK, class DUALMASK>
  void
  RBasis<PRIMALMASK, DUALMASK>::reconstruct(const InfiniteVector<double, Index>& c,
					    const int j,
					    InfiniteVector<double, Index>& d) const
  {
    for (InfiniteVector<double, Index>::const_iterator it(c.begin()), itend(c.end());
	 it != itend; ++it)
      {
	InfiniteVector<double, Index> help;
	reconstruct_1(it.index(), j, help);
	d.add(*it, help);
      }
  }

  template <class PRIMALMASK, class DUALMASK>
  void
  RBasis<PRIMALMASK, DUALMASK>::reconstruct_t(const InfiniteVector<double, Index>& c,
					      const int j,
					      InfiniteVector<double, Index>& d) const
  {
    for (InfiniteVector<double, Index>::const_iterator it(c.begin()), itend(c.end());
	 it != itend; ++it)
      {
	InfiniteVector<double, Index> help;
	reconstruct_t_1(it.index(), j, help);
	d.add(*it, help);
      }
  }

  template <class PRIMALMASK, class DUALMASK>
  void
  RBasis<PRIMALMASK, DUALMASK>::decompose_1(const Index& lambda,
					    const int j0,
					    InfiniteVector<double, Index>& c) const
  {
    assert(lambda.j() >= j0);
    c.clear();
    if (lambda.e() == 1)
      {
	// the true wavelet coefficients don't have to be modified
	c[lambda] = 1.0;
      }
    else
      {
	// a generator on a (possibly) fine level
	if (lambda.j() == j0)
	  {
	    // generators from the coarsest level can be copied
	    c[lambda] = 1.0;
	  }
	else
	  {
	    // j>j0, perform multiscale decomposition
	    
	    const int aTbegin = aT().abegin();
	    const int aTend   = aT().aend();
	    const int bTbegin = 1-a().aend();
	    const int bTend   = 1-a().abegin();
	    
	    // compute d_{j-1}
  	    for (int l((int)ceil((lambda.k()-bTend)/2.0));
 		 l <= (int)floor((lambda.k()-bTbegin)/2.0); l++)
  	      c[Index(lambda.j()-1, 1, l)] =
		M_SQRT1_2 * bT(lambda.k()-2*l);
	    
 	    // compute c_{j_0} via recursion
	    if (lambda.j() == j0-1) {
	      for (int l((int)ceil((lambda.k()-aTend)/2.0));
		   l <= (int)floor((lambda.k()-aTbegin)/2.0); l++) {
		c.add_coefficient(Index(lambda.j()-1, 0, l),
				  M_SQRT1_2 * aT(lambda.k()-2*l));
	      }
	    } else {
	      for (int l((int)ceil((lambda.k()-aTend)/2.0));
		   l <= (int)floor((lambda.k()-aTbegin)/2.0); l++) {
		InfiniteVector<double, Index> d;
		decompose_1(Index(lambda.j()-1, 0, l), j0, d);
		c.add(M_SQRT1_2 * aT(lambda.k()-2*l), d);
	      }
	    }
	  }
      }
  }

  template <class PRIMALMASK, class DUALMASK>
  void
  RBasis<PRIMALMASK, DUALMASK>::decompose_t_1(const Index& lambda,
					      const int j0,
					      InfiniteVector<double, Index>& c) const
  {
    assert(lambda.j() >= j0);
    c.clear();
    if (lambda.e() == 1)
      {
	// the true wavelet coefficients don't have to be modified
	c[lambda] = 1.0;
      }
    else
      {
	// a generator on a (possibly) fine level
	if (lambda.j() == j0)
	  {
	    // generators from the coarsest level can be copied
	    c[lambda] = 1.0;
	  }
	else
	  {
	    // j>j0, perform multiscale decomposition
	    
	    const int abegin = a().abegin();
	    const int aend   = a().aend();
	    const int bbegin = 1-aT().aend();
	    const int bend   = 1-aT().abegin();
	    
	    // compute d_{j-1}
 	    for (int l((int)ceil((lambda.k()-bend)/2.0));
		 l <= (int)floor((lambda.k()-bbegin)/2.0); l++)
 	      c[Index(lambda.j()-1, 1, l)] =
		M_SQRT1_2 * b(lambda.k()-2*l);
	    
	    // compute c_{j_0} via recursion
	    if (lambda.j() == j0-1) {
	      for (int l((int)ceil((lambda.k()-aend)/2.0));
		   l <= (int)floor((lambda.k()-abegin)/2.0); l++) {
		c.add_coefficient(Index(lambda.j()-1, 0, l),
				  M_SQRT1_2 * a(lambda.k()-2*l));
	      }
	    } else {
	      for (int l((int)ceil((lambda.k()-aend)/2.0));
		   l <= (int)floor((lambda.k()-abegin)/2.0); l++) {
		InfiniteVector<double, Index> d;
		decompose_t_1(Index(lambda.j()-1, 0, l), j0, d);
		c.add(M_SQRT1_2 * a(lambda.k()-2*l), d);
	      }
	    }
	  }
      }
  }

  template <class PRIMALMASK, class DUALMASK>
  void
  RBasis<PRIMALMASK, DUALMASK>::reconstruct_1(const Index& lambda,
					      const int j,
					      InfiniteVector<double, Index>& c) const
  {
    if (lambda.j() >= j)
      {
	// then we can just copy \psi_\lambda
	c[lambda] += 1.0;
      }
    else
      {
	// reconstruct by recursion

	const int abegin = a().abegin();
	const int aend   = a().aend();
	const int bbegin = 1-aT().aend();
	const int bend   = 1-aT().abegin();

	if (lambda.e() == 0)
	  {
	    if (lambda.j()+1>=j) {
	      for (int l(2*lambda.k()+abegin); l <= 2*lambda.k()+aend; l++)
		c[Index(lambda.j()+1, 0, l)]
		  += M_SQRT1_2 * a(l-2*lambda.k());
	    } else {
	      for (int l(2*lambda.k()+abegin); l <= 2*lambda.k()+aend; l++) {
		InfiniteVector<double, Index> d;
		reconstruct_1(Index(lambda.j()+1, 0, l), j, d);
		c.add(M_SQRT1_2 * a(l-2*lambda.k()), d);
	      }
	    }
 	  }
 	else
 	  {
	    if (lambda.j()+1>=j) {
	      for (int l(2*lambda.k()+bbegin); l <= 2*lambda.k()+bend; l++)
		c[Index(lambda.j()+1, 0, l)]
		  += M_SQRT1_2 * b(l-2*lambda.k());
	    } else {
	      for (int l(2*lambda.k()+bbegin); l <= 2*lambda.k()+bend; l++) {
		InfiniteVector<double, Index> d;
		reconstruct_1(Index(lambda.j()+1, 0, l), j, d);
		c.add(M_SQRT1_2 * b(l-2*lambda.k()), d);
 	      }
	    }
	  }
      }
  }

  template <class PRIMALMASK, class DUALMASK>
  void
  RBasis<PRIMALMASK, DUALMASK>::reconstruct_t_1(const Index& lambda,
						const int j,
						InfiniteVector<double, Index>& c) const
  {
    if (lambda.j() >= j)
      {
	// then we can just copy \psi_\lambda
	c[lambda] += 1.0;
      }
    else
      {
	// reconstruct by recursion

	const int aTbegin = aT().abegin();
	const int aTend   = aT().aend();
	const int bTbegin = 1-a().aend();
	const int bTend   = 1-a().abegin();
	    
	if (lambda.e() == 0)
	  {
	    if (lambda.j()+1>=j) {
	      for (int l(2*lambda.k()+aTbegin); l <= 2*lambda.k()+aTend; l++) {
		c[Index(lambda.j()+1, 0, l)]
		  += M_SQRT1_2 * aT(l-2*lambda.k());
	      }
	    } else {
	      for (int l(2*lambda.k()+aTbegin); l <= 2*lambda.k()+aTend; l++) {
 		InfiniteVector<double, Index> d;
 		reconstruct_t_1(Index(lambda.j()+1, 0, l), j, d);
 		c.add(M_SQRT1_2 * aT(l-2*lambda.k()), d);
 	      }
	    }
 	  }
 	else
 	  {
	    if (lambda.j()+1>=j) {
	      for (int l(2*lambda.k()+bTbegin); l <= 2*lambda.k()+bTend; l++) {
		c[Index(lambda.j()+1, 0, l)]
		  += M_SQRT1_2 * bT(l-2*lambda.k());
	      }
	    } else {
	      for (int l(2*lambda.k()+bTbegin); l <= 2*lambda.k()+bTend; l++) {
 		InfiniteVector<double, Index> d;
 		reconstruct_t_1(Index(lambda.j()+1, 0, l), j, d);
 		c.add(M_SQRT1_2 * bT(l-2*lambda.k()), d);
 	      }
	    }
	  }
      }
  }
  
  template <class PRIMALMASK, class DUALMASK>
  SampledMapping<1>
  RBasis<PRIMALMASK, DUALMASK>::evaluate(const unsigned int derivative,
					 const RIndex& lambda,
					 const bool primal,
					 const int A, const int B,
					 const int resolution) const
  {
    if (lambda.e() != 0)
      {
 	InfiniteVector<double, Index> coeffs, gcoeffs;
 	coeffs[lambda] = 1.0;
 	return evaluate(derivative, coeffs, primal, A, B, resolution);
      }
    
    return (primal
	    ? a().evaluate(derivative,
			   lambda.j(),
			   lambda.k(),
			   A,
			   B,
			   resolution)
 	    : aT().evaluate(derivative,
			    lambda.j(),
			    lambda.k(),
			    A,
			    B,
			    resolution));
  }
  
  template <class PRIMALMASK, class DUALMASK>
  SampledMapping<1>
  RBasis<PRIMALMASK, DUALMASK>::evaluate(const unsigned int derivative,
					 const InfiniteVector<double, Index>& coeffs,
					 const bool primal,
					 const int A, const int B,
					 const int resolution) const
  {
    Grid<1> grid(A, B, (B-A)*(1<<resolution));
    Array1D<double> values((B-A)*(1<<resolution)+1);
    for (unsigned int i(0); i < values.size(); i++) values[i] = 0;
    
    if (coeffs.size() > 0)
      {
	// determine maximal level
	int jmax(0);
	for (InfiniteVector<double, Index>::const_iterator it(coeffs.begin()), itend(coeffs.end());
	     it != itend; ++it)
	  jmax = std::max(it.index().j()+it.index().e(), jmax);
	
	// switch to generator representation
	InfiniteVector<double, Index> gcoeffs;
	if (primal)
	  reconstruct(coeffs,jmax,gcoeffs);
	else
	  reconstruct_t(coeffs,jmax,gcoeffs);
	
	for (InfiniteVector<double,RIndex>::const_iterator it(gcoeffs.begin()), itend(gcoeffs.end());
	     it != itend; ++it)
	  {
	    SampledMapping<1> help(evaluate(derivative, it.index(), primal, A, B, resolution));
	    for (unsigned int i(0); i < values.size(); i++)
	      values[i] += *it * help.values()[i];
	  }
      }
    
    return SampledMapping<1>(grid, values);
  }
}
