// implementation of r_basis.h

#include <cmath>
#include <cassert>

namespace WaveletTL
{
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
	d += *it * help;
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
	d += *it * help;
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
	d += *it * help;
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
	d += *it * help;
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
	    
	    const int abegin(a().begin().index()),
	      aend(a().rbegin().index()),
	      bbegin(1-at().rbegin().index()),
	      bend(1-at().begin().index());
	    
	    // compute d_{j-1}
 	    for (int l((int)ceil((lambda.k()-bend)/2.0));
		 l <= (int)floor((lambda.k()-bbegin)/2.0); l++)
 	      c[Index(lambda.j()-1, 1, l)] = M_SQRT1_2 * b(lambda.k()-2*l);

	    // compute c_{j_0} via recursion
	    for (int l((int)ceil((lambda.k()-aend)/2.0));
		 l <= (int)floor((lambda.k()-abegin)/2.0); l++)
	      {
		InfiniteVector<double, Index> d;
		decompose_1(Index(lambda.j()-1, 0, l), j0, d);
		c += M_SQRT1_2 * a().get_coefficient(lambda.k()-2*l) * d;
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
	    
	    const int atbegin(at().begin().index()),
	      atend(at().rbegin().index()),
	      btbegin(1-a().rbegin().index()),
	      btend(1-a().begin().index());
	    
	    // compute d_{j-1}
 	    for (int l((int)ceil((lambda.k()-btend)/2.0));
		 l <= (int)floor((lambda.k()-btbegin)/2.0); l++)
 	      c[Index(lambda.j()-1, 1, l)] = M_SQRT1_2 * bt(lambda.k()-2*l);
	    
	    // compute c_{j_0} via recursion
	    for (int l((int)ceil((lambda.k()-atend)/2.0));
		 l <= (int)floor((lambda.k()-atbegin)/2.0); l++)
	      {
		InfiniteVector<double, Index> d;
		decompose_t_1(Index(lambda.j()-1, 0, l), j0, d);
		c += M_SQRT1_2 * at().get_coefficient(lambda.k()-2*l) * d;
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

	const int atbegin(at().begin().index()),
	  atend(at().rbegin().index()),
	  btbegin(1-a().rbegin().index()),
	  btend(1-a().begin().index());

	if (lambda.e() == 0)
	  {
 	    for (int l(2*lambda.k()+atbegin); l <= 2*lambda.k()+atend; l++)
 	      {
 		InfiniteVector<double, Index> d;
 		reconstruct_1(Index(lambda.j()+1, 0, l), j, d);
 		c += M_SQRT1_2 * at().get_coefficient(l-2*lambda.k()) * d;
 	      }
 	  }
 	else
 	  {
 	    for (int l(2*lambda.k()+btbegin); l <= 2*lambda.k()+btend; l++)
 	      {
 		InfiniteVector<double, Index> d;
 		reconstruct_1(Index(lambda.j()+1, 0, l), j, d);
 		c += M_SQRT1_2 * bt(l-2*lambda.k()) * d;
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

	const int abegin(a().begin().index()),
	  aend(a().rbegin().index()),
	  bbegin(1-at().rbegin().index()),
	  bend(1-at().begin().index());
	    
	if (lambda.e() == 0)
	  {
 	    for (int l(2*lambda.k()+abegin); l <= 2*lambda.k()+aend; l++)
 	      {
 		InfiniteVector<double, Index> d;
 		reconstruct_t_1(Index(lambda.j()+1, 0, l), j, d);
 		c += M_SQRT1_2 * a().get_coefficient(l-2*lambda.k()) * d;
 	      }
 	  }
 	else
 	  {
 	    for (int l(2*lambda.k()+bbegin); l <= 2*lambda.k()+bend; l++)
 	      {
 		InfiniteVector<double, Index> d;
 		reconstruct_t_1(Index(lambda.j()+1, 0, l), j, d);
 		c += M_SQRT1_2 * b(l-2*lambda.k()) * d;
 	      }
	  }
      }
  }

  template <class PRIMALMASK, class DUALMASK>
  template <unsigned int DERIVATIVE>
  SampledMapping<1>
  RBasis<PRIMALMASK, DUALMASK>::evaluate(const RIndex& lambda,
					 const bool primal,
					 const int A, const int B,
					 const int resolution) const
  {
    if (lambda.e() != 0)
      {
 	InfiniteVector<double, Index> coeffs, gcoeffs;
 	coeffs[lambda] = 1.0;
 	if (primal)
 	  reconstruct_1(lambda, lambda.j()+1, gcoeffs);
 	else
 	  reconstruct_t_1(lambda, lambda.j()+1, gcoeffs);
 	return evaluate<DERIVATIVE>(gcoeffs, primal, A, B, resolution);
      }
    
    return (primal
	    ? a().evaluate<DERIVATIVE>(lambda.j(), lambda.k(), A, B, resolution)
 	    : at().evaluate<DERIVATIVE>(lambda.j(), lambda.k(), A, B, resolution));
  }

  template <class PRIMALMASK, class DUALMASK>
  template <unsigned int DERIVATIVE>
  SampledMapping<1>
  RBasis<PRIMALMASK, DUALMASK>::evaluate(const InfiniteVector<double, Index>& coeffs,
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
	    SampledMapping<1> help(evaluate<DERIVATIVE>(it.index(), primal, A, B, resolution));
	    for (unsigned int i(0); i < values.size(); i++)
	      values[i] += *it * help.values()[i];
	  }
      }
    
    return SampledMapping<1>(grid, values);
  }
}
