// implementation for tp_basis.h

namespace WaveletTL
{
  template <class BASIS1, class BASIS2>
  TensorProductBasis<BASIS1,BASIS2>::TensorProductBasis()
  {
    j0_ = std::max(basis1_.j0(), basis2_.j0());
  }

  template <class BASIS1, class BASIS2>
  void
  TensorProductBasis<BASIS1,BASIS2>::decompose(const InfiniteVector<double, Index>& c,
					       const int j0,
					       InfiniteVector<double, Index>& d) const {
    for (typename InfiniteVector<double, Index>::const_iterator it(c.begin()), itend(c.end());
	 it != itend; ++it) {
      InfiniteVector<double, Index> help;
      decompose_1(it.index(), j0, help);
      d.add(*it, help);
    }
  }
  
  template <class BASIS1, class BASIS2>
  void
  TensorProductBasis<BASIS1,BASIS2>::decompose_t(const InfiniteVector<double, Index>& c,
						 const int j0,
						 InfiniteVector<double, Index>& d) const {
    for (typename InfiniteVector<double, Index>::const_iterator it(c.begin()), itend(c.end());
	 it != itend; ++it) {
      InfiniteVector<double, Index> help;
      decompose_t_1(it.index(), j0, help);
      d.add(*it, help);
    }
  }

  template <class BASIS1, class BASIS2>
  void
  TensorProductBasis<BASIS1,BASIS2>::reconstruct(const InfiniteVector<double, Index>& c,
						 const int j,
						 InfiniteVector<double, Index>& d) const {
    for (typename InfiniteVector<double, Index>::const_iterator it(c.begin()), itend(c.end());
	 it != itend; ++it) {
      InfiniteVector<double, Index> help;
      reconstruct_1(it.index(), j, help);
      d.add(*it, help);
    }
  }
  
  template <class BASIS1, class BASIS2>
  void
  TensorProductBasis<BASIS1,BASIS2>::reconstruct_t(const InfiniteVector<double, Index>& c,
						   const int j,
						   InfiniteVector<double, Index>& d) const {
    for (typename InfiniteVector<double, Index>::const_iterator it(c.begin()), itend(c.end());
	 it != itend; ++it) {
      InfiniteVector<double, Index> help;
      reconstruct_t_1(it.index(), j, help);
      d.add(*it, help);
    }
  }

  template <class BASIS1, class BASIS2>
  void
  TensorProductBasis<BASIS1,BASIS2>::decompose_1(const Index& lambda,
						 const int j0,
						 InfiniteVector<double, Index>& c) const {
    assert(lambda.j() >= j0);
    c.clear();
    if (lambda.index1().e() != 0 || lambda.index2().e() != 0) {
      // the true wavelet coefficients don't have to be modified
      c.set_coefficient(lambda, 1.0);
    } else {
      // a generator on a (possibly) fine level
      if (lambda.j() == j0) {
 	// generators from the coarsest level can be copied
 	c.set_coefficient(lambda, 1.0);
      }	else {
 	// j>j0, perform multiscale decomposition

 	typedef typename BASIS1::Index Index1;
 	typedef typename BASIS2::Index Index2;
 	InfiniteVector<double,Index1> c1;
 	InfiniteVector<double,Index2> c2;
	basis1().decompose_1(lambda.index1(), lambda.j()-1, c1);
  	basis2().decompose_1(lambda.index2(), lambda.j()-1, c2);

 	for (typename InfiniteVector<double,Index1>::const_iterator it1(c1.begin()), it1end(c1.end());
  	     it1 != it1end; ++it1)
  	  for (typename InfiniteVector<double,Index2>::const_iterator it2(c2.begin()), it2end(c2.end());
  	       it2 != it2end; ++it2) {
	    if (it1.index().e() == 0 && it2.index().e() == 0) { // generators have to be refined further
	      InfiniteVector<double,Index> d;
	      decompose_1(Index(this, it1.index(), it2.index()), j0, d);
	      c.add(*it1 * *it2, d);
	    } else
	      c.set_coefficient(Index(this, it1.index(), it2.index()), *it1 * *it2);
	  }
      }
    }
  }
  
  template <class BASIS1, class BASIS2>
  void
  TensorProductBasis<BASIS1,BASIS2>::decompose_t_1(const Index& lambda,
						   const int j0,
						   InfiniteVector<double, Index>& c) const {
  }

  template <class BASIS1, class BASIS2>
  void
  TensorProductBasis<BASIS1,BASIS2>::reconstruct_1(const Index& lambda,
						   const int j,
						   InfiniteVector<double, Index>& c) const {
    if (lambda.j() >= j) {
      // then we can just copy \psi_\lambda
      c.add_coefficient(lambda, 1.0);
    } else {
      // reconstruct by recursion
      
      typedef typename BASIS1::Index Index1;
      typedef typename BASIS2::Index Index2;
      InfiniteVector<double,Index1> c1;
      InfiniteVector<double,Index2> c2;
      basis1().reconstruct_1(lambda.index1(), lambda.j()+1, c1);
      basis2().reconstruct_1(lambda.index2(), lambda.j()+1, c2);

      for (typename InfiniteVector<double,Index1>::const_iterator it1(c1.begin()), it1end(c1.end());
	   it1 != it1end; ++it1)
	for (typename InfiniteVector<double,Index2>::const_iterator it2(c2.begin()), it2end(c2.end());
	     it2 != it2end; ++it2) {
	  InfiniteVector<double,Index> d;
	  reconstruct_1(Index(this, it1.index(), it2.index()), j, d);
	  c.add(*it1 * *it2, d);
	}
    }
  }
  
  template <class BASIS1, class BASIS2>
  void
  TensorProductBasis<BASIS1,BASIS2>::reconstruct_t_1(const Index& lambda,
						     const int j,
						     InfiniteVector<double, Index>& c) const {
  }

}
