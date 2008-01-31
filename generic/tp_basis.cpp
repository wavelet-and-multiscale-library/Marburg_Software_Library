// implementation for tp_basis.h

namespace WaveletTL
{
  template <class BASIS0, class BASIS1>
  TensorProductBasis<BASIS0,BASIS1>::TensorProductBasis()
  {
  }

  template <class BASIS0, class BASIS1>
  TensorProductIndex<BASIS0,BASIS1>
  TensorProductBasis<BASIS0,BASIS1>::first_generator(const int j)
  {
    assert(j >= j0());
    return TensorProductIndex<BASIS0,BASIS1>(BASIS0::first_generator(j),
 					     BASIS1::first_generator(j));
  }

  template <class BASIS0, class BASIS1>
  TensorProductIndex<BASIS0,BASIS1>
  TensorProductBasis<BASIS0,BASIS1>::last_generator(const int j)
  {
    assert(j >= j0());
    return TensorProductIndex<BASIS0,BASIS1>(BASIS0::last_generator(j),
 					     BASIS1::last_generator(j));
  }

  template <class BASIS0, class BASIS1>
  TensorProductIndex<BASIS0,BASIS1>
  TensorProductBasis<BASIS0,BASIS1>::first_wavelet(const int j)
  {
    assert(j >= j0());
    return TensorProductIndex<BASIS0,BASIS1>(BASIS0::first_generator(j),
 					     BASIS1::first_wavelet(j));
  }

  template <class BASIS0, class BASIS1>
  TensorProductIndex<BASIS0,BASIS1>
  TensorProductBasis<BASIS0,BASIS1>::last_wavelet(const int j)
  {
    assert(j >= j0());
    return TensorProductIndex<BASIS0,BASIS1>(BASIS0::last_wavelet(j),
 					     BASIS1::last_wavelet(j));
  }

  template <class BASIS0, class BASIS1>
  void
  TensorProductBasis<BASIS0,BASIS1>::decompose(const InfiniteVector<double, Index>& c,
					       const int j0,
					       InfiniteVector<double, Index>& d) const {
    InfiniteVector<double, Index> help;
    for (typename InfiniteVector<double, Index>::const_iterator it(c.begin()), itend(c.end());
	 it != itend; ++it) {
      decompose_1(it.index(), j0, help); // calls help.clear() first
      d.add(*it, help);
    }
  }
  
//   template <class BASIS0, class BASIS1>
//   void
//   TensorProductBasis<BASIS0,BASIS1>::decompose_t(const InfiniteVector<double, Index>& c,
// 						 const int j0,
// 						 InfiniteVector<double, Index>& d) const {
//     InfiniteVector<double, Index> help;
//     for (typename InfiniteVector<double, Index>::const_iterator it(c.begin()), itend(c.end());
// 	 it != itend; ++it) {
//       decompose_t_1(it.index(), j0, help); // calls help.clear() first
//       d.add(*it, help);
//     }
//   }

  template <class BASIS0, class BASIS1>
  void
  TensorProductBasis<BASIS0,BASIS1>::reconstruct(const InfiniteVector<double, Index>& c,
						 const int j,
						 InfiniteVector<double, Index>& d) const {
    for (typename InfiniteVector<double, Index>::const_iterator it(c.begin()), itend(c.end());
	 it != itend; ++it) {
      InfiniteVector<double, Index> help;
      reconstruct_1(it.index(), j, help);
      d.add(*it, help);
    }
  }
  
//   template <class BASIS0, class BASIS1>
//   void
//   TensorProductBasis<BASIS0,BASIS1>::reconstruct_t(const InfiniteVector<double, Index>& c,
// 						   const int j,
// 						   InfiniteVector<double, Index>& d) const {
//     for (typename InfiniteVector<double, Index>::const_iterator it(c.begin()), itend(c.end());
// 	 it != itend; ++it) {
//       InfiniteVector<double, Index> help;
//       reconstruct_t_1(it.index(), j, help);
//       d.add(*it, help);
//     }
//   }

  template <class BASIS0, class BASIS1>
  void
  TensorProductBasis<BASIS0,BASIS1>::decompose_1(const Index& lambda,
						 const int j0,
						 InfiniteVector<double, Index>& c) const {
    assert(lambda.j() >= j0);
    c.clear();
    if (lambda.index0().e() != 0 || lambda.index1().e() != 0) {
      // the true wavelet coefficients don't have to be modified
      c.set_coefficient(lambda, 1.0);
    } else {
      // a generator on a (possibly) fine level
      if (lambda.j() == j0) {
 	// generators from the coarsest level can be copied
 	c.set_coefficient(lambda, 1.0);
      }	else {
 	// j>j0, perform multiscale decomposition

 	typedef typename BASIS0::Index Index0;
 	typedef typename BASIS1::Index Index1;
 	InfiniteVector<double,Index0> c1;
 	InfiniteVector<double,Index1> c2;
	basis0().decompose_1(lambda.index0(), lambda.j()-1, c1);
  	basis1().decompose_1(lambda.index1(), lambda.j()-1, c2);

 	for (typename InfiniteVector<double,Index0>::const_iterator it1(c1.begin()), it1end(c1.end());
  	     it1 != it1end; ++it1)
  	  for (typename InfiniteVector<double,Index1>::const_iterator it2(c2.begin()), it2end(c2.end());
  	       it2 != it2end; ++it2) {
// 	    if (it1.index().e() == 0 && it2.index().e() == 0) { // generators have to be refined further
	    InfiniteVector<double,Index> d;
	    decompose_1(Index(it1.index(), it2.index()), j0, d);
	    c.add(*it1 * *it2, d);
// 	    } else
// 	      c.set_coefficient(Index(this, it1.index(), it2.index()), *it1 * *it2);
	  }
      }
    }
  }
  
//   template <class BASIS0, class BASIS1>
//   void
//   TensorProductBasis<BASIS0,BASIS1>::decompose_t_1(const Index& lambda,
// 						   const int j0,
// 						   InfiniteVector<double, Index>& c) const {
//     assert(lambda.j() >= j0);
//     c.clear();
//     if (lambda.index0().e() != 0 || lambda.index1().e() != 0) {
//       // the true wavelet coefficients don't have to be modified
//       c.set_coefficient(lambda, 1.0);
//     } else {
//       // a generator on a (possibly) fine level
//       if (lambda.j() == j0) {
//  	// generators from the coarsest level can be copied
//  	c.set_coefficient(lambda, 1.0);
//       }	else {
//  	// j>j0, perform multiscale decomposition

//  	typedef typename BASIS0::Index Index0;
//  	typedef typename BASIS1::Index Index1;
//  	InfiniteVector<double,Index0> c1;
//  	InfiniteVector<double,Index1> c2;
// 	basis0().decompose_t_1(lambda.index0(), lambda.j()-1, c1);
//   	basis1().decompose_t_1(lambda.index1(), lambda.j()-1, c2);

//  	for (typename InfiniteVector<double,Index0>::const_iterator it1(c1.begin()), it1end(c1.end());
//   	     it1 != it1end; ++it1)
//   	  for (typename InfiniteVector<double,Index1>::const_iterator it2(c2.begin()), it2end(c2.end());
//   	       it2 != it2end; ++it2) {
// // 	    if (it1.index().e() == 0 && it2.index().e() == 0) { // generators have to be refined further
// 	      InfiniteVector<double,Index> d;
// 	      decompose_t_1(Index(it1.index(), it2.index()), j0, d);
// 	      c.add(*it1 * *it2, d);
// // 	    } else
// // 	      c.set_coefficient(Index(this, it1.index(), it2.index()), *it1 * *it2);
// 	  }
//       }
//     }
//   }

  template <class BASIS0, class BASIS1>
  void
  TensorProductBasis<BASIS0,BASIS1>::reconstruct_1(const Index& lambda,
						   const int j,
						   InfiniteVector<double, Index>& c) const {
    if (lambda.j() >= j) {
      // then we can just copy \psi_\lambda
      c.add_coefficient(lambda, 1.0);
    } else {
      // reconstruct by recursion
      
      typedef typename BASIS0::Index Index0;
      typedef typename BASIS1::Index Index1;
      InfiniteVector<double,Index0> c1;
      InfiniteVector<double,Index1> c2;
      basis0().reconstruct_1(lambda.index0(), lambda.j()+1, c1);
      basis1().reconstruct_1(lambda.index1(), lambda.j()+1, c2);

      for (typename InfiniteVector<double,Index0>::const_iterator it1(c1.begin()), it1end(c1.end());
	   it1 != it1end; ++it1)
	for (typename InfiniteVector<double,Index1>::const_iterator it2(c2.begin()), it2end(c2.end());
	     it2 != it2end; ++it2) {
	  InfiniteVector<double,Index> d;
	  reconstruct_1(Index(it1.index(), it2.index()), j, d);
	  c.add(*it1 * *it2, d);
	}
    }
  }
  
//   template <class BASIS0, class BASIS1>
//   void
//   TensorProductBasis<BASIS0,BASIS1>::reconstruct_t_1(const Index& lambda,
// 						     const int j,
// 						     InfiniteVector<double, Index>& c) const {
//     if (lambda.j() >= j) {
//       // then we can just copy \psi_\lambda
//       c.add_coefficient(lambda, 1.0);
//     } else {
//       // reconstruct by recursion
      
//       typedef typename BASIS0::Index Index0;
//       typedef typename BASIS1::Index Index1;
//       InfiniteVector<double,Index0> c1;
//       InfiniteVector<double,Index1> c2;
//       basis0().reconstruct_t_1(lambda.index0(), lambda.j()+1, c1);
//       basis1().reconstruct_t_1(lambda.index1(), lambda.j()+1, c2);

//       for (typename InfiniteVector<double,Index0>::const_iterator it1(c1.begin()), it1end(c1.end());
// 	   it1 != it1end; ++it1)
// 	for (typename InfiniteVector<double,Index1>::const_iterator it2(c2.begin()), it2end(c2.end());
// 	     it2 != it2end; ++it2) {
// 	  InfiniteVector<double,Index> d;
// 	  reconstruct_t_1(Index(it1.index(), it2.index()), j, d);
// 	  c.add(*it1 * *it2, d);
// 	}
//     }
//   }

}
