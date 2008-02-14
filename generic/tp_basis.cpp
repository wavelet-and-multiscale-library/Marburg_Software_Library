// implementation for tp_basis.h

#include <algebra/kronecker_matrix.h>

using MathTL::KroneckerMatrix;

namespace WaveletTL
{
  template <class BASIS0, class BASIS1>
  TensorProductBasis<BASIS0,BASIS1>::TensorProductBasis()
  {
  }

  template <class BASIS0, class BASIS1>
  inline int
  TensorProductBasis<BASIS0,BASIS1>::Deltasize(const int j)
  {
    assert(j >= j0());
    return BASIS0::Deltasize(j) * BASIS1::Deltasize(j);
  }
  
  template <class BASIS0, class BASIS1>
  inline int
  TensorProductBasis<BASIS0,BASIS1>::Nabla01size(const int j)
  {
    assert(j >= j0());
    return BASIS0::Deltasize(j) * (1<<j);
  }

  template <class BASIS0, class BASIS1>
  inline int
  TensorProductBasis<BASIS0,BASIS1>::Nabla10size(const int j)
  {
    assert(j >= j0());
    return BASIS1::Deltasize(j) * (1<<j);
  }
  
  template <class BASIS0, class BASIS1>
  inline int
  TensorProductBasis<BASIS0,BASIS1>::Nabla11size(const int j)
  {
    assert(j >= j0());
    return 1<<(2*j);
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
  template <class V>
  void
  TensorProductBasis<BASIS0,BASIS1>::apply_Mj0(const int j, const V& x, V& y,
					       const size_type x_offset, const size_type y_offset,
					       const bool add_to) const
  {
    basis0_.Mj0_.set_level(j);
    basis1_.Mj0_.set_level(j);

    KroneckerHelper<double, typename BASIS0::QuasiStationaryMatrixType, typename BASIS1::QuasiStationaryMatrixType>
      K(basis0_.Mj0_, basis1_.Mj0_);
    K.apply(x, y, x_offset, y_offset, add_to);
  }

  template <class BASIS0, class BASIS1>
  template <class V>
  void
  TensorProductBasis<BASIS0,BASIS1>::apply_Mj1_01(const int j, const V& x, V& y,
						  const size_type x_offset, const size_type y_offset,
						  const bool add_to) const
  {
    basis0_.Mj0_.set_level(j);
    basis1_.Mj1_.set_level(j);

    KroneckerHelper<double, typename BASIS0::QuasiStationaryMatrixType, typename BASIS1::QuasiStationaryMatrixType>
      K(basis0_.Mj0_, basis1_.Mj1_);
    K.apply(x, y, x_offset, y_offset, add_to);
  }

  template <class BASIS0, class BASIS1>
  template <class V>
  void
  TensorProductBasis<BASIS0,BASIS1>::apply_Mj1_10(const int j, const V& x, V& y,
						  const size_type x_offset, const size_type y_offset,
						  const bool add_to) const
  {
    basis0_.Mj1_.set_level(j);
    basis1_.Mj0_.set_level(j);

    KroneckerHelper<double, typename BASIS0::QuasiStationaryMatrixType, typename BASIS1::QuasiStationaryMatrixType>
      K(basis0_.Mj1_, basis1_.Mj0_);
    K.apply(x, y, x_offset, y_offset, add_to);
  }

  template <class BASIS0, class BASIS1>
  template <class V>
  void
  TensorProductBasis<BASIS0,BASIS1>::apply_Mj1_11(const int j, const V& x, V& y,
						  const size_type x_offset, const size_type y_offset,
						  const bool add_to) const
  {
    basis0_.Mj1_.set_level(j);
    basis1_.Mj1_.set_level(j);

    KroneckerHelper<double, typename BASIS0::QuasiStationaryMatrixType, typename BASIS1::QuasiStationaryMatrixType>
      K(basis0_.Mj1_, basis1_.Mj1_);
    K.apply(x, y, x_offset, y_offset, add_to);
  }
  
  template <class BASIS0, class BASIS1>
  template <class V>
  void
  TensorProductBasis<BASIS0,BASIS1>::apply_Mj0T_transposed(const int j, const V& x, V& y,
							   const size_type x_offset, const size_type y_offset,
							   const bool add_to) const
  {
    basis0_.Mj0T_.set_level(j);
    basis1_.Mj0T_.set_level(j);

    KroneckerHelper<double, typename BASIS0::QuasiStationaryMatrixType, typename BASIS1::QuasiStationaryMatrixType>
      K(basis0_.Mj0T_, basis1_.Mj0T_);
    K.apply_transposed(x, y, x_offset, y_offset, add_to);
  }

  template <class BASIS0, class BASIS1>
  template <class V>
  void
  TensorProductBasis<BASIS0,BASIS1>::apply_Mj1T_01_transposed(const int j, const V& x, V& y,
							      const size_type x_offset, const size_type y_offset,
							      const bool add_to) const
  {
    basis0_.Mj0T_.set_level(j);
    basis1_.Mj1T_.set_level(j);

    KroneckerHelper<double, typename BASIS0::QuasiStationaryMatrixType, typename BASIS1::QuasiStationaryMatrixType>
      K(basis0_.Mj0T_, basis1_.Mj1T_);
    K.apply_transposed(x, y, x_offset, y_offset, add_to);
  }

  template <class BASIS0, class BASIS1>
  template <class V>
  void
  TensorProductBasis<BASIS0,BASIS1>::apply_Mj1T_10_transposed(const int j, const V& x, V& y,
							      const size_type x_offset, const size_type y_offset,
							      const bool add_to) const
  {
    basis0_.Mj1T_.set_level(j);
    basis1_.Mj0T_.set_level(j);

    KroneckerHelper<double, typename BASIS0::QuasiStationaryMatrixType, typename BASIS1::QuasiStationaryMatrixType>
      K(basis0_.Mj1T_, basis1_.Mj0T_);
    K.apply_transposed(x, y, x_offset, y_offset, add_to);
  }
  
  template <class BASIS0, class BASIS1>
  template <class V>
  void
  TensorProductBasis<BASIS0,BASIS1>::apply_Mj1T_11_transposed(const int j, const V& x, V& y,
							      const size_type x_offset, const size_type y_offset,
							      const bool add_to) const
  {
    basis0_.Mj1T_.set_level(j);
    basis1_.Mj1T_.set_level(j);

    KroneckerHelper<double, typename BASIS0::QuasiStationaryMatrixType, typename BASIS1::QuasiStationaryMatrixType>
      K(basis0_.Mj1T_, basis1_.Mj1T_);
    K.apply_transposed(x, y, x_offset, y_offset, add_to);
  }
  
  template <class BASIS0, class BASIS1>
  template <class V>
  void
  TensorProductBasis<BASIS0,BASIS1>::apply_Mj(const int j, const V& x, V& y) const
  {
    // decompose x appropriately
    apply_Mj0   (j, x, y, 0, 0, false);                                         // apply Mj0
    apply_Mj1_01(j, x, y, Deltasize(j), 0, true);                               // apply Mj1, e=(0,1)
    apply_Mj1_10(j, x, y, Deltasize(j)+Nabla01size(j), 0, true);                // apply Mj1, e=(1,0)
    apply_Mj1_11(j, x, y, Deltasize(j)+Nabla01size(j)+Nabla10size(j), 0, true); // apply Mj1, e=(1,1)
  }

  template <class BASIS0, class BASIS1>
  template <class V>
  void
  TensorProductBasis<BASIS0,BASIS1>::apply_Gj(const int j, const V& x, V& y) const
  {
    // write into the block vector y
    apply_Mj0T_transposed   (j, x, y, 0, 0, false);                                          // write Mj0T block
    apply_Mj1T_01_transposed(j, x, y, 0, Deltasize(j), false);                               // write Mj1T block, e=(0,1)
    apply_Mj1T_10_transposed(j, x, y, 0, Deltasize(j)+Nabla01size(j), false);                // write Mj1T block, e=(1,0)
    apply_Mj1T_11_transposed(j, x, y, 0, Deltasize(j)+Nabla01size(j)+Nabla10size(j), false); // write Mj1T block, e=(1,1)
  }

  template <class BASIS0, class BASIS1>
  template <class V>
  void
  TensorProductBasis<BASIS0,BASIS1>::apply_Tj(const int j, const V& x, V& y) const
  { 
    y = x;
    V z(x);
    apply_Mj(j0(), z, y);
    for (int k = j0()+1; k <= j; k++) {
      apply_Mj(k, y, z);
      y.swap(z);
    }
  }

  template <class BASIS0, class BASIS1>
  template <class V>
  void
  TensorProductBasis<BASIS0,BASIS1>::apply_Tjinv(const int j, const V& x, V& y) const
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
