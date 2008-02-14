// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2008                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_TP_BASIS_H
#define _WAVELETTL_TP_BASIS_H

#include <cassert>
#include <algebra/vector.h>
#include <algebra/infinite_vector.h>

#include <generic/tp_index.h>

namespace WaveletTL
{
  /*!
    Template class for a tensor product wavelet basis from two wavelet bases
    Psi, Xi over bounded domains.
    The template parameters BASISi may or may not allow the specification of
    boundary conditions.
  */
  template <class BASIS0, class BASIS1>
  class TensorProductBasis
  {
  public:
    //! size_type, for convenience
    typedef Vector<double>::size_type size_type;

    //! default constructor
    TensorProductBasis();

    //! space dimension of the underlying domain
    static const int space_dimension = BASIS0::space_dimension+BASIS1::space_dimension;

    //! coarsest possible level j0
    static const int j0()
    {
      assert(BASIS0::j0() == BASIS1::j0());
      return BASIS0::j0();
    }    

    //! wavelet index class
    typedef TensorProductIndex<BASIS0,BASIS1> Index;

    //! read access to the first basis
    const BASIS0& basis0() const { return basis0_; }

    //! read access to the second basis
    const BASIS1& basis1() const { return basis1_; }

    //! size of Delta_j
    static int Deltasize(const int j);

    //! sizes of the different wavelet index sets
    static int Nabla01size(const int j);
    static int Nabla10size(const int j);
    static int Nabla11size(const int j);

    //! index of first generator on level j >= j0
    static Index first_generator(const int j);
      
    //! index of last generator on level j >= j0
    static Index last_generator(const int j);
      
    //! index of first wavelet on level j >= j0
    static Index first_wavelet(const int j);
      
    //! index of last wavelet on level j >= j0
    static Index last_wavelet(const int j);

    /*!
      apply Mj0 to some vector x (partial "reconstruct");
      the routine writes only into the first part of y, i.e,
      y might be larger than necessary, which is helpful for other routines;
      offsets and an add_to flag can be specified also
    */
    template <class V>
    void apply_Mj0(const int j, const V& x, V& y,
		   const size_type x_offset, const size_type y_offset,
		   const bool add_to) const;

    /*!
      an analogous routine for Mj1, e=(0,1)
    */
    template <class V>
    void apply_Mj1_01(const int j, const V& x, V& y,
		      const size_type x_offset, const size_type y_offset,
		      const bool add_to) const;

    /*!
      an analogous routine for Mj1, e=(1,0)
    */
    template <class V>
    void apply_Mj1_10(const int j, const V& x, V& y,
		      const size_type x_offset, const size_type y_offset,
		      const bool add_to) const;

    /*!
      an analogous routine for Mj1, e=(1,1)
    */
    template <class V>
    void apply_Mj1_11(const int j, const V& x, V& y,
		      const size_type x_offset, const size_type y_offset,
		      const bool add_to) const;

    /*!
      apply Mj=(Mj0 Mj1) to some vector x ("reconstruct");
      the routine writes only into the first part of y, i.e,
      y might be larger than necessary, which is helpful for apply_Tj
    */
    template <class V>
    void apply_Mj(const int j, const V& x, V& y) const;

    //! DECOMPOSE routine, simple version
    /*!
      Constructs for a given single wavelet index lambda a coefficient set c,
      such that
        \psi_lambda = \sum_{\lambda'}c_{\lambda'}\psi_{\lambda'}
      where the multiscale decomposition starts with the coarsest
      generator level jmin.
     */
    void decompose_1(const Index& lambda, const int jmin,
		     InfiniteVector<double, Index>& c) const;

    //! DECOMPOSE routine, full version
    /*!
      constructs for a given coefficient set c another one v with level >= jmin,
      such that
        \sum_{\lambda}c_\lambda\psi_lambda = \sum_{\lambda'}v_{\lambda'}\psi_{\lambda'}
    */
    void decompose(const InfiniteVector<double, Index>& c, const int jmin,
		   InfiniteVector<double, Index>& v) const;

//     //! dual DECOMPOSE routine, simple version
//     /*!
//       Constructs for a given single wavelet index lambda a coefficient set c,
//       such that
//         \tilde\psi_lambda = \sum_{\lambda'}c_{\lambda'}\tilde\psi_{\lambda'}
//       where the multiscale decomposition starts with the coarsest
//       generator level jmin.
//      */
//     void decompose_t_1(const Index& lambda, const int jmin,
// 		       InfiniteVector<double, Index>& c) const;

//     //! dual DECOMPOSE routine, full version
//     /*!
//       constructs for a given coefficient set c another one v with level >= jmin,
//       such that
//         \sum_{\lambda}c_\lambda\tilde\psi_lambda = \sum_{\lambda'}d_{\lambda'}\tilde\psi_{\lambda'}
//     */
//     void decompose_t(const InfiniteVector<double, Index>& c, const int jmin,
// 		     InfiniteVector<double, Index>& v) const;

    //! RECONSTRUCT routine, simple version
    /*!
      Constructs for a given single wavelet index lambda a coefficient set c,
      such that
        \psi_lambda = \sum_{\lambda'}c_{\lambda'}\psi_{\lambda'}
      where always |\lambda'|>=j
     */
    void reconstruct_1(const Index& lambda, const int j,
		       InfiniteVector<double, Index>& c) const;

    //! RECONSTRUCT routine, full version
    /*!
      Constructs for a given coefficient set c another one v,
      such that
        \sum_{\lambda}c_\lambda\psi_lambda = \sum_{\lambda'}v_{\lambda'}\psi_{\lambda'}
      where always |\lambda'|>=j
    */
    void reconstruct(const InfiniteVector<double, Index>& c, const int j,
		     InfiniteVector<double, Index>& v) const;

//     //! dual RECONSTRUCT routine, simple version
//     /*!
//       Constructs for a given single wavelet index lambda a coefficient set c,
//       such that
//         \tilde\psi_lambda = \sum_{\lambda'}c_{\lambda'}\tilde\psi_{\lambda'}
//       where always |\lambda'|>=j
//      */
//     void reconstruct_t_1(const Index& lambda, const int j,
// 			 InfiniteVector<double, Index>& c) const;

//     //! dual RECONSTRUCT routine, full version
//     /*!
//       Constructs for a given coefficient set c another one v,
//       such that
//         \sum_{\lambda}c_\lambda\tilde\psi_\lambda = \sum_{\lambda'}v_{\lambda'}\tilde\psi_{\lambda'}
//       where always |\lambda'|>=j
//     */
//     void reconstruct_t(const InfiniteVector<double, Index>& c, const int j,
// 		       InfiniteVector<double, Index>& v) const;

  protected:
    //! instances of the two 1D bases
    BASIS0 basis0_;
    BASIS1 basis1_;
  };
}

#include <generic/tp_basis.cpp>

#endif
