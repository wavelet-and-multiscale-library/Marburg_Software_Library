// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_SPLINE_BASIS_H
#define _WAVELETTL_SPLINE_BASIS_H

#include <interval/spline_basis_data.h>

namespace WaveletTL
{
  /*!
    Template class for spline wavelet bases on the interval.
  */
  template <int d, int dT>
  class SplineBasis
    : public SplineBasisData<d,dT>
  {
  public:
    /*!
      constructor from the primal and dual boundary condition orders,
      cf. SplineBasisData
    */
    SplineBasis(const char* flavor,
		const char* options,
		const int s0, const int s1, const int sT0, const int sT1);

    //! coarsest possible level
    inline const int j0() const { return SplineBasisData<d,dT>::j0_; }

    //! size_type, for convenience
    typedef Vector<double>::size_type size_type;

    //! read access to the primal b.c. order at x=0
    const int get_s0() const { return SplineBasisData<d,dT>::s0_; }

    //! read access to the primal b.c. order at x=1
    const int get_s1() const { return SplineBasisData<d,dT>::s1_; }

    //! read access to the primal b.c. order at x=0
    const int get_sT0() const { return SplineBasisData<d,dT>::sT0_; }

    //! read access to the primal b.c. order at x=1
    const int get_sT1() const { return SplineBasisData<d,dT>::sT1_; }

    //! extremal generator indices
    inline const int DeltaLmin() const { return DeltaLmin_; }
    inline const int DeltaRmax(const int j) const { return (1<<j)+DeltaRmax_offset; }

    //! size of Delta_j
    inline const int Deltasize(const int j) const { return DeltaRmax(j)-DeltaLmin()+1; }
    
    //! boundary indices in \nabla_j
    inline const int Nablamin() const { return 0; }
    inline const int Nablamax(const int j) const { return (1<<j)-1; }
    
    //! size of Nabla_j
    inline const int Nablasize(const int j) const { return 1<<j; }
    
    /*!
      apply Mj=(Mj0 Mj1) to some vector x ("reconstruct");
      the routine writes only into the first part of y, i.e,
      y might be larger than necessary, which is helpful for apply_Tj
    */
    template <class V>
    void apply_Mj(const int j, const V& x, V& y) const;

    //! apply Mj^T to some vector x
    template <class V>
    void apply_Mj_transposed(const int j, const V& x, V& y) const;

    //! apply Gj=(Mj0T Mj1T)^T to some vector x ("decompose")
    template <class V>
    void apply_Gj(const int j, const V& x, V& y) const;

    //! apply G_j^T to some vector x
    template <class V>
    void apply_Gj_transposed(const int j, const V& x, V& y) const;

    //! apply Tj=Mj*diag(M_{j-1},I)*...*diag(M_{j_0},I), i.e., several "reconstructions" at once
    template <class V>
    void apply_Tj(const int j, const V& x, V& y) const;

    /*!
      apply Tj^T,
      this generic variant works with "infinite" vector classes like
      std::map<size_type,double> or InfiniteVector<size_type,double>
    */
    template <class V>
    void apply_Tj_transposed(const int j, const V& x, V& y) const;

    //! ... dito, specialization for V=Vector<double>
    void apply_Tj_transposed(const int j, const Vector<double>& x, Vector<double>& y) const;

    //! apply Tj^{-1}, several "decompositions" at once
    void apply_Tjinv(const int j, const Vector<double>& x, Vector<double>& y) const;

  protected:
    int DeltaLmin_, DeltaRmax_offset;
  };
}

#include <interval/spline_basis.cpp>

#endif
