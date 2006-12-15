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
#include <interval/i_index.h>

#include <interval/spline_support.h>
#include <interval/spline_evaluate.h>

namespace WaveletTL
{
  /*!
    Template class for spline wavelet bases on the interval.
  */
  template <int d, int dT, SplineBasisFlavor flavor>
  class SplineBasis
    : public SplineBasisData<d,dT,flavor>
  {
  public:
    /*!
      constructor from the primal and dual boundary condition orders,
      cf. SplineBasisData
    */
    SplineBasis(const char* options,
		const int s0, const int s1, const int sT0, const int sT1);

    //! coarsest possible level
    inline const int j0() const { return SplineBasisData<d,dT,flavor>::j0_; }

    //! wavelet index class
    typedef IntervalIndex<SplineBasis<d,dT,flavor> > Index;
    
    //! size_type, for convenience
    typedef Vector<double>::size_type size_type;

    //! geometric type of the support sets
    typedef struct {
      int j;
      int k1;
      int k2;
    } Support;

    /*!
      Compute an interval 2^{-j}[k1,k2] which contains the support of a
      single primal [P] generator or wavelet \psi_\lambda.
      (j == lambda.j()+lambda.e() is neglected for performance reasons)
    */
    void support(const Index& lambda, int& k1, int& k2) const;

    //! space dimension of the underlying domain
    static const int space_dimension = 1;

    //! critical Sobolev regularity for the primal generators/wavelets
    static double primal_regularity() { return d - 0.5; }

    //! degree of polynomial reproduction for the primal generators/wavelets
    static unsigned int primal_polynomial_degree() { return d; }

    //! number of vanishing moments for the primal wavelets
    static unsigned int primal_vanishing_moments() { return dT; }

    //! read access to the primal b.c. order at x=0
    const int get_s0() const { return SplineBasisData<d,dT,flavor>::s0_; }

    //! read access to the primal b.c. order at x=1
    const int get_s1() const { return SplineBasisData<d,dT,flavor>::s1_; }

    //! read access to the primal b.c. order at x=0
    const int get_sT0() const { return SplineBasisData<d,dT,flavor>::sT0_; }

    //! read access to the primal b.c. order at x=1
    const int get_sT1() const { return SplineBasisData<d,dT,flavor>::sT1_; }

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
    
    //! index of first (leftmost) generator on level j >= j0
    Index first_generator(const int j) const;

    //! index of last (rightmost) generator on level j >= j0
    Index last_generator(const int j) const;

    //! index of first (leftmost) wavelet on level j >= j0
    Index first_wavelet(const int j) const;

    //! index of last (rightmost) wavelet on level j >= j0
    Index last_wavelet(const int j) const;

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
