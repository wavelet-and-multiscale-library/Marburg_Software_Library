// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2007                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_SPLINE_BASIS_H
#define _WAVELETTL_SPLINE_BASIS_H

#include <geometry/sampled_mapping.h>
#include <algebra/infinite_vector.h>

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
      single primal spline generator or wavelet \psi_\lambda.
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

    
    /*!
      Evaluate an arbitrary linear combination of primal wavelets
      on a dyadic subgrid of [0,1].
    */
    SampledMapping<1>
    evaluate
    (const InfiniteVector<double, typename SplineBasis<d,dT,flavor>::Index>& coeffs,
     const int resolution) const;
    
    /*!
      point evaluation of (derivatives) of a single primal generator
      or wavelet \psi_\lambda
    */
    double
    evaluate
    (const unsigned int derivative,
     const typename SplineBasis<d,dT,flavor>::Index& lambda,
     const double x) const;
    
    /*!
      point evaluation of (derivatives) of a single primal generator
      or wavelet \psi_\lambda at several points simultaneously
    */
    void
    evaluate
    (const unsigned int derivative,
     const typename SplineBasis<d,dT,flavor>::Index& lambda,
     const Array1D<double>& points, Array1D<double>& values) const;

    /*!
      point evaluation of 0-th and first derivative of a single primal generator
      or wavelet \psi_\lambda at several points simultaneously
    */
    void
    evaluate
    (const typename SplineBasis<d,dT,flavor>::Index& lambda,
     const Array1D<double>& points, Array1D<double>& funcvalues, Array1D<double>& dervalues) const;


    /*!
      For a given function, compute all integrals w.r.t. the primal
      or dual generators/wavelets \psi_\lambda with |\lambda|\le jmax.
      - When integrating against the primal functions, the integrand has to be smooth
        to be accurately reproduced by the dual basis.
      - When integration against dual functions is specified,
        we integrate against the primal ones instead and multiply the resulting
        coefficients with the inverse of the primal gramian.

      Maybe a thresholding of the returned coefficients is helpful (e.g. for
      expansions of spline functions).
    */
    void
    expand
    (const Function<1>* f,
     const bool primal,
     const int jmax,
     InfiniteVector<double, typename SplineBasis<d,dT,flavor>::Index>& coeffs) const;
    
    /*!
      analogous routine for Vector<double> output
    */
    void
    expand
    (const Function<1>* f,
     const bool primal,
     const int jmax,
     Vector<double>& coeffs) const;
    

  protected:
    int DeltaLmin_, DeltaRmax_offset;
  };

  /*!
    template specialization to flavor==DS_construction
  */
  template <int d, int dT>
  class SplineBasis<d,dT,DS_construction>
    : public SplineBasisData<d,dT,DS_construction>
  {
  public:
    /*!
      constructor from the primal and dual boundary condition orders,
      cf. SplineBasisData
    */
    SplineBasis(const char* options,
		const int s0, const int s1, const int sT0, const int sT1);

    //! coarsest possible level
    inline const int j0() const { return SplineBasisData<d,dT,DS_construction>::j0_; }
    
    //! wavelet index class
    typedef IntervalIndex<SplineBasis<d,dT,DS_construction> > Index;
    
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
      single primal spline generator or wavelet \psi_\lambda.
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
    const int get_s0() const { return SplineBasisData<d,dT,DS_construction>::s0_; }

    //! read access to the primal b.c. order at x=1
    const int get_s1() const { return SplineBasisData<d,dT,DS_construction>::s1_; }

    //! read access to the primal b.c. order at x=0
    const int get_sT0() const { return SplineBasisData<d,dT,DS_construction>::sT0_; }

    //! read access to the primal b.c. order at x=1
    const int get_sT1() const { return SplineBasisData<d,dT,DS_construction>::sT1_; }

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

    
    /*!
      Evaluate a single primal/dual generator or wavelet \psi_\lambda
      on a dyadic subgrid of [0,1].
    */
    SampledMapping<1>
    evaluate
    (const typename SplineBasis<d,dT,DS_construction>::Index& lambda,
     const int resolution) const;

    /*!
      Evaluate an arbitrary linear combination of primal wavelets
      on a dyadic subgrid of [0,1].
    */
    SampledMapping<1>
    evaluate
    (const InfiniteVector<double, typename SplineBasis<d,dT,DS_construction>::Index>& coeffs,
     const int resolution) const;
    
    /*!
      point evaluation of (derivatives) of a single primal generator
      or wavelet \psi_\lambda
    */
    double evaluate(const unsigned int derivative,
		    const typename SplineBasis<d,dT,DS_construction>::Index& lambda,
		    const double x) const;
    
    /*!
      point evaluation of (derivatives) of a single primal generator
      or wavelet \psi_\lambda at several points simultaneously
    */
    void evaluate(const unsigned int derivative,
		  const typename SplineBasis<d,dT,DS_construction>::Index& lambda,
		  const Array1D<double>& points, Array1D<double>& values) const;

    //! a version without an Index object
    void evaluate(const unsigned int derivative,
		  const int j, const int e, const int k,
		  const Array1D<double>& points, Array1D<double>& values) const;
    
    /*!
      point evaluation of 0-th and first derivative of a single primal generator
      or wavelet \psi_\lambda at several points simultaneously
    */
    void evaluate(const typename SplineBasis<d,dT,DS_construction>::Index& lambda,
		  const Array1D<double>& points, Array1D<double>& funcvalues, Array1D<double>& dervalues) const;


    /*!
      For a given function, compute all integrals w.r.t. the primal
      or dual generators/wavelets \psi_\lambda with |\lambda|\le jmax.
      - When integrating against the primal functions, the integrand has to be smooth
        to be accurately reproduced by the dual basis.
      - When integration against dual functions is specified,
        we integrate against the primal ones instead and multiply the resulting
        coefficients with the inverse of the primal gramian.

      Maybe a thresholding of the returned coefficients is helpful (e.g. for
      expansions of spline functions).
    */
    void
    expand
    (const Function<1>* f,
     const bool primal,
     const int jmax,
     InfiniteVector<double, typename SplineBasis<d,dT,DS_construction>::Index>& coeffs) const;
    
    /*!
      analogous routine for Vector<double> output
    */
    void
    expand
    (const Function<1>* f,
     const bool primal,
     const int jmax,
     Vector<double>& coeffs) const;
    
    /*!
      helper function, integrate a smooth function f against a
      primal DS generator or wavelet
    */
    double
    integrate
    (const Function<1>* f,
     const typename SplineBasis<d,dT,DS_construction>::Index& lambda) const;

    /*!
      helper function, integrate two primal DS generators or wavelets
      against each other (for the Gramian)
    */
    double
    integrate
    (const typename SplineBasis<d,dT,DS_construction>::Index& lambda,
     const typename SplineBasis<d,dT,DS_construction>::Index& mu) const;
    
  protected:
    int DeltaLmin_, DeltaRmax_offset;
  };

}

#include <interval/spline_basis.cpp>

#endif
