// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2008                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_SPLINE_BASIS_H
#define _WAVELETTL_SPLINE_BASIS_H

#include <geometry/sampled_mapping.h>
#include <algebra/infinite_vector.h>

#include <Rd/cdf_utils.h>
#include <interval/spline_basis_data.h>
#include <interval/i_index.h>

#include <interval/spline_support.h>
#include <interval/spline_evaluate.h>

namespace WaveletTL
{
  /*!
    Template class for spline wavelet bases on the interval.

    If unsure about the suitable value of J0, choose
      J0 = SplineBasisData_j0<d,dT,flavor,s0,s1,sT0,sT1>::j0>
  */
  template <int d, int dT, SplineBasisFlavor flavor, int s0, int s1, int sT0, int sT1, int J0>
  class SplineBasis
    : public SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>
  {
  public:
    using SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::check;

    //! default constructor
    SplineBasis();

    //! wavelet index class
    typedef IntervalIndex2<SplineBasis<d,dT,flavor,s0,s1,sT0,sT1,J0> > Index;
    
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

    //! coarsest level
    static const int j0() { return J0; }

    //! critical Sobolev regularity for the primal generators/wavelets
    static double primal_regularity() { return d - 0.5; }

    //! degree of polynomial reproduction for the primal generators/wavelets
    static unsigned int primal_polynomial_degree() { return d; }

    //! number of vanishing moments for the primal wavelets
    static unsigned int primal_vanishing_moments() { return dT; }

    //! read access to the primal b.c. order at x=0
    static const int get_s0() { return s0; }

    //! read access to the primal b.c. order at x=1
    static const int get_s1() { return s1; }

    //! read access to the primal b.c. order at x=0
    static const int get_sT0() { return sT0; }

    //! read access to the primal b.c. order at x=1
    static const int get_sT1() { return sT1; }

    //! extremal generator indices
    static const int DeltaLmin();
    static const int DeltaRmax_offset();
    static const int DeltaRmax(const int j) { return (1<<j)+DeltaRmax_offset(); }
    
    //! size of Delta_j
    static const int Deltasize(const int j) { return DeltaRmax(j)-DeltaLmin()+1; }
    
    //! boundary indices in \nabla_j
    static const int Nablamin() { return 0; }
    static const int Nablamax(const int j) { return (1<<j)-1; }
    
    //! size of Nabla_j
    static const int Nablasize(const int j) { return 1<<j; }
    
    //! index of first (leftmost) generator on level j >= j0
    static Index first_generator(const int j);

    //! index of last (rightmost) generator on level j >= j0
    static Index last_generator(const int j);

    //! index of first (leftmost) wavelet on level j >= j0
    static Index first_wavelet(const int j);

    //! index of last (rightmost) wavelet on level j >= j0
    static Index last_wavelet(const int j);

    /*!
      index of first function with type e
      (mainly for TensorProductBasis)
    */
    static Index first_index(const int j, const int e);

    /*!
      index of last function with type e
      (mainly for TensorProductBasis)
    */
    static Index last_index(const int j, const int e);

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
    (const Index& lambda,
     const int resolution) const;

    /*!
      Evaluate an arbitrary linear combination of primal wavelets
      on a dyadic subgrid of [0,1].
    */
    SampledMapping<1>
    evaluate
    (const InfiniteVector<double, Index>& coeffs,
     const int resolution) const;
    
    /*!
      point evaluation of (derivatives) of a single primal generator
      or wavelet \psi_\lambda
    */
    double
    evaluate
    (const unsigned int derivative,
     const Index& lambda,
     const double x) const;
    
    /*!
      point evaluation of (derivatives) of a single primal generator
      or wavelet \psi_\lambda at several points simultaneously
    */
    void
    evaluate
    (const unsigned int derivative,
     const Index& lambda,
     const Array1D<double>& points, Array1D<double>& values) const;

    /*!
      point evaluation of 0-th and first derivative of a single primal generator
      or wavelet \psi_\lambda at several points simultaneously
    */
    void
    evaluate
    (const Index& lambda,
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
     InfiniteVector<double, Index>& coeffs) const;
    
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
     const Index& lambda) const;

    /*!
      helper function, integrate two primal DS generators or wavelets
      against each other (for the Gramian)
    */
    double
    integrate
    (const Index& lambda,
     const Index& mu) const;

    void set_jmax(const int jmax) {
      jmax_ = jmax;
      setup_full_collection();
    }

    //! get the wavelet index corresponding to a specified number
    const inline Index* get_wavelet (const int number) const {
      return &full_collection[number];
    }

    //! number of wavelets between coarsest and finest level
    const int degrees_of_freedom() const { return full_collection.size(); };

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

  protected:
    //! finest possible level
    int jmax_;

    //! setup full collectin of wavelets between j0 and jmax_ as long as a jmax_ has been specified
    void setup_full_collection();

    //! collection of all wavelets between coarsest and finest level
    Array1D<Index> full_collection;

  };

  /*!
    template specialization to flavor==P_construction

    If unsure about the suitable value of J0, choose
      J0 = SplineBasisData_j0<d,dT,P_construction,s0,s1,sT0,sT1>::j0>
  */
  template <int d, int dT, int s0, int s1, int sT0, int sT1, int J0>
  class SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1,J0>
    : public SplineBasisData<d,dT,P_construction,s0,s1,sT0,sT1>
  {
  public:
    using SplineBasisData<d,dT,P_construction,s0,s1,sT0,sT1>::check;

    //! default constructor
    SplineBasis();

    //! wavelet index class
    typedef IntervalIndex2<SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1,J0> > Index;
    
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

    /*!
      Decide whether the supports of two generators/wavelets \psi_\lambda and
      \psi_\nu have an intersection of positive measure and compute it
      in the form 2^{-j}[k1,k2]. If the return value is false, the computed
      support intersection will have no meaningful values, for performance reasons.
    */
    bool intersect_supports(const Index& lambda,
			    const Index& nu,
			    Support& supp) const;

    //! space dimension of the underlying domain
    static const int space_dimension = 1;

    //! coarsest level
    static const int j0() { return J0; }

    //! critical Sobolev regularity for the primal generators/wavelets
    static double primal_regularity() { return d - 0.5; }

    //! degree of polynomial reproduction for the primal generators/wavelets
    static unsigned int primal_polynomial_degree() { return d; }

    //! number of vanishing moments for the primal wavelets
    static unsigned int primal_vanishing_moments() { return dT; }

    //! read access to the primal b.c. order at x=0
    static const int get_s0() { return s0; }

    //! read access to the primal b.c. order at x=1
    static const int get_s1() { return s1; }

    //! read access to the primal b.c. order at x=0
    static const int get_sT0() { return sT0; }

    //! read access to the primal b.c. order at x=1
    static const int get_sT1() { return sT1; }

    //! extremal generator indices
    static const int DeltaLmin();
    static const int DeltaRmax_offset();
    static const int DeltaRmax(const int j) { return (1<<j)+DeltaRmax_offset(); }
    
    //! size of Delta_j
    static const int Deltasize(const int j) { return DeltaRmax(j)-DeltaLmin()+1; }
    
    //! boundary indices in \nabla_j
    static const int Nablamin() { return 0; }
    static const int Nablamax(const int j) { return (1<<j)-1; }
    
    //! size of Nabla_j
    static const int Nablasize(const int j) { return 1<<j; }
    
    //! index of first (leftmost) generator on level j >= j0
    static Index first_generator(const int j);

    //! index of last (rightmost) generator on level j >= j0
    static Index last_generator(const int j);

    //! index of first (leftmost) wavelet on level j >= j0
    static Index first_wavelet(const int j);

    //! index of last (rightmost) wavelet on level j >= j0
    static Index last_wavelet(const int j);
 
    /*!
      index of first function with type e
      (mainly for TensorProductBasis)
    */
    static Index first_index(const int j, const int e);

    /*!
      index of last function with type e
      (mainly for TensorProductBasis)
    */
    static Index last_index(const int j, const int e);

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
    (const Index& lambda,
     const int resolution) const;

    /*!
      Evaluate an arbitrary linear combination of primal wavelets
      on a dyadic subgrid of [0,1].
    */
    SampledMapping<1>
    evaluate
    (const InfiniteVector<double, Index>& coeffs,
     const int resolution) const;
    
    /*!
      point evaluation of (derivatives) of a single primal generator
      or wavelet \psi_\lambda
    */
    double
    evaluate
    (const unsigned int derivative,
     const Index& lambda,
     const double x) const;
    
    /*!
      point evaluation of (derivatives) of a single primal generator
      or wavelet \psi_\lambda at several points simultaneously
    */
    void
    evaluate
    (const unsigned int derivative,
     const Index& lambda,
     const Array1D<double>& points, Array1D<double>& values) const;

    /*!
      point evaluation of 0-th and first derivative of a single primal generator
      or wavelet \psi_\lambda at several points simultaneously
    */
    void
    evaluate
    (const Index& lambda,
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
     InfiniteVector<double, Index>& coeffs) const;
    
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
      helper function, integrate two primal generators or wavelets
      against each other (for the Gramian)
    */
    double
    integrate
    (const Index& lambda,
     const Index& mu) const;

    //! set maximal level
    void set_jmax(const int jmax) {
      jmax_ = jmax;
      setup_full_collection();
    }

    //! get the wavelet index corresponding to a specified number
    const inline Index* get_wavelet (const int number) const {
      return &full_collection[number];
    }

    //! number of wavelets between coarsest and finest level
    const int degrees_of_freedom() const { return full_collection.size(); };

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

  protected:
    //! finest possible level
    int jmax_;
    
    //! setup full collectin of wavelets between j0 and jmax_ as long as a jmax_ has been specified
    void setup_full_collection();
    
    //! collection of all wavelets between coarsest and finest level
    Array1D<Index> full_collection;
  };

}

#include <interval/spline_basis.cpp>

#endif
