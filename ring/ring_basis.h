// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2008                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_RING_BASIS_H
#define _WAVELETTL_RING_BASIS_H

#include <iostream>

#include <algebra/vector.h>
#include <algebra/infinite_vector.h>
#include <utils/array1d.h>
#include <utils/function.h>
#include <geometry/sampled_mapping.h>
#include <geometry/ring_chart.h>

#include <Rd/cdf_basis.h>
#include <interval/periodic.h>
#include <interval/spline_basis.h>
#include <ring/ring_index.h>

using MathTL::Vector;
using MathTL::InfiniteVector;
using MathTL::Array1D;
using MathTL::Function;

namespace WaveletTL
{
  /*!
    Template class for a wavelet basis on the ring-shaped domain
      R = {(x,y) : r_0 <= ||(x,y)|| <= r_1 }
    with complementary boundary conditions of order s0 at the interior boundary
    and of order s1 at the outer boundary.

    The basis functions psi_lambda are constructed in such a way that we use the chart
      kappa: (0,1)^2 -> R,
      kappa(s,phi) = r(s)(cos(2*pi*phi),sin(2*pi*phi))
    with
      r(s) = r_0+s*(r_1-r_0).
    The wavelets are then defined as
      psi_lambda(x,y) = psi^0_lambda(kappa^{-1}(x,y))/|det Dkappa(kappa^{-1}(x,y))|^{1/2},
    i.e.,
      psi_lambda(kappa(s,phi)) = psi^0_lambda(s,phi)/(r(s))^{1/2}.
  */
  template <int d, int dt, int s0, int s1>
  class RingBasis
  {
  public:
    //! size_type, for convenience
    typedef Vector<double>::size_type size_type;

    //! wavelet index class
    typedef RingIndex<d,dt,s0,s1> Index;

    //! type of 1D basis in angular direction
    typedef PeriodicBasis<CDFBasis<d,dt> > Basis0;
    
    //! type of 1D indices in angular direction
    typedef typename Basis0::Index Index0;
    
    //! type of 1D basis in radial direction
    typedef SplineBasis<d,dt,P_construction,s0,s1,0,0,
			SplineBasisData_j0<d,dt,P_construction,s0,s1,0,0>::j0> Basis1;

    //! type of 1D indices in radial direction
    typedef typename Basis1::Index Index1;

    //! default constructor, also sets r0 and r1
    RingBasis(const double r0 = 0.5, const double r1 = 2.0);
    
    //! critical Sobolev regularity for the primal generators/wavelets
    static double primal_regularity() { return Basis0::primal_regularity(); }
    
    //! degree of polynomial reproduction for the primal generators/wavelets
    static unsigned int primal_polynomial_degree() { return d; }

    //! number of vanishing moments for the primal wavelets
    static unsigned int primal_vanishing_moments() { return dt; }

    /*!
      coarsest possible level j0; we may assume that the interval basis in radial direction
      has a larger minimal level than the corresponding periodic basis
    */
    static const int j0() {
      return Basis1::j0();
    }

    //! index of first generator on level j >= j0
    static Index first_generator(const int j);

    //! index of last generator on level j >= j0
    static Index last_generator(const int j);

    //! index of first wavelet on level j >= j0
    static Index first_wavelet(const int j);

    //! index of last wavelet on level j >= j0
    static Index last_wavelet(const int j);

    //! size of Delta_j
    static int Deltasize(const int j);

    //! sizes of the different wavelet index sets
    static int Nabla01size(const int j);
    static int Nabla10size(const int j);
    static int Nabla11size(const int j);

    //! the chart
    const RingChart& chart() const { return chart_; }

    //! the radii
    const double r0() const { return r0_; }
    const double r1() const { return r1_; }

    //! the 1D bases
    const Basis0& basis0() const { return basis0_; }
    const Basis1& basis1() const { return basis1_; }

    /*!
      Evaluate a single primal/dual generator or wavelet \psi_\lambda
      on a "dyadic" subgrid of the ring, pulled back to the unit square
    */
    SampledMapping<2> evaluate(const Index& lambda,
			       const int resolution) const;
    
    /*!
      Evaluate an arbitrary linear combination of primal/dual wavelets
    */
    SampledMapping<2> evaluate(const InfiniteVector<double,Index>& coeffs,
 			       const int resolution) const;
    
    /*!
      apply Mj0 to some vector x (partial "reconstruct")
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
      apply Mj0T^T to some vector x (partial "decompose")
    */
    template <class V>
    void apply_Mj0T_transposed(const int j, const V& x, V& y,
			       const size_type x_offset, const size_type y_offset,
			       const bool add_to) const;

    /*!
      an analogous routine for Mj1T^T, e=(0,1)
    */
    template <class V>
    void apply_Mj1T_01_transposed(const int j, const V& x, V& y,
				  const size_type x_offset, const size_type y_offset,
				  const bool add_to) const;

    /*!
      an analogous routine for Mj1T^T, e=(1,0)
    */
    template <class V>
    void apply_Mj1T_10_transposed(const int j, const V& x, V& y,
				  const size_type x_offset, const size_type y_offset,
				  const bool add_to) const;

    /*!
      an analogous routine for Mj1T^T, e=(1,1)
    */
    template <class V>
    void apply_Mj1T_11_transposed(const int j, const V& x, V& y,
				  const size_type x_offset, const size_type y_offset,
				  const bool add_to) const;

    /*!
      apply Mj=(Mj0 Mj1) to some vector x ("reconstruct");
      the routine writes only into the first part of y, i.e,
      y might be larger than necessary, which is helpful for apply_Tj
    */
    template <class V>
    void apply_Mj(const int j, const V& x, V& y) const;

    //! apply Gj=(Mj0T Mj1T)^T to some vector x ("decompose")
    template <class V>
    void apply_Gj(const int j, const V& x, V& y) const;

    //! apply Tj=Mj*diag(M_{j-1},I)*...*diag(M_{j_0},I), i.e., several "reconstructions" at once
    template <class V>
    void apply_Tj(const int j, const V& x, V& y) const;

    //! apply Tj^{-1}, several "decompositions" at once
    template <class V>
    void apply_Tjinv(const int j, const V& x, V& y) const;

    /*!
      For a given function, compute all integrals w.r.t. the primal
      or dual generators/wavelets \psi_\lambda with |\lambda|\le jmax.
    */
    void
    expand
    (const Function<2>* f,
     const bool primal,
     const int jmax,
     InfiniteVector<double, Index>& coeffs) const;

    /*!
      For a given function, compute all integrals w.r.t. the primal
      or dual generators/wavelets \psi_\lambda with |\lambda|\le jmax.
    */
    void
    expand
    (const Function<2>* f,
     const bool primal,
     const int jmax,
     Vector<double>& coeffs) const;

    /*!
      helper function, integrate a smooth function f against a
      primal generator or wavelet
    */
    double
    integrate
    (const Function<2>* f,
     const Index& lambda) const;

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

  protected:
    //! an instance of the periodic basis (angular direction)
    Basis0 basis0_;

    //! an instance of the nonperiodic basis (radial direction)
    Basis1 basis1_;

    //! radii
    double r0_, r1_;

    //! an instance of the corresponding chart
    RingChart chart_;

    //! finest possible level
    int jmax_;
    
    //! setup full collectin of wavelets between j0 and jmax_ as long as a jmax_ has been specified
    void setup_full_collection();
    
    //! collection of all wavelets between coarsest and finest level
    Array1D<Index> full_collection;
  };
}

// include implementation
#include <ring/ring_basis.cpp>

#endif
