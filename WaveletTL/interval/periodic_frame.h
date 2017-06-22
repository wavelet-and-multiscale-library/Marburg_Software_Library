// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner, Philipp Keding                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_PERIODIC_FRAME_H
#define _WAVELETTL_PERIODIC_FRAME_H

#include <iostream>

#include <algebra/vector.h>
#include <algebra/infinite_vector.h>
#include <algebra/qs_matrix.h>
#include <utils/array1d.h>
#include <utils/function.h>
#include <geometry/sampled_mapping.h>

#include <interval/i_q_index.h>
#include <interval/periodic_frame_support.h>

using MathTL::Vector;
using MathTL::InfiniteVector;
using MathTL::Array1D;
using MathTL::Function;
using MathTL::SampledMapping;
using MathTL::PeriodicQuasiStationaryMatrix;

namespace WaveletTL
{
  /*!
    Template class for a periodic, quarklet frame on the unit interval [0,1],
    derived from a quarklet frame on R (which is specified as a
    template parameter). The quarklet frame on R should allow for point evaluations.

    The periodized scaling functions (and analogously the quarklets) look like

      phi^per_{p,j,k}(x) = \sum_{l\in\mathbb Z} phi_{p,j,k}(x+l)

    References:
    [D] Daubechies,
        Ten Lectures On Wavelets, pp. 304ff
  */
  template <class RFRAME>
  class PeriodicFrame
  {
  public:
    //! quarklet index class
    typedef IntervalQIndex2<PeriodicFrame<RFRAME> > Index;

    //! type of the matrices Mj0, Mj1, Mj0T, Mj1T
    typedef PeriodicQuasiStationaryMatrix<double> QuasiStationaryMatrixType;
    

    
    //! default constructor
    PeriodicFrame();
    
    
    /*!
      coarsest possible level j0;
      j0 is chosen such large that the primal and dual wavelets on level j0
      do not "overlap themselves". Note that
        |\supp\psi| = |\supp\tilde\psi| = (L+Lt)/2 - 1
      so that it suffices to choose 2^{-j0}((L+Lt)/2-1) <= 1
    */
    static const int j0() {
      return (int) ceil(log((RFRAME::primal_mask::length + RFRAME::dual_mask::length)/2-1.)/M_LN2);
    }
    
    //! size_type, for convenience
    typedef Vector<double>::size_type size_type;
    
    /*!
      geometric type of the support sets
      (note that if k1 > k2, the support consists of two intervals instead of one)
    */
    typedef struct {
      int j;
      int k1;
      int k2;
    } Support;

    /*!
      Compute an interval 2^{-j}[k1,k2] which contains the support of a
      single primal spline generator or wavelet \psi_\lambda.
      Note that due to the periodization, it may happen that k2<k1.
      (j == lambda.j()+lambda.e() is implicitly assumed for performance reasons)
    */
    static void support(const Index& lambda, int& k1, int& k2);

    /*!
      check whether the support sets of psi_lambda and psi_mu intersect
    */
    static bool intersect_supports(const Index& lambda, const Index& mu);
    
    //! space dimension of the underlying domain
    static const int space_dimension = 1;

    /*!
      critical Sobolev regularity for the primal generators/wavelets
    */
    static inline double primal_regularity() { return RFRAME::primal_regularity(); }

    /*!
      degree of polynomial reproduction for the primal generators/wavelets
    */
    static unsigned int primal_polynomial_degree() { return RFRAME::primal_polynomial_degree(); }

    /*!
      number of vanishing moments for the primal wavelets
    */
    static inline unsigned int primal_vanishing_moments() { return RFRAME::primal_vanishing_moments(); }

    //! bounds for the generator indices
    static const int DeltaLmin() { return 0; }
    static const int DeltaRmax(const int j) { return (1<<j) - 1; }

    //! bounds for the wavelet indices
    static const int Nablamin() { return 0; }
    static const int Nablamax(const int j) { return (1<<j) - 1; }
    
    //! size of Delta_j
    static const int Deltasize(const int j) { return 1<<j; }
    
    //! size of Nabla_j
    static const int Nablasize(const int j) { return 1<<j; }

    //! index of first (leftmost) generator on level j >= j0
    static Index first_generator(const int j, const int p = 0);

    //! index of last (rightmost) generator on level j >= j0
    static Index last_generator(const int j, const int p = 0);

    //! index of first (leftmost) wavelet on level j >= j0
    static Index first_wavelet(const int j, const int p = 0);

    //! index of last (rightmost) wavelet on level j >= j0
    static Index last_wavelet(const int j, const int p = 0);

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
    
    
    // set the indices jmax_ and pmax_
    
    void set_jpmax(const int jmax, const int pmax){
       
       jmax_=jmax;
       pmax_=pmax;       
       setup_full_collection();
    }
    
    int get_jmax_() const{
        return jmax_;
    }
    
    int get_pmax_() const{
        return pmax_;
    }
    
    
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
      an analogous routine for Mj1
    */
    template <class V>
    void apply_Mj1(const int j, const V& x, V& y,
		   const size_type x_offset, const size_type y_offset,
		   const bool add_to) const;
    
    /*!
      apply Mj=(Mj0 Mj1) to some vector x ("reconstruct");
      the routine writes only into the first part of y, i.e,
      y might be larger than necessary, which is helpful for apply_Tj
    */
    template <class V>
    void apply_Mj(const int j, const V& x, V& y) const;
    
    
    //! apply Tj=Mj*diag(M_{j-1},I)*...*diag(M_{j_0},I), i.e., several "reconstructions" at once
    template <class V>
    void apply_Tj(const int j, const V& x, V& y) const;

/*!
      Evaluate a single primal/dual generator or wavelet \psi_\lambda
      on a dyadic subgrid of [0,1].
    */
    SampledMapping<1>
    evaluate
    (const Index& lambda,
     const int resolution, const bool normalization = 1) const;

    /*!
      Evaluate an arbitrary linear combination of primal wavelets
      on a dyadic subgrid of [0,1].
    */
    SampledMapping<1>
    evaluate
    (const InfiniteVector<double, Index>& coeffs,
     const int resolution,
     const int derivative = 0, const bool normalization = 1) const;
    
    /*!
      point evaluation of (derivatives) of a single primal generator
      or wavelet \psi_\lambda
    */
    double evaluate(const unsigned int derivative,
		    const Index& lambda,
		    const double x, const bool normalization = 1) const;
    
    /*!
      point evaluation of (derivatives) of a single primal generator
      or wavelet \psi_\lambda at several points simultaneously
    */
    void evaluate(const unsigned int derivative,
		  const Index& lambda,
		  const Array1D<double>& points,
		  Array1D<double>& values, const bool normalization = 1) const;
    
    /*!
      point evaluation of 0-th and first derivative of a single primal generator
      or wavelet \psi_\lambda at several points simultaneously
    */
    void evaluate(const Index& lambda,
		  const Array1D<double>& points,
		  Array1D<double>& funcvalues,
		  Array1D<double>& dervalues, const bool normalization = 0) const;

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
     const int pmax,
     InfiniteVector<double, Index>& coeffs, const bool normalization = 0) const;

    /*!
      analogous routine for Vector<double> output
    */
    /*void
    expand
    (const Function<1>* f,
     const bool primal,
     const int jmax,
     Vector<double>& coeffs) const;*/

    /*!
      helper function, integrate a smooth function f against a
      periodic primal generator or wavelet
    */
    double
    integrate
    (const Function<1>* f,
     const Index& lambda, const bool& normalization = 0) const;
    
    /*
      helper function, integrate two primal generators or wavelets
      against each other (for the Gramian or Laplacian) /PROBLEMS WITH NORMALIZATION =1 AND DERIVATIVE=0 @PHK look at routines in cachedquarkletproblem f(lambda)
      there are also problems. everything need to be fixed for the gramian.
      bei zwei generatoren muss beachtet werden normalissierter wert= alter wert-integral(phi_lambda)*integral(phi_mu) 
     still errors: if two generators, we need to subtract 1 time the multiple of the integrals  */
    double
    integrate
    (const unsigned int& derivative,
     const Index& lambda,
     const Index& mu, const bool& normalization = 1) const;
    
    const int degrees_of_freedom() const { return full_collection.size(); };
    
    //! get the wavelet index corresponding to a specified number
    const inline Index* get_quarklet (const int number) const {
      return &full_collection[number];
    }

    
  private:
    //! an instance of RFRAME
    RFRAME r_frame;
    
  protected:
      
    //! coarsest possible level
    int j0_;

    //! finest possible level
    int jmax_;
    
    //! biggest polynomial degree
    int pmax_;
    
     //! setup full collectin of wavelets between j0_ and jmax_ as long as a jmax_ has been specified
    void setup_full_collection();

    //! collection of all wavelets between coarsest and finest level
    Array1D<Index> full_collection;
    
  public:
    //! refinement matrices
    QuasiStationaryMatrixType Mj0_, Mj1_, Mj0T_, Mj1T_;
        


  
  };
}

// include implementation
#include <interval/periodic_frame.cpp>

#endif
