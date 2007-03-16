// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2007                                            |
// | Thorsten Raasch, Manuel Werner, Andreas Schneider                  |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_A_BASIS_H
#define _WAVELETTL_A_BASIS_H

#include <algebra/vector.h>
#include <algebra/matrix.h>
#include <algebra/infinite_vector.h>
#include <algebra/piecewise.h>
#include <utils/array1d.h>
#include <interval/i_multi_index.h>

using MathTL::Vector;
using MathTL::Matrix;
using MathTL::InfiniteVector;
using MathTL::Array1D;
using MathTL::Polynomial;
using MathTL::Piecewise;

namespace WaveletTL
{
  /*!
    Template class for the (multi)wavelet bases on the interval as introduced in [A].

    The template parameter n corresponds to the order of the splines used
    in the basis.

    References:
    [A] Alpert:
        A Class of Bases in $L^2$ for the Sparse Representation of Integral Operators
  */
  template <unsigned int n>
  class ABasis
  {
  public:
    /*!
      constructor
    */
    ABasis();

    //! coarsest possible level
    inline const int j0() const { return j0_; }

    //! wavelet index class
    typedef IntervalMultiIndex< ABasis<n> > Index;

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
      single (primal) generator or wavelet with index lambda.
      Note that wavelets have granularity j = lambda.j()+1, so the
      returned k values are respective to lambda.j()+lambda.e().
    */
    void primal_support(const Index& lambda, int& k1, int& k2) const;

    //! space dimension of the underlying domain
    static const int space_dimension = 1;

    //! number of compononents of the multi-wavelet basis
    static const unsigned int number_of_components = n;

    //! critical Sobolev regularity for the (primal) generators/wavelets
    static double primal_regularity() { return 0.5; }

    //! degree of polynomial reproduction for the primal generators/wavelets
    static unsigned int primal_polynomial_degree() { return n-1; }

    //! number of vanishing moments for the primal wavelets
    static unsigned int primal_vanishing_moments() { return n; }

    //! size of Delta_j
    inline const int Deltasize(const int j) const { return 1<<j; }

    //! size of Nabla_j
    inline const int Nablasize(const int j) const { return 1<<j; }
    
    //! extremal generator indices in Delta_j
    inline const int DeltaLmin() const { return 0; }
    inline const int DeltaRmax(const int j) const { return (1<<j)-1; }
    
    //! extremal wavelet indices in Nabla_j
    inline const int Nablamin() const { return 0; }
    inline const int Nablamax(const int j) const { return (1<<j)-1; }

    //! index of first (leftmost) generator on level j >= j0
    Index first_generator(const int j) const;

    //! index of last (rightmost) generator on level j >= j0
    Index last_generator(const int j) const;

    //! index of first (leftmost) wavelet on level j >= j0
    Index first_wavelet(const int j) const;

    //! index of last (rightmost) wavelet on level j >= j0
    Index last_wavelet(const int j) const;

    //! return (derivatives of) a wavelet or generator function with given index
    Piecewise<double> get_function(const Index index, const unsigned int derivative = 0) const;

    /*!
      point evaluation of (derivatives) of a single generator or wavelet
      \psi_\lambda
     */
    double primal_evaluate(const unsigned int derivative, const Index& lambda, const double x) const;

    /*!
      point evaluation of (derivatives) of a single generator or wavelet
      \psi_\lambda at several points simultaneously
    */
    void primal_evaluate(const unsigned int derivative, const Index& lambda,
                         const Array1D<double>& points, Array1D<double>& values) const;

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
    //! coarsest possible level
    static const int j0_ = 0;

    //! The generator functions
    Array1D<Piecewise<double> > generators;

    //! The wavelet functions
    Array1D<Piecewise<double> > wavelets;

    //! general setup routine which is shared by the different constructors
    void setup();
  };
}

#include <interval/a_basis.cpp>

#endif  // _WAVELETTL_A_BASIS_H
