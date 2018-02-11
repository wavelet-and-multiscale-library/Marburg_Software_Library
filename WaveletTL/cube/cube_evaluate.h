// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_CUBE_EVALUATE_H
#define _WAVELETTL_CUBE_EVALUATE_H

#include <geometry/sampled_mapping.h>
#include <algebra/infinite_vector.h>
#include <cube/cube_basis.h>
#include <geometry/point.h>

using MathTL::SampledMapping;
using MathTL::InfiniteVector;

extern double time_consumption_of_calc_deriv;
extern double time_consumption_of_calc_deriv_outer;

namespace WaveletTL
{
  template <class IBASIS, unsigned int DIM> class CubeBasis;

  /*!
    Evaluate a single primal/dual generator or wavelet \psi_\lambda
    on a dyadic subgrid of [0,1]^d.
  */
  template <class IBASIS, unsigned int DIM>
  SampledMapping<DIM> evaluate(const CubeBasis<IBASIS,DIM>& basis,
			       const typename CubeBasis<IBASIS,DIM>::Index& lambda,
			       const bool primal,
			       const int resolution);

  /*!
    Evaluate an arbitrary linear combination of primal/dual wavelets
    on a dyadic subgrid of [0,1]^d.
  */
  template <class IBASIS, unsigned int DIM>
  SampledMapping<DIM> evaluate(const CubeBasis<IBASIS,DIM>& basis,
			       const InfiniteVector<double, typename CubeBasis<IBASIS,DIM>::Index>& coeffs,
			       const bool primal,
			       const int resolution);

  
#ifdef P_POISSON
  /*!
    Evaluate a single primal generator or wavelet \psi_\lambda
    in a single point in [0,1]^d.
  */
  template <class IBASIS, unsigned int DIM>
  double evaluate(const CubeBasis<IBASIS,DIM>& basis,
	              const typename CubeBasis<IBASIS,DIM>::Index& lambda,
	              const Point<DIM> p);

  /*!
    Evaluate an arbitrary linear combination of primal wavelets
    in a single point in [0,1]^d.
  */
  template <class IBASIS, unsigned int DIM>
  double evaluate(const CubeBasis<IBASIS,DIM>& basis,
                  const InfiniteVector<double, typename CubeBasis<IBASIS,DIM>::Index>& coeffs,
	              const Point<DIM> p);

  /*!
    Evaluate the derivative of a single primal generator or wavelet \psi_\lambda
    in a single point in [0,1]^d.
  */
  template <class IBASIS, unsigned int DIM>
  void evaluate(const CubeBasis<IBASIS,DIM>& basis,
	              const typename CubeBasis<IBASIS,DIM>::Index& lambda,
	              const Point<DIM> p,
	              FixedArray1D<double,DIM>& deriv_values);

  /*!
    Evaluate the derivative of an arbitrary linear combination of primal wavelets
    in a single point in [0,1]^d.
  */
  template <class IBASIS, unsigned int DIM>
  void evaluate(const CubeBasis<IBASIS,DIM>& basis,
                  const InfiniteVector<double, typename CubeBasis<IBASIS,DIM>::Index>& coeffs,
	              const Point<DIM> p,
	              FixedArray1D<double,DIM>& deriv_values);


  /*! C.Hartmann: added on 25.09.15
    Evaluate a single primal generator or wavelet \psi_\lambda
    on an arbitrary grid of [0,1]^d. The grid is given by a tensor product
    of 1-dimensional grids on [0,1].
  */
  template <class IBASIS, unsigned int DIM>
  Matrix<double> evaluate(const CubeBasis<IBASIS,DIM>& basis,
                          const typename CubeBasis<IBASIS,DIM>::Index& lambda,
                          const FixedArray1D<Array1D<double>,DIM>& grid);


  /*! C.Hartmann: added on 25.09.15
    Evaluate an arbitrary linear combination of primal wavelets
    on an arbitrary grid of [0,1]^d. The grid is given by a tensor product
    of 1-dimensional grids on [0,1].
  */
  template <class IBASIS, unsigned int DIM>
  Matrix<double> evaluate(const CubeBasis<IBASIS,DIM>& basis,
                          const InfiniteVector<double, typename CubeBasis<IBASIS,DIM>::Index>& coeffs,
                          const FixedArray1D<Array1D<double>,DIM>& grid);

  /*! C.Hartmann: added on 26.09.15
    Evaluate the first derivative of a single primal generator or wavelet \psi_\lambda
    on an arbitrary grid of [0,1]^d. The grid is given by a tensor product of 1-dimensional grids on [0,1].
    Function values are stored in the matrices deriv_x and deriv_y. Matrix dimensions must match grid dimension!
    ATTENTION: Only dimension = 2 and evaluation of primal generators/wavelets supported a.t.m.
  */
  template <class IBASIS, unsigned int DIM>
  void
  evaluate(const CubeBasis<IBASIS,DIM>& basis,
           const typename CubeBasis<IBASIS,DIM>::Index& lambda,
           const FixedArray1D<Array1D<double>,DIM>& grid,
           Matrix<double>& deriv_x, Matrix<double>& deriv_y);


  /*! C.Hartmann: added on 26.09.15
    Evaluate the first derivative of an arbitrary linear combination of primal wavelets
    on an arbitrary grid of [0,1]^d. The grid is given by a tensor product of 1-dimensional grids on [0,1].
    Function values are stored in the matrices deriv_x and deriv_y. Matrix dimensions must match grid dimension!
    ATTENTION: Only dimension = 2 and evaluation of primal generators/wavelets supported a.t.m.
  */
  template <class IBASIS, unsigned int DIM>
  void
  evaluate(const CubeBasis<IBASIS,DIM>& basis,
           const InfiniteVector<double, typename CubeBasis<IBASIS,DIM>::Index>& coeffs,
           const FixedArray1D<Array1D<double>,DIM>& grid,
           Matrix<double>& deriv_x, Matrix<double>& deriv_y);


//! tuned Versions (Christoph)

  template <class IBASIS, unsigned int DIM>
  void evaluate_tuned(const CubeBasis<IBASIS,DIM>& basis,
	              const typename CubeBasis<IBASIS,DIM>::Index& lambda,
	              const Point<DIM> p,
	              FixedArray1D<double,DIM>& deriv_values);

  template <class IBASIS, unsigned int DIM>
  void evaluate_tuned(const CubeBasis<IBASIS,DIM>& basis,
                  const InfiniteVector<double, typename CubeBasis<IBASIS,DIM>::Index>& coeffs,
	              const Point<DIM> p,
	              FixedArray1D<double,DIM>& deriv_values);

  //! (only for dimension=2)
  template <class IBASIS>
  void evaluate_tuned_dim2(const CubeBasis<IBASIS,2>& basis,
                  const InfiniteVector<double, typename CubeBasis<IBASIS,2>::Index>& coeffs,
	              const Point<2> p,
	              FixedArray1D<double,2>& deriv_values);

  //! (only for dimension=2, linear (2,2) Primbs basis with homog. bc's)
  template <class IBASIS>
  void evaluate_tuned_dim2_lin_primbs(const CubeBasis<IBASIS,2>& basis,
                  const InfiniteVector<double, typename CubeBasis<IBASIS,2>::Index>& coeffs,
	              const Point<2> p,
	              FixedArray1D<double,2>& deriv_values);

  //! (only for dimension=2, linear (2,2) Primbs basis with homog. bc's)
  template <class IBASIS, unsigned int DIM>
  double evaluate_tuned_dim2_lin_primbs(const CubeBasis<IBASIS,DIM>& basis,
                  const InfiniteVector<double, typename CubeBasis<IBASIS,DIM>::Index>& coeffs,
	              const Point<DIM> p);
#endif

}

#include <cube/cube_evaluate.cpp>

#endif
