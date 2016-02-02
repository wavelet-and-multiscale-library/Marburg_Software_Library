// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_PQ_EVALUATE_H
#define _WAVELETTL_PQ_EVALUATE_H

#include <geometry/sampled_mapping.h>
#include <algebra/infinite_vector.h>
#include <interval/pq_frame.h>

using MathTL::SampledMapping;
using MathTL::InfiniteVector;

namespace WaveletTL
{
  template <int d, int dT> class PQFrame;

  /*!
    Evaluate a single primal/dual generator or wavelet \psi_\lambda
    on a dyadic subgrid of [0,1].
  */
  template <int d, int dT>
  SampledMapping<1> evaluate(const PQFrame<d,dT>& basis,
			     const typename PQFrame<d,dT>::Index& lambda,
			     const bool primal,
			     const int resolution);
  
  /*!
    Evaluate an arbitrary linear combination of primal or dual
    wavelets on a dyadic subgrid of [0,1].
  */
  template <int d, int dT>
  SampledMapping<1> evaluate(const PQFrame<d,dT>& basis,
			     const InfiniteVector<double, typename PQFrame<d,dT>::Index>& coeffs,
			     const bool primal,
			     const int resolution);

  /*!
    point evaluation of (derivatives) of a single primal [P] generator
    or wavelet \psi_\lambda
    The Index version is slow. Usage is disadvised!
  */
  template <int d, int dT>
  double evaluate(const PQFrame<d,dT>& basis, const unsigned int derivative,
		  const typename PQFrame<d,dT>::Index& lambda,
		  const double x);

  template <int d, int dT>
  double evaluate(const PQFrame<d,dT>& basis, const unsigned int derivative,
		  const int p, const int j, const int e, const int k,
		  const double x);

  /*!
    Expand of a single primal [P] generator or wavelet \psi_\lambda as a PP Funktion
   * not changed to quarklet-setting. Use is not clear @PHK
  */
  template <int d, int dT>
  Piecewise<double> expandAsPP(const PQFrame<d,dT>& basis, const typename PQFrame<d,dT>::Index& lambda);


  /*!
    point evaluation of (derivatives) of a single primal [P] generator
    or wavelet \psi_\lambda at several points simultaneously
   * 
   * Function using Index is slower than p,j,e,k version. Usage is disadvised
  */
  template <int d, int dT>
  void evaluate(const PQFrame<d,dT>& basis, const unsigned int derivative,
		const typename PQFrame<d,dT>::Index& lambda,
		const Array1D<double>& points, Array1D<double>& values);
  template <int d, int dT>
  void evaluate(const PQFrame<d,dT>& basis, const unsigned int derivative,
		const int p, const int j, const int e, const int k,
		const Array1D<double>& points, Array1D<double>& values);

  /*!
    point evaluation of 0-th and first derivative of a single primal [P] generator
    or wavelet \psi_\lambda at several points simultaneously
  */
  template <int d, int dT>
  void evaluate(const PQFrame<d,dT>& basis,
		const typename PQFrame<d,dT>::Index& lambda,
		const Array1D<double>& points, Array1D<double>& funcvalues, Array1D<double>& dervalues);
  
  //not implemented @PHK
  template <int d, int dT>
  void evaluate2(const PQFrame<d,dT>& basis,
		const typename PQFrame<d,dT>::Index& lambda,
		const Array1D<double>& points, Array1D<double>& funcvalues, Array1D<double>& dervalues);
  
  template <int d, int dT>
  void evaluate(const PQFrame<d,dT>& basis,
		const int p, const int j, const int e, const int k,
		const Array1D<double>& points, Array1D<double>& funcvalues, Array1D<double>& dervalues);


}

#include <interval/pq_evaluate.cpp>

#endif
