// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2007                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_SPLINE_EVALUATE_H
#define _WAVELETTL_SPLINE_EVALUATE_H

#include <geometry/sampled_mapping.h>
#include <algebra/infinite_vector.h>
#include <interval/spline_basis.h>

using MathTL::SampledMapping;
using MathTL::InfiniteVector;

namespace WaveletTL
{
  template <int d, int dT, SplineBasisFlavor flavor> class SplineBasis;
  template <int d, int dT> class SplineBasis<d,dT,DS_construction>;

  /*!
    Evaluate a single primal/dual generator or wavelet \psi_\lambda
    on a dyadic subgrid of [0,1].
  */
  template <int d, int dT, SplineBasisFlavor flavor>
  SampledMapping<1> evaluate(const SplineBasis<d,dT,flavor>& basis,
			     const typename SplineBasis<d,dT,flavor>::Index& lambda,
			     const int resolution);

  /*!
    Evaluate an arbitrary linear combination of primal wavelets
    on a dyadic subgrid of [0,1].
  */
  template <int d, int dT, SplineBasisFlavor flavor>
  SampledMapping<1> evaluate(const SplineBasis<d,dT,flavor>& basis,
			     const InfiniteVector<double, typename SplineBasis<d,dT,flavor>::Index>& coeffs,
			     const int resolution);

  //! for compatibility also with useless "primal" parameter
  template <int d, int dT, SplineBasisFlavor flavor>
  SampledMapping<1> evaluate(const SplineBasis<d,dT,flavor>& basis,
			     const InfiniteVector<double, typename SplineBasis<d,dT,flavor>::Index>& coeffs,
			     const bool primal,
			     const int resolution)
  {
    return evaluate(basis, coeffs, resolution);
  }

  /*!
    point evaluation of (derivatives) of a single primal generator
    or wavelet \psi_\lambda
  */
  template <int d, int dT, SplineBasisFlavor flavor>
  double evaluate(const SplineBasis<d,dT,flavor>& basis, const unsigned int derivative,
		  const typename SplineBasis<d,dT,flavor>::Index& lambda,
		  const double x);

  /*!
    point evaluation of (derivatives) of a single primal generator
    or wavelet \psi_\lambda at several points simultaneously
  */
  template <int d, int dT, SplineBasisFlavor flavor>
  void evaluate(const SplineBasis<d,dT,flavor>& basis, const unsigned int derivative,
		const typename SplineBasis<d,dT,flavor>::Index& lambda,
		const Array1D<double>& points, Array1D<double>& values);

  /*!
    point evaluation of 0-th and first derivative of a single primal generator
    or wavelet \psi_\lambda at several points simultaneously
  */
  template <int d, int dT, SplineBasisFlavor flavor>
  void evaluate(const SplineBasis<d,dT,flavor>& basis,
		const typename SplineBasis<d,dT,flavor>::Index& lambda,
		const Array1D<double>& points, Array1D<double>& funcvalues, Array1D<double>& dervalues);

}

#include <interval/spline_evaluate.cpp>

#endif
