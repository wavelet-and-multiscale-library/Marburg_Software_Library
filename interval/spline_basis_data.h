// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_SPLINE_BASIS_DATA_H
#define _WAVELETTL_SPLINE_BASIS_DATA_H

#include <string.h>
#include <algebra/qs_matrix.h>

using namespace MathTL;

namespace WaveletTL
{
  /*!
    The following structure holds all the pieces of information
    needed for the construction of a spline wavelet basis on the interval.

   */
  template <int d, int dT>
  class SplineBasisData
  {
  public:
    /*!
      constructor from the primal and dual boundary condition orders
      plus additional information
    */
    SplineBasisData(const char* flavor,
		    const char* options,
		    const int s0, const int s1, const int sT0, const int sT1);

    //! destructor
    virtual ~SplineBasisData();

    //! check integrity of the internal data
    void check() const;

  protected:
    // flavor (e.g., "DS", "P" etc.), determines also the shape of the primal generators
    std::string flavor_;

    // boundary condition orders
    const int s0_, s1_, sT0_, sT1_;

    // coarsest level
    int j0_;
    
    // refinement matrices (without prefactor 1/sqrt(2)!)
    QuasiStationaryMatrix<double> *Mj0_, *Mj1_, *Mj0T_, *Mj1T_;
  };
}

#include <interval/spline_basis_data.cpp>

#endif
