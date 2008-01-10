// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2008                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_SPLINE_BASIS_DATA_H
#define _WAVELETTL_SPLINE_BASIS_DATA_H

#include <string.h>
#include <algebra/matrix.h>
#include <algebra/qs_matrix.h>
#include <interval/spline_basis_flavor.h>

using namespace MathTL;

namespace WaveletTL
{
  /*!
    A class to hold all the pieces of information
    needed for the construction of a spline wavelet basis on the interval.
  */
  template <int d, int dT, SplineBasisFlavor flavor, int s0, int s1, int sT0, int sT1>
  class SplineBasisData
  {
  public:
    //! default constructor
    SplineBasisData();

    //! destructor
    virtual ~SplineBasisData();

    //! coarsest level
    static const int j0();

    //! check integrity of the internal data
    void check() const;

  protected:
    //! refinement matrices
    QuasiStationaryMatrix<double> Mj0_, Mj1_, Mj0T_, Mj1T_;

    //! expansion coefficients of the primal generators w.r.t. restricted cardinal B-splines
    Matrix<double> CLA_, CRA_, CLAT_, CRAT_;

    /*!
      initial stable completion with homogeneous b.c.'s (for composite basis),
      available for bio5(-energy)

      the first and last row is zero, we don't store it
    */
    QuasiStationaryMatrix<double> Mj1c_;
  };

  /*!
    template specialization to SplineBasisFlavor==P_construction,
    here we do not have to store neither the CLA* matrices nor Mj1c
  */
  template <int d, int dT, int s0, int s1, int sT0, int sT1>
  class SplineBasisData<d,dT,P_construction,s0,s1,sT0,sT1>
  {
  public:
    //! default constructor
    SplineBasisData();
    
    //! destructor
    virtual ~SplineBasisData();
    
    //! coarsest level
    static const int j0();

    //! check integrity of the internal data
    void check() const;
 
  protected:
    //! refinement matrices
    QuasiStationaryMatrix<double> Mj0_, Mj1_, Mj0T_, Mj1T_;
  };

}

#include <interval/spline_basis_data.cpp>

#endif
