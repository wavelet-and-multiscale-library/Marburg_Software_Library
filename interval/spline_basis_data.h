// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2007                                            |
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
    The following structure holds all the pieces of information
    needed for the construction of a spline wavelet basis on the interval.
    The primal functions are either B-splines or linear combinations thereof.
   */
  template <int d, int dT, SplineBasisFlavor flavor>
  class SplineBasisData
  {
  public:
    //! constructor from the primal and dual boundary condition orders
    SplineBasisData(const char* options,
		    const int s0, const int s1, const int sT0, const int sT1);

    //! destructor
    virtual ~SplineBasisData();

    //! check integrity of the internal data
    void check() const;

  protected:
    //! options
    std::string options_;

    //! boundary condition orders
    const int s0_, s1_, sT0_, sT1_;

    //! coarsest level
    int j0_;
    
    //! refinement matrices
    QuasiStationaryMatrix<double> Mj0_, Mj1_, Mj0T_, Mj1T_;
  };

  /*!
    template specialization to SplineBasisFlavor==DS_construction,
    here we have to store further matrices for the spline expansion of the generators
  */
  template <int d, int dT>
  class SplineBasisData<d,dT,DS_construction>
  {
  public:
    /*!
      constructor from the primal and dual boundary condition orders;
      the variable "options" should be set to the appropriate
      biorthogonalization method:
        "bio5"        <-> BernsteinSVD
        "bio5-energy" <-> BernsteinSVD + energy inner product orthogonalization from [Ba01]
    */
    SplineBasisData(const char* options,
		    const int s0, const int s1, const int sT0, const int sT1);
 
    //! destructor
    virtual ~SplineBasisData();
 
    //! check integrity of the internal data
    void check() const;

  protected:
    //! options
    std::string options_;

    //! boundary condition orders
    const int s0_, s1_, sT0_, sT1_;

    //! coarsest level
    int j0_;
 
    //! refinement matrices
    QuasiStationaryMatrix<double> Mj0_, Mj1_, Mj0T_, Mj1T_;

    //! expansion coefficients of the primal generators w.r.t. restricted cardinal B-splines
    Matrix<double> CLA_, CRA_, CLAT_, CRAT_;

    //! initial stable completion with homogeneous b.c.'s (for composite basis), available for bio5(-energy)
    QuasiStationaryMatrix<double> Mj1c_;
  };

}

#include <interval/spline_basis_data.cpp>

#endif
