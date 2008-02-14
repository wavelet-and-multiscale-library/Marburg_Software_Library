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
    A helper class to provide the coarsest possible level j0
    of a spline wavelet basis on the interval
  */
  template <int d, int dT, SplineBasisFlavor flavor, int s0, int s1, int sT0, int sT1>
  struct SplineBasisData_j0
  {
  public:
    static const int j0;
  };
  
  /*!
    A helper class to hold all the pieces of information (apart from j0)
    needed for the construction of a spline wavelet basis on the interval.

  */
  template <int d, int dT, SplineBasisFlavor flavor, int s0, int s1, int sT0, int sT1>
  class SplineBasisData
  {
  public:
    /*!
      default constructor
      (yields matrix for j0==SplineBasisData_j0<d,dT,flavor,s0,s1,sT0,sT1>)
    */
    SplineBasisData();

    //! destructor
    virtual ~SplineBasisData();

    //! check integrity of the internal data
    void check() const;

    //! refinement matrices
    QuasiStationaryMatrix<double> Mj0_, Mj1_, Mj0T_, Mj1T_;

  protected:
    /*!
      The interior wavelets are CDF wavelets, multiplied by this factor, see also
      [Pr06, Abb. 5.1 & (5.12)].
      Note that the interior wavelets from the right half of the interval are
      also reflected at 0.5 whenever d%2 != 0
    */
    int CDF_factor;
    
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
    /*!
      default constructor
      (yields matrix for j0==SplineBasisData_j0<d,dT,flavor,s0,s1,sT0,sT1>)
    */
    SplineBasisData();
    
    //! destructor
    virtual ~SplineBasisData();
    
    //! check integrity of the internal data
    void check() const;
 
    //! refinement matrices
    QuasiStationaryMatrix<double> Mj0_, Mj1_, Mj0T_, Mj1T_;

  protected:
    /*!
      The interior wavelets are CDF wavelets, multiplied by the following factor.
      Note that the interior wavelets from the right half of the interval are
      also reflected at 0.5 whenever d%2 != 0
    */
    double CDF_factor;    
  };

}

#include <interval/spline_basis_data.cpp>

#endif
