// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_CUBE_BASIS_H
#define _WAVELETTL_CUBE_BASIS_H

#include <algebra/vector.h>
#include <algebra/infinite_vector.h>
#include <utils/fixed_array1d.h>

#include <cube/cube_index.h>

namespace WaveletTL
{
  /*!
    Template class for wavelet bases on the d-dimensional (hyper)cube [0,1]^d.
    For each of the 2*d facets, you can specify the orders s_i, sT_i
    of the Dirichlet boundary conditions of the primal and dual basis
    in the normal direction (this is tailored for the DSBasis constructor...).
  */
  template <class IBASIS, unsigned int DIM = 2>
  class CubeBasis
  {
  public:
    //! default constructor (no b.c.'s)
    CubeBasis();

    //! destructor
    ~CubeBasis();

    //! coarsest possible level j0
    inline const int j0() const { return j0_; }

    //! wavelet index class
    typedef CubeIndex<IBASIS,DIM> Index;

    //! read access to the bases
    const FixedArray1D<IBASIS*,DIM> bases() const { return bases_; }

  protected:
    //! coarsest possible level j0
    int j0_;

    //! the instances of the 1D bases
    FixedArray1D<IBASIS*,DIM> bases_;
  };
}

#include <cube/cube_basis.cpp>

#endif
