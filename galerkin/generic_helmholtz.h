// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2008                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_GENERIC_HELMHOLTZ_H
#define _WAVELETTL_GENERIC_HELMHOLTZ_H

#include <iostream>
#include <algebra/vector.h>
#include <galerkin/full_laplacian.h>

using namespace MathTL;

namespace WaveletTL
{
  /*!
    Generic class for the Helmholtz equation in a given wavelet basis,
    where the matrix components (Gramian+Laplacian) are read from a file.
  */
  template <class WBASIS>
  class GenericFullHelmholtz
  {
  public:
    /*!
      type of indexes and size type (cf. STL containers)
     */
    typedef typename Vector<double>::size_type size_type;

    /*!
      constructor taking an information object on the underlying wavelet
      basis and the filenames to read the matrix components from.
    */
    GenericFullHelmholtz(const WBASIS& basis,
			 const char* filenameG,
			 const char* filenameA,
			 const int jmax,
			 const double alpha,
			 const PreconditioningType precond = dyadic);

    /*!
      row dimension
    */
    const size_type row_dimension() const;

    /*!
      column dimension
    */
    const size_type column_dimension() const;

    /*!
      evaluate the diagonal preconditioner D
    */
    double D(const size_type k) const;

    /*!
      read-only access to a single matrix entry
    */
    const double get_entry(const size_type row, const size_type column) const;

    /*!
      read-only access to a diagonal entry of the unpreconditioned system
    */
    const double diagonal(const size_type row) const;

    /*!
      matrix-vector multiplication Mx = (*this) * x;
      we assume that the vector Mx has the correct size and
      is not identical to x
    */
    template <class VECTOR>
    void apply(const VECTOR& x, VECTOR& Mx,
	       const bool preconditioning = true) const;
    
    /*!
      set reaction coefficient alpha
    */
    void set_alpha(const double alpha) const;

  protected:
    const WBASIS& basis_;
    mutable double alpha_;
    const int jmax_;
    PreconditioningType precond_;
    SparseMatrix<double> G_, A_;
    mutable Vector<double> D_;
    void setup_D() const;
  };
}

#include <galerkin/generic_helmholtz.cpp>

#endif
