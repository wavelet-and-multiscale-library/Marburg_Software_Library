// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_FULL_LAPLACIAN_H
#define _WAVELETTL_FULL_LAPLACIAN_H

#include <iostream>
#include <map>
#include <algebra/vector.h>
#include <interval/spline_basis.h>

using namespace MathTL;

namespace WaveletTL
{
  /*!
    The following class models the discretization of the Laplacian operator
    in a primal spline wavelet basis on a level j.
    We use the decomposition

      <A Psi_j,Psi_j>^T = T_{j-1} * <A Phi_j,Phi_j>^T * T_{j-1}^T

    where T_j models the multiscale transformation within V_{j+1}
    (T_j=I for j==j0-1).
    The system is diagonally preconditioned with the inverse of
    (D)_{\lambda,\lambda}=2^{|\lambda|}
  */
  template <int d, int dT>
  class FullLaplacian
  {
  public:
    /*!
      type of indexes and size type (cf. STL containers)
     */
    typedef typename Vector<double>::size_type size_type;

    /*!
      constructor taking an information object on some spline wavelet basis
    */
    FullLaplacian(const SplineBasis<d,dT>& sb);

    /*!
      row dimension
    */
    const size_type row_dimension() const;

    /*!
      column dimension
    */
    const size_type column_dimension() const;

    /*!
      set level j (stiffness matrix in V_j)
    */
    void set_level(const int j) const;

    /*!
      read-only access to a single matrix entry
    */
    const double get_entry(const size_type row, const size_type column) const;

    /*!
      matrix-vector multiplication Mx = (*this) * x;
      we assume that the vector Mx has the correct size and
      is not identical to x
    */
    template <class VECTOR>
    void apply(const VECTOR& x, VECTOR& Mx) const;

    //! specialization to std::map<size_type,double>
    void apply(const std::map<size_type,double>& x,
	       std::map<size_type,double>& Mx) const;
    
    /*!
      stream output with user-defined tabwidth and precision
      (cf. deal.II)
    */
    void print(std::ostream& os,
	       const unsigned int tabwidth = 8,
	       const unsigned int precision = 3) const;

  protected:
    const SplineBasis<d,dT>& sb_;
    mutable int j_;
  };

  /*!
    Matlab-style stream output
  */
  template <int d, int dT>
  std::ostream& operator << (std::ostream& os, const FullLaplacian<d,dT>& M);
}

#include <galerkin/full_laplacian.cpp>

#endif
