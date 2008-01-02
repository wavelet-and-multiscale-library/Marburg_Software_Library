// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2008                                            |
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
    three diagonal preconditioning possibilities
   */
  enum PreconditioningType {
    no_precond,
    dyadic, // 2^{-js}
    energy  // sqrt(a(\psi_\lambda,\psi_\lambda))
  };

  /*!
    The following class models the discretization of the Laplacian operator
    in a primal spline wavelet basis on a level j.
    We use the decomposition

      <A Psi_j,Psi_j>^T = T_{j-1} * <A Phi_j,Phi_j>^T * T_{j-1}^T

    where T_j models the multiscale transformation within V_{j+1}
    (T_j=I for j==j0-1).
    The system is diagonally preconditioned with the inverse of either
      (D)_{\lambda,\lambda}=2^{|\lambda|}
    (in the case dyadic==true) or
      (D)_{\lambda,\lambda}=a(\psi_\lambda,\psi_\lambda)
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
    FullLaplacian(const SplineBasis<d,dT,P_construction>& sb,
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
      set level j (stiffness matrix in V_j)
    */
    void set_level(const int j) const;

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
    
    //! specialization to std::map<size_type,double>
    void apply(const std::map<size_type,double>& x,
	       std::map<size_type,double>& Mx,
	       const bool preconditioning = true) const;
    
    //! conversion to a sparse matrix
    void to_sparse(SparseMatrix<double>& S) const;
    
    /*!
      stream output with user-defined tabwidth and precision
      (cf. deal.II)
    */
    void print(std::ostream& os,
	       const unsigned int tabwidth = 8,
	       const unsigned int precision = 3) const;

  protected:
    const SplineBasis<d,dT,P_construction>& sb_;
    PreconditioningType precond_;
    mutable int j_;

    mutable Vector<double> D_;
    void setup_D() const;
  };

  /*!
    Matlab-style stream output
  */
  template <int d, int dT>
  std::ostream& operator << (std::ostream& os, const FullLaplacian<d,dT>& M);
}

#include <galerkin/full_laplacian.cpp>

#endif
