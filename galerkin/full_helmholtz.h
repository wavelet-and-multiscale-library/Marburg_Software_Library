// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2008                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_FULL_HELMHOLTZ_H
#define _WAVELETTL_FULL_HELMHOLTZ_H

#include <iostream>
#include <map>
#include <algebra/vector.h>
#include <interval/spline_basis.h>
#include <galerkin/full_laplacian.h>
#include <galerkin/full_gramian.h>

using namespace MathTL;

namespace WaveletTL
{
  /*!
    The following class models the discretization of the
    1D Helmholtz operator alpha*I-Delta
    in a primal spline wavelet basis on a level j.

    Same as in FullLaplacian, there are two possible preconditioning
    strategies available.

    We assume that the basis has homogeneous b.c. at least at one interval end,
    i.e., that s0+s1>=1.
  */
  template <int d, int dT, int s0, int s1, int J0>
  class FullHelmholtz
  {
  public:
    /*!
      type of indexes and size type (cf. STL containers)
     */
    typedef typename Vector<double>::size_type size_type;

    /*!
      constructor taking an information object on some spline wavelet basis
    */
    FullHelmholtz(const SplineBasis<d,dT,P_construction,s0,s1,0,0,J0>& sb,
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
      set level j (stiffness matrix in V_j)
    */
    void set_level(const int j) const;

    /*!
      set reaction coefficient alpha
    */
    void set_alpha(const double alpha) const;

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

    FullGramian<d,dT,s0,s1,0,0,J0> G_;
    FullLaplacian<d,dT,s0,s1,J0> A_;
    const SplineBasis<d,dT,P_construction,s0,s1,0,0,J0>& sb_;

  protected:
    mutable double alpha_;
    PreconditioningType precond_;

    mutable int j_;

    mutable Vector<double> D_;
    void setup_D() const;
  };

  /*!
    Matlab-style stream output
  */
  template <int d, int dT, int s0, int s1, int J0>
  std::ostream& operator << (std::ostream& os, const FullHelmholtz<d,dT,s0,s1,J0>& M);
}

#include <galerkin/full_helmholtz.cpp>

#endif
