// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_PRECOND_H
#define _WAVELETTL_PRECOND_H

#include <algebra/infinite_vector.h>
#include <algebra/infinite_matrix.h>

using MathTL::InfiniteVector;
using MathTL::InfiniteDiagonalMatrix;

namespace WaveletTL
{
  /*!
    Base class for preconditioners of infinite linear systems Ax=b, i.e.,
    for those arising in wavelet-Galerkin discretizations of PDEs

    The class will be used to model both preconditioning
    - from the left  (solve P^{-1}Ax=P^{-1}b instead of Ax=b), and
    - from the right (solve AQ^{-1}y=b and x=Q^{-1}y instead of Ax=b)

    Preconditioners in wavelet-Galerkin schemes should have the property
    that P^{-1} and Q^{-1} preserve the active support sets, at least for
    the proper wavelet coefficients (e!=0). In other words, the preconditioner
    should behave like a diagonal matrix for rows corresponding to
    wavelet coefficients. For fully diagonal preconditioners,
    this is of course no problem. However, in the case of a special treatment of
    the "generator blocks" in the stiffness matrix, an application
    of P^{-1} or Q^{-1} might introduce additional generator coefficients
    which were not yet present in the argument vector.

    [B] M. Benzi:
        Preconditioning Techniques for Large Linear Systems: A Survey
        J. Comput. Phys. 182(2002), 418-477
   */
  template <class INDEX>
  class InfinitePreconditioner
  {
  public:
    /*!
      purely virtual destructor
    */
    virtual ~InfinitePreconditioner() = 0;

    /*!
      Apply the left preconditioner to a given vector y, i.e., solve Px = y.
      At least the active proper wavelet coefficients (e!=0) should be preserved.
    */
    virtual void apply_left_preconditioner
    (const InfiniteVector<double,INDEX>& y,
     InfiniteVector<double,INDEX>& x) const = 0;
    
    /*!
      Apply the right preconditioner to a given vector y, i.e., solve Qx = y.
      At least the active proper wavelet coefficients (e!=0) should be preserved.
    */
    virtual void apply_right_preconditioner
    (const InfiniteVector<double,INDEX>& y,
     InfiniteVector<double,INDEX>& x) const = 0;
  };
  
  /*!
    Base class for symmetric preconditioners, i.e., P=Q
  */
  template <class INDEX>
  class InfiniteSymmetricPreconditioner
    : public InfinitePreconditioner<INDEX>
  {
  public:
    /*!
      Apply the left preconditioner to a given vector y, i.e., solve Px = y.
    */
    void apply_left_preconditioner(const InfiniteVector<double,INDEX>& y,
				   InfiniteVector<double,INDEX>& x) const;
    
    /*!
      Apply the right preconditioner to a given vector y, i.e., solve Px = y.
      At least the active proper wavelet coefficients (e!=0) should be preserved.
    */
    void apply_right_preconditioner(const InfiniteVector<double,INDEX>& y,
				    InfiniteVector<double,INDEX>& x) const;
    
    /*!
      Apply the preconditioner to a given vector y, i.e., solve Px = y.
      At least the active proper wavelet coefficients (e!=0) should be preserved.
    */
    virtual void apply_preconditioner(const InfiniteVector<double,INDEX>& y,
				      InfiniteVector<double,INDEX>& x) const = 0;
  };
  
  /*!
    Base class for fully diagonal symmetric preconditioners P=Q=D=diag(d_i).
    In fact, this class can be derived from InfiniteDiagonalMatrix,
    making the fast routine InfiniteVector::scale() accessible.
  */
  template <class INDEX>
  class FullyDiagonalPreconditioner
    : public InfiniteSymmetricPreconditioner<INDEX>,
      public InfiniteDiagonalMatrix<double,INDEX>
  {
  public:
    /*!
      evaluate the diagonal preconditioner D
    */
    virtual double diag(const INDEX& lambda) const = 0;
    
    /*!
      apply preconditioner, x=P^{-1}y=Q^{-1}y=D^{-1}y
    */
    void apply_preconditioner(const InfiniteVector<double,INDEX>& y,
			      InfiniteVector<double,INDEX>& x) const;
  };
  
  /*!
    When working with wavelet-Galerkin schemes, any differential or integral operator
    of order 2t induces a natural diagonal symmetric preconditioner
    D=diag(2^{t|lambda|}), which is modeled by this class.
  */
  template <class INDEX>
  class WaveletNEPreconditioner
    : public FullyDiagonalPreconditioner<INDEX>
  {
    /*!
      (half) operator order t
     */
    virtual double operator_order() const = 0;

    /*!
      evaluate the diagonal preconditioner D
    */
    double diag(const INDEX& lambda) const;
  };
}

#include <galerkin/precond.cpp>

#endif
