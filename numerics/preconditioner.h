// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2007                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_PRECONDITIONER_H
#define _MATHTL_PRECONDITIONER_H

#include <algebra/vector.h>

namespace MathTL
{
  /*!
    Abstract base class for preconditioners P of (finite) linear systems Ax=b.

    In applications, the class will be used to model both preconditioning
    - from the left  (solve P^{-1}Ax=P^{-1}b instead of Ax=b), and
    - from the right (solve AQ^{-1}y=b and x=Q^{-1}y instead of Ax=b)

    Preconditioners P behave like square matrices, the solution of linear
    systems Px=y should be cheap.

    References:
    
    [B] M. Benzi:
        Preconditioning Techniques for Large Linear Systems: A Survey
        J. Comput. Phys. 182(2002), 418-477
   */
  template <class VECTOR>
  class Preconditioner
  {
  public:
    /*!
      type of indices and size type (cf. STL containers)
    */
    typedef typename VECTOR::size_type size_type;

    /*!
      purely virtual destructor
    */
    virtual ~Preconditioner() = 0;
    
    /*!
      row dimension
    */
    virtual const size_type row_dimension() const = 0;

    /*!
      column dimension (== row dimension)
    */
    const size_type column_dimension() const { return row_dimension(); }

    /*!
      apply P, i.e., reverse the preconditioning
    */
    virtual void apply(const VECTOR& x, VECTOR& Px) const = 0;

    /*!
      apply P^{-1}, i.e., perform the preconditioning
    */
    virtual void apply_preconditioner(const VECTOR& Px, VECTOR& x) const = 0;
  };

  /*!
    Identity matrix as a preconditioner
  */
  template <class MATRIX, class VECTOR>
  class IdentityPreconditioner
    : public Preconditioner<VECTOR>
  {
  public:
    /*!
      type of indices and size type (cf. STL containers)
    */
    typedef typename VECTOR::size_type size_type;

    /*!
      default constructor, takes the matrix A as input parameter
    */
    IdentityPreconditioner(const MATRIX& A);

    /*!
      row dimension
    */
    const size_type row_dimension() const { return A.row_dimension(); }

    /*!
      apply P, i.e., reverse the preconditioning
    */
    void apply(const VECTOR& x, VECTOR& Px) const;

    /*!
      apply P^{-1}, i.e., perform the preconditioning
    */
    void apply_preconditioner(const VECTOR& Px, VECTOR& x) const;

  protected:
    /*!
      pointer to the matrix class under consideration
    */
    const MATRIX& A;
  };

  /*!
    As a nontrivial example class, provide the Jacobi preconditioner, i.e., P=diag(A).
  */
  template <class MATRIX, class VECTOR>
  class JacobiPreconditioner
    : public Preconditioner<VECTOR>
  {
  public:
    /*!
      type of indices and size type (cf. STL containers)
    */
    typedef typename VECTOR::size_type size_type;

    /*!
      default constructor, takes the matrix A as input parameter
    */
    JacobiPreconditioner(const MATRIX& A);

    /*!
      row dimension
    */
    const size_type row_dimension() const { return A.row_dimension(); }

    /*!
      apply P, i.e., reverse the preconditioning
    */
    void apply(const VECTOR& x, VECTOR& Px) const;

    /*!
      apply P^{-1}, i.e., perform the preconditioning
    */
    void apply_preconditioner(const VECTOR& Px, VECTOR& x) const;

  protected:
    /*!
      pointer to the matrix class under consideration
    */
    const MATRIX& A;
  };
}

#include <numerics/preconditioner.cpp>

#endif
