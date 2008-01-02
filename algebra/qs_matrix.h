// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2008                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_QS_MATRIX_H
#define _WAVELETTL_QS_MATRIX_H

#include <cmath>
#include <map>
#include <algebra/matrix.h>
#include <algebra/sparse_matrix.h>

namespace MathTL
{
  /*!
    This class models quasi-stationary matrices like M_{j,0}, M_{j,1} or their dual
    counterparts. A quasi-stationary matrix is given by an upper left
    and a lower right corner block, and two bands for the left and right half,
    respectively. The bands are given by two stencil column vectors.
    For an odd number of columns, by convention,
    the central band stems from the _right_ half.
   */
  template <class C>
  class QuasiStationaryMatrix
  {
  public:
    /*!
      type of indexes and size type (cf. STL containers)
     */
    typedef typename Vector<C>::size_type size_type;

    //! default constructor, yields (2^j)x(2^j) zero matrix
    QuasiStationaryMatrix();

    /*!
      constructor from the matrix ingredients
      j0             : coarsest level
      mj0,nj0        : number of rows and columns on coarsest level
      ML,MR          : upper left and lower right corner block (may also be empty)
      bandL,bandR    : left and right band
      offsetL,offsetR: left and right row offset of the bands
      factor         : constant pre-factor
     */
    QuasiStationaryMatrix(const int j0,
			  const size_type mj0, const size_type nj0,
			  const Matrix<C>& ML, const Matrix<C>& MR,
			  const Vector<C>& bandL, const Vector<C>& bandR,
			  const size_type offsetL, const size_type offsetR,
			  const double factor = 1.0);

    //! assignment
    QuasiStationaryMatrix<C>& operator = (const QuasiStationaryMatrix<C>& M);
    
    /*!
      row dimension
    */
    const size_type row_dimension() const;

    /*!
      column dimension
    */
    const size_type column_dimension() const;

    /*!
      set level j
    */
    void set_level(const int j) const;

    /*!
      read-only access to a single matrix entry
    */
    const C get_entry(const size_type row, const size_type column) const;

    /*!
      construct a sparse matrix
    */
    void to_sparse(SparseMatrix<C>& S) const;

    /*!
      matrix-vector multiplication Mx = (*this) * x;
      it is possible to specify an offset which parts
      of the vectors shall be used, thus enabling in-place algorithms,
      the flag "add_to" toggles whether the result
      is added to Mx (true) or Mx is overwritten (false)
    */
    template <class VECTOR>
    void apply(const VECTOR& x, VECTOR& Mx,
	       const size_type x_offset = 0,
	       const size_type Mx_offset = 0,
	       const bool add_to = false) const;

    /*!
      matrix-vector multiplication Mx = (*this) * x,
      where the vector is modeled by std::map
    */
    void apply(const std::map<size_type, C>& x, std::map<size_type, C>& Mx,
	       const size_type x_offset = 0,
	       const size_type Mx_offset = 0,
	       const bool add_to = false) const;
    
    /*!
      transposed matrix-vector multiplication Mtx = (*this)^T * x;
      again potentially with offsets for both input and output vector
      and with an "add_to" flag
    */
    template <class VECTOR>
    void apply_transposed(const VECTOR& x, VECTOR& Mtx,
			  const size_type x_offset = 0,
			  const size_type Mtx_offset = 0,
			  const bool add_to = false) const;
    
    /*!
      transposed matrix-vector multiplication Mtx = (*this)^T * x;
      where the vector is modeled by std::map
    */
    void apply_transposed(const std::map<size_type, C>& x, std::map<size_type, C>& Mtx,
			  const size_type x_offset = 0,
			  const size_type Mtx_offset = 0,
			  const bool add_to = false) const;
    
    /*!
      apply the central block of the matrix (i.e. without first/last row/column) to a vector
    */
    void apply_central_block
    (const std::map<size_type, C>& x, std::map<size_type, C>& Mx,
     const size_type x_offset = 0,
     const size_type Mx_offset = 0,
     const bool add_to = false) const;

    /*!
      apply the transpose of the central block of the matrix (i.e. without first/last row/column) to a vector
    */
    void apply_central_block_transposed
    (const std::map<size_type, C>& x, std::map<size_type, C>& Mtx,
     const size_type x_offset = 0,
     const size_type Mtx_offset = 0,
     const bool add_to = false) const;

    /*!
      apply the central columns of the matrix (i.e. without first/last column) to a vector
    */
    void apply_central_columns
    (const std::map<size_type, C>& x, std::map<size_type, C>& Mx,
     const size_type x_offset = 0,
     const size_type Mx_offset = 0,
     const bool add_to = false) const;

    /*!
      apply the transpose of the central columns of the matrix (i.e. without first/last column) to a vector
    */
    void apply_central_columns_transposed
    (const std::map<size_type, C>& x, std::map<size_type, C>& Mtx,
     const size_type x_offset = 0,
     const size_type Mtx_offset = 0,
     const bool add_to = false) const;
    
    /*!
      stream output with user-defined tabwidth and precision
      (cf. deal.II)
    */
    void print(std::ostream& os,
	       const unsigned int tabwidth = 8,
	       const unsigned int precision = 3) const;

  protected:
    int j0_;
    mutable int j_;
    size_type mj0_, nj0_;
    Matrix<C> ML_, MR_;
    Vector<C> bandL_,bandR_;
    size_type offsetL_, offsetR_;
    double factor_;
  };

  /*!
    Matlab-style stream output for quasi-stationary matrices
  */
  template <class C>
  std::ostream& operator << (std::ostream& os, const QuasiStationaryMatrix<C>& M);
}

#include <algebra/qs_matrix.cpp>

#endif
