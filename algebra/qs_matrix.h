// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_QS_MATRIX_H
#define _WAVELETTL_QS_MATRIX_H

#include <algebra/matrix.h>

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

    /*!
      constructor from the matrix ingredients
      j0             : coarsest level
      mj0,nj0        : number of rows and columns on coarsest level
      ML,MR          : upper left and lower right corner block
      bandL,bandR    : left and right band
      offsetL,offsetR: left and right row offset of the bands
     */
    QuasiStationaryMatrix(const int j0,
			  const unsigned int mj0, const unsigned int nj0,
			  const Matrix<C>& ML, const Matrix<C>& MR,
			  const Vector<C>& bandL, const Vector<C>& bandR,
			  const int offsetL, const int offsetR);
    
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
    void set_level(const int j);

    /*!
      read-only access to a single matrix entry
    */
    const C get_entry(const size_type row, const size_type column) const;

    /*!
      matrix-vector multiplication Mx = (*this) * x;
      we assume that the vector Mx has the correct size and
      is not identical to x
    */
    template <class VECTOR>
    void apply(const VECTOR& x, VECTOR& Mx) const;

    /*!
      transposed matrix-vector multiplication Mtx = (*this)^T * x;
      we assume that the vector Mtx has the correct size and
      is not identical to x
    */
    template <class VECTOR>
    void apply_transposed(const VECTOR& x, VECTOR& Mtx) const;

    /*!
      stream output with user-defined tabwidth and precision
      (cf. deal.II)
    */
    void print(std::ostream& os,
	       const unsigned int tabwidth = 8,
	       const unsigned int precision = 3) const;

  protected:
    int j0_, j_;
    unsigned int mj0_, nj0_;
    Matrix<C> ML_, MR_;
    Vector<C> bandL_,bandR_;
    int offsetL_, offsetR_;
  };

  /*!
    Matlab-style stream output for quasi-stationary matrices
  */
  template <class C>
  std::ostream& operator << (std::ostream& os, const QuasiStationaryMatrix<C>& M);
}

#include <algebra/qs_matrix.cpp>

#endif
