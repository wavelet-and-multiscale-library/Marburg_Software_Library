// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_BLOCK_MATRIX_H
#define _MATHTL_BLOCK_MATRIX_H

#include <iostream>
#include <utils/array1d.h>
#include <algebra/vector.h>
#include <algebra/matrix_block.h>

// matrix norms, for convenience
#include <algebra/matrix_norms.h>

namespace MathTL
{
  /*!
    A class for block matrices with almost arbitrary subblocks.
  */
  template <class C>
  class BlockMatrix
    : public MatrixBlock<C>
  {
  public:
    /*!
      type of indexes and size type (cf. STL containers)
     */
    typedef typename Vector<C>::size_type size_type;
    
    //! default constructor, yields a zero (nxn)-block (each 1x1) matrix, empty per default
    explicit BlockMatrix(const size_type n = 0);

    //! row dimension
    const size_type row_dimension() const { return rowdim_; }
    
    //! column dimension
    const size_type column_dimension() const { return coldim_; }

    //! number of blocks in row direction (horizontal)
    const size_type row_blocks() const;

    //! number of blocks in column direction (vertical)
    const size_type column_blocks() const;

    //! get matrix block
    const MatrixBlock<C>* get_block(const size_type row, const size_type column) const;

    //! return true if matrix is empty (cf. STL containers)
    bool empty() const;

    /*!
      stream output with user-defined tabwidth and precision
      (cf. deal.II)
    */
    void print(std::ostream& os,
	       const unsigned int tabwidth = 10,
	       const unsigned int precision = 3) const;

  protected:
    //! pointers to the blocks
    Array1D<MatrixBlock<C>*> blocks;

    //! column dimensions of the blocks in row direction (horizontal)
    Array1D<size_type> row_blocks_columns;

    //! row dimensions of the blocks in column direction (vertical)
    Array1D<size_type> column_blocks_rows;
    
    //! (overall) row dimension
    size_type rowdim_;

    //! (overall) column dimension
    size_type coldim_;
  };

  /*!
    stream output for block matrices
  */
  template <class C>
  std::ostream& operator << (std::ostream& os, const BlockMatrix<C>& M);
  
}

// include implementation of inline functions
#include <algebra/block_matrix.cpp>

#endif
