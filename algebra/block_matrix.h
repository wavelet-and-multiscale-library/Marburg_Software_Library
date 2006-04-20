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

    //! constructor for (mxn) block (each 1x1) matrices, empty per default
    BlockMatrix(const size_type block_rows, const size_type block_columns);

    //! destructor (all subblocks will be deleted)
    ~BlockMatrix();
    
    //! (overall) row dimension
    const size_type row_dimension() const { return rowdim_; }
    
    //! (overall) column dimension
    const size_type column_dimension() const { return coldim_; }

    //! number of "block rows", i.e., blocks in column direction (vertical)
    const size_type block_rows() const;

    //! number of "block columns", i.e., blocks in row direction (horizontal)
    const size_type block_columns() const;

    //! get a single matrix block
    const MatrixBlock<C>* get_block(const size_type block_row, const size_type block_column) const;

    /*!
      Set a single matrix block. Note that the old block will be deleted, and ~BlockMatrix() will
      delete this block!
      Be sure to resize the corresponding row/column block properly before or after this call.
      (this is necessary to have zero row/column blocks without having to allocate them)
    */
    void set_block(const size_type row, const size_type column, MatrixBlock<C>* block);
    
    //! resize a row block (i.e., its number of rows)
    void resize_block_row(const size_type block_row, const size_type rows);

    //! resize a block column (i.e., its number of columns)
    void resize_block_column(const size_type block_row, const size_type columns);

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

    //! row dimensions of the block rows
    Array1D<size_type> block_rows_rows;
    
    //! column dimensions of the block columns
    Array1D<size_type> block_columns_columns;

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
