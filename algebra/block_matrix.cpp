// implementation for block_matrix.h

namespace MathTL
{
  template <class C>
  BlockMatrix<C>::BlockMatrix(const size_type n)
    : blocks(n*n),
      block_rows_rows(n), block_columns_columns(n),
      rowdim_(n), coldim_(n)
  {
    for (size_type i = 0; i < n; i++) {
      block_rows_rows[i] = 1;
      block_columns_columns[i] = 1;
      for (size_type j = 0; j < n; j++)
	blocks[i*n+j] = 0;
    }
  }

  template <class C>
  BlockMatrix<C>::BlockMatrix(const size_type block_rows, const size_type block_columns)
    : blocks(block_rows*block_columns),
      block_rows_rows(block_rows), block_columns_columns(block_columns),
      rowdim_(block_rows), coldim_(block_columns)
  {
    for (size_type i = 0; i < block_rows; i++) {
      block_rows_rows[i] = 1;
      for (size_type j = 0; j < block_columns; j++)
	blocks[i*block_columns+j] = 0;
    }
    for (size_type j = 0; j < block_columns; j++)
      block_columns_columns[j] = 1;
  }

  template <class C>
  BlockMatrix<C>::~BlockMatrix()
  {
    for (size_type i = 0; i < block_rows()*block_columns(); i++) {
      if (blocks[i])
	delete blocks[i];
    } 
  }
  
  template <class C>
  inline
  bool BlockMatrix<C>::empty() const
  {
    return row_dimension()*column_dimension()==0;
  }

  template <class C>
  inline
  const typename BlockMatrix<C>::size_type
  BlockMatrix<C>::block_rows() const
  {
    return block_rows_rows.size();
  }

  template <class C>
  inline
  const typename BlockMatrix<C>::size_type
  BlockMatrix<C>::block_columns() const
  {
    return block_columns_columns.size();
  }
  
  template <class C>
  inline
  const MatrixBlock<C>*
  BlockMatrix<C>::get_block(const size_type block_row, const size_type block_column) const
  {
    return blocks[block_row*block_columns()+block_column];
  }

  template <class C>
  inline
  void
  BlockMatrix<C>::set_block(const size_type row, const size_type column, MatrixBlock<C>* block)
  {
    if (blocks[row*block_columns()+column])
      delete blocks[row*block_columns()+column];
    blocks[row*block_columns()+column] = block;
  }

  template <class C>
  inline
  void
  BlockMatrix<C>::resize_block_row(const size_type block_row, const size_type rows)
  {
    block_rows_rows[block_row] = rows;
  }

  template <class C>
  inline
  void
  BlockMatrix<C>::resize_block_column(const size_type block_column, const size_type columns)
  {
    block_columns_columns[block_column] = columns;
  }

  template <class C>
  void BlockMatrix<C>::print(std::ostream &os,
			     const unsigned int tabwidth,
			     const unsigned int precision) const
  {
    // quick hack, for the moment
    os << *this;
  }

  template <class C>
  inline
  std::ostream& operator << (std::ostream& os, const BlockMatrix<C>& M)
  {
    if (M.empty())
      os << "[]" << std::endl; // Matlab style
    else
      {
 	for (typename BlockMatrix<C>::size_type i(0); i < M.block_rows(); ++i)
	  for (typename BlockMatrix<C>::size_type j(0); j < M.block_columns(); ++j) {
	    os << "block (" << i << "," << j << "):" << std::endl;
	    if (M.get_block(i, j)==0)
	      os << "ZERO" << std::endl;
	    else
	      M.get_block(i, j)->print(os);
	  }
      }
    
    return os;
  }
  
}
