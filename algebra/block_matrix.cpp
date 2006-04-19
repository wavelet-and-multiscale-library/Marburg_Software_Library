// implementation for block_matrix.h

namespace MathTL
{
  template <class C>
  BlockMatrix<C>::BlockMatrix(const size_type n)
    : blocks(n*n),
      row_blocks_columns(n), column_blocks_rows(n),
      rowdim_(n), coldim_(n)
  {
    for (size_type i = 0; i < n; i++) {
      row_blocks_columns[i] = 1;
      column_blocks_rows[i] = 1;
      for (size_type j = 0; j < n; j++)
	blocks[i*n+j] = 0;
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
  BlockMatrix<C>::row_blocks() const
  {
    return row_blocks_columns.size();
  }

  template <class C>
  inline
  const typename BlockMatrix<C>::size_type
  BlockMatrix<C>::column_blocks() const
  {
    return column_blocks_rows.size();
  }
  
  template <class C>
  inline
  const MatrixBlock<C>*
  BlockMatrix<C>::get_block(const size_type row, const size_type column) const
  {
    return blocks[row*row_blocks()+column];
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
 	for (typename BlockMatrix<C>::size_type i(0); i < M.column_blocks(); ++i)
	  for (typename BlockMatrix<C>::size_type j(0); j < M.row_blocks(); ++j) {
	    os << "block (" << i << "," << j << "):" << endl;
	    if (M.get_block(i, j)==0)
	      os << "ZERO" << endl;
	    else
	      M.get_block(i, j)->print(os);
	  }
      }
    
    return os;
  }
  
}
