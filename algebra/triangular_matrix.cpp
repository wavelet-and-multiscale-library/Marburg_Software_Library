// implementation of MathTL::*TriangularMatrix inline functions

#include <algorithm>
#include <cassert>
#include <iomanip>
#include <sstream>

namespace MathTL
{
  template <class C>
  inline
  typename LowerTriangularMatrix<C>::size_type
  LowerTriangularMatrix<C>::triangle_size(const size_type rows, const size_type columns)
  {
    typename LowerTriangularMatrix<C>::size_type m(std::min(rows,columns));
    return (rows <= columns ? m*(m+1)/2 : m*(m+1)/2 + (rows-columns)*columns);
  }

  template <class C>
  inline
  typename LowerTriangularMatrix<C>::size_type
  LowerTriangularMatrix<C>::triangle_index(const size_type row,
			    const size_type column,
			    const size_type rowdim,
			    const size_type coldim)
  {
    typename LowerTriangularMatrix<C>::size_type index(0);
    if (column > row)
      index = triangle_index(column, row, coldim, rowdim);
    else
      {
	if (rowdim <= coldim || row < coldim)
	  index = column+row*(row+1)/2;
	else
	  index = coldim*(coldim+1)/2 + (row-coldim)*coldim + column;
      }
    return index;
  }

  template <class C>
  LowerTriangularMatrix<C>::LowerTriangularMatrix(const size_type n)
    : entries_(triangle_size(n,n)), rowdim_(n), coldim_(n)
  {
  }

  template <class C>
  inline
  LowerTriangularMatrix<C>::LowerTriangularMatrix(const LowerTriangularMatrix<C>& M)
    : entries_(M.entries_), rowdim_(M.rowdim_), coldim_(M.coldim_)
  {
  }

  template <class C>
  LowerTriangularMatrix<C>::LowerTriangularMatrix(const size_type rows,
						  const size_type columns)
    : entries_(triangle_size(rows,columns)), rowdim_(rows), coldim_(columns)
  {
  }

  template <class C>
  LowerTriangularMatrix<C>::LowerTriangularMatrix(const size_type rows,
						  const size_type columns,
						  const char* str,
						  const bool byrow)
    : entries_(triangle_size(rows, columns), false),
      rowdim_(rows), coldim_(columns)
  {
    std::istringstream ins(str);
    if (byrow)
      {
	for (size_type i(0); i < rowdim_ && ins.good(); i++)
	  for (size_type j(0); j <= std::min(i, coldim_-1) && ins.good(); j++)
	    ins >> (*this).operator () (i,j);
      }
    else
      {
	for (size_type j(0); j < coldim_ && ins.good(); j++)
	  for (size_type i(j); i < rowdim_ && ins.good(); i++)
	    ins >> (*this).operator () (i,j);
      }
  }

  template <class C>
  inline
  const typename LowerTriangularMatrix<C>::size_type
  LowerTriangularMatrix<C>::row_dimension() const
  {
    return rowdim_;
  }

  template <class C>
  inline
  const typename LowerTriangularMatrix<C>::size_type
  LowerTriangularMatrix<C>::column_dimension() const
  {
    return coldim_;
  }

  template <class C>
  inline
  const typename LowerTriangularMatrix<C>::size_type
  LowerTriangularMatrix<C>::size() const
  {
    return triangle_size(row_dimension(),column_dimension());
  }

  template <class C>
  inline
  const typename LowerTriangularMatrix<C>::size_type
  LowerTriangularMatrix<C>::memory_consumption() const
  {
    return sizeof(*this) + triangle_size(row_dimension(),column_dimension())*sizeof(C);
  }

  template <class C>
  inline
  bool LowerTriangularMatrix<C>::empty() const
  {
    return size()==0;
  }

  template <class C>
  inline
  const C LowerTriangularMatrix<C>::operator () (const size_type row,
						 const size_type column) const
  {
    return row >= column ? entries_[triangle_index(row,column,rowdim_,coldim_)] : 0;
  }

  template <class C>
  inline
  C& LowerTriangularMatrix<C>::operator () (const size_type row,
					    const size_type column)
  {
    assert(row >= column);
    return entries_[triangle_index(row,column,rowdim_,coldim_)];
  }

  template <class C>
  template <class C2>
  bool LowerTriangularMatrix<C>::operator == (const LowerTriangularMatrix<C2>& M) const
  {
    if (rowdim_ != M.rowdim_ || coldim_ != M.coldim_) return false;
    return std::equal(entries_.begin(), entries_.end(), M.entries_.begin());
  }
  
  template <class C>
  template <class C2>
  inline
  bool LowerTriangularMatrix<C>::operator != (const LowerTriangularMatrix<C2>& M) const
  {
    return !((*this) == M);
  }

  template <class C>
  LowerTriangularMatrix<C>&
  LowerTriangularMatrix<C>::operator = (const LowerTriangularMatrix<C>& M)
  {
    rowdim_ = M.rowdim_;
    coldim_ = M.coldim_;
    entries_ = M.entries_;

    return *this;
  }

  template <class C>
  template <class VECTOR>
  void LowerTriangularMatrix<C>::apply(const VECTOR& x, VECTOR& Mx) const
  {
    assert(Mx.size() == rowdim_);
    
    for (typename LowerTriangularMatrix<C>::size_type i(0); i < rowdim_; i++)
      {
	Mx[i] = 0;
	for (typename LowerTriangularMatrix<C>::size_type j(0);
	     j <= std::min(i, coldim_-1); j++)
	  Mx[i] += this->operator () (i, j) * x[j];
      }
  }

  template <class C>
  template <class VECTOR>
  void LowerTriangularMatrix<C>::apply_transposed(const VECTOR& x, VECTOR& Mtx) const
  {
    assert(Mtx.size() == coldim_);
    
    for (typename Matrix<C>::size_type i(0); i < coldim_; i++)
      {
	Mtx[i] = 0;
	for (typename Matrix<C>::size_type j(i);
	     j < rowdim_; j++)
	  Mtx[i] += this->operator () (j, i) * x[j];
      }
  }

  template <class C>
  void LowerTriangularMatrix<C>::print(std::ostream &os,
				       const unsigned int tabwidth,
				       const unsigned int precision) const
  {
    if (empty())
      os << "[]" << std::endl; // Matlab style
    else
      {
	for (typename LowerTriangularMatrix<C>::size_type i(0);
	     i < row_dimension(); ++i)
	  {
	    for (typename LowerTriangularMatrix<C>::size_type j(0);
		 j < column_dimension(); ++j)
	      os << std::setw(tabwidth) << std::setprecision(precision)
		 << this->operator () (i, j);
	    os << std::endl;
	  }
      }
  }

  template <class C>
  inline
  std::ostream& operator << (std::ostream& os, const LowerTriangularMatrix<C>& M)
  {
    M.print(os);
    return os;
  }
}
