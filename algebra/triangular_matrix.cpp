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
  void LowerTriangularMatrix<C>::resize(const size_type rows, const size_type columns)
  {
    entries_.resize(triangle_size(rows,columns));
    rowdim_ = rows;
    coldim_ = columns;
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
    return row >= column ? entries_[triangle_index(row,column,rowdim_,coldim_)] : C(0);
  }

  template <class C>
  inline
  const C LowerTriangularMatrix<C>::get_entry(const size_type row,
					      const size_type column) const
  {
    return this->operator () (row, column);
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
  inline
  void LowerTriangularMatrix<C>::set_entry(const size_type row,
					   const size_type column,
					   const C value)
  {
    this->operator () (row, column) = value;
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
  inline
  void
  LowerTriangularMatrix<C>::scale(const C s)
  {
    entries_.scale(s);
  }
  
  template <class C>
  void
  LowerTriangularMatrix<C>::inverse(LowerTriangularMatrix<C>& MInv) const
  {
    assert(rowdim_ == coldim_);
    MInv.resize(rowdim_, coldim_);
    for (typename LowerTriangularMatrix<C>::size_type j(0); j < rowdim_; j++) {
      // j-th column of M^{-1} solves Mx=e_j
      for (typename LowerTriangularMatrix<C>::size_type i(j); i < rowdim_; i++) {
	double help = i == j ? 1.0 : 0.0;
	for (typename LowerTriangularMatrix<C>::size_type k(0); k < i; k++)
	  help -= get_entry(i, k) * MInv.get_entry(k, j);
	MInv.set_entry(i, j, help/get_entry(i, i));
      }
    }
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
    
    for (typename LowerTriangularMatrix<C>::size_type i(0); i < coldim_; i++)
      {
	Mtx[i] = 0;
	for (typename LowerTriangularMatrix<C>::size_type j(i);
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
	unsigned int old_precision = os.precision(precision);
	for (typename LowerTriangularMatrix<C>::size_type i(0);
	     i < row_dimension(); ++i)
	  {
	    for (typename LowerTriangularMatrix<C>::size_type j(0);
		 j < column_dimension(); ++j)
	      os << std::setw(tabwidth) << std::setprecision(precision)
		 << this->operator () (i, j);
	    os << std::endl;
	  }
	os.precision(old_precision);
      }
  }

  template <class C>
  inline
  std::ostream& operator << (std::ostream& os, const LowerTriangularMatrix<C>& M)
  {
    M.print(os);
    return os;
  }

  template <class C>
  inline
  UpperTriangularMatrix<C>::UpperTriangularMatrix(const size_type n)
    : LowerTriangularMatrix<C>(n)
  {
  }

  template <class C>
  inline
  UpperTriangularMatrix<C>::UpperTriangularMatrix(const UpperTriangularMatrix<C>& M)
    : LowerTriangularMatrix<C>(M)
  {
  }

  template <class C>
  UpperTriangularMatrix<C>::UpperTriangularMatrix(const size_type rows,
						  const size_type columns)
    : LowerTriangularMatrix<C>(columns, rows)
  {
  }

  template <class C>
  UpperTriangularMatrix<C>::UpperTriangularMatrix(const size_type rows,
						  const size_type columns,
						  const char* str,
						  const bool byrow)
    : LowerTriangularMatrix<C>(columns,rows,str,!byrow)
  {
  }

  template <class C>
  inline
  const typename UpperTriangularMatrix<C>::size_type
  UpperTriangularMatrix<C>::row_dimension() const
  {
    return LowerTriangularMatrix<C>::column_dimension();
  }

  template <class C>
  inline
  const typename UpperTriangularMatrix<C>::size_type
  UpperTriangularMatrix<C>::column_dimension() const
  {
    return LowerTriangularMatrix<C>::row_dimension();
  }

  template <class C>
  inline
  void UpperTriangularMatrix<C>::resize(const size_type rows, const size_type columns)
  {
    LowerTriangularMatrix<C>::resize(columns, rows);
  }

  template <class C>
  inline
  const C UpperTriangularMatrix<C>::operator () (const size_type row,
						 const size_type column) const
  {
    return LowerTriangularMatrix<C>::operator () (column, row);
  }

  template <class C>
  inline
  const C UpperTriangularMatrix<C>::get_entry(const size_type row,
					      const size_type column) const
  {
    return this->operator () (row, column);
  }

  template <class C>
  inline
  C& UpperTriangularMatrix<C>::operator () (const size_type row,
					    const size_type column)
  {
    return LowerTriangularMatrix<C>::operator () (column, row);
  }

  template <class C>
  inline
  void UpperTriangularMatrix<C>::set_entry(const size_type row,
					   const size_type column,
					   const C value)
  {
    this->operator () (row, column) = value;
  }

  template <class C>
  template <class C2>
  inline
  bool UpperTriangularMatrix<C>::operator == (const UpperTriangularMatrix<C2>& M) const
  {
    return LowerTriangularMatrix<C>::operator == (M);
  }
  
  template <class C>
  template <class C2>
  inline
  bool UpperTriangularMatrix<C>::operator != (const UpperTriangularMatrix<C2>& M) const
  {
    return !((*this) == M);
  }

  template <class C>
  UpperTriangularMatrix<C>&
  UpperTriangularMatrix<C>::operator = (const UpperTriangularMatrix<C>& M)
  {
    LowerTriangularMatrix<C>::operator = (M);
    return *this;
  }

  template <class C>
  template <class VECTOR>
  void UpperTriangularMatrix<C>::apply(const VECTOR& x, VECTOR& Mx) const
  {
    LowerTriangularMatrix<C>::apply_transposed(x, Mx);
  }

  template <class C>
  template <class VECTOR>
  void UpperTriangularMatrix<C>::apply_transposed(const VECTOR& x, VECTOR& Mtx) const
  {
    LowerTriangularMatrix<C>::apply(x, Mtx);
  }

  template <class C>
  void UpperTriangularMatrix<C>::print(std::ostream &os,
				       const unsigned int tabwidth,
				       const unsigned int precision) const
  {
    if (LowerTriangularMatrix<C>::empty())
      os << "[]" << std::endl; // Matlab style
    else
      {
	unsigned int old_precision = os.precision(precision);
	for (typename UpperTriangularMatrix<C>::size_type i(0);
	     i < row_dimension(); ++i)
	  {
	    for (typename UpperTriangularMatrix<C>::size_type j(0);
		 j < column_dimension(); ++j)
	      os << std::setw(tabwidth) << std::setprecision(precision)
		 << this->operator () (i, j);
	    os << std::endl;
	  }
	os.precision(old_precision);
      }
  }

  template <class C>
  inline
  std::ostream& operator << (std::ostream& os, const UpperTriangularMatrix<C>& M)
  {
    M.print(os);
    return os;
  }
}
