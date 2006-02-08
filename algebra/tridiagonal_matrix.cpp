// implementation for tridiagonal_matrix.h

#include <iomanip>

namespace MathTL
{
  template <class C>
  TridiagonalMatrix<C>::TridiagonalMatrix(const size_type n)
    : a_(n-1), b_(n), c_(n-1), rowdim_(n)
  {
  }

  template <class C>
  TridiagonalMatrix<C>::TridiagonalMatrix(const TridiagonalMatrix<C>& M)
    : a_(M.a_), b_(M.b_), c_(M.c_), rowdim_(M.rowdim_)
  {
  }

  template <class C>
  inline
  const typename TridiagonalMatrix<C>::size_type
  TridiagonalMatrix<C>::row_dimension() const
  {
    return rowdim_;
  }

  template <class C>
  inline
  const typename TridiagonalMatrix<C>::size_type
  TridiagonalMatrix<C>::column_dimension() const
  {
    return row_dimension(); // only square matrices!
  }

  template <class C>
  inline
  const typename TridiagonalMatrix<C>::size_type
  TridiagonalMatrix<C>::size() const
  {
    return 3*row_dimension()-2;
  }

  template <class C>
  void TridiagonalMatrix<C>::resize(const size_type rows)
  {
    a_.resize(rows-1);
    b_.resize(rows);
    c_.resize(rows-1);
    rowdim_ = rows;
  }

  template <class C>
  inline
  const typename TridiagonalMatrix<C>::size_type
  TridiagonalMatrix<C>::memory_consumption() const
  {
    return sizeof(*this) + size()*sizeof(C);
  }
  
  template <class C>
  inline
  bool TridiagonalMatrix<C>::empty() const
  {
    return size()==0;
  }

  template <class C>
  inline
  const C TridiagonalMatrix<C>::operator () (const size_type row,
					     const size_type column) const
  {
    if (row == column)   return b_[row];
    if (row+1 == column) return c_[row];
    if (row == column+1) return a_[column];
    
    return 0;
  }

  template <class C>
  inline
  const C TridiagonalMatrix<C>::get_entry(const size_type row,
					  const size_type column) const
  {
    return this->operator () (row, column);
  }

  template <class C>
  inline
  C& TridiagonalMatrix<C>::operator () (const size_type row,
					const size_type column)
  {
    assert(row+1 >= column && column+1 >= row);

    if (row+1 == column) return c_[row];
    if (row == column+1) return a_[column];
    
    return b_[row];
  }

  template <class C>
  inline
  void TridiagonalMatrix<C>::set_entry(const size_type row,
				       const size_type column,
				       const C value)
  {
    this->operator () (row, column) = value;
  }

  template <class C>
  template <class C2>
  bool TridiagonalMatrix<C>::operator == (const TridiagonalMatrix<C2>& M) const
  {
    if (rowdim_ != M.rowdim_) return false;
    if (!(std::equal(a_.begin(), a_.end(), M.a_.begin()))) return false;
    if (!(std::equal(b_.begin(), b_.end(), M.b_.begin()))) return false;
    if (!(std::equal(c_.begin(), c_.end(), M.c_.begin()))) return false;
    return true;
  }
  
  template <class C>
  template <class C2>
  inline
  bool TridiagonalMatrix<C>::operator != (const TridiagonalMatrix<C2>& M) const
  {
    return !((*this) == M);
  }

  template <class C>
  TridiagonalMatrix<C>&
  TridiagonalMatrix<C>::operator = (const TridiagonalMatrix<C>& M)
  {
    rowdim_ = M.rowdim_;
    a_ = M.a_;
    b_ = M.b_;
    c_ = M.c_;

    return *this;
  }

  template <class C>
  template <class VECTOR>
  void TridiagonalMatrix<C>::apply(const VECTOR& x, VECTOR& Mx) const
  {
    assert(Mx.size() == rowdim_);
    
    for (typename TridiagonalMatrix<C>::size_type i(0); i < rowdim_; i++)
      {
 	Mx[i] = b_[i] * x[i];
	
 	if (i > 0)
 	  Mx[i] += a_[i-1] * x[i-1];

 	if (i+1 < rowdim_)
 	  Mx[i] += c_[i] * x[i+1];
      }
  }

  template <class C>
  template <class VECTOR>
  void TridiagonalMatrix<C>::apply_transposed(const VECTOR& x, VECTOR& Mtx) const
  {
    assert(Mtx.size() == rowdim_);
    
    for (typename TridiagonalMatrix<C>::size_type i(0); i < rowdim_; i++)
      {
 	Mtx[i] = b_[i] * x[i];
	
	if (i > 0)
	  Mtx[i] += c_[i-1] * x[i-1];

	if (i+1 < rowdim_)
	  Mtx[i] += a_[i] * x[i+1];
      }
  }

  template <class C>
  void TridiagonalMatrix<C>::print(std::ostream &os,
				   const unsigned int tabwidth,
				   const unsigned int precision) const
  {
    if (empty())
      os << "[]" << std::endl; // Matlab style
    else
      {
	unsigned int old_precision = os.precision(precision);
	for (typename TridiagonalMatrix<C>::size_type i(0);
	     i < row_dimension(); ++i)
	  {
	    for (typename TridiagonalMatrix<C>::size_type j(0);
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
  std::ostream& operator << (std::ostream& os, const TridiagonalMatrix<C>& M)
  {
    M.print(os);
    return os;
  }
}
