// implementation of MathTL::Matrix inline functions

#include <iomanip>

namespace MathTL
{
  template <class C>
  inline
  Matrix<C>::Matrix(const size_type n)
    : entries_(n*n), rowdim_(n), coldim_(n)
  {
  }

  template <class C>
  inline
  Matrix<C>::Matrix(const Matrix<C>& M)
    : entries_(M.entries_), rowdim_(M.rowdim_), coldim_(M.coldim_)
  {
  }

  template <class C>
  inline
  Matrix<C>::Matrix(const size_type row_dimension, const size_type column_dimension)
    : entries_(row_dimension*column_dimension), rowdim_(row_dimension),
      coldim_(column_dimension)
  {
  }

  template <class C>
  inline
  const typename Matrix<C>::size_type
  Matrix<C>::row_dimension() const
  {
    return rowdim_;
  }

  template <class C>
  inline
  const typename Matrix<C>::size_type
  Matrix<C>::column_dimension() const
  {
    return coldim_;
  }

  template <class C>
  inline
  const typename Matrix<C>::size_type
  Matrix<C>::size() const
  {
    return row_dimension()*column_dimension();
  }

  template <class C>
  inline
  const typename Matrix<C>::size_type
  Matrix<C>::memory_consumption() const
  {
    return sizeof(*this) + row_dimension()*column_dimension()*sizeof(C);
  }

  template <class C>
  inline
  bool Matrix<C>::empty() const
  {
    return size()==0;
  }

  template <class C>
  template <class C2>
  bool Matrix<C>::operator == (const Matrix<C2>& M) const
  {
    if (rowdim_ != M.rowdim_ || coldim_ != M.coldim_) return false;
    return std::equal(entries_.begin(), entries_.end(), M.entries_.begin());
  }
  
  template <class C>
  template <class C2>
  inline
  bool Matrix<C>::operator != (const Matrix<C2>& M) const
  {
    return !((*this) == M);
  }

  template <class C>
  inline
  const C Matrix<C>::operator () (const size_type row, const size_type column) const
  {
    return entries_[row+column*rowdim_];
  }

  template <class C>
  inline
  C& Matrix<C>::operator () (const size_type row, const size_type column)
  {
    return entries_[row+column*rowdim_];
  }

  template <class C>
  Matrix<C>& Matrix<C>::operator = (const Matrix<C>& M)
  {
    assert(rowdim_ = M.rowdim_);
    assert(coldim_ = M.coldim_);

    entries_ = M.entries_;

    return *this;
  }

  template <class C>
  inline
  void Matrix<C>::swap(Matrix<C>& M)
  {
    std::swap(rowdim_, M.rowdim_);
    std::swap(coldim_, M.coldim_);
    entries_.swap(M.entries_);
  }

  template <class C>
  inline
  void swap(Matrix<C>& M1, Matrix<C>& M2)
  {
    M1.swap(M2);
  }

  template <class C>
  void Matrix<C>::print(std::ostream &os,
			const unsigned int tabwidth,
			const unsigned int precision) const
  {
    if (empty())
      os << "[]" << std::endl; // Matlab style
    else
      {
	for (typename Matrix<C>::size_type i(0); i < row_dimension(); ++i)
	  {
	    for (typename Matrix<C>::size_type j(0); j < column_dimension(); ++j)
	      os << std::setw(tabwidth) << std::setprecision(precision)
		 << this->operator () (i, j);
	    os << std::endl;
	  }
      }
  }
  
  template <class C>
  inline
  std::ostream& operator << (std::ostream& os, const Matrix<C>& M)
  {
    M.print(os);
    return os;
  }
}
