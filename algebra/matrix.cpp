// implementation of MathTL::Matrix inline functions

#include <iomanip>
#include <sstream>

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
  Matrix<C>::Matrix(const SymmetricMatrix<C>& M)
    : entries_(M.row_dimension()*M.column_dimension()),
      rowdim_(M.row_dimension()), coldim_(M.column_dimension())
  {
    for (size_type i(0); i < rowdim_; i++)
      for (size_type j(0); j <= i; j++)
	{
	  this->operator () (i, j) = this->operator () (j, i) = M(i, j);
	}
  }

  template <class C>
  inline
  Matrix<C>::Matrix(const LowerTriangularMatrix<C>& M)
    : entries_(M.row_dimension()*M.column_dimension()),
      rowdim_(M.row_dimension()), coldim_(M.column_dimension())
  {
    for (size_type i(0); i < rowdim_; i++)
      for (size_type j(0); j <= i; j++)
	{
	  this->operator () (i, j) = M(i, j);
	}
  }

  template <class C>
  inline
  Matrix<C>::Matrix(const size_type row_dimension,
		    const size_type column_dimension)
    : entries_(row_dimension*column_dimension), rowdim_(row_dimension),
      coldim_(column_dimension)
  {
  }

  template <class C>
  Matrix<C>::Matrix(const size_type row_dimension,
		    const size_type column_dimension,
		    const char* str,
		    const bool byrow)
    : entries_(row_dimension*column_dimension, false),
      rowdim_(row_dimension), coldim_(column_dimension)
  {
    std::istringstream ins(str);
    if (byrow)
      {
	for (size_type i(0); i < rowdim_ && ins.good(); i++)
	  for (size_type j(0); j < coldim_ && ins.good(); j++)
	    ins >> (*this).operator () (i,j);
      }
    else
      {
	for (size_type j(0); j < coldim_ && ins.good(); j++)
	  for (size_type i(0); i < rowdim_ && ins.good(); i++)
	    ins >> (*this).operator () (i,j);
      }
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
  void Matrix<C>::resize(const size_type rows, const size_type columns)
  {
    entries_.resize(rows*columns);
    rowdim_ = rows;
    coldim_ = columns;
  }

  template <class C>
  void Matrix<C>::diagonal(const size_type n, const C diag)
  {
    resize(n, n);

    for (size_type row(0); row < n; row++)
      set_entry(row, row, diag);
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
  inline
  const C Matrix<C>::operator () (const size_type row, const size_type column) const
  {
    return entries_[row+column*rowdim_];
  }

  template <class C>
  inline
  const C Matrix<C>::get_entry(const size_type row, const size_type column) const
  {
    return this->operator () (row, column);
  }

  template <class C>
  inline
  C& Matrix<C>::operator () (const size_type row, const size_type column)
  {
    return entries_[row+column*rowdim_];
  }

  template <class C>
  inline
  void Matrix<C>::set_entry(const size_type row, const size_type column, const C value)
  {
    this->operator () (row, column) = value;
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
  Matrix<C>& Matrix<C>::operator = (const Matrix<C>& M)
  {
    rowdim_ = M.rowdim_;
    coldim_ = M.coldim_;
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
  void Matrix<C>::mirror(Matrix<C>& M) const
  {
    M.resize(rowdim_, coldim_);
    for (typename Matrix<C>::size_type j(0); j < rowdim_; j++)
      for (typename Matrix<C>::size_type k(0); k < coldim_; k++)
	M(j, k) = this->operator () (rowdim_-j-1, coldim_-k-1);
  }

  template <class C>
  template <class VECTOR>
  void Matrix<C>::apply(const VECTOR& x, VECTOR& Mx) const
  {
    assert(Mx.size() == rowdim_);
    
    for (typename Matrix<C>::size_type i(0); i < rowdim_; i++)
      {
	Mx[i] = 0;
	for (typename Matrix<C>::size_type j(0);
	     j < coldim_; j++)
	  Mx[i] += this->operator () (i, j) * x[j];
      }
  }

  template <class C>
  template <class VECTOR>
  void Matrix<C>::apply_transposed(const VECTOR& x, VECTOR& Mtx) const
  {
    assert(Mtx.size() == coldim_);
    
    for (typename Matrix<C>::size_type i(0); i < coldim_; i++)
      {
	Mtx[i] = 0;
	for (typename Matrix<C>::size_type j(0);
	     j < rowdim_; j++)
	  Mtx[i] += this->operator () (j, i) * x[j];
      }
  }

  template <class C>
  inline
  void Matrix<C>::compress(const double eta)
  {
    entries_.compress(eta);
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
	unsigned int old_precision = os.precision(precision);
	for (typename Matrix<C>::size_type i(0); i < row_dimension(); ++i)
	  {
	    for (typename Matrix<C>::size_type j(0); j < column_dimension(); ++j)
	      os << std::setw(tabwidth) << std::setprecision(precision)
		 << this->operator () (i, j);
	    os << std::endl;
	  }
	os.precision(old_precision);
      }
  }
  
  template <class C>
  Matrix<C> operator - (const Matrix<C>& M, const Matrix<C>& N)
  {
    assert(M.column_dimension() == N.column_dimension());
    assert(M.row_dimension() == N.row_dimension());
    typedef typename Matrix<C>::size_type size_type;

    Matrix<C> R(M.row_dimension(), M.column_dimension());
    for (size_type i(0); i < M.row_dimension(); i++)
      for (size_type j(0); j < M.column_dimension(); j++)
	{
	  R(i, j) = M(i, j) - N(i, j);
	}

    return R;
  }

  template <class C>
  Matrix<C> operator * (const Matrix<C>& M, const Matrix<C>& N)
  {
    assert(M.column_dimension() == N.row_dimension());
    typedef typename Matrix<C>::size_type size_type;

    Matrix<C> R(M.row_dimension(), N.column_dimension());
    for (size_type i(0); i < M.row_dimension(); i++)
      for (size_type j(0); j < N.column_dimension(); j++)
	{
	  double help(0);
	  for (size_type k(0); k < N.row_dimension(); k++)
	    help += M(i, k) * N(k, j);
	  R(i, j) = help;
	}

    return R;
  }

  template <class C>
  Matrix<C> transpose(const Matrix<C>& M)
  {
    typedef typename Matrix<C>::size_type size_type;

    Matrix<C> R(M.column_dimension(), M.row_dimension());
    for (size_type i(0); i < M.row_dimension(); i++)
      for (size_type j(0); j < M.column_dimension(); j++)
	R(j, i) = M(i, j);

    return R;
  }

  template <class C>
  inline
  std::ostream& operator << (std::ostream& os, const Matrix<C>& M)
  {
    M.print(os);
    return os;
  }
}
