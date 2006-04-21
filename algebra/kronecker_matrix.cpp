// implementation for kronecker_matrix.h

#include <iomanip>

namespace MathTL
{
  template <class C, class MATRIX1, class MATRIX2>
  KroneckerMatrix<C,MATRIX1,MATRIX2>::KroneckerMatrix(const MATRIX1& A_, const MATRIX2& B_)
    : A(A_), B(B_)
  {
  }

  template <class C, class MATRIX1, class MATRIX2>
  KroneckerMatrix<C,MATRIX1,MATRIX2>::KroneckerMatrix
  (const KroneckerMatrix<C,MATRIX1,MATRIX2>& M)
    : A(M.A), B(M.B)
  {
  }

  template <class C, class MATRIX1, class MATRIX2>
  inline
  const typename KroneckerMatrix<C,MATRIX1,MATRIX2>::size_type
  KroneckerMatrix<C,MATRIX1,MATRIX2>::row_dimension() const
  {
    return A.row_dimension() * B.row_dimension();
  }

  template <class C, class MATRIX1, class MATRIX2>
  inline
  const typename KroneckerMatrix<C,MATRIX1,MATRIX2>::size_type
  KroneckerMatrix<C,MATRIX1,MATRIX2>::column_dimension() const
  {
    return A.column_dimension() * B.column_dimension();
  }

  template <class C, class MATRIX1, class MATRIX2>
  inline
  bool KroneckerMatrix<C,MATRIX1,MATRIX2>::empty() const
  {
    return A.empty() || B.empty();
  }

  template <class C, class MATRIX1, class MATRIX2>
  inline
  const C
  KroneckerMatrix<C,MATRIX1,MATRIX2>::operator () (const size_type row,
						   const size_type column) const
  {
    return A(row / B.row_dimension(), column / B.column_dimension())
      * B(row % B.row_dimension(), column % B.column_dimension());
  }

  template <class C, class MATRIX1, class MATRIX2>
  inline
  const C
  KroneckerMatrix<C,MATRIX1,MATRIX2>::get_entry(const size_type row,
						const size_type column) const
  {
    return this->operator () (row, column);
  }

  template <class C, class MATRIX1, class MATRIX2>
  template <class VECTOR>
  void
  KroneckerMatrix<C,MATRIX1,MATRIX2>::apply(const VECTOR& x, VECTOR& Mx) const
  {
    assert(Mx.size() == row_dimension());
    
    for (size_type i(0); i < row_dimension(); i++) {
      Mx[i] = 0;
      for (size_type j(0); j < column_dimension(); j++)
	Mx[i] += this->operator () (i, j) * x[j];
    }
  }
  
  template <class C, class MATRIX1, class MATRIX2>
  void
  KroneckerMatrix<C,MATRIX1,MATRIX2>::apply(const Vector<C>& x, Vector<C>& Mx) const
  {
    assert(Mx.size() == row_dimension());
    
    for (size_type i(0); i < row_dimension(); i++) {
      Mx[i] = 0;
      for (size_type j(0); j < column_dimension(); j++)
	Mx[i] += this->operator () (i, j) * x[j];
    }
  }
  
  template <class C, class MATRIX1, class MATRIX2>
  template <class VECTOR>
  void
  KroneckerMatrix<C,MATRIX1,MATRIX2>::apply_transposed(const VECTOR& x, VECTOR& Mtx) const
  {
    assert(Mtx.size() == column_dimension());
    
    for (size_type i(0); i < column_dimension(); i++) {
      Mtx[i] = 0;
      for (size_type j(0); j < row_dimension(); j++)
	Mtx[i] += this->operator () (j, i) * x[j];
    }
  }
  
  template <class C, class MATRIX1, class MATRIX2>
  void
  KroneckerMatrix<C,MATRIX1,MATRIX2>::apply_transposed(const Vector<C>& x, Vector<C>& Mtx) const
  {
    assert(Mtx.size() == column_dimension());
    
    for (size_type i(0); i < column_dimension(); i++) {
      Mtx[i] = 0;
      for (size_type j(0); j < row_dimension(); j++)
	Mtx[i] += this->operator () (j, i) * x[j];
    }
  }
  
  template <class C, class MATRIX1, class MATRIX2>
  void KroneckerMatrix<C,MATRIX1,MATRIX2>::print(std::ostream &os,
						 const unsigned int tabwidth,
						 const unsigned int precision) const
  {
    if (empty())
      os << "[]" << std::endl; // Matlab style
    else
      {
 	unsigned int old_precision = os.precision(precision);
 	for (size_type i(0); i < row_dimension(); ++i) {
	  for (size_type j(0); j < column_dimension(); ++j)
	    os << std::setw(tabwidth) << std::setprecision(precision)
	       << this->operator () (i, j);
	  os << std::endl;
	}
 	os.precision(old_precision);
      }
  }

  template <class C, class MATRIX1, class MATRIX2>
  inline
  std::ostream& operator << (std::ostream& os, const KroneckerMatrix<C,MATRIX1,MATRIX2>& M)
  {
    M.print(os);
    return os;
  }
}
