// implementation for kronecker_matrix.h

#include <iomanip>

namespace MathTL
{
  template <class C, class MATRIX1, class MATRIX2>
  KroneckerMatrix<C,MATRIX1,MATRIX2>::KroneckerMatrix(const MATRIX1& A_, const MATRIX2& B_,
						      const double factor)
    : A(A_), B(B_)
  {
    A.scale(factor);
  }

  template <class C, class MATRIX1, class MATRIX2>
  KroneckerMatrix<C,MATRIX1,MATRIX2>::KroneckerMatrix
  (const KroneckerMatrix<C,MATRIX1,MATRIX2>& M)
    : A(M.A), B(M.B)
  {
  }

  template <class C, class MATRIX1, class MATRIX2>
  MatrixBlock<C>*
  KroneckerMatrix<C,MATRIX1,MATRIX2>::clone() const
  {
    return new KroneckerMatrix<C,MATRIX1,MATRIX2>(*this);
  }

  template <class C, class MATRIX1, class MATRIX2>
  MatrixBlock<C>*
  KroneckerMatrix<C,MATRIX1,MATRIX2>::clone_transposed() const
  {
    return new KroneckerMatrix<C,MATRIX1,MATRIX2>(transpose(A),
						  transpose(B));
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
    return A.get_entry(row / B.row_dimension(), column / B.column_dimension())
      * B.get_entry(row % B.row_dimension(), column % B.column_dimension());
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
    
    for (size_type Arow(0), Mxrow(0); Arow < A.row_dimension(); Arow++)
      for (size_type Brow(0); Brow < B.row_dimension(); Brow++, Mxrow++) {
	C help(0);
	
	for (size_type Acol(0), xrow(0); Acol < A.column_dimension(); Acol++) {
	  C helpA(A.get_entry(Arow, Acol));
	  if (helpA != 0) {
	    for (size_type Bcol(0); Bcol < B.column_dimension(); Bcol++, xrow++)
	      help += helpA * B.get_entry(Brow, Bcol) * x[xrow];
	  } else {
	    xrow += B.column_dimension(); // !
	  }
	}
	
	Mx[Mxrow] = help;
      }
  }
  
  template <class C, class MATRIX1, class MATRIX2>
  void
  KroneckerMatrix<C,MATRIX1,MATRIX2>::apply(const Vector<C>& x, Vector<C>& Mx) const
  {
    assert(Mx.size() == row_dimension());

    for (size_type Arow(0), Mxrow(0); Arow < A.row_dimension(); Arow++)
      for (size_type Brow(0); Brow < B.row_dimension(); Brow++, Mxrow++) {
	C help(0);
	
	for (size_type Acol(0), xrow(0); Acol < A.column_dimension(); Acol++) {
	  C helpA(A.get_entry(Arow, Acol));
	  if (helpA != 0) {
	    for (size_type Bcol(0); Bcol < B.column_dimension(); Bcol++, xrow++)
	      help += helpA * B.get_entry(Brow, Bcol) * x[xrow];
	  } else {
	    xrow += B.column_dimension(); // !
	  }
	}
	
	Mx[Mxrow] = help;
      }
  }
  
  template <class C, class MATRIX1, class MATRIX2>
  template <class VECTOR>
  void
  KroneckerMatrix<C,MATRIX1,MATRIX2>::apply_transposed(const VECTOR& x, VECTOR& Mtx) const
  {
    assert(Mtx.size() == column_dimension());
    
    for (size_type Acol(0), Mtxrow(0); Acol < A.column_dimension(); Acol++)
      for (size_type Bcol(0); Bcol < B.column_dimension(); Bcol++, Mtxrow++) {
	C help(0);
	
	for (size_type Arow(0), xrow(0); Arow < A.row_dimension(); Arow++) {
	  C helpA(A.get_entry(Arow, Acol));
	  if (helpA != 0) {
	    for (size_type Brow(0); Brow < B.row_dimension(); Brow++, xrow++)
	      help += helpA * B.get_entry(Brow, Bcol) * x[xrow];
	  } else {
	    xrow += B.row_dimension(); // !
	  }
	}
	
	Mtx[Mtxrow] = help;
      }
  }
  
  template <class C, class MATRIX1, class MATRIX2>
  void
  KroneckerMatrix<C,MATRIX1,MATRIX2>::apply_transposed(const Vector<C>& x, Vector<C>& Mtx) const
  {
    assert(Mtx.size() == column_dimension());
    
    for (size_type Acol(0), Mtxrow(0); Acol < A.column_dimension(); Acol++)
      for (size_type Bcol(0); Bcol < B.column_dimension(); Bcol++, Mtxrow++) {
	C help(0);
	
	for (size_type Arow(0), xrow(0); Arow < A.row_dimension(); Arow++) {
	  C helpA(A.get_entry(Arow, Acol));
	  if (helpA != 0) {
	    for (size_type Brow(0); Brow < B.row_dimension(); Brow++, xrow++)
	      help += helpA * B.get_entry(Brow, Bcol) * x[xrow];
	  } else {
	    xrow += B.row_dimension(); // !
	  }
	}
	
	Mtx[Mtxrow] = help;
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
