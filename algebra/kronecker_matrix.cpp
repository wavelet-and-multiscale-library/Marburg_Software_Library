// implementation for kronecker_matrix.h

#include <iomanip>
#include <utils/array1d.h>

namespace MathTL
{
  template <class C, class MATRIX1, class MATRIX2>
  KroneckerMatrix<C,MATRIX1,MATRIX2>::KroneckerMatrix(const MATRIX1& A_, const MATRIX2& B_,
						      const double factor)
    : A(A_), B(B_)
  {
    if (factor != 1.0)
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
    
#if 1
    // new variant, using only the apply() methods

    // for readability:
    const size_type m = A.row_dimension();
    const size_type n = A.column_dimension();
    const size_type p = B.row_dimension();
    const size_type q = B.column_dimension();

    // first decompose x into the corresponding blocks
    Array1D<Vector<C> > x_blocks(n);
    for (size_type j = 0; j < n; j++) {
      x_blocks[j].resize(q);
      for (size_type l = 0; l < q; l++)
	x_blocks[j][l] = x[j*q+l];
    }

    // apply B to these blocks
    Array1D<Vector<C> > Bx_blocks(n);
    for (size_type j = 0; j < n; j++) {
      Bx_blocks[j].resize(p);
      B.apply(x_blocks[j], Bx_blocks[j]);
    }
    
    // resort the Bx[i] blocks
    Array1D<Vector<C> > w_blocks(p);
    for (size_type k = 0; k < p; k++) {
      w_blocks[k].resize(n);
      for (size_type j = 0; j < n; j++)
	w_blocks[k][j] = Bx_blocks[j][k];
    }
    
    // apply A to the w[k] blocks
    Array1D<Vector<C> > Aw_blocks(p);
    for (size_type k = 0; k < p; k++) {
      Aw_blocks[k].resize(m);
      A.apply(w_blocks[k], Aw_blocks[k]);
    }
    
    // resort the Aw[k] blocks into the result
    for (size_type i = 0; i < m; i++)
      for (size_type k = 0; k < p; k++)
	Mx[i*p+k] = Aw_blocks[k][i];
#else
    // old code, resorting to single entries
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
#endif
  }
  
  template <class C, class MATRIX1, class MATRIX2>
  void
  KroneckerMatrix<C,MATRIX1,MATRIX2>::apply(const Vector<C>& x, Vector<C>& Mx) const
  {
    assert(Mx.size() == row_dimension());

#if 1
    // new variant, using only the apply() methods

    // for readability:
    const size_type m = A.row_dimension();
    const size_type n = A.column_dimension();
    const size_type p = B.row_dimension();
    const size_type q = B.column_dimension();

    // first decompose x into the corresponding blocks
    Array1D<Vector<C> > x_blocks(n);
    for (size_type j = 0; j < n; j++) {
      x_blocks[j].resize(q);
      for (size_type l = 0; l < q; l++)
	x_blocks[j][l] = x[j*q+l];
    }

    // apply B to these blocks
    Array1D<Vector<C> > Bx_blocks(n);
    for (size_type j = 0; j < n; j++) {
      Bx_blocks[j].resize(p);
      B.apply(x_blocks[j], Bx_blocks[j]);
    }
    
    // resort the Bx[i] blocks
    Array1D<Vector<C> > w_blocks(p);
    for (size_type k = 0; k < p; k++) {
      w_blocks[k].resize(n);
      for (size_type j = 0; j < n; j++)
	w_blocks[k][j] = Bx_blocks[j][k];
    }
    
    // apply A to the w[k] blocks
    Array1D<Vector<C> > Aw_blocks(p);
    for (size_type k = 0; k < p; k++) {
      Aw_blocks[k].resize(m);
      A.apply(w_blocks[k], Aw_blocks[k]);
    }
    
    // resort the Aw[k] blocks into the result
    for (size_type i = 0; i < m; i++)
      for (size_type k = 0; k < p; k++)
	Mx[i*p+k] = Aw_blocks[k][i];
#else
    // old code, resorting to single entries
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
#endif
  }
  
  template <class C, class MATRIX1, class MATRIX2>
  template <class VECTOR>
  void
  KroneckerMatrix<C,MATRIX1,MATRIX2>::apply_transposed(const VECTOR& x, VECTOR& Mtx) const
  {
    assert(Mtx.size() == column_dimension());
    
#if 1
    // new variant, using only the apply() methods

    // for readability:
    const size_type m = A.row_dimension();
    const size_type n = A.column_dimension();
    const size_type p = B.row_dimension();
    const size_type q = B.column_dimension();

    // first decompose x into the corresponding blocks
    Array1D<Vector<C> > x_blocks(m);
    for (size_type i = 0; i < m; i++) {
      x_blocks[i].resize(p);
      for (size_type k = 0; k < p; k++)
	x_blocks[i][k] = x[i*p+k];
    }

    // apply B^T to these blocks
    Array1D<Vector<C> > BTx_blocks(m);
    for (size_type i = 0; i < m; i++) {
      BTx_blocks[i].resize(q);
      B.apply_transposed(x_blocks[i], BTx_blocks[i]);
    }
    
    // resort the BTx[i] blocks
    Array1D<Vector<C> > w_blocks(q);
    for (size_type l = 0; l < q; l++) {
      w_blocks[l].resize(m);
      for (size_type i = 0; i < m; i++)
 	w_blocks[l][i] = BTx_blocks[i][l];
    }
    
    // apply A^T to the w[l] blocks
    Array1D<Vector<C> > ATw_blocks(q);
    for (size_type l = 0; l < q; l++) {
      ATw_blocks[l].resize(n);
      A.apply_transposed(w_blocks[l], ATw_blocks[l]);
    }
    
    // resort the ATw[l] blocks into the result
    for (size_type j = 0; j < n; j++)
      for (size_type l = 0; l < q; l++)
 	Mtx[j*q+l] = ATw_blocks[l][j];
#else
    // old code, resorting to single entries
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
#endif
  }
  
  template <class C, class MATRIX1, class MATRIX2>
  void
  KroneckerMatrix<C,MATRIX1,MATRIX2>::apply_transposed(const Vector<C>& x, Vector<C>& Mtx) const
  {
    assert(Mtx.size() == column_dimension());
    
#if 1
    // new variant, using only the apply() methods

    // for readability:
    const size_type m = A.row_dimension();
    const size_type n = A.column_dimension();
    const size_type p = B.row_dimension();
    const size_type q = B.column_dimension();

    // first decompose x into the corresponding blocks
    Array1D<Vector<C> > x_blocks(m);
    for (size_type i = 0; i < m; i++) {
      x_blocks[i].resize(p);
      for (size_type k = 0; k < p; k++)
	x_blocks[i][k] = x[i*p+k];
    }

    // apply B^T to these blocks
    Array1D<Vector<C> > BTx_blocks(m);
    for (size_type i = 0; i < m; i++) {
      BTx_blocks[i].resize(q);
      B.apply_transposed(x_blocks[i], BTx_blocks[i]);
    }
    
    // resort the BTx[i] blocks
    Array1D<Vector<C> > w_blocks(q);
    for (size_type l = 0; l < q; l++) {
      w_blocks[l].resize(m);
      for (size_type i = 0; i < m; i++)
 	w_blocks[l][i] = BTx_blocks[i][l];
    }
    
    // apply A^T to the w[l] blocks
    Array1D<Vector<C> > ATw_blocks(q);
    for (size_type l = 0; l < q; l++) {
      ATw_blocks[l].resize(n);
      A.apply_transposed(w_blocks[l], ATw_blocks[l]);
    }
    
    // resort the ATw[l] blocks into the result
    for (size_type j = 0; j < n; j++)
      for (size_type l = 0; l < q; l++)
 	Mtx[j*q+l] = ATw_blocks[l][j];
#else
    // old code, resorting to single entries
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
#endif
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
