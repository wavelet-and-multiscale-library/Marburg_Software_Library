// implementation of MathTL::*TriangularMatrix inline functions

#include <algorithm>
#include <cassert>

namespace MathTL
{
  template <class C>
  inline
  SymmetricMatrix<C>::SymmetricMatrix(const size_type n)
    : LowerTriangularMatrix<C>(n)
  {
  }
  
  template <class C>
  inline
  SymmetricMatrix<C>::SymmetricMatrix(const SymmetricMatrix<C>& M)
    : LowerTriangularMatrix<C>(M)
  {
  }

  template <class C>
  SymmetricMatrix<C>::SymmetricMatrix(const size_type n,
				      const char* str,
				      const bool byrow)
    : LowerTriangularMatrix<C>(n,n,str,byrow)
  {
  }

  template <class C>
  inline
  void SymmetricMatrix<C>::resize(const size_type n)
  {
    LowerTriangularMatrix<C>::resize(n);
  }

  template <class C>
  inline
  const C SymmetricMatrix<C>::operator () (const size_type row,
					   const size_type column) const
  {
    return LowerTriangularMatrix<C>::entries_[triangle_index(row,column,
							     LowerTriangularMatrix<C>::rowdim_,
							     LowerTriangularMatrix<C>::coldim_)];
  }
  
  template <class C>
  inline
  C& SymmetricMatrix<C>::operator () (const size_type row,
				      const size_type column)
  {
    return LowerTriangularMatrix<C>::entries_[triangle_index(row,column,
							     LowerTriangularMatrix<C>::rowdim_,
							     LowerTriangularMatrix<C>::coldim_)];    
  }

  template <class C>
  template <class C2>
  inline
  bool SymmetricMatrix<C>::operator == (const SymmetricMatrix<C2>& M) const
  {
    return LowerTriangularMatrix<C>::operator == (M);
  }
  
  template <class C>
  template <class C2>
  inline
  bool SymmetricMatrix<C>::operator != (const SymmetricMatrix<C2>& M) const
  {
    return !((*this) == M);
  }
  
  template <class C>
  SymmetricMatrix<C>&
  SymmetricMatrix<C>::operator = (const SymmetricMatrix<C>& M)
  {
    LowerTriangularMatrix<C>::operator = (M);
    return *this;
  }

  template <class C>
  template <class VECTOR>
  void SymmetricMatrix<C>::apply(const VECTOR& x, VECTOR& Mx) const
  {
    assert(Mx.size() == LowerTriangularMatrix<C>::rowdim_);
    
    for (typename SymmetricMatrix<C>::size_type i(0);
	 i < LowerTriangularMatrix<C>::rowdim_; i++)
      {
	Mx[i] = 0;
	for (typename SymmetricMatrix<C>::size_type j(0);
	     j < LowerTriangularMatrix<C>::coldim_; j++)
	  Mx[i] += this->operator () (i, j) * x[j];
      }
  }

  template <class C>
  template <class VECTOR>
  void SymmetricMatrix<C>::apply_transposed(const VECTOR& x, VECTOR& Mtx) const
  {
    apply(x, Mtx);
  }

  template <class C>
  void SymmetricMatrix<C>::print(std::ostream &os,
				 const unsigned int tabwidth,
				 const unsigned int precision) const
  {
    if (LowerTriangularMatrix<C>::empty())
      os << "[]" << std::endl; // Matlab style
    else
      {
	unsigned int old_precision = os.precision(precision);
	for (typename SymmetricMatrix<C>::size_type i(0);
	     i < SymmetricMatrix<C>::row_dimension(); ++i)
	  {
	    for (typename SymmetricMatrix<C>::size_type j(0);
		 j < SymmetricMatrix<C>::column_dimension(); ++j)
	      os << std::setw(tabwidth) << std::setprecision(precision)
		 << this->operator () (i, j);
	    os << std::endl;
	  }
	os.precision(old_precision);
      }
  }

  template <class C>
  inline
  std::ostream& operator << (std::ostream& os, const SymmetricMatrix<C>& M)
  {
    M.print(os);
    return os;
  }
}
