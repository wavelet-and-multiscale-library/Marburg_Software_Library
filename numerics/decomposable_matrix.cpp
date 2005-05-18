// implementation of MathTL::DecomposableMatrix inline functions

#include <iomanip>
#include <sstream>

namespace MathTL
{
  template <class C>
  inline
  DecomposableMatrix<C>::DecomposableMatrix(const size_type n)
    : Matrix<C>(n), decomposed(none)
  {
  }

  template <class C>
  inline
  DecomposableMatrix<C>::DecomposableMatrix(const DecomposableMatrix<C>& M)
    : Matrix<C>(M), decomposed(M.decomposed)
  {
    switch(decomposed)
      {
      case QR:
	break;
      case LU:
	P = M.P;
	D = M.D;
	break;
      case none:
      default:
	break;
      }
  }

  template <class C>
  inline
  const typename DecomposableMatrix<C>::size_type
  DecomposableMatrix<C>::row_dimension() const
  {
    return Matrix<C>::row_dimension();
  }

  template <class C>
  inline
  const typename DecomposableMatrix<C>::size_type
  DecomposableMatrix<C>::column_dimension() const
  {
    return Matrix<C>::column_dimension();
  }

  template <class C>
  inline
  const typename DecomposableMatrix<C>::size_type
  DecomposableMatrix<C>::size() const
  {
    return Matrix<C>::size();
  }

  template <class C>
  inline
  void DecomposableMatrix<C>::print(std::ostream &os,
			const unsigned int tabwidth,
			const unsigned int precision) const
  {
    switch(decomposed)
      {
      case LU:
	break;
      case QR:
	break;
      case none:
      default:
	Matrix<C>::print(os, tabwidth, precision);
	break;
      }
  }

  template <class C>
  inline
  std::ostream& operator << (std::ostream& os, const DecomposableMatrix<C>& M)
  {
    M.print(os);
    return os;
  }
}
