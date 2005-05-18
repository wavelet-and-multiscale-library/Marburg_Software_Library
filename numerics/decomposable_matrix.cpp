// implementation of MathTL::DecomposableMatrix inline functions

#include <iomanip>
#include <sstream>

using std::cout;
using std::endl;

namespace MathTL
{
  template <class C>
  inline
  DecomposableMatrix<C>::DecomposableMatrix(const size_type n)
    : Matrix<C>(n), decomposition(none)
  {
  }

  template <class C>
  inline
  DecomposableMatrix<C>::DecomposableMatrix(const DecomposableMatrix<C>& M)
    : Matrix<C>(M), decomposition(M.decomposition)
  {
    switch(decomposition)
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
  DecomposableMatrix<C>::DecomposableMatrix(const Matrix<C>& M)
    : Matrix<C>(M), decomposition(none)
  {
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
  void DecomposableMatrix<C>::decompose(DecompositionType d)
  {
    if (decomposition != none)
      revert_decomposition();

    switch(d)
      {
      case LU:
	LU_decomposition();
	break;
      case QR:
	QR_decomposition();
	break;
      case none:
      default:
	break;
      }
  }

  template <class C>
  void DecomposableMatrix<C>::LU_decomposition()
  {
    // row equilibration
    D.resize(row_dimension());
    for (size_type i(0); i < row_dimension(); i++)
      {
 	double rowsum(0);
 	for (size_type j(0); j < column_dimension(); j++)
 	  rowsum += fabs(Matrix<C>::get_entry(i, j));

  	for (size_type j(0); j < column_dimension(); j++)
  	  Matrix<C>::set_entry(i, j, Matrix<C>::get_entry(i, j) / rowsum);

  	D[i] = rowsum;
      }

    P.resize(row_dimension());
    for (size_type i(0); i < row_dimension(); i++)
      P[i] = i;

    for (size_type k(0); k < column_dimension(); k++)
      {
	// row pivoting
	size_type prow(k);
	C pivot(Matrix<C>::get_entry(P[k], k));
	for (size_type i(k+1); i < row_dimension(); i++)
	  {
	    C help(Matrix<C>::get_entry(P[i], k));
	    if (fabs(help) > fabs(pivot))
	      {
		prow = i;
		pivot = help;
	      }   
	  }
	std::swap(P[k], P[prow]);

	// forward elimination
	for (size_type i(k+1); i < row_dimension(); i++)
	  {
	    double lik(Matrix<C>::get_entry(P[i], k) / pivot);
	    Matrix<C>::set_entry(P[i], k, lik);
	    for (size_type j(k+1); j < column_dimension(); j++)
	      Matrix<C>::set_entry(P[i], j, Matrix<C>::get_entry(P[i], j) - lik * Matrix<C>::get_entry(P[k], j));
	  }
      }

    decomposition = LU;
  }

  template <class C>
  void DecomposableMatrix<C>::QR_decomposition()
  {
  }

  template <class C>
  void DecomposableMatrix<C>::revert_decomposition()
  {
    Matrix<C> M(row_dimension(), column_dimension());

    switch(decomposition)
      {
      case LU:
	for (size_type j(0); j < column_dimension(); j++)
	  {
	    for (size_type i(0); i < row_dimension(); i++)
	      {
		C mij(0);
		for (size_type k(0); k <= i && k <= j; k++)
		  {
		    mij += (k == i && i < column_dimension())
		      ?                                 Matrix<C>::get_entry(P[k], j)
		      : Matrix<C>::get_entry(P[i], k) * Matrix<C>::get_entry(P[k], j);
		  }
		M.set_entry(P[i], j, mij * D[P[i]]);
	      }
	  }
	Matrix<C>::swap(M);
	break;
      case QR:
	break;
      case none:
      default:
	break;
      }

    decomposition = none;
  }

  template <class C>
  inline
  void DecomposableMatrix<C>::print(std::ostream &os,
			const unsigned int tabwidth,
			const unsigned int precision) const
  {
    switch(decomposition)
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
