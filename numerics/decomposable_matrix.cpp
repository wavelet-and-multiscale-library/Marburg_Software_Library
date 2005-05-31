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
      case QU:
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
    if (decomposition != d || decomposition != none)
      revert_decomposition();

    switch(d)
      {
      case LU:
	LU_decomposition();
	break;
      case QU:
	QU_decomposition();
	break;
      case none:
      default:
	break;
      }
  }

  template <class C>
  void DecomposableMatrix<C>::solve(const Vector<C>& b, Vector<C>& x,
				    DecompositionType d)
  {
    assert(d != none);

    if (d != decomposition)
      decompose(d);
    
    Vector<C> z(column_dimension());

    switch(d)
      {
      case LU:
	assert(row_dimension() == column_dimension());
	
	// PDAx=PDb=LUx, so we solve Lz=PDb first:
	for (size_type i(0); i < column_dimension(); i++)
	  {
	    z[i] = b[P[i]] / D[P[i]];
	    for (size_type j(0); j < i; j++)
	      z[i] -= Matrix<C>::get_entry(P[i], j) * z[j];
	  }

	// solve Ux=z
 	x.resize(column_dimension());
 	for (size_type i(column_dimension()-1);;)
 	  {
	    double xi = z[i];
	    for (size_type j(i+1); j < column_dimension(); j++)
	      xi -= Matrix<C>::get_entry(P[i], j) * x[j];
	    x[i] = xi / Matrix<C>::get_entry(P[i], i);

	    if (i>0)
	      i--;
	    else
	      break;
 	  }

	break;
      case QU:
	break;
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
  void DecomposableMatrix<C>::QU_decomposition()
  {
    for (size_type k(0); k < column_dimension(); k++)
      {
	double sqrnorm(0);
	for (size_type i(k); i < row_dimension(); i++)
	  sqrnorm += Matrix<C>::get_entry(i, k) * Matrix<C>::get_entry(i, k);

	double akk(Matrix<C>::get_entry(k, k));
	double alphak = (akk == 0 ? sqrt(sqrnorm) : -akk/fabs(akk)*sqrt(sqrnorm));

	Matrix<C>::set_entry(k, k, alphak);
	for (size_type i(k+1); i < row_dimension(); i++)
	  Matrix<C>::set_entry(i, k, Matrix<C>::get_entry(i, k) / (akk-alphak));

	double vTv(1); // v_1=1 by construction
	for (size_type i(k+1); i < row_dimension(); i++)
	  vTv += Matrix<C>::get_entry(i, k) * Matrix<C>::get_entry(i, k);

	// apply Q_v = I - 2*v*v^T/(v^Tv) to the other columns
	for (size_type j(k+1); j < column_dimension(); j++)
	  {
	    double vTx(Matrix<C>::get_entry(k, j)); // v_1=1
	    for (size_type i(k+1); i < row_dimension(); i++)
	      vTx += Matrix<C>::get_entry(i, k) * Matrix<C>::get_entry(i, j);

	    Matrix<C>::set_entry(k, j, Matrix<C>::get_entry(k, j) - 2*vTx/vTv); // v_1=1
	    for (size_type i(k+1); i < row_dimension(); i++)
	      Matrix<C>::set_entry(i, j, Matrix<C>::get_entry(i, j) - 2*vTx/vTv*Matrix<C>::get_entry(i, k));
	  }
      }

    decomposition = QU;
  }

  template <class C>
  void DecomposableMatrix<C>::revert_decomposition()
  {
    Matrix<C> M(row_dimension(), column_dimension());
    Vector<C> v(row_dimension(), false);

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
      case QU:
	for (size_type k(column_dimension()-1);;)
	  {
 	    // copy Householder vector v
 	    v[0] = 1;
 	    for (size_type i(k+1); i < row_dimension(); i++)
 	      v[i-k] = Matrix<C>::get_entry(i, k);

 	    double vTv(1);
 	    for (size_type i(k+1); i < row_dimension(); i++)
 	      vTv += v[i-k]*v[i-k];

	    for (size_type j(k); j < column_dimension(); j++)
	      {
		double vTx(Matrix<C>::get_entry(k, j)); // v_1=1
		if (j > k)
		  for (size_type i(k+1); i < row_dimension(); i++)
		    vTx += v[i-k] * Matrix<C>::get_entry(i, j);

		for (size_type i(k); i < row_dimension(); i++)
		  {
		    if (j == k && i > j ) // a_{i,j}=0
		      Matrix<C>::set_entry(i, j, -2*vTx/vTv*v[i-k]);
		    else // j > k || i <= j
		      Matrix<C>::set_entry(i, j, Matrix<C>::get_entry(i, j) - 2*vTx/vTv*v[i-k]);
		  }
	      }

	    if (k > 0)
	      k--;
	    else
	      break;
	  }
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
      case QU:
// 	break;
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
