#include <cmath>
#include <cassert>
#include <iostream>
#include <utils/tiny_tools.h>

using std::cout;
using std::endl;

namespace MathTL
{
  template <class C>
  bool CholeskyDecomposition(const SymmetricMatrix<C>& A,
			     LowerTriangularMatrix<C>& L)
  {
    typedef typename SymmetricMatrix<C>::size_type size_type;
    size_type n(A.row_dimension());
    bool r(true);

    L.resize(n, n);
    for (size_type j(0); j < n; j++)
      {
	double d(0.0);
	for (size_type k(0); k < j; k++)
	  {
	    C s(0);
	    for (size_type i(0); i < k; i++)
	      s += L(k, i) * L(j, i);
	    L(j, k) = s = (A(j, k) - s) / L(k, k);
	    d += s*s;
	    r = (r && (A(k, j) == A(j, k)));
	  }
	d = A(j, j) - d;
	r = (r && (d > 0.0));
	L(j, j) = sqrt(d > 0.0 ? d : 0.0);
      }

    return r;
  }

  //
  //
  // QU decomposition

  template <class C>
  QUDecomposition<C>::QUDecomposition(const Matrix<C>& A)
  {
    // initialization
    rowdim_ = A.row_dimension();
    coldim_ = A.column_dimension();
    QU_.resize(rowdim_,coldim_);
    typedef typename Matrix<C>::size_type size_type;
    for (size_type row(0); row < rowdim_; row++)
      for (size_type col(0); col < coldim_; col++)
	QU_(row,col) = A(row,col);
    Udiag_.resize(coldim_, false);

    // main loop
    for (size_type k(0); k < coldim_; k++)
      {
	// 2-norm of k-th column
	C nrm(0);
	for (size_type i(k); i < rowdim_; i++) nrm = hypot(nrm,QU_(i,k));
	if (nrm != 0.0)
	  {
	    // construct k-th Householder vector
            if (QU_(k,k) < 0) nrm = -nrm;
            for (size_type i(k); i < rowdim_; i++) QU_(i,k) /= nrm;
            QU_(k,k) += 1.0;

            // transformation of the remaining columns
            for (size_type j(k+1); j < coldim_; j++)
	      {
               C s(0);
               for (size_type i(k); i < rowdim_; i++)
		 s += QU_(i,k)*QU_(i,j);
               s = -s/QU_(k,k);
               for (size_type i(k); i < rowdim_; i++)
		 QU_(i,j) += s*QU_(i,k);
	      }
	  }
	Udiag_[k] = -nrm;
      }
  }

  template <class C>
  bool QUDecomposition<C>::hasFullRank() const
  {
    for (typename Matrix<C>::size_type j(0); j < coldim_; j++) 
      {
	if (Udiag_[j] == 0)
	  return false;
      }
    return true;
  }

  template <class C>
  void QUDecomposition<C>::getU(UpperTriangularMatrix<C>& U) const
  {
    U.resize(coldim_, coldim_);
    typedef typename Matrix<C>::size_type size_type;
    for (size_type i(0); i < coldim_; i++)
      for (size_type j(i); j < coldim_; j++)
	{
	  if (i < j)
	    U(i,j) = QU_(i,j);
	  else
	    U(i,j) = Udiag_[i];
	}
  }

  template <class C>
  void QUDecomposition<C>::getQ(Matrix<C>& Q) const
  {
    Q.resize(rowdim_,coldim_);
    typedef typename Matrix<C>::size_type size_type;
    for (size_type k(coldim_-1); k >= 0;)
      {
	for (size_type i(0); i < rowdim_; i++) Q(i,k) = 0.0;
	Q(k,k) = 1.0;
	for (size_type j(k); j < coldim_; j++)
	  {
	    if (QU_(k,k) != 0)
	      {
		C s(0);
		for (size_type i(k); i < rowdim_; i++)
		  s += QU_(i,k)*Q(i,j);
		s = -s/QU_(k,k);
		for (size_type i(k); i < rowdim_; i++)
		  Q(i,j) += s*QU_(i,k);
	      }
	  }

	if (k == 0) // unsigned!!!
	  break;
	else
	  --k;
      }
  }

  template <class C>
  void QUDecomposition<C>::solve(const Vector<C>& b, Vector<C>& x) const
  {
    x.resize(coldim_);
    Vector<C> Qtb(coldim_);

    // maybe one can speed up the following code by using QU_ directly
    UpperTriangularMatrix<double> U;
    getU(U);
    Matrix<double> Q;
    getQ(Q);

    typedef typename Matrix<C>::size_type size_type;
    for (size_type i(0); i < coldim_; i++)
      {
 	Qtb[i] = 0;
 	for (size_type j(0); j < rowdim_; j++)
	  Qtb[i] += Q(j, i) * b[j];
      }

    // backsubstitution
    for (size_type i(coldim_-1); i >= 0;)
      {
	x[i] = Qtb[i];
	for (size_type j(i+1); j < coldim_; j++)
	  x[i] -= U(i, j) * x[j];
	x[i] /= U(i, i);

	if (i == 0) // unsigned!!!
	  break;
	else
	  --i;
      }   
  }

  template <class C>
  void QUDecomposition<C>::inverse(Matrix<C>& AInv) const
  {
    assert(coldim_ == rowdim_);
    AInv.resize(coldim_, coldim_);

    // maybe one can speed up the following code by using QU_ directly
    UpperTriangularMatrix<double> U;
    getU(U);
    Matrix<double> Q;
    getQ(Q);

    // special backsubstitution
    typedef typename Matrix<C>::size_type size_type;
    for (size_type i(0); i < coldim_; i++)
      {
	AInv(coldim_-1, i) = Q(i, coldim_-1) / U(coldim_-1, coldim_-1);
	for (size_type j(coldim_-2); j >= 0;)
	  {
	    double c(Q(i, j));
	    for (size_type k(j+1); k <= coldim_-1; k++)
	      c -= U(j, k) * AInv(k, i);
	    AInv(j, i) = c / U(j, j);
	    
	    if (j == 0)
	      break;
	    else
	      --j;
	  }
      }

    AInv.compress();
  }

  //
  //
  // SVD

  template <class C>
  SVD<C>::SVD(const Matrix<C>& A)
  {
    // the following code is essentially taken from the SYNAPS svd routine

    assert(rowdim_>1 && coldim_>1);

    typedef typename Matrix<C>::size_type size_type;

    size_type i, j, k, EstColRank, RotCount, SweepCount, slimit;
    C eps(1e-10), e2, tol, vt, p, x0, y0, q, r, c0, s0=0, d1, d2;
    rowdim_ = A.row_dimension();
    coldim_ = A.column_dimension();
    assert(rowdim_ >= coldim_);
    slimit = coldim_/4;

    // working copy of A
    U_.resize(rowdim_,coldim_);
    for (size_type row(0); row < rowdim_; row++)
      for (size_type col(0); col < coldim_; col++)
 	U_(row,col) = A(row,col);
    
    Sdiag_.resize(coldim_);
    V_.resize(coldim_,coldim_);

    if (slimit < 6.0)
      slimit = 6;
    SweepCount = 0;
    e2 = C(10.0)*C(rowdim_)*eps*eps;
    tol = eps*C(.1);
    EstColRank = coldim_;
    for (i=0; i<coldim_; i++)
      {
	for (j=0; j<coldim_; j++)
	  V_(j,i) = C(0.0);
	V_(i,i) = C(1.0);
      }

    RotCount = EstColRank*(EstColRank-1)/2;
    while (RotCount != 0 && SweepCount <= slimit)
      {
	RotCount = EstColRank*(EstColRank-1)/2;
	SweepCount++;
	for (j=0; j<EstColRank-1; j++)
	  {
	    for (k=j+1; k<EstColRank; k++)
	      {
		p = q = r = C(0.0);
		for (i=0; i<rowdim_; i++)
		  {
		    x0 = U_(i,j); y0 = U_(i,k);
		    p += x0*y0; q += x0*x0; r += y0*y0;
		  }
		Sdiag_(j) = q; Sdiag_(k) = r;
		if (q >= r)
		  {
		    if (q<=e2*Sdiag_(0) || (abs(p)<=(tol*q))) RotCount--;
		    else
		      {
			p /= q; r = -r/q+1; vt = sqrt(p*p*4+r*r);
			c0 = sqrt(abs((r/vt+1)*C(.5))); 
			s0 = p/(vt*c0);
			for (i=0; i<rowdim_; i++)
			  {
			    d1 = U_(i,j); d2 = U_(i,k);
			    U_(i,j) = d1*c0+d2*s0; U_(i,k) = -d1*s0+d2*c0; 
			  }
			for (i=0; i<coldim_; i++)        
			  {
			    d1 = V_(j,i); d2 = V_(k,i);
			    V_(j,i) = d1*c0+d2*s0; 
			    V_(k,i) = -d1*s0+d2*c0; 
			  }
		      }
		  }
		else
		  {
		    p /= r; q = q/r-C(1); vt = sqrt(p*p*4+q*q);
		    s0 = sqrt(abs(C(.5)*(-q/vt+1)));
		    if (p<0) s0 = -s0;
		    c0 = p/(vt*s0);
		    for (i=0; i<rowdim_; i++)
		      {
			d1 = U_(i,j); d2 = U_(i,k);
			U_(i,j) = d1*c0+d2*s0; U_(i,k) = -d1*s0+d2*c0;
		      }
		    for (i=0; i<coldim_; i++)
		      {
			d1 = V_(j,i); d2 = V_(k,i);
			V_(j,i) = d1*c0+d2*s0; V_(k,i) = -d1*s0+d2*c0;
		      }
		  }
	      }
	  }
	while (EstColRank>=3 && Sdiag_(EstColRank-1)<=Sdiag_(0)*tol+tol*tol)
	  EstColRank--;
      }

    for(i=0; i<Sdiag_.size(); i++) Sdiag_(i) = sqrt(Sdiag_(i));
  }

  template <class C>
  void SVD<C>::getUS(Matrix<C>& US) const
  {
    US.resize(rowdim_,coldim_);
    US = U_;
  }

  template <class C>
  void SVD<C>::getU(Matrix<C>& U) const
  {
    U.resize(rowdim_,coldim_);
    typedef typename Matrix<C>::size_type size_type;
    for (size_type i(0); i < rowdim_; i++)
      for (size_type j(0); j < coldim_; j++)
	U(i,j) = U_(i,j) / Sdiag_[j];
  }

  template <class C>
  void SVD<C>::getV(Matrix<C>& V) const
  {
    V.resize(rowdim_,coldim_);
    V = V_;
  }

  template <class C>
  Hessenberg<C>::Hessenberg(const Matrix<C>& A)
    : Matrix<C>(A)
  {
    assert(A.row_dimension() == A.column_dimension());

    // TODO: implement this!

    cout << "after Hessenberg():" << endl << *this;
  }
}
