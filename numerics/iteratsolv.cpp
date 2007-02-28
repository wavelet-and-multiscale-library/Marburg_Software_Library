// implementation of iterative solvers

#include <numerics/preconditioner.h>

namespace MathTL
{
  template <class VECTOR, class MATRIX>
  void Richardson(const MATRIX &A, const VECTOR &b, VECTOR &xk, const double omega,
		  const double tol, const unsigned int maxiter, unsigned int& iterations)
  {
    // first residual
    VECTOR rk(A.row_dimension(), false);
    A.apply(xk, rk);
    rk -= b; // rk=A*xk-b
    double error = l2_norm(rk);

    for (iterations = 1; iterations <= maxiter && error > tol; iterations++)
      {
	xk.add(-omega, rk);

	A.apply(xk, rk);
	rk -= b;
	error = l2_norm(rk);
      }
  }

  template <class VECTOR, class MATRIX>
  void Jacobi(const MATRIX &A, const VECTOR &b, VECTOR &xk,
	      const double tol, const unsigned int maxiter, unsigned int& iterations)
  {
    const typename MATRIX::size_type n(A.column_dimension());
    
    // first residual
    VECTOR rk(A.row_dimension(), false);
    A.apply(xk, rk);
    rk -= b; // rk=A*xk-b
    double error = l2_norm(rk);

    for (iterations = 1; iterations <= maxiter && error > tol; iterations++)
      {
	for (typename VECTOR::size_type i(0); i < n; i++)
	  xk[i] -= rk[i] / A(i, i);

	A.apply(xk, rk);
	rk -= b;

	error = l2_norm(rk);
      }
  }

  template <class VECTOR, class MATRIX>
  void GaussSeidel(const MATRIX &A, const VECTOR &b, VECTOR &xk,
		   const double tol, const unsigned int maxiter, unsigned int& iterations)
  {
    const typename MATRIX::size_type n(A.column_dimension());

    VECTOR rk(A.row_dimension(), false);
    double error = 2*tol;
    
    for (iterations = 0; iterations <= maxiter && error > tol; iterations++)
      {
	for (typename MATRIX::size_type i(0); i < n; i++)
	  {
	    xk[i] = b[i];
	    for (typename MATRIX::size_type j(0); j < n; j++)
	      if (j != i)
		xk[i] -= A(i, j) * xk[j];
	    xk[i] /= A(i, i);
	  }
	
	A.apply(xk, rk);
	rk -= b;

	error = l2_norm(rk);
      }
  }

  template <class VECTOR, class MATRIX>
  void SOR(const MATRIX &A, const VECTOR &b, VECTOR &xk, const double omega,
	   const double tol, const unsigned int maxiter, unsigned int& iterations)
  {
    const typename MATRIX::size_type n(A.column_dimension());

    VECTOR rk(A.row_dimension(), false);
    double error = 2*tol;
    
    for (iterations = 0; iterations <= maxiter && error > tol; iterations++)
      {
	for (typename MATRIX::size_type i(0); i < n; i++)
	  {
	    double sigma(b[i]);
	    for (typename MATRIX::size_type j(0); j < n; j++)
	      if (j != i)
		sigma -= A(i, j) * xk[j];
	    sigma /= A(i, i);
	    xk[i] += omega * (sigma - xk[i]);
	  }
	
	A.apply(xk, rk);
	rk -= b;

	error = l2_norm(rk);
      }
  }

  template <class VECTOR, class MATRIX>
  void SSOR(const MATRIX &A, const VECTOR &b, VECTOR &xk, const double omega,
	   const double tol, const unsigned int maxiter, unsigned int& iterations)
  {
    const typename MATRIX::size_type n(A.column_dimension());
    
    VECTOR rk(A.row_dimension(), false), xkpuf(xk);
    double error = 2*tol;
    
    for (iterations = 0; iterations <= maxiter && error > tol; iterations++)
      {
	for (typename MATRIX::size_type i(0); i < n; i++)
	  {
	    double sigma(0);
	    for (typename MATRIX::size_type j(0); j+1 <= i; j++)
	      sigma += A(i, j) * xkpuf[j];
	    for (typename MATRIX::size_type j(i+1); j < n; j++)
	      sigma += A(i, j) * xk[j];
	    sigma = (b[i]-sigma) / A(i, i);
	    xkpuf[i] = xk[i] + omega * (sigma - xk[i]);
	  }

	for (typename MATRIX::size_type i(n); i >= 1; i--)
	  {
	    double sigma(0);
	    for (typename MATRIX::size_type j(0); j+1 <= i-1; j++)
	      sigma += A(i-1, j) * xkpuf[j];
	    for (typename MATRIX::size_type j(i); j < n; j++)
	      sigma += A(i-1, j) * xk[j];
	    sigma = (b[i-1]-sigma) / A(i-1, i-1);
	    xk[i-1] = xkpuf[i-1] + omega * (sigma - xkpuf[i-1]);
	  }
	
	A.apply(xk, rk);
	rk -= b;

	error = l2_norm(rk);
      }
  }

  template <class VECTOR, class MATRIX>
  inline
  bool CG(const MATRIX &A, const VECTOR &b, VECTOR &xk,
	  const double tol, const unsigned int maxiter, unsigned int& iterations)
  {
    IdentityPreconditioner<MATRIX,VECTOR> I(A);
    return PCG(A, b, I, xk, tol, maxiter, iterations);

//     iterations = 0;
    
//     VECTOR gk(A.row_dimension(), false),
//       dk(A.row_dimension(), false),
//       Adk(A.row_dimension(), false);
    
//     // first residual
//     A.apply(xk, gk);
//     gk -= b; // gk = A*xk-b
    
//     // first direction
//     dk = gk;
//     dk.scale(-1.0); // dk = -gk
    
//     double error = l2_norm(gk);
//     for (iterations = 1; error > tol && iterations <= maxiter; iterations++)
//       {
// 	const double denom = gk * gk;
// 	A.apply(dk, Adk);
// 	const double alpha = denom / (dk * Adk);
// 	xk.add(alpha, dk);
// 	gk.add(alpha, Adk);
// 	error = l2_norm(gk);
// 	const double beta = (gk * gk) / denom;
// 	dk *= beta; 
// 	dk -= gk;
//       }
    
//     return (iterations <= maxiter);
  }

  template <class VECTOR, class MATRIX, class PREC>
  bool PCG(const MATRIX &A, const VECTOR &b, const PREC& P, VECTOR &xk,
	   const double tol, const unsigned int maxiter, unsigned int& iterations)
  {
    // see: "Templates for the Solution of Linear Systems: Building Blocks for Iterative Methods"

    VECTOR rk(A.row_dimension(), false),
      zk(A.row_dimension(), false),
      pk(A.row_dimension(), false),
      Apk(A.row_dimension(), false);

    // first (negative) residual
    A.apply(xk, rk);
    rk.subtract(b);
    const double normr0 = l2_norm_sqr(rk);
    double normrk = normr0, rhok = 0, oldrhok = 0;
    for (iterations = 1; normrk/normr0 > tol*tol && iterations <= maxiter; iterations++)
    //for (iterations = 1; normrk > tol*tol && iterations <= maxiter; iterations++)
      {
	P.apply_preconditioner(rk, zk);
	rhok = rk * zk;

	if (iterations == 1) // TODO: shift this case upwards!
	  pk = zk;
	else
	  pk.sadd(rhok/oldrhok, zk);

	A.apply(pk, Apk);
	const double alpha = rhok/(pk*Apk);
	xk.add(-alpha,  pk);
	rk.add(-alpha, Apk);
	normrk = l2_norm_sqr(rk);

	//cout << "normrk = " << sqrt(normrk) << endl;

	oldrhok = rhok;
      }

    return (iterations <= maxiter);
  }
}
