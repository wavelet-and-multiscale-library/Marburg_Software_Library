// implementation of iterative solvers

namespace MathTL
{
  template <class VECTOR, class MATRIX>
  void Richardson(const MATRIX &A, const VECTOR &b, VECTOR &xk, const double omega,
		  const double tol, const int maxiter, int& iterations)
  {
    // first residual
    VECTOR rk;
    A.APPLY(xk, rk);
    rk -= b;
    double error = norm(rk);

    for (iterations = 1; iterations <= maxiter && error > tol; iterations++)
      {
	xk -= omega * rk;
	A.APPLY(xk, rk);
	rk -= b;
	error = norm(rk);
      }
  }

  template <class VECTOR, class MATRIX>
  void Jacobi(const MATRIX &A, const VECTOR &b, VECTOR &xk,
	      const double tol, const int maxiter, int& iterations)
  {
    const int n(A.coldim());
    
    // first residual
    VECTOR rk;
    A.APPLY(xk, rk);
    rk -= b;
    double error = norm(rk);

    for (iterations = 1; iterations <= maxiter && error > tol; iterations++)
      {
	for (int i(0); i < n; i++)
	  xk(i) -= rk(i) / A(i, i);

	A.APPLY(xk, rk);
	rk -= b;
	error = norm(rk);
      }
  }

  template <class VECTOR, class MATRIX>
  void GaussSeidel(const MATRIX &A, const VECTOR &b, VECTOR &xk,
		   const double tol, const int maxiter, int& iterations)
  {
    const int n(A.coldim());

    VECTOR rk;
    double error = 2*tol;
    
    for (iterations = 0; iterations <= maxiter && error > tol; iterations++)
      {
	for (int i(0); i < n; i++)
	  {
	    xk(i) = b(i);
	    for (int j(0); j < n; j++)
	      if (j != i)
		xk(i) -= A(i, j) * xk(j);
	    xk(i) /= A(i, i);
	  }
	
	A.APPLY(xk, rk);
	rk -= b;
	error = norm(rk);
      }
  }

  template <class VECTOR, class MATRIX>
  void SOR(const MATRIX &A, const VECTOR &b, VECTOR &xk, const double omega,
	   const double tol, const int maxiter, int& iterations)
  {
    const int n(A.coldim());

    VECTOR rk;
    double error = 2*tol;
    
    for (iterations = 0; iterations <= maxiter && error > tol; iterations++)
      {
	for (int i(0); i < n; i++)
	  {
	    double sigma(b(i));
	    for (int j(0); j < n; j++)
	      if (j != i)
		sigma -= A(i, j) * xk(j);
	    sigma /= A(i, i);
	    xk(i) += omega * (sigma - xk(i));
	  }
	
	A.APPLY(xk, rk);
	rk -= b;
	error = norm(rk);
      }
  }

  template <class VECTOR, class MATRIX>
  void SSOR(const MATRIX &A, const VECTOR &b, VECTOR &xk, const double omega,
	   const double tol, const int maxiter, int& iterations)
  {
    const int n(A.coldim());
    
    VECTOR rk, xkpuf(xk);
    double error = 2*tol;
    
    for (iterations = 0; iterations <= maxiter && error > tol; iterations++)
      {
	for (int i(0); i < n; i++)
	  {
	    double sigma(0);
	    for (int j(0); j <= i-1; j++)
	      sigma += A(i, j) * xkpuf(j);
	    for (int j(i+1); j < n; j++)
	      sigma += A(i, j) * xk(j);
	    sigma = (b(i)-sigma) / A(i, i);
	    xkpuf(i) = xk(i) + omega * (sigma - xk(i));
	  }
	
	for (int i(n-1); i >= 0; i--)
	  {
	    double sigma(0);
	    for (int j(0); j <= i-1; j++)
	      sigma += A(i, j) * xkpuf(j);
	    for (int j(i+1); j < n; j++)
	      sigma += A(i, j) * xk(j);
	    sigma = (b(i)-sigma) / A(i, i);
	    xk(i) = xkpuf(i) + omega * (sigma - xkpuf(i));
	  }
	
	A.APPLY(xk, rk);
	rk -= b;
	error = norm(rk);
      }
  }

  template <class VECTOR, class MATRIX>
  bool CG(const MATRIX &A, const VECTOR &b, VECTOR &xk,
	  const double tol, const int maxiter, int& iterations)
  {
    iterations = 0;
    
    VECTOR gk, dk, Adk;
    
    // first residual
    A.APPLY(xk, gk);
    gk -= b;
    
    // first direction
    dk = -gk;
    
    double error = norm(gk);
    for (iterations = 1; error > tol && iterations <= maxiter; iterations++)
      {
	const double denom = gk * gk;
	A.APPLY(dk, Adk);
	const double alpha = denom / (dk * Adk);
	xk += alpha * dk;
	gk += alpha * Adk;
	error = norm(gk);
	const double beta = (gk * gk) / denom;
	dk *= beta; 
	dk -= gk;
      }
    
    return (iterations <= maxiter);
  }
}
