#include <cassert>
#include <cmath>
#include <numerics/iteratsolv.h>
#include <algebra/atra.h>

namespace MathTL
{
  template <class VECTOR, class MATRIX>
  double PowerIteration(const MATRIX &A, VECTOR &xk,
			const double tol, const int maxit, int &iterations)
  {
    double lambdak(0), error = 2.0 * tol;
    
    // start with xk=(1,...,1)^T
    xk.resize(A.rowdim());
    xk = 1;

    VECTOR yk, diff;
    A.APPLY(xk, yk);

    for (iterations = 1; iterations < maxit && error > tol; iterations++)
      {
 	xk = 1./norm(yk)*yk;
 	A.APPLY(xk, yk);
 	lambdak = xk * yk;
 	diff = yk - lambdak * xk;
 	error = norm(diff) / fabs(lambdak);
      }

    xk /= norm(xk, 1.0);

    return lambdak;
  }
  
  template <class VECTOR, class MATRIX>
  double InversePowerIteration(const MATRIX &A, VECTOR &xk,
 			       const double tol, const int maxit, int &iterations)
  {
    double lambdak(0), error = 2.0 * tol;
    
    xk = VECTOR(A.rowdim());
    xk = 1;
    VECTOR yk(A.coldim());
    yk = 1; // initial value for CG iteration
    
    VECTOR zk, diff;
    int it, CGits;
    
    CG(A, xk, yk, tol, 200, CGits);
    
    for (it = 1; it < maxit && (error > tol); it++)
      {
 	xk = 1./norm(yk)*yk;
 	CG(A, xk, yk, tol/100.0, 200, CGits);
 	lambdak = xk * yk;
 	diff = yk - lambdak * xk;
 	error = norm(diff) / fabs(lambdak);
      }
    
    xk /= norm(xk, 1.0);

    return lambdak;
  }

  template <class VECTOR, class MATRIX, class MATRIX2>
  void SymmEigenvalues(const MATRIX& A, VECTOR& evals, MATRIX2& evecs)
  {
    // the following code stems (at the moment) from JAMA
    assert(A.rowdim() == A.coldim());
    int n = A.rowdim();

    evals.resize(n);
    VECTOR ehelp(n);

    // use evecs as working copy of A
    evecs.resize(n,n);
    for (int i(0); i < n; i++)
      for (int j(0); j < n; j++)
	evecs(i,j) = A(i,j);
    
    // transform A to tridiagonal form via symmetric Householder reduction
    for (int j = 0; j < n; j++)
      evals(j) = evecs(n-1,j);
    
    for (int i = n-1; i > 0; i--)
      {
	// Scale to avoid under/overflow.
	double scale = 0.0;
	double h = 0.0;
	for (int k = 0; k < i; k++) {
	  scale = scale + abs(evals(k));
	}
	if (scale == 0.0)
	  {
	    ehelp(i) = evals(i-1);
	    for (int j = 0; j < i; j++) {
	      evals(j) = evecs(i-1,j);
	      evecs(i,j) = 0.0;
	      evecs(j,i) = 0.0;
	    }
	  }
	else
	  {
	    // Generate Householder vector.
	    for (int k = 0; k < i; k++) {
	      evals(k) /= scale;
	      h += evals(k) * evals(k);
	    }
	    double f = evals(i-1);
	    double g = sqrt(h);
	    if (f > 0) {
	      g = -g;
	    }
	    ehelp(i) = scale * g;
	    h = h - f * g;
	    evals(i-1) = f - g;
	    for (int j = 0; j < i; j++) {
	      ehelp(j) = 0.0;
	    }

	// Apply similarity transformation to remaining columns.
	    for (int j = 0; j < i; j++) {
	      f = evals(j);
	      evecs(j,i) = f;
	      g = ehelp(j) + evecs(j,j) * f;
	      for (int k = j+1; k <= i-1; k++) {
		g += evecs(k,j) * evals(k);
		ehelp(k) += evecs(k,j) * f;
	      }
	      ehelp(j) = g;
	    }
	    
	    f = 0.0;
	    for (int j = 0; j < i; j++) {
	      ehelp(j) /= h;
	      f += ehelp(j) * evals(j);
	    }
	    
	    double hh = f / (h+h);
	    for (int j = 0; j < i; j++) {
	      ehelp(j) -= hh*evals(j);
	    }
	    for (int j = 0; j < i; j++) {
	      f = evals(j);
	      g = ehelp(j);
	      for (int k = j; k <= i-1; k++) {
		evecs(k,j) -= (f*ehelp(k) + g*evals(k));
	      }
	      evals(j) = evecs(i-1,j);
	      evecs(i,j) = 0.0;
	    }
	  }
	evals(i) = h;
      }
    
    // Accumulate transformations.
    for (int i = 0; i < n-1; i++) {
      evecs(n-1,i) = evecs(i,i);
      evecs(i,i) = 1.0;
      double h = evals(i+1);
      if (h != 0.0) {
	for (int k = 0; k <= i; k++) {
	  evals(k) = evecs(k,i+1)/h;
	}
	for (int j = 0; j <= i; j++) {
	  double g = 0.0;
	  for (int k = 0; k <= i; k++) {
	    g += evecs(k,i+1) * evecs(k,j);
	  }
	  for (int k = 0; k <= i; k++) {
	    evecs(k,j) -= g * evals(k);
	  }
	}
      }
      for (int k = 0; k <= i; k++) {
	evecs(k,i+1) = 0.0;
      }
    }
    for (int j = 0; j < n; j++) {
      evals(j) = evecs(n-1,j);
      evecs(n-1,j) = 0.0;
    }
    evecs(n-1,n-1) = 1.0;
    ehelp(0) = 0.0;
    
    // diagonalization
    for (int i = 1; i < n; i++) {
      ehelp(i-1) = ehelp(i);
    }
    ehelp(n-1) = 0.0;
    
    double f = 0.0;
    double tst1 = 0.0;
    double eps = pow(2.0,-52.0);
    for (int l = 0; l < n; l++) {
      // Find small subdiagonal element
      
      tst1 = std::max(tst1,abs(evals(l)) + abs(ehelp(l)));
      int m = l;
      
      // Original while-loop from Java code
      while (m < n) {
	if (abs(ehelp(m)) <= eps*tst1) {
	  break;
	}
	m++;
      }
      
      // If m == l, d[l] is an eigenvalue,
      // otherwise, iterate.
      
      if (m > l) {
	int iter = 0;
	do {
	  iter = iter + 1;  // (Could check iteration count here.)
	  
	  // Compute implicit shift
	  
	  double g = evals(l);
	  double p = (evals(l+1) - g) / (2.0 * ehelp(l));
	  double r = hypot(p,1.0);
	  if (p < 0) {
	    r = -r;
	  }
	  evals(l) = ehelp(l) / (p + r);
	  evals(l+1) = ehelp(l) * (p + r);
	  double dl1 = evals(l+1);
	  double h = g - evals(l);
	  for (int i = l+2; i < n; i++) {
	    evals(i) -= h;
	  }
	  f = f + h;
	  
	  // Implicit QL transformation.
	  p = evals(m);
	  double c = 1.0;
	  double c2 = c;
	  double c3 = c;
	  double el1 = ehelp(l+1);
	  double s = 0.0;
	  double s2 = 0.0;
	  for (int i = m-1; i >= l; i--) {
	    c3 = c2;
	    c2 = c;
	    s2 = s;
	    g = c * ehelp(i);
	    h = c * p;
	    r = hypot(p,ehelp(i));
	    ehelp(i+1) = s * r;
	    s = ehelp(i) / r;
	    c = p / r;
	    p = c * evals(i) - s * g;
	    evals(i+1) = h + s * (c * g + s * evals(i));
	    
	    // Accumulate transformation.
	    for (int k = 0; k < n; k++) {
	      h = evecs(k,i+1);
	      evecs(k,i+1) = s * evecs(k,i) + c * h;
	      evecs(k,i) = c * evecs(k,i) - s * h;
	    }
	  }
	  p = -s * s2 * c3 * el1 * ehelp(l) / dl1;
	  ehelp(l) = s * p;
	  evals(l) = c * p;
	  
	  // Check for convergence.
	  
	} while (abs(ehelp(l)) > eps*tst1);
      }
      evals(l) += f;
      ehelp(l) = 0.0;
    }

    // Sort eigenvalues and corresponding vectors.
    for (int i = 0; i < n-1; i++) {
      int k = i;
      double p = evals(i);
      for (int j = i+1; j < n; j++) {
 	if (evals(j) < p) {
 	  k = j;
 	  p = evals(j);
 	}
      }
      if (k != i) {
 	evals(k) = evals(i);
 	evals(i) = p;
 	for (int j = 0; j < n; j++) {
 	  p = evecs(j,i);
 	  evecs(j,i) = evecs(j,k);
 	  evecs(j,k) = p;
 	}
      }
    }
  } 

  template <class MATRIX>
  double CondSymm(const MATRIX& A, double tol, int maxit)
  {
    RawVector<double> x;
    int iterations;
    
    return fabs(PowerIteration(A, x, tol, maxit, iterations)
		* InversePowerIteration(A, x, tol, maxit, iterations));
  }

  template <class MATRIX>
  double CondNonSymm(const MATRIX& A, double tol, int maxit)
  {
    
    return sqrt(CondSymm(AtrA<MATRIX>(A), tol, maxit));
  }

}
