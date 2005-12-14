#include <cassert>
#include <cmath>
#include <algebra/vector.h>
#include <algebra/atra.h>
#include <algebra/shifted_matrix.h>
#include <algebra/symmetric_matrix.h>
#include <numerics/iteratsolv.h>
#include <utils/tiny_tools.h>
#include <utils/random.h>

namespace MathTL
{
  template <class VECTOR, class MATRIX>
  double PowerIteration(const MATRIX &A, VECTOR &xk,
			const double tol, const unsigned int maxit, unsigned int &iterations)
  {
    double lambdak(0), error = 2.0 * tol;

    assert(xk.size() == A.row_dimension());

    // start with xk=(1,...,1)^T
//     xk.resize(A.row_dimension(), false);
//     xk = 1;

    VECTOR yk(xk.size(), false), diff(xk.size(), false);
    A.apply(xk, yk);

    for (iterations = 1; iterations < maxit && error > tol; iterations++)
      {
 	xk = yk; xk.scale(1./l2_norm(yk));
 	A.apply(xk, yk);
 	lambdak = xk * yk;
 	diff = yk; diff.add(-lambdak, xk);

 	error = l2_norm(diff) / abs(lambdak);
      }

    xk.scale(1.0/l1_norm(xk));

    return lambdak;
  }
  
  template <class VECTOR, class MATRIX>
  inline
  double InversePowerIteration(const MATRIX &A, VECTOR &xk,
 			       const double tol, const unsigned int maxit,
			       unsigned int &iterations)
  {
    return InversePowerIteration(A, xk, 0, tol, maxit, iterations);
  }

  template <class VECTOR, class MATRIX>
  double InversePowerIteration(const MATRIX& A, VECTOR& zk, const double lambda,
			       const double tol, const unsigned int maxit, unsigned int &iterations)
  {
    double lambdak(lambda), error = 2.0 * tol, diff;
    
    assert(zk.size() == A.row_dimension());

    ShiftedMatrix<MATRIX> ATilde(A, lambda);
    VECTOR qk(zk.size(), false);

    for (iterations = 1; iterations < maxit && error > tol; iterations++)
      {
	unsigned int ik, CGits;
	diff = lambdak;
	qk = zk; qk.scale(1./linfty_norm(zk, ik)); // ||qk||_infty=1
//  	if (iterations > 3) ATilde.set_lambda(lambdak);
 	CG(ATilde, qk, zk, tol, 200, CGits); // be cautious in the non-symmetric case...
//   	lambdak += qk[ik]/zk[ik];
  	lambdak = lambda + qk[ik]/zk[ik];
	diff -= lambdak;
	error = fabs(diff);
// 	cout << "# inverse power iteration, after iteration " << iterations
// 	     << ", lambdak=" << lambdak << ", difference " << error << endl;
      }

    zk.scale(1./linfty_norm(zk));

    return lambdak;
  }

  template <class VECTOR, class MATRIX, class MATRIX2>
  void SymmEigenvalues(const MATRIX& A, VECTOR& evals, MATRIX2& evecs)
  {
    // the following code stems (at the moment) from JAMA
    assert(A.row_dimension() == A.column_dimension());
    typedef typename MATRIX::size_type size_type;
    size_type n = A.row_dimension();

    evals.resize(n, false);
    VECTOR ehelp(n);

    // use evecs as working copy of A
    evecs.resize(n,n);
    for (size_type i(0); i < n; i++)
      for (size_type j(0); j < n; j++)
	evecs(i,j) = A(i,j);
    
    // transform A to tridiagonal form via symmetric Householder reduction
    for (size_type j(0); j < n; j++)
      evals[j] = evecs(n-1,j);
    
    for (size_type i(n-1); i > 0;)
      {
	// Scale to avoid under/overflow.
	double scale = 0.0;
	double h = 0.0;
	for (size_type k(0); k < i; k++) {
	  scale = scale + abs(evals[k]);
	}
	if (scale == 0.0)
	  {
	    ehelp[i] = evals[i-1];
	    for (size_type j(0); j < i; j++) {
	      evals[j] = evecs(i-1,j);
	      evecs(i,j) = 0.0;
	      evecs(j,i) = 0.0;
	    }
	  }
	else
	  {
	    // Generate Householder vector.
	    for (size_type k(0); k < i; k++) {
	      evals[k] /= scale;
	      h += evals[k] * evals[k];
	    }
	    double f = evals[i-1];
	    double g = sqrt(h);
	    if (f > 0) {
	      g = -g;
	    }
	    ehelp[i] = scale * g;
	    h = h - f * g;
	    evals[i-1] = f - g;
	    for (size_type j = 0; j < i; j++) {
	      ehelp[j] = 0.0;
	    }
	    
	    // Apply similarity transformation to remaining columns.
	    for (size_type j(0); j < i; j++) {
	      f = evals[j];
	      evecs(j,i) = f;
	      g = ehelp[j] + evecs(j,j) * f;
	      for (size_type k(j+1); k <= i-1; k++) {
		g += evecs(k,j) * evals[k];
		ehelp[k] += evecs(k,j) * f;
	      }
	      ehelp[j] = g;
	    }
	    
	    f = 0.0;
	    for (size_type j(0); j < i; j++) {
	      ehelp[j] /= h;
	      f += ehelp[j] * evals[j];
	    }
	    
	    double hh = f / (h+h);
	    for (size_type j(0); j < i; j++) {
	      ehelp[j] -= hh*evals[j];
	    }
	    for (size_type j(0); j < i; j++) {
	      f = evals[j];
	      g = ehelp[j];
	      for (size_type k(j); k <= i-1; k++) {
		evecs(k,j) -= (f*ehelp[k] + g*evals[k]);
	      }
	      evals[j] = evecs(i-1,j);
	      evecs(i,j) = 0.0;
	    }
	  }
	evals[i] = h;

	if (i > 0) --i;
      }
    
    // Accumulate transformations.
    for (size_type i(0); i < n-1; i++) {
      evecs(n-1,i) = evecs(i,i);
      evecs(i,i) = 1.0;
      double h = evals(i+1);
      if (h != 0.0) {
	for (size_type k(0); k <= i; k++) {
	  evals[k] = evecs(k,i+1)/h;
	}
	for (size_type j(0); j <= i; j++) {
	  double g = 0.0;
	  for (size_type k(0); k <= i; k++) {
	    g += evecs(k,i+1) * evecs(k,j);
	  }
	  for (size_type k(0); k <= i; k++) {
	    evecs(k,j) -= g * evals[k];
	  }
	}
      }
      for (size_type k(0); k <= i; k++) {
	evecs(k,i+1) = 0.0;
      }
    }
    for (size_type j(0); j < n; j++) {
      evals[j] = evecs(n-1,j);
      evecs(n-1,j) = 0.0;
    }
    evecs(n-1,n-1) = 1.0;
    ehelp[0] = 0.0;
    
    // diagonalization
    for (size_type i(1); i < n; i++) {
      ehelp[i-1] = ehelp[i];
    }
    ehelp[n-1] = 0.0;
    
    double f = 0.0;
    double tst1 = 0.0;
    double eps = pow(2.0,-52.0);
    for (size_type l(0); l < n; l++) {
      // Find small subdiagonal element
      
      tst1 = std::max(tst1, abs(evals[l])+abs(ehelp[l]));
      size_type m(l);
      
      // Original while-loop from Java code
      while (m < n) {
 	if (abs(ehelp[m]) <= eps*tst1) {
 	  break;
 	}
 	m++;
      }
      
      // If m == l, d[l] is an eigenvalue,
      // otherwise, iterate.
      
      if (m > l) {
 	int iter = 0;
 	do {
 	  iter++;  // (Could check iteration count here.)
	  
 	  // Compute implicit shift
	  
 	  double g = evals[l];
 	  double p = (evals[l+1] - g) / (2.0 * ehelp[l]);
 	  double r = hypot(p,1.0);
 	  if (p < 0) {
 	    r = -r;
 	  }
	  evals[l] = ehelp[l] / (p + r);
	  evals[l+1] = ehelp[l] * (p + r);
	  double dl1 = evals[l+1];
	  double h = g - evals[l];
	  for (size_type i(l+2); i < n; i++) {
	    evals[i] -= h;
	  }
	  f = f + h;
	  
	  // Implicit QL transformation.
	  p = evals[m];
	  double c = 1.0;
	  double c2 = c;
	  double c3 = c;
	  double el1 = ehelp[l+1];
	  double s = 0.0;
	  double s2 = 0.0;
	  for (size_type i(m-1); i >= l;) {
	    c3 = c2;
	    c2 = c;
	    s2 = s;
	    g = c * ehelp[i];
	    h = c * p;
	    r = hypot(p,ehelp[i]);
	    ehelp[i+1] = s * r;
	    s = ehelp[i] / r;
	    c = p / r;
	    p = c * evals[i] - s * g;
	    evals[i+1] = h + s * (c * g + s * evals[i]);
	    
	    // Accumulate transformation.
	    for (size_type k(0); k < n; k++) {
	      h = evecs(k,i+1);
	      evecs(k,i+1) = s * evecs(k,i) + c * h;
	      evecs(k,i) = c * evecs(k,i) - s * h;
	    }

	    if (i > l)
	      --i;
	    else
	      break;
	  }
	  p = -s * s2 * c3 * el1 * ehelp[l] / dl1;
	  ehelp[l] = s * p;
	  evals[l] = c * p;
	  
	  // Check for convergence.
	  
	} while (abs(ehelp[l]) > eps*tst1);
      }
      evals[l] += f;
      ehelp[l] = 0.0;
    }

    // Sort eigenvalues and corresponding vectors.
    for (size_type i(0); i < n-1; i++) {
      size_type k(i);
      double p = evals[i];
      for (size_type j(i+1); j < n; j++) {
 	if (evals[j] < p) {
 	  k = j;
 	  p = evals[j];
 	}
      }
      if (k != i) {
 	evals[k] = evals[i];
 	evals[i] = p;
 	for (size_type j(0); j < n; j++) {
 	  p = evecs(j,i);
 	  evecs(j,i) = evecs(j,k);
 	  evecs(j,k) = p;
 	}
      }
    }
  } 

  template <class MATRIX>
  double CondSymm(const MATRIX& A, double tol, unsigned int maxit)
  {
    Vector<double> x(A.row_dimension(), false);
    x = 1;
    unsigned int iterations;
    
    return abs(PowerIteration(A, x, tol, maxit, iterations)
	       * InversePowerIteration(A, x, tol, maxit, iterations));
  }

  template <class MATRIX>
  double CondNonSymm(const MATRIX& A, double tol, unsigned int maxit)
  {
    return sqrt(CondSymm(AtrA<MATRIX>(A), tol, maxit));
  }

  template <class MATRIX>
  void LanczosIteration(const MATRIX& A, const double tol,
			double& lambdamin, double& lambdamax,
			const unsigned int maxit, unsigned int& k)
  {
    assert(A.row_dimension() == A.column_dimension());
    const unsigned int n = A.row_dimension();

    Array1D<double> alpha(maxit), gamma(maxit+1);

    Vector<double> dk(n); // start vector d^{(0)}
    for (unsigned int i(0); i < dk.size(); i++) dk[i] = random_double();
//     cout << "d^{(0)}=" << dk << endl;
    gamma[0] = l2_norm(dk);
//     cout << "gamma_{0}=" << gamma[0] << endl;
    Vector<double> qk(n), qkold(n);
//     cout << "q^{(0)}=" << qk << endl;

    double change = 2*tol;

    for (k = 1; k <= maxit && fabs(gamma[k-1])>1e-14 && change > tol; k++) {
#if _MATHTL_LANCZOS_VERBOSITY >= 1
      cout << "Lanczos iteration, k=" << k << endl;
#endif
      qkold = qk; qk = dk; qk.scale(1./gamma[k-1]);
//       cout << "q^{(" << k << ")}=" << qk << endl;
      A.apply(qk, dk); // d^{(k)}=Aq^{(k)}
      alpha[k-1] = qk * dk;
//       cout << "alpha_{" << k-1 << "}=" << alpha[k-1] << endl;
      dk.add(-alpha[k-1], qk);
      dk.add(-gamma[k-1], qkold);
//       cout << "d^{(" << k << ")}=" << dk << endl;
      gamma[k] = l2_norm(dk);
//       cout << "gamma_{" << k << "}=" << gamma[k] << endl;

      SymmetricMatrix<double> M(k);
      for (unsigned int i = 0; i < k; i++) {
 	M.set_entry(i, i, alpha[i]);
 	if (i > 0)   M.set_entry(i, i-1, gamma[i]);
 	if (i < k-1) M.set_entry(i, i+1, gamma[i+1]);
      }

//       cout << "Mk=" << endl;
//       cout << M;

      Matrix<double> evecs;
      Vector<double> evals;
      SymmEigenvalues(M, evals, evecs);

//       cout << "eigenvalues of M_{" << k << "}: " << evals << endl;

      if (k > 2)
	change = std::max(fabs((lambdamin-evals[0])/lambdamin),
			  fabs((lambdamax-evals[k-1])/lambdamax));
      
      if (k >= 2) {
	lambdamin = evals[0];
	lambdamax = evals[k-1];
      }
    }
  }
}
