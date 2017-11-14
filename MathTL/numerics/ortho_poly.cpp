// implementation for ortho_poly.h

#include <cassert>
#include <cmath>
#include <iostream>
#include <algebra/triangular_matrix.h>
#include <algebra/symmetric_matrix.h>
#include <numerics/matrix_decomp.h>

namespace MathTL
{
  inline
  double OrthogonalPolynomial::operator () (const unsigned int n, const double x) const
  {
    double pk(1.0); // p_0(x)
    
    if (n > 0) {
      double pkminus1(pk), pkminus2(0.0);
      for (unsigned int k(1); k <= n; k++) {
	pk = (x-a(k))*pkminus1 - b(k)*pkminus2;
	
	if (k < n) {
	  pkminus2 = pkminus1;
	  pkminus1 = pk;
	}
      }
    }
    
    return pk;
  }
  
  inline
  Polynomial<double> OrthogonalPolynomial::assemble(const unsigned int n) const
  {
    Polynomial<double> pk(1.0); // p_0

    if (n > 0) {
      Polynomial<double> pkminus1(pk), pkminus2, q;
      q.set_coefficient(1, 1.0);
      for (unsigned int k(1); k <= n; k++) {
	q.set_coefficient(0, -a(k)); // q(x)=x-a(k)
	pk = q*pkminus1 - b(k)*pkminus2;
	
	if (k < n) {
	  pkminus2 = pkminus1;
	  pkminus1 = pk;
	}
      }
    }
    
    return pk;
  }

  inline
  double OrthogonalPolynomial::forwardSummation(const Vector<double>& coeffs, const double x) const
  {
    double S(0.0);

    if (coeffs.size() >= 1) {
      S = coeffs(0); // p_0(x)=1
      
      double pk(1.0), pkminus1(1.0), pkminus2(0.0);
      for (unsigned int k(1); k < coeffs.size(); k++) {
	pk = (x-a(k))*pkminus1 - b(k)*pkminus2;
	S += coeffs[k] * pk;
	
	if (k < coeffs.size()-1) {
	  pkminus2 = pkminus1;
	  pkminus1 = pk;
	}
      }
    }
    
    return S;
  }

  inline
  double OrthogonalPolynomial::adjointSummation(const Vector<double>& coeffs, const double x) const
  {
    assert(coeffs.size() >= 1);
    
    double uk(0), ukplus1(0), ukplus2(0);
    
    // calculate u_N,...,u_1
    for (unsigned int k(coeffs.size()-1); k >= 1; --k) {
      uk = (x-a(k+1))*ukplus1 - b(k+1)*ukplus2 + coeffs[k];
      if (k > 1) {
	ukplus2 = ukplus1;
	ukplus1 = uk;
      }
    }
    
    // calculate u_0 = -b_2*u_2 + \alpha_0
    double u0(-b(2)*ukplus1 + coeffs(0));
	
    return u0 + uk*(x-a(1));
  }

  inline
  double Monomial::a(const unsigned int k) const
  {
    return 0.0;
  }
  
  inline
  double Monomial::b(const unsigned int k) const
  {
    return 0.0;
  }

  inline
  double ChebyshevPolynomial::a(const unsigned int k) const
  {
    return 0.0;
  }

  inline
  double ChebyshevPolynomial::b(const unsigned int k) const
  {
    double r(0);

    if (k >= 2) {
      if (k == 2)
	r = 0.5;
      else
	r = 0.25;
    }
    
    return r;
  }

  inline
  double LegendrePolynomial::a(const unsigned int k) const
  {
    return 0.0;
  }

  inline
  double LegendrePolynomial::b(const unsigned int k) const
  {
    double r(0);

    if (k >= 2) {
      if (k == 2)
	r = 1.0/3.0;
      else
	r = (k-1.0)*(k-1.0)/((2*k-1.0)*(2*k-3.0));
    }
    
    return r;
  }

  inline
  GenMomentsPolynomial::GenMomentsPolynomial(const Array1D<double>& moments,
					     const double a, const double b,
					     const unsigned int N)
    : ak_(N), bk_(N)
  {
    assert(N >= 1);
    assert(moments.size() >= 2*N+1);
    
    // setup Hankel matrix of the moments
    SymmetricMatrix<double> GammaN(N+1);
    for (unsigned int alpha(0); alpha <= N; alpha++)
      for (unsigned int beta(alpha); beta <= N; beta++) // only setup upper triangular part
	GammaN(alpha, beta) = moments[alpha+beta]; // requires moments[0],...,moments[2*N]

    // perform Cholesky decomposition GammaN=L*L^T=R^T*R
    LowerTriangularMatrix<double> L;
    bool result = CholeskyDecomposition(GammaN, L);
    if (result) {
      // Due to Mysovskih, the polynomials
      //   q_j(x):=\sum_{l=0}^js_{l,j}x^l, 0<=j<=N
      // where the s_{m,n} are the entries of R^{-1}=(L^T)^{-1}
      // are already orthogonal w.r.t. the weighted scalar product.
      // By comparison of the coefficients at x^{k-1} and x^{k-2}, one can directly
      // compute the three-term recursion coefficients alpha_k, beta_k in
      //   p_k(t) = (t-alpha_k)p_{k-1}(t) - beta_kp_{k-2}(t),
      // where k=1..N.
      // (see documentation)

      ak_[0] = L(1,0)/L(0,0); // a_1
      bk_[0] = moments[0]; // b_1 (by convention, but it is almost never needed)
      
      for (unsigned int n(1); n <= N-1; n++) {
	ak_[n] = L(n+1,n)/L(n,n) - L(n,n-1)/L(n-1,n-1);
	if (n==1)
	  bk_[n] = L(2,0)/L(0,0) - pow(L(1,0)/L(0,0), 2);
	else
	  bk_[n] = (L(n-1,n-2)*L(n,n-1))/(L(n-2,n-2)*L(n-1,n-1))
	    - L(n,n-2)/L(n-2,n-2)
	    - pow(L(n,n-1)/L(n-1,n-1), 2)
	    + L(n+1,n-1)/L(n-1,n-1);
      }
    }
  }

  inline
  GenMomentsPolynomial::GenMomentsPolynomial(const Array1D<double>& moments,
					     const OrthogonalPolynomial& T,
					     const double a, const double b,
					     const unsigned int N)
    : ak_(N), bk_(N) 
  {
    assert(N >= 1);
    assert(moments.size() >= 2*N);
    
    // set up the matrix S of the products <T_m,T_n>_w columnwise,
    // thereby also computing the a_k and b_k
    LowerTriangularMatrix<double> S(2*N, 2*N);
    for (unsigned int m(0); m < 2*N; m++)
      S(m, 0) = moments[m];
    
    ak_[0] = T.a(1) + S(1,0)/S(0,0); // a_1
    bk_[0] = moments[0]; // b_1 (by convention, but it is almost never needed)
    
    for (unsigned int n(1); n <= N-1; n++) {
      for (unsigned int m(n); m <= 2*N-n-1; m++) {
	S(m,n) = T.b(m+1)*S(m-1,n-1)+(T.a(m+1)-ak_[n-1])*S(m,n-1)+S(m+1,n-1);
	if (n > 1)
	  S(m,n) -= bk_[n-1]*S(m,n-2);
      }
      
      bk_[n] = S(n,n)/S(n-1,n-1); // b_{n+1}
      ak_[n] = T.a(n+1) + S(n+1,n)/S(n,n)-S(n,n-1)/S(n-1,n-1); // a_{n+1}
    }
  }
  
  inline
  double GenMomentsPolynomial::a(const unsigned int k) const
  {
    assert(k >= 1 && k <= ak_.size());
    return ak_[k-1];
  }
  
  inline
  double GenMomentsPolynomial::b(const unsigned int k) const
  {
    assert(k >= 1 && k <= bk_.size());
    return bk_[k-1];
  }
}
