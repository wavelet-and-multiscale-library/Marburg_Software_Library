// implementation for ortho_poly.h

#include <cassert>
#include <iostream>
#include <algebra/triangular_matrix.h>

namespace MathTL
{
  double OrthogonalPolynomial::operator () (const unsigned int n, const double x) const
  {
    double pk(1.0); // p_0(x)

    if (n > 0)
      {
	double pkminus1(pk), pkminus2(0.0);
	for (unsigned int k(1); k <= n; k++)
	  {
	    pk = (x-a(k))*pkminus1 - b(k)*pkminus2;

	    if (k < n)
	      {
		pkminus2 = pkminus1;
		pkminus1 = pk;
	      }
	  }
      }

    return pk;
  }

  Polynomial<double> OrthogonalPolynomial::assemble(const unsigned int n) const
  {
    Polynomial<double> pk(1.0); // p_0

    if (n > 0)
      {
	Polynomial<double> pkminus1(pk), pkminus2, q;
	q.set_coefficient(1, 1.0);
	for (unsigned int k(1); k <= n; k++)
	  {
	    q.set_coefficient(0, -a(k)); // q(x)=x-a(k)
	    pk = q*pkminus1 - b(k)*pkminus2;

	    if (k < n)
	      {
		pkminus2 = pkminus1;
		pkminus1 = pk;
	      }
	  }
      }

    return pk;
  }

  double OrthogonalPolynomial::forwardSummation(const Vector<double>& coeffs, const double x) const
  {
    double S(0.0);

    if (coeffs.size() >= 1)
      {
	S = coeffs(0); // p_0(x)=1

	double pk(1.0), pkminus1(1.0), pkminus2(0.0);
	for (unsigned int k(1); k < coeffs.size(); k++)
	  {
	    pk = (x-a(k))*pkminus1 - b(k)*pkminus2;
	    S += coeffs[k] * pk;

	    if (k < coeffs.size()-1)
	      {
		pkminus2 = pkminus1;
		pkminus1 = pk;
	      }
	  }
      }
    
    return S;
  }

  double OrthogonalPolynomial::adjointSummation(const Vector<double>& coeffs, const double x) const
  {
    assert(coeffs.size() >= 1);
    
    double uk(0), ukplus1(0), ukplus2(0);
    
    // calculate u_N,...,u_1
    for (unsigned int k(coeffs.size()-1); k >= 1; --k)
      {
	uk = (x-a(k+1))*ukplus1 - b(k+1)*ukplus2 + coeffs[k];
	if (k > 1)
	  {
	    ukplus2 = ukplus1;
	    ukplus1 = uk;
	  }
      }

    // calculate u_0 = -b_2*u_2 + \alpha_0
    double u0(-b(2)*ukplus1 + coeffs(0));
	
    return u0 + uk*(x-a(1));
  }

  double Monomial::a(const unsigned int k) const
  {
    return 0.0;
  }
  
  double Monomial::b(const unsigned int k) const
  {
    return 0.0;
  }

  double ChebyshevPolynomial::a(const unsigned int k) const
  {
    return 0.0;
  }

  double ChebyshevPolynomial::b(const unsigned int k) const
  {
    double r(0);

    if (k >= 2)
      {
	if (k == 2)
	  r = 0.5;
	else
	  r = 0.25;
      }

    return r;
  }

  double LegendrePolynomial::a(const unsigned int k) const
  {
    return 0.0;
  }

  double LegendrePolynomial::b(const unsigned int k) const
  {
    double r(0);

    if (k >= 2)
      {
	if (k == 2)
	  r = 1.0/3.0;
	else
	  r = (k-1.0)*(k-1.0)/((2*k-1.0)*(2*k-3.0));
      }
    
    return r;
  }

  GenMomentsPolynomial::GenMomentsPolynomial(const Array1D<double>& moments,
					     const OrthogonalPolynomial& T,
					     const double a, const double b,
					     const unsigned int N,
					     GenMomentsPolynomial::GMPMethod method)
    : ak_(N), bk_(N)
  {
    switch(method)
      {
      case SackDonovan:
	{
	  assert(N >= 1);
	  assert(moments.size() >= 2*N);
	  
	  // set up the matrix S of the products <T_m,T_n>_w columnwise,
	  // thereby also computing the a_k and b_k
	  LowerTriangularMatrix<double> S(2*N, 2*N);
	  for (unsigned m(0); m < 2*N; m++)
	    S(m, 0) = moments[m];
	  
	  ak_[0] = T.a(1) + S(1,0)/S(0,0); // a_1
	  bk_[0] = 0.0; // b_1
	  
	  for (unsigned int n(1); n <= N-1; n++)
	    {
	      for (unsigned int m(n); m <= 2*N-n-1; m++)
		{
		  S(m,n) = T.b(m+1)*S(m-1,n-1)+(T.a(m+1)-ak_[n-1])*S(m,n-1)+S(m+1,n-1);
		  if (n > 1)
		    S(m,n) -= bk_[n-1]*S(m,n-2);
		}
	      
	      bk_[n] = S(n,n)/S(n-1,n-1); // b_{n+1}
	      ak_[n] = T.a(n+1) + S(n+1,n)/S(n,n)-S(n,n-1)/S(n-1,n-1); // a_{n+1}
	    }
	}
      default:
	break;
      }
  }

  double GenMomentsPolynomial::a(const unsigned int k) const
  {
    assert(k >= 1 && k <= ak_.size());
    return ak_[k-1];
  }

  double GenMomentsPolynomial::b(const unsigned int k) const
  {
    assert(k >= 1 && k <= bk_.size());
    return bk_[k-1];
  }
}
