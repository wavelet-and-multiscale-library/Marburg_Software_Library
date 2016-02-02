// implementation of recursion.h inline functions

#include <cassert>

namespace MathTL
{
  double Recursion::operator () (const unsigned int n) const
  {
    double pk(p0_);

    if (n > 0)
      {
	if (n == 1)
	  pk = p1_;
	else
	  {
	    double pkminus1(p1_), pkminus2(p0_);
	    for (unsigned int k(2); k <= n; k++)
	      {
		pk = a(k)*pkminus1 + b(k)*pkminus2;
		if (k < n)
		  {
		    pkminus2 = pkminus1;
		    pkminus1 = pk;
		  }
	      }
	  }
      }

    return pk;
  }

  double Recursion::forwardSummation(const Vector<double>& coeffs) const
  {
    assert(coeffs.size() >= 2);

    double S(coeffs[0] * p0_ + coeffs[1] * p1_);

    double pk, pkminus1(p1_), pkminus2(p0_);
    for (unsigned int k(2); k < coeffs.size(); k++)
      {
	pk = a(k)*pkminus1 + b(k)*pkminus2;
	S += coeffs[k] * pk;
	if (k < coeffs.size()-1)
	  {
	    pkminus2 = pkminus1;
	    pkminus1 = pk;
	  }
      }
    
    return S;
  }

  double Recursion::adjointSummation(const Vector<double>& coeffs) const
  {
    assert(coeffs.size() >= 1);

    double uk(0), ukplus1(0), ukplus2(0);
    
    // calculate u_N,...,u_1
    for (unsigned int k(coeffs.size()-1); k >= 1; --k)
      {
	uk = a(k+1)*ukplus1 + b(k+1)*ukplus2 + coeffs[k];
	if (k > 1)
	  {
	    ukplus2 = ukplus1;
	    ukplus1 = uk;
	  }
      }

    // calculate u_0 = b_2*u_2 + \alpha_0
    double u0(b(2)*ukplus1 + coeffs[0]);
	
    return u0*p0_ + uk*p1_;
  }
  
  JacobiRecursion::JacobiRecursion(const double a, const double b, const double x)
    : a_(a), b_(b), x_(x)
  {
    p0_ = 1.0;
    p1_ = (a_-b_)/2.0 + (1.0+(a_+b_)/2.0)*x_;
  }

  double JacobiRecursion::a(const unsigned int k) const
  {
    return (2*k+a_+b_-1)*(a_*a_-b_*b_+(2*k+a_+b_-2)*(2*k+a_+b_)*x_)/(2*k*(k+a_+b_)*(2*k+a_+b_-2));
  }
  
  double JacobiRecursion::b(const unsigned int k) const
  {
    return -(k+a_-1)*(k+b_-1)*(2*k+a_+b_)/(k*(k+a_+b_)*(2*k+a_+b_-2));
  }
}
