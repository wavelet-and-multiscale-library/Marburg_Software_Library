// implementation of MathTL::Polynomial inline functions

#include <cassert>
#include <numerics/gauss_data.h>
// #include <iostream>
// #include <iomanip>
// #include <vector>
// #include <cmath>
// #include <algorithm>

// using namespace std;

namespace MathTL
{
  template <class C>
  Polynomial<C>::Polynomial()
    : Vector<C>(1)
  {
  }

  template <class C>
  inline
  unsigned int Polynomial<C>::degree() const
  {
    return Vector<C>::size()-1;
  }

  template <class C>
  inline
  C Polynomial<C>::get_coefficient(const unsigned int k) const
  {
    return Vector<C>::operator [] (k);
  }
  
  template <class C>
  void Polynomial<C>::set_coefficient(const unsigned int k,
				      const C coeff)
  {
    if (k > degree())
      {
	Vector<C> help(*this);
	Vector<C>::resize(k+1);
	std::copy(help.begin(), help.end(), Vector<C>::begin());
      }

    Vector<C>::operator [] (k) = coeff;
  }


//   Polynomial::Polynomial(const Polynomial &p)
//   {
//     a.resize(p.getCoeffs()->size(), 0.0);
//     a = *p.getCoeffs();
//   }

//   Polynomial::Polynomial(const CoeffsType *coeffs)
//   {
//     a.resize(coeffs->size(), 0.0);
//     a = *coeffs;
//   }

//   Polynomial::Polynomial(const double c)
//   {
//     a.resize(1, c);
//   }

//   //
//   //
//   // default destructor

//   Polynomial::~Polynomial()
//   {
//   }


//   //
//   //
//   // Polynomial methods


//   const CoeffsType *Polynomial::getCoeffs() const
//   {
//     return &a;
//   }

//   void Polynomial::setCoeff(const int k, const double coeff)
//   {
//     assert(k >= 0);

//     if (k > deg() || deg() == 0)
//       {
// 	CoeffsType help(a);
// 	a.resize(k+1, 0.0);
// 	std::copy(&help[0], &help[help.size()], &a[0]);
//       }
  
//     a[k] = coeff;
//   }

//   void Polynomial::setCoeffs(const CoeffsType *coeffs)
//   {
//     a.resize(coeffs->size(), 0.0);
//     a = *coeffs;
//   }


//   Polynomial Polynomial::sdifferentiate() const
//   {
//     CoeffsType c(0.0, deg());

//     for (int n = 1; n <= deg(); n++)
//       c[n-1] = n * a[n];
  
//     Polynomial r(&c);

//     return r;
//   }

//   Polynomial Polynomial::sintegrate() const
//   {
//     CoeffsType c(0.0, deg()+2);

//     for (int n = 0; n <= deg(); n++)
//       c[n+1] = a[n] / (n + 1.0);

//     Polynomial r(&c);

//     return r;
//   }

//   double Polynomial::integrate(const double a,
// 			       const double b,
// 			       const bool quadrature) const
//   {
//     assert(a <= b);

//     double r = 0;

//     if(quadrature)
//       {
// 	int N = (int)ceil((deg()+1)/2.0); // we must ensure 2*N-1>=deg()

// 	for(int i = 1; i <= N; i++)
// 	  r += GaussWeights[N-1][i-1] * this->operator()((a+b + (b-a)*GaussPoints[N-1][i-1])/2.0);
      
// 	r *= b-a;
//       }
//     else
//       {
// 	Polynomial P(sintegrate());
// 	r =  P(b)-P(a);
//       }

//     return r;
//   }

//   void Polynomial::substituteIntoMe(const Polynomial &p)
//   {
//     if (deg() > 0)
//       {
// 	Polynomial q;
// 	double a0 = a[0]; // not changed by substitution

// 	Polynomial r;
// 	for (int expo = 1; expo <= deg(); expo++)
// 	  {
// 	    if (a[expo] != 0)
// 	      {
// 		q = p;
// 		for (int l = 2; l <= expo; l++)
// 		  q *= p;
// 		q *= a[expo];

// 		r += q;
// 	      }
// 	  }

// 	setCoeffs(r.getCoeffs());
// 	a[0] += a0;
//       }
//   }

//   Polynomial Polynomial::substituteInto(const Polynomial &p) const
//   {
//     Polynomial r(&a);
//     r.substituteIntoMe(p);

//     return r;
//   }

//   Polynomial Polynomial::power(const int k) const
//   {
//     assert(k >= 0);

//     Polynomial r;

//     if (k == 0)
//       r = 1;
//     else
//       {
// 	r = Polynomial(*this);
// 	for (int l=2; l <= k; l++)
// 	  r *= (*this);
//       }
  
//     return r;
//   }

//   //
//   //
//   // Polynomial operators

//   double Polynomial::operator () (const double x) const
//   {
//     double value = a[deg()];
//     for (int k = deg(); k > 0; k--)
//       value = value * x + a[k-1];
  
//     return value;
//   }

//   double Polynomial::derivative(const double x) const
//   {
//     double r = 0;

//     if (deg() > 0)
//       {
// 	CoeffsType help(0.0, deg());

// 	help[deg()-1] = a[deg()];
// 	for (int k = deg()-1; k > 0; k--)
// 	  help[k-1] = help[k] * x + a[k];

// 	r = help[deg()-1];
// 	for (int k = deg()-1; k > 0; k--)
// 	  r = r * x + help[k-1];
//       }

//     return r;
//   }

//   Polynomial& Polynomial::operator = (const Polynomial &p)
//   {
//     setCoeffs(p.getCoeffs());

//     return *this;
//   }

//   Polynomial& Polynomial::operator = (const double c)
//   {
//     a.resize(1, c);

//     return *this;
//   }

//   Polynomial& Polynomial::operator += (const Polynomial &p)
//   {
//     if (p.deg() > deg() || deg() == 0)
//       {
// 	CoeffsType help(0.0, p.deg()+1);

// 	for (int n(0); n <= deg(); n++)
// 	  help[n] = getCoeff(n);

// 	a.resize(help.size(), 0.0);
// 	a = help;
//       }

//     for (int n = 0; n <= p.deg(); n++)
//       a[n] += p.getCoeff(n);
  
//     return *this;
//   }

//   Polynomial& Polynomial::operator -= (const Polynomial &p)
//   {
//     if (p.deg() > deg() || deg() == 0)
//       {
// 	CoeffsType help(0.0, p.deg()+1);

// 	for (int n(0); n <= deg(); n++)
// 	  help[n] = getCoeff(n);

// 	a.resize(help.size(), 0.0);
// 	a = help;
//       }

//     for (int n = 0; n <= p.deg(); n++)
//       a[n] -= p.getCoeff(n);
  
//     return *this;
//   }

//   Polynomial& Polynomial::operator *= (const double c)
//   {
//     if (c == 0)
//       {
// 	a.resize(0);
//       }
//     else
//       for (int n = 0; n <= deg(); n++)
// 	a[n] *= c;
  
//     return *this;
//   }

//   Polynomial& Polynomial::operator *= (const Polynomial &p)
//   {
//     CoeffsType b(0.0, deg()+p.deg()+1);
//     const CoeffsType *coeffs = p.getCoeffs();

//     for (int n = 0; n <= deg(); n++)
//       for (int m = 0; m <= p.deg(); m++)
// 	b[n+m] += a[n] * (*coeffs)[m];

//     a.resize(b.size(), 0.0);
//     std::copy(&b[0], &b[b.size()], &a[0]);
  
//     return *this;
//   }

//   //
//   //
//   // (friend) operators

//   Polynomial operator + (const Polynomial &p, const Polynomial &q)
//   {
//     Polynomial r(p);
//     r += q;

//     return r;
//   }

//   Polynomial operator - (const Polynomial &p, const Polynomial &q)
//   {
//     Polynomial r(p);
//     r -= q;

//     return r;
//   }

//   Polynomial operator - (const Polynomial &p)
//   {
//     Polynomial r(p);
//     r *= -1;

//     return r;
//   }

//   Polynomial operator * (const double c, const Polynomial &p)
//   {
//     Polynomial r(p);
//     r *= c;

//     return r;
//   }

//   Polynomial operator * (const Polynomial &p, const Polynomial &q)
//   {
//     Polynomial r(p);
//     r *= q;

//     return r;
//   }

  template <class C>
  std::ostream& operator << (std::ostream &s, const Polynomial<C> &p)
  {
    double c;
    int oldprecision = s.precision();
    std::ios::fmtflags oldflags = s.flags();
    //  s.setf(ios::scientific, ios::floatfield);
    s.precision(12);

    bool first = true;
    for (unsigned int k = p.degree(); k >= 1; k--)
      {
	c = p.get_coefficient(k);
	if (c != 0)
	  {
	    if (!first)
	      {
		if (c >= 0)
		  s << "+";
	      }
	    if (fabs(c) != 1)
	      s << c;
	    else
	      if (c == -1)
		s << "-";
	    s << "x";
	    if (k > 1)
	      s << "^" << k;
	    first = false;
	  }
      }
  
    c = p.get_coefficient(0);

    if ((c!=0)||(p.degree()==0))
      {
	if (!first)
	  {
	    if (c >= 0)
	      s << "+";
	  }
	s << c;
      }
  
    s.setf(oldflags);
    s.precision(oldprecision);

    return s;
  }
}
