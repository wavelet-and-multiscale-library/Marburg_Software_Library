// implementation of MathTL::Polynomial inline functions

#include <cassert>
#include <cmath>
#include <numerics/gauss_data.h>
#include <utils/tiny_tools.h>

namespace MathTL
{
  template <class C>
  Polynomial<C>::Polynomial()
    : Vector<C>(1), Function<1,C>()
  {
  }
  
  template <class C>
  Polynomial<C>::Polynomial(const Polynomial<C>& p)
    : Vector<C>(p)
  {
  }

  template <class C>
  Polynomial<C>::Polynomial(const Vector<C>& coeffs)
    : Vector<C>(coeffs)
  {
  }

  template <class C>
  Polynomial<C>::Polynomial(const C c)
    : Vector<C>(1)
  {
    set_coefficient(0, c);
  }

  template <class C>
  Polynomial<C>::~Polynomial()
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
  inline
  void Polynomial<C>::get_coefficients(Vector<C>& coeffs) const
  {
    coeffs = *this;
  }
  
  template <class C>
  void Polynomial<C>::set_coefficient(const unsigned int k,
				      const C coeff)
  {
    if (k > degree())
      {
	Vector<C> help;
	help.swap(*this);
	Vector<C>::resize(k+1);
	std::copy(help.begin(), help.end(), Vector<C>::begin());
      }

    Vector<C>::operator [] (k) = coeff;

    trim();
  }

  template <class C>
  void Polynomial<C>::set_coefficients(const Vector<C>& coeffs)
  {
    assert(coeffs.size() > 0);
    Vector<C>::operator = (coeffs);

    trim();
  }

  template <class C>
  C Polynomial<C>::value(const C x) const
  {
    C value(Vector<C>::operator [] (degree()));

    for (unsigned int k(degree()); k > 0; k--)
      value = value * x + get_coefficient(k-1);

    return value;
  }

  template <class C>
  C Polynomial<C>::value(const C x, const unsigned int derivative) const
  {
    // storage for full Horner scheme
    Vector<C> horner(*this);

    for (unsigned row(0); row <= derivative; row++)
      {
	for (unsigned n(degree()); n > row; n--)
	  horner[n-1] += x*horner[n];
      }

    return horner[derivative] * faculty(derivative);
  }

  template <class C>
  inline
  C Polynomial<C>::value(const Point<1>& p,
			 const unsigned int component) const
  {
    return value(p[0]);
  }

  template <class C>
  inline
  void Polynomial<C>::vector_value(const Point<1> &p,
				   Vector<C>& values) const
  {
    values.resize(1, false);
    values[0] = value(p[0]);
  }

  template <class C>
  void Polynomial<C>::scale(const C s)
  {
    C factor(1.0);
    for (typename Vector<C>::iterator it(Vector<C>::begin()), itend(Vector<C>::end());
	 it != itend; ++it)
      {
	*it *= factor;
	factor *= s;
      }

    trim(); // TODO: check s=0
  }
  
  template <class C>
  void Polynomial<C>::shift(const C s)
  {
    Vector<C> new_coeffs(*this);
    
    for (unsigned int d(1); d < new_coeffs.size(); d++)
      {
	unsigned int n(d);
	unsigned int binomial(1);
	C s_power(s);
	for (unsigned int k(0); k < d; k++)
	  {
	    binomial = (binomial*(n-k))/(k+1);
	    new_coeffs[d-k-1] +=
	      new_coeffs[d] * binomial * s_power;
	    s_power *= s;
	  }
      }
    
    swap(new_coeffs);

    // no trimming necessary
  }

  template <class C>
  void Polynomial<C>::chain(const Polynomial<C>& p)
  {
    if (degree() > 0)
      {
	// maybe the following can be optimized

  	Polynomial<C> q;
  	C a0(get_coefficient(0)); // not changed by substitution

  	Polynomial<C> r;
  	for (unsigned int expo = 1; expo <= degree(); expo++)
  	  {
  	    if (get_coefficient(expo) != 0)
  	      {
  		q = p;
  		for (unsigned int l = 2; l <= expo; l++)
  		  q *= p;
  		q *= get_coefficient(expo);

  		r += q;
  	      }
  	  }

 	Vector<C>::swap(r);
 	set_coefficient(0, get_coefficient(0) + a0);
      }

    // no trimming necessary
  }

  template <class C>
  Polynomial<C>& Polynomial<C>::operator = (const Polynomial<C>& p)
  {
    Vector<C>::operator = (p);
    trim();
    return *this;
  }

  template <class C>
  Polynomial<C>& Polynomial<C>::operator = (const C c)
  {
    Vector<C>::resize(1, false);
    set_coefficient(0, c);
    return *this; // no trimming necessary
  }

  template <class C>
  void Polynomial<C>::add(const Polynomial<C>& p)
  {
    Vector<C> help(std::max(degree(), p.degree())+1);
    if (degree()+1 == help.size())
      {
	std::copy(p.begin(), p.end(), help.begin());
	help.add(*this);
      }
    else
      {
	std::copy(Vector<C>::begin(), Vector<C>::end(), help.begin());
	help.add(p);
      }
    swap(help);

    trim();
  }

  template <class C>
  void Polynomial<C>::add(const C s, const Polynomial<C>& p)
  {
    Vector<C> help(std::max(degree(), p.degree())+1);
    if (degree()+1 == help.size())
      {
	std::copy(p.begin(), p.end(), help.begin());
	help.sadd(s, *this); // help <- s*help + *this
      }
    else
      {
	std::copy(Vector<C>::begin(), Vector<C>::end(), help.begin());
	help.add(s, p); // help <- help + s*p
      }
    swap(help);

    trim();
  }

  template <class C>
  void Polynomial<C>::sadd(const C s, const Polynomial<C>& p)
  {
    Vector<C> help(std::max(degree(), p.degree())+1);
    if (degree()+1 == help.size())
      {
	std::copy(p.begin(), p.end(), help.begin());
	help.add(s, *this); // help <- help + s*(*this)
      }
    else
      {
	std::copy(Vector<C>::begin(), Vector<C>::end(), help.begin());
	help.sadd(s, p); // help <- s*help + p
      }
    swap(help);

    trim();
  }

  template <class C>
  inline
  Polynomial<C>& Polynomial<C>::operator += (const Polynomial<C>& p)
  {
    add(p); // trims
    return *this;
  }

  template <class C>
  inline
  Polynomial<C> Polynomial<C>::operator + (const Polynomial<C>& p) const
  {
    return (Polynomial<C>(*this) += p); // trims
  }

  template <class C>
  void Polynomial<C>::subtract(const Polynomial<C>& p)
  {
    Vector<C> help(std::max(degree(), p.degree())+1);
    if (degree()+1 == help.size())
      {
	std::copy(p.begin(), p.end(), help.begin());
	help.sadd(C(-1), *this);
      }
    else
      {
	std::copy(Vector<C>::begin(), Vector<C>::end(), help.begin());
	help.add(C(-1), p);
      }
    swap(help);

    trim();
  }

  template <class C>
  inline
  Polynomial<C>& Polynomial<C>::operator -= (const Polynomial<C>& p)
  {
    subtract(p); // trims
    return *this;
  }

  template <class C>
  inline
  Polynomial<C> Polynomial<C>::operator - () const
  {
    return (Polynomial<C>() -= *this); // trims
  }

  template <class C>
  inline
  Polynomial<C> Polynomial<C>::operator - (const Polynomial<C>& p) const
  {
    return (Polynomial<C>(*this) -= p);
  }

  template <class C>
  Polynomial<C>& Polynomial<C>::operator *= (const C c)
  {
    if (c == C(0))
      Vector<C>::resize(1); // p(x)=0, degree changes
    else
      Vector<C>::operator *= (c);

    return *this;
  }

  template <class C>
  inline
  Polynomial<C> Polynomial<C>::operator * (const C c) const
  {
    return (Polynomial(*this) *= c);
  }

  template <class C>
  void
  Polynomial<C>::multiply(const Polynomial<C>& p)
  {
    Vector<C> coeffs(degree()+p.degree()+1);
    for (unsigned int n(0); n <= degree(); n++)
      for (unsigned int m(0); m <= p.degree(); m++)
	coeffs[n+m] += Vector<C>::operator [] (n) * p.get_coefficient(m);

    swap(coeffs);

    trim(); // for safety
  }

  template <class C>
  inline
  Polynomial<C>& Polynomial<C>::operator *= (const Polynomial<C>& p)
  {
    multiply(p); // trims
    return *this;
  }

  template <class C>
  inline
  Polynomial<C> Polynomial<C>::operator * (const Polynomial<C>& p)
  {
    return (Polynomial(*this) *= p); // trims
  }

  template <class C>
  Polynomial<C> Polynomial<C>::power(const unsigned int k) const
  {
    Polynomial<C> r(1);

    for (unsigned int i(0); i < k; i++)
      r.multiply(*this);

    return r;
  }

  template <class C>
  void Polynomial<C>::divide(const Polynomial<C>& q,
			     Polynomial<C>& p, Polynomial<C>& r) const
  {
    assert(degree() >= q.degree());

    r = *this;
    p = 0;
    for (int k(r.degree() - q.degree());
	 k >= 0; k = r.degree() - q.degree())
      {
 	double factor = r.get_coefficient(r.degree()) / q.get_coefficient(q.degree());
 	p.set_coefficient(k, factor);
 	Polynomial monom;
 	monom.set_coefficient(k, factor);
 	r.subtract(monom * q);
      }
  }

  template <class C>
  void Polynomial<C>::divide(const Polynomial<C>& q, Polynomial<C>& p) const
  {
    Polynomial<C> r;
    divide(q, p, r);
  }

  template <class C>
  Polynomial<C> Polynomial<C>::differentiate() const
  {
    Vector<C> coeffs(std::max(degree(),(unsigned int)1));

    for (unsigned int n(1); n <= degree(); n++)
      coeffs[n-1] = n * get_coefficient(n);

    return Polynomial<C>(coeffs);
  }

  template <class C>
  Polynomial<C> Polynomial<C>::integrate() const
  {
    Vector<C> coeffs(degree()+2);
    
    for (unsigned int n(0); n <= degree(); n++)
      coeffs[n+1] = get_coefficient(n) / (n + 1.0);

    return Polynomial<C>(coeffs);
  }

  template <class C>
  double Polynomial<C>::integrate(const double a,
				  const double b,
				  const bool quadrature) const
  {
    assert(a <= b);

    double r = 0;

    if (quadrature)
      {
	const unsigned int N = (unsigned int)ceil((degree()+1)/2.0); // ensures 2*N-1>=degree()

 	for(unsigned int i(0); i < N; i++)
    	  r += GaussWeights[N-1][i] * value((a+b + (b-a)*GaussPoints[N-1][i])/2.0);
	
	r *= b-a;
      }
    else
      {
	Polynomial<C> P(integrate());
	r = P.value(b)-P.value(a);
      }

    return r;
  }

  template <class C>
  Polynomial<C> operator * (const C c, const Polynomial<C>& p)
  {
    return (Polynomial<C>(p) *= c);
  }

  template <class C>
  void Polynomial<C>::trim()
  {
    // determine "true" degree
    unsigned int deg(0);
    for(unsigned int i(1); i <= degree(); i++)
      {
	if (get_coefficient(i) != 0)
	  deg = i;
      } 

    if (deg < degree())
      {
	Vector<C> trimmed(std::max(1u, deg+1));
	for (unsigned int i(0); i < trimmed.size(); i++)
	  trimmed[i] = get_coefficient(i);
	swap(trimmed);
      }
  }

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
