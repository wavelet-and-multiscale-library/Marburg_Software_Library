// implementation of LaurentPolynomial inline functions

#include <algorithm>

namespace MathTL
{
  template <class R>
  LaurentPolynomial<R>::LaurentPolynomial()
    : std::map<int,R>(), Function<1,R>()
  {
  }

  template <class R>
  LaurentPolynomial<R>::LaurentPolynomial(const LaurentPolynomial<R>& p)
    : std::map<int,R>(p)
  {
  }

  template <class R>
  LaurentPolynomial<R>::LaurentPolynomial(const R value)
    : std::map<int,R>(), Function<1,R>()
  {
    this->operator [] (0) = value;
  }

  template <class R>
  LaurentPolynomial<R>::~LaurentPolynomial()
  {
  }

  template <class R>
  LaurentPolynomial<R>& LaurentPolynomial<R>::operator = (const LaurentPolynomial<R>& p)
  {
    std::map<int,R>::operator = (p);

    return *this;
  }

  template <class R>
  inline
  R LaurentPolynomial<R>::get_coefficient(const int k) const
  {
    // we must not use std::map<int,R>::operator [] for reading,
    // since it may add unwanted zero elements!
    typename std::map<int,R>::const_iterator ir(lower_bound(k));
    if (it != end() && !key_comp()(k,it->first))
      return it->second;

    return R(0);
  }

  template <class R>
  inline
  void LaurentPolynomial<R>::set_coefficient(const int k,
					     const R coeff)
  {
    std::map<int,R>::operator [] (k) = coeff;
  }

  template <class R>
  R LaurentPolynomial<R>::value(const R x) const
  {
//     C value(Vector<C>::operator [] (degree()));

//     for (unsigned int k(degree()); k > 0; k--)
//       value = value * x + get_coefficient(k-1);

//     return value;
    return 0.0;
  }

  template <class R>
  inline
  R LaurentPolynomial<R>::value(const Point<1>& p,
				const unsigned int component) const
  {
    return value(p[0]);
  }

  template <class R>
  inline
  void LaurentPolynomial<R>::vector_value(const Point<1> &p,
					  Vector<R>& values) const
  {
    values.resize(1, false);
    values[0] = value(p[0]);
  }

//   template <class C>
//   std::ostream& operator << (std::ostream &s, const Polynomial<C> &p)
//   {
//     double c;
//     int oldprecision = s.precision();
//     std::ios::fmtflags oldflags = s.flags();
//     //  s.setf(ios::scientific, ios::floatfield);
//     s.precision(12);

//     bool first = true;
//     for (unsigned int k = p.degree(); k >= 1; k--)
//       {
// 	c = p.get_coefficient(k);
// 	if (c != 0)
// 	  {
// 	    if (!first)
// 	      {
// 		if (c >= 0)
// 		  s << "+";
// 	      }
// 	    if (fabs(c) != 1)
// 	      s << c;
// 	    else
// 	      if (c == -1)
// 		s << "-";
// 	    s << "x";
// 	    if (k > 1)
// 	      s << "^" << k;
// 	    first = false;
// 	  }
//       }
  
//     c = p.get_coefficient(0);

//     if ((c!=0)||(p.degree()==0))
//       {
// 	if (!first)
// 	  {
// 	    if (c >= 0)
// 	      s << "+";
// 	  }
// 	s << c;
//       }
  
//     s.setf(oldflags);
//     s.precision(oldprecision);

//     return s;
//   }

  template <class R>
  std::ostream& operator << (std::ostream& s, const LaurentPolynomial<R>& p)
  {
    for (typename LaurentPolynomial<R>::const_iterator it(p.begin());
	 it != p.end(); ++it)
      {
 	if (it != p.begin())
 	  s << "+";
 	s << "(" << it->second << ")";
 	if (it->first != 0)
 	  {
 	    s << "z";
 	    if (it->first != 1)
 	      {
 		s << "^";
 		if (it->first < 0)
 		  s << "{" << it->first << "}";
 		else
 		  s << it->first;
 	      }
 	  }
      }
    
    return s;
  }
}
