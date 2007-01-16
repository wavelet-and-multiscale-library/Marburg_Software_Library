// implementation of MathTL::Piecewise inline functions

#include <algebra/piecewise.h>
#include <algebra/polynomial.h>
#include <cmath>
#include <iostream>
//#include <iomanip>

namespace MathTL
{
  /*! default constructor */
  template <class C>
  Piecewise<C>::Piecewise()
    :  Function<1,C>()
  {
    expansion.clear();
    granularity = 0;
  }


  /*! copy constructor */
  template <class C>
  Piecewise<C>::Piecewise(const Piecewise<C> &p)
  {
    expansion.clear();
    granularity = p.get_granularity();
    const typename Piecewise<C>::PiecesType* pieces = p.get_expansion();
    typename Piecewise<C>::PiecesType::const_iterator iter;
    for (iter=pieces->begin(); iter != pieces->end(); ++iter)
      expansion.insert(typename Piecewise<C>::PiecesType::value_type(iter->first, iter->second));
  }


  /*! constructor for predefined granularity */
  template <class C>
  Piecewise<C>::Piecewise(const int j)
//    :  Function<1,C>()
  {
    granularity = j;
  }


  /*! destructor */
  template <class C>
  Piecewise<C>::~Piecewise()
  {
    expansion.clear();
  }


  /*! reading access, local expansion on one subinterval */
  template <class C>
  Polynomial<C> Piecewise<C>::get_local_expansion(const int k) const
  {
    Polynomial<C> p;
    if (expansion.find(k)!=expansion.end())
      p = expansion.find(k)->second;

    return p;
  }


  /*! reading access, all polynomials */
  template <class C>
  const typename Piecewise<C>::PiecesType *Piecewise<C>::get_expansion() const
  {
    return &expansion;
  }


  /*! writing access, local expansion */
  template <class C>
  void Piecewise<C>::set_local_expansion(const int k, const Polynomial<C> &p)
  {
    expansion[k] = p;
  }


  /*! clip */
  template <class C>
  void Piecewise<C>::clip_me(const int k1, const int k2)
  {
    typename Piecewise<C>::PiecesType help(expansion);
    typename Piecewise<C>::PiecesType::const_iterator iter;
    for (iter=help.begin(); iter!=help.end(); ++iter)
      if ((iter->first<k1)||(iter->first>=k2))
        expansion.erase(iter->first);
  }

  template <class C>
  Piecewise<C> Piecewise<C>::clip(const int k1, const int k2) const
{
  Piecewise<C> r(*this);
  r.clip_me(k1, k2);

  return r;
}



  /*! shift */
  template <class C>
  void Piecewise<C>::shift_me(const int k)
  {
    Polynomial<C> p;
    p.set_coefficient(0, ldexp((double) -k, -granularity));
    p.set_coefficient(1, 1);
    typename Piecewise<C>::PiecesType help(expansion);
    expansion.clear();
    typename Piecewise<C>::PiecesType::const_iterator iter;
    for (iter=help.begin(); iter!=help.end(); ++iter)
      expansion.insert(typename Piecewise<C>::PiecesType::value_type(iter->first+k, iter->second.substitute_into(p)));

//   if (k>0)
//     {
//       map<int, Polynomial<C> >::reverse_iterator iter;
//       for (iter=expansion.rbegin(); iter!=expansion.rend(); ++iter)
// 	{
// 	  expansion[iter->first + k] = iter->second;
// 	  expansion.erase(iter->first);
//  	}
//     }
//   else
//     if (k<0)
//       {
// 	map<int, Polynomial<C> >::iterator iter;
// 	for (iter=expansion.begin(); iter!=expansion.end(); ++iter)
// 	  {
// 	    expansion[iter->first + k] = iter->second;
// 	    expansion.erase(iter->first);
// 	  }
//       }
  
//   map<int, Polynomial<C> > help;
//   map<int, Polynomial<C> >::const_iterator iter;
//   for (iter=expansion.begin(); iter!=expansion.end(); ++iter)
//     help[iter->first + k] = iter->second;
//   expansion = help;
  }

  template <class C>
  Piecewise<C> Piecewise<C>::shift(const int k) const
  {
    Piecewise r(*this);
    r.shift_me(k);

    return r;
  }


  /*! dilate */
  template <class C>
  void Piecewise<C>::dilate_me(const int j)
  {
    granularity += j;
    scale( sqrt((double)(1<<j))); // scale by 2^{j/2}
  }

  template <class C>
  Piecewise<C> Piecewise<C>::dilate(const int j) const
  {
    Piecewise r(*this);
    r.dilate_me(j);

    return r;
  }


  /*! increase granularity */
  template <class C>
  void Piecewise<C>::split_me(const int jnew)
  {
    if (jnew<granularity)
    {
      std::cout << std::endl << "Piecewise::split_me() error; granularity jnew=" << jnew
                << " coarser than current granularity j=" << granularity << "!" << std::endl;
      abort();
    }
    else
    {
      if (jnew>granularity)
      {
        const int d = 1<<(jnew-granularity);
          typename Piecewise<C>::PiecesType help;
          typename Piecewise<C>::PiecesType::const_iterator iter;
          for (iter=expansion.begin(); iter!=expansion.end(); ++iter)
            for (int k=d*iter->first; k<d*(iter->first+1); k++)
              help.insert(typename Piecewise<C>::PiecesType::value_type(k, iter->second));

          expansion.clear();
          for (iter=help.begin(); iter!=help.end(); ++iter)
            expansion.insert(typename Piecewise<C>::PiecesType::value_type(iter->first, iter->second));
          granularity = jnew;
      }
    }
  }

  template <class C>
  Piecewise<C> Piecewise<C>::split(const int jnew) const
  {
    Piecewise r(*this);
    r.split_me(jnew);

    return r;
  }


  /*! symbolic differentiation */
  template <class C>
  Piecewise<C> Piecewise<C>::differentiate() const
  {
    Piecewise<C> r(granularity);
    typename Piecewise<C>::PiecesType::const_iterator iter;
    for (iter=expansion.begin(); iter!=expansion.end(); ++iter)
      if(iter->second.degree() > 0)
        r.set_local_expansion(iter->first, iter->second.differentiate());

    return r;
  }


  /*! integration, entire support */
  template <class C>
  double Piecewise<C>::integrate(const bool quadrature) const
  {
    double help = ldexp(1.0, -granularity);
    double r = 0;
    typename Piecewise<C>::PiecesType::const_iterator iter;
    for (iter=expansion.begin(); iter!=expansion.end(); ++iter)
      r += iter->second.integrate(help*iter->first, help*(iter->first+1), quadrature);

    return r;
  }


  /*! integration, specific support */
  template <class C>
  double Piecewise<C>::integrate(const int k1, const int k2,
                                 const bool quadrature) const
  {
    double help = ldexp(1.0, -granularity);
    double r = 0;
    typename Piecewise<C>::PiecesType::const_iterator iter;
    for (iter=expansion.begin(); iter!=expansion.end(); ++iter)
      if ((iter->first>=k1)&&(iter->first<k2))
        r += iter->second.integrate(help*iter->first, help*(iter->first+1), quadrature);
      
    return r;
  }


  /*
    Piecewise operators ...
  */

  /*! assignment of another piecewise */
  template <class C>
  Piecewise<C>& Piecewise<C>::operator = (const Piecewise<C>& p)
  {
    expansion.clear();
    granularity = p.get_granularity();
    const typename Piecewise<C>::PiecesType* pieces = p.get_expansion();
    typename Piecewise<C>::PiecesType::const_iterator iter;
    for (iter=pieces->begin(); iter != pieces->end(); ++iter)
      expansion.insert(typename Piecewise<C>::PiecesType::value_type(iter->first, iter->second));

    return *this;
  }

  /*! point evaluation */
  template <class C>
  C Piecewise<C>::value (const C x) const
  {
    int k = (int) floor(ldexp(1.0, granularity)*x);
    typename Piecewise<C>::PiecesType::const_iterator iter = expansion.find(k);
    if (iter != expansion.end())
      return iter->second.value(x);
    else
      return C(0);
  }

  /*! point evaluation (calls above value(const C)) */
  template <class C>
  inline
  C Piecewise<C>::value(const Point<1>& p,
                        const unsigned int component) const
  {
    return value(p[0]);
  }

  template <class C>
  inline
  void Piecewise<C>::vector_value(const Point<1>& p,
                                  Vector<C>& values) const
  {
    values.resize(1, false);
    values[0] = value(p[0]);
  }

  /*! point evaluation operator */
  template <class C>
  C Piecewise<C>::operator () (const C x) const
  {
    return value(x);
  }

  /*! point evaluation of first derivative */
  template <class C>
  C Piecewise<C>::derivative(const C x) const
  {
    int k = (int) floor(ldexp(1.0, granularity)*x);
    typename Piecewise<C>::PiecesType::const_iterator iter = expansion.find(k);
    if (iter != expansion.end())
      return (iter->second).value(x,1); // evaluate first derivative
    else
      return C(0);
  }


  /*
    sums
  */
  /*! add a polynomial to this piecewise */
  template <class C>
  Piecewise<C>& Piecewise<C>::add (const Polynomial<C> &p)
  {
    // note: this operation is independent from current granularity
    typename Piecewise<C>::PiecesType::iterator iter;
    for (iter=expansion.begin(); iter!=expansion.end(); ++iter)
      iter->second += p;
    return *this;
  }

  /*! add a polynomial to this piecewise */
  template <class C>
  Piecewise<C>& Piecewise<C>::operator += (const Polynomial<C> &p)
  {
    return add(p);
  }

  /*! add an other piecewise to this piecewise */
  template <class C>
  Piecewise<C>& Piecewise<C>::add (const Piecewise<C> &p)
  {
    typename Piecewise<C>::PiecesType::const_iterator iter;
    if (p.get_granularity()>=granularity)
    {
      split_me(p.get_granularity());
      for (iter=p.expansion.begin(); iter!=p.expansion.end(); ++iter)
        if (expansion.find(iter->first)!=expansion.end())
          expansion[iter->first] = expansion[iter->first] + iter->second;
        else
          expansion[iter->first] = iter->second;
    }
    else
    {
      Piecewise<C> q(p);
      q.split_me(granularity);
      for (iter=q.expansion.begin(); iter!=q.expansion.end(); ++iter)
      if (expansion.find(iter->first)!=expansion.end())
        expansion[iter->first] = expansion[iter->first] + iter->second;
      else
        expansion[iter->first] = iter->second;
    }

    return *this;
  }

  /*! add an other piecewise to this piecewise */
  template <class C>
  Piecewise<C>& Piecewise<C>::operator += (const Piecewise<C> &p)
  {
    return add(p);
  }


  /*
    differences
  */
  /*! subtract a polynomial from this piecewise */
  template <class C>
  Piecewise<C>& Piecewise<C>::operator -= (const Polynomial<C> &p)
  {
    // note: this operation is independent from current granularity
    typename Piecewise<C>::PiecesType::iterator iter;
    for (iter=expansion.begin(); iter!=expansion.end(); ++iter)
      iter->second -= p;

    return *this;
  }

  /*! subtract an other piecewise from this piecewise */
  template <class C>
  Piecewise<C>& Piecewise<C>::operator -= (const Piecewise<C> &p)
  {
    typename Piecewise<C>::PiecesType::const_iterator iter;
    if (p.get_granularity()>=granularity)
    {
      split_me(p.get_granularity());
      for (iter=p.expansion.begin(); iter!=p.expansion.end(); ++iter)
        if (expansion.find(iter->first)!=expansion.end())
          expansion[iter->first] = expansion[iter->first] - iter->second;
        else
          expansion[iter->first] = -iter->second;
    }
    else
    {
      Piecewise<C> q(p);
      q.split_me(granularity);
      for (iter=q.expansion.begin(); iter!=q.expansion.end(); ++iter)
        if (expansion.find(iter->first)!=expansion.end())
          expansion[iter->first] = expansion[iter->first] - iter->second;
        else
          expansion[iter->first] = -iter->second;
    }

    return *this;
  }


  /*
    products
  */
  /*! in-place multiplication with a constant */
  template <class C>
  Piecewise<C>& Piecewise<C>::scale (const C c)
  {
    // independent from granularity
    typename Piecewise<C>::PiecesType::iterator iter;
    for (iter=expansion.begin(); iter!=expansion.end(); ++iter)
      iter->second *= c;

    return *this;
  }

  /*! in-place multiplication with a constant */
  template <class C>
  Piecewise<C>& Piecewise<C>::operator *= (const C c)
  {
    return scale(c);
  }

  template <class C>
  Piecewise<C>& Piecewise<C>::operator *= (const Polynomial<C>& p)
  {
    // independent from granularity
    typename Piecewise<C>::PiecesType::iterator iter;
    for (iter=expansion.begin(); iter!=expansion.end(); ++iter)
      iter->second *= p;

    return *this;
  }

  template <class C>
  Piecewise<C>& Piecewise<C>::operator *= (const Piecewise &p)
  {
    Piecewise<C> q(p);
    typename Piecewise<C>::PiecesType::const_iterator iter;

    if (q.get_granularity()>=granularity)
      split_me(q.get_granularity());
    else
      q.split_me(granularity);

    Piecewise<C> r(*this);

    for (iter = r.expansion.begin(); iter != r.expansion.end(); ++iter)
      if (q.expansion.find(iter->first) != q.expansion.end())
        expansion[iter->first] = q.expansion[iter->first] * iter->second;
      else
        expansion.erase(iter->first);

    return *this;
  }


  /*! inner product with another piecewise */
  template <class C>
  C Piecewise<C>::inner_product(const Piecewise<C> &p) const
  {
    Piecewise<C> r(*this); // make a copy
    r *= p;
    return r.integrate(true);
  }

  /*! inner product with a polynomial */
  template <class C>
  C Piecewise<C>::inner_product(const Polynomial<C> &p) const
  {
    Piecewise<C> r(*this); // make a copy
    r *= p;
    return r.integrate(true);
  }


  /*! Error output */
  template <class C>
  void Piecewise<C>::MatError(char* str) const
  {
    std::cerr << std::endl << "Piecewise error : " << str << std::endl;
    abort();
  }


  /*
    non-class-member operators ...
  */

  /*
    sums
  */
  template <class C>
  Piecewise<C> operator + (const Polynomial<C> &p, const Piecewise<C> &q)
  {
    Piecewise<C> r(q);
    r += p;

    return r;
  }

  template <class C>
  Piecewise<C> operator + (const Piecewise<C> &p, const Polynomial<C> &q)
  {
    Piecewise<C> r(p);
    r += q;

    return r;
  }

  template <class C>
  Piecewise<C> operator + (const Piecewise<C> &p, const Piecewise<C> &q)
  {
    Piecewise<C> r(p);
    r += q;

    return r;
  }


  /*
    differences
  */
  template <class C>
  Piecewise<C> operator - (const Polynomial<C> &p, const Piecewise<C> &q)
  {
    Piecewise<C> r(q);
    r *= -1.0;
    r += p;

    return r;
  }

  template <class C>
  Piecewise<C> operator - (const Piecewise<C> &p, const Polynomial<C> &q)
  {
    Piecewise<C> r(p);
    r -= q;

    return r;
  }

  template <class C>
  Piecewise<C> operator - (const Piecewise<C> &p, const Piecewise<C> &q)
  {
    Piecewise<C> r(p);
    r -= q;

    return r;
  }


  /*
    products
  */
  template <class C>
  Piecewise<C> operator * (const C c, const Piecewise<C> &q)
  {
    Piecewise<C> r(q);
    r *= c;

    return r;
  }

  template <class C>
  Piecewise<C> operator * (const Polynomial<C> &p, const Piecewise<C> &q)
  {
    Piecewise<C> r(q);
    r *= p;

    return r;
  }

  template <class C>
  Piecewise<C> operator * (const Piecewise<C> &p, const Piecewise<C> &q)
  {
    Piecewise<C> r(q);
    r *= p;

    return r;
  }


  /*! output into a stream */
  template <class C>
  std::ostream& operator << (std::ostream& s, const Piecewise<C>& p)
  {
    int oldprecision = s.precision();
    std::ios::fmtflags oldflags = s.flags();
    s.precision(12);
    const typename Piecewise<C>::PiecesType *pieces = p.get_expansion();
    typename Piecewise<C>::PiecesType::const_iterator iter;
    for (iter=pieces->begin(); iter!=pieces->end(); ++iter)
    {
      s << "[" << ldexp(1.0, -p.get_granularity())* iter->first << ","
        << ldexp(1.0, -p.get_granularity())* (iter->first+1) << "]: "
        << iter->second << std::endl;
    }

    s.setf(oldflags);
    s.precision(oldprecision);

    return s;
  }
}
