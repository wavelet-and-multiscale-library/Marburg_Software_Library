// implementation for multi_laurent_polynomial.h

namespace MathTL
{
  template <class R, unsigned int DIMENSION>
  MultivariateLaurentPolynomial<R, DIMENSION>::MultivariateLaurentPolynomial()
    : InfiniteVector<R, MultiIndex<int, DIMENSION> >(), Function<DIMENSION, R>()
  {
  }

  template <class R, unsigned int DIMENSION>
  MultivariateLaurentPolynomial<R, DIMENSION>::MultivariateLaurentPolynomial
  (const MultivariateLaurentPolynomial<R, DIMENSION>& p)
    : InfiniteVector<R, MultiIndex<int, DIMENSION> >(p)
  {
  }
  
  template <class R, unsigned int DIMENSION>
  MultivariateLaurentPolynomial<R, DIMENSION>::MultivariateLaurentPolynomial
  (const R value)
    : InfiniteVector<R, MultiIndex<int, DIMENSION> >(), Function<DIMENSION, R>()
  {
    set_coefficient(MultiIndex<int, DIMENSION>(), value);
  }

  template <class R, unsigned int DIMENSION>
  inline
  R MultivariateLaurentPolynomial<R, DIMENSION>::get_coefficient
  (const MultiIndex<int, DIMENSION>& k) const
  {
    return InfiniteVector<R, MultiIndex<int, DIMENSION> >::operator [] (k);
  }

  template <class R, unsigned int DIMENSION>
  inline
  void MultivariateLaurentPolynomial<R, DIMENSION>::set_coefficient
  (const MultiIndex<int, DIMENSION>& k, const R coeff)
  {
    InfiniteVector<R, MultiIndex<int, DIMENSION> >::operator [] (k) = coeff;
  }

  template <class R, unsigned int DIMENSION>
  R MultivariateLaurentPolynomial<R, DIMENSION>::value
  (const Point<DIMENSION>& p,
   const unsigned int component) const
  {
    R r(0);

    // we don't care about complexity here:
    for (typename MultivariateLaurentPolynomial<R, DIMENSION>::const_iterator it(begin());
  	 it != end(); ++it)
      {
	R help(1);
	for (unsigned int i(0); i < DIMENSION; i++)
	  {
	    assert(p[i] != 0);
	    double xi = (it.index()[i] < 0 ? 1.0/p[i] : p[i]);
	    for (int m(1); m <= abs(it.index()[i]); m++)
	      help *= xi;
	  }

	r += *it * help;
      }

    return r;
  }

  template <class R, unsigned int DIMENSION>
  inline
  void MultivariateLaurentPolynomial<R, DIMENSION>::vector_value
  (const Point<DIMENSION> &p,
   Vector<R>& values) const
  {
    values.resize(1, false);
    values[0] = value(p);
  }

  template <class R, unsigned int DIMENSION>
  inline
  typename MultivariateLaurentPolynomial<R, DIMENSION>::const_iterator
  MultivariateLaurentPolynomial<R, DIMENSION>::begin() const
  {
    return InfiniteVector<R, MultiIndex<int, DIMENSION> >::begin();
  }

  template <class R, unsigned int DIMENSION>
  inline
  typename MultivariateLaurentPolynomial<R, DIMENSION>::const_iterator
  MultivariateLaurentPolynomial<R, DIMENSION>::end() const
  {
    return InfiniteVector<R, MultiIndex<int, DIMENSION> >::end();
  }  

  template <class R, unsigned int DIMENSION>
  inline
  typename MultivariateLaurentPolynomial<R, DIMENSION>::const_reverse_iterator
  MultivariateLaurentPolynomial<R, DIMENSION>::rbegin() const
  {
    return InfiniteVector<R, MultiIndex<int, DIMENSION> >::rbegin();
  }
  
  template <class R, unsigned int DIMENSION>
  inline
  typename MultivariateLaurentPolynomial<R, DIMENSION>::const_reverse_iterator
  MultivariateLaurentPolynomial<R, DIMENSION>::rend() const
  {
    return InfiniteVector<R, MultiIndex<int, DIMENSION> >::rend();
  }

  template <class R, unsigned int DIMENSION>
  std::ostream& operator <<
    (std::ostream& s, const MultivariateLaurentPolynomial<R, DIMENSION>& p)
  {
    for (typename MultivariateLaurentPolynomial<R, DIMENSION>::const_iterator it(p.begin());
  	 it != p.end(); ++it)
      {
  	if (it != p.begin())
  	  s << "+";
  	s << "(" << *it << ")";
	for (unsigned int i(0); i < DIMENSION; i++)
	  {
	    if (it.index()[i] != 0)
	      {
		s << "z_" << i+1;
		if (it.index()[i] != 1)
		  {
		    s << "^";
		    if (it.index()[i] < 0)
		      s << "{" << it.index()[i] << "}";
		    else
		      s << it.index()[i];
		  }
	      }
  	  }
      }
    
    return s;
  }
}
