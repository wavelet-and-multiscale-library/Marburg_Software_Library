// implementation of multi_differences.h template functions

namespace MathTL
{
  template <unsigned int DIMENSION, unsigned int DIRECTION>
  inline
  InfiniteVector<double, MultiIndex<int, DIMENSION> >
  multivariate_forward_difference
  (const InfiniteVector<double, MultiIndex<int, DIMENSION> >& a)
  {
    typedef InfiniteVector<double, MultiIndex<int, DIMENSION> > vtype;
    vtype r;

    if (!a.empty())
      {
	for (typename vtype::const_iterator it(a.begin()); it != a.end(); ++it)
	  {
	    MultiIndex<int, DIMENSION> index1(it.index());
	    index1[DIRECTION] -= 1;
	    r.set_coefficient(index1, a.get_coefficient(it.index()));
	  }

	r -= a;
      }

    return r;
  }

  template <unsigned int DIMENSION, unsigned int DIRECTION>
  inline
  InfiniteVector<double, MultiIndex<int, DIMENSION> >
  multivariate_backward_difference
  (const InfiniteVector<double, MultiIndex<int, DIMENSION> >& a)
  {
    typedef InfiniteVector<double, MultiIndex<int, DIMENSION> > vtype;
    vtype r;
    
    if (!a.empty())
      {
	for (typename vtype::const_iterator it(a.begin()); it != a.end(); ++it)
	  {
	    MultiIndex<int, DIMENSION> index1(it.index());
	    index1[DIRECTION] += 1;
	    r.set_coefficient(index1, -a.get_coefficient(it.index()));
	  }

	r += a;
      }

    return r;
  }
}
