// implementation of differences.h template functions

// #ifndef _MATHTL_DIFFERENCES_CPP
// #define _MATHTL_DIFFERENCES_CPP

// #include <numerics/differences.h>

namespace MathTL
{
  template <unsigned int K>
  inline
  InfiniteVector<double, int> forward_difference(const InfiniteVector<double, int>& a)
  {
    return forward_difference<1>(forward_difference<K-1>(a));
  }

  template <>
  inline
  InfiniteVector<double, int> forward_difference<0>(const InfiniteVector<double, int>& a)
  {
    return InfiniteVector<double, int>(a);
  }

  template <>
  InfiniteVector<double, int> forward_difference<1>(const InfiniteVector<double, int>& a)
  {
    InfiniteVector<double, int> r;

    if (!a.empty())
      {
 	r[a.begin().index()-1] = *(a.begin());
	
	for (int k = a.begin().index(); k <= a.rbegin().index(); k++)
	  {
	    double help = a[k+1]-a[k];
	    if (help != 0)
	      r[k] = help;
	  }
      }
	
    return r;
  }

  template <unsigned int K>
  inline
  InfiniteVector<double, int> backward_difference(const InfiniteVector<double, int>& a)
  {
    return backward_difference<1>(backward_difference<K-1>(a));
  }

  template <>
  inline
  InfiniteVector<double, int> backward_difference<0>(const InfiniteVector<double, int>& a)
  {
    return InfiniteVector<double, int>(a);
  }

  template <>
  InfiniteVector<double, int> backward_difference<1>(const InfiniteVector<double, int>& a)
  {
    InfiniteVector<double, int> r;

    if (!a.empty())
      {
 	r[a.rbegin().index()+1] = -*(a.rbegin());
	
	for (int k = a.begin().index(); k <= a.rbegin().index(); k++)
	  {
	    double help = a[k]-a[k-1];
	    if (help != 0)
	      r[k] = help;
	  }
      }
	
    return r;
  }
}

// #endif
