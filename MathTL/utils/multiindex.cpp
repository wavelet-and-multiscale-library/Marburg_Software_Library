// implementation for multiindex.h
#include <cassert>
#include <algorithm>
#include <utils/tiny_tools.h>

namespace MathTL
{
  template <class I, unsigned int DIMENSION>
  MultiIndex<I, DIMENSION>::MultiIndex()
    : FixedArray1D<I, DIMENSION>()
  {
    for (unsigned int i(0); i < DIMENSION; i++)
      FixedArray1D<I, DIMENSION>::operator [] (i) = I(0);
  }

  template <class I, unsigned int DIMENSION>
  MultiIndex<I, DIMENSION>::MultiIndex(const MultiIndex<I, DIMENSION>& lambda)
    : FixedArray1D<I, DIMENSION>(lambda)
  {
  }

  template <class I, unsigned int DIMENSION>
  MultiIndex<I, DIMENSION>::MultiIndex(const I& i0)
    : FixedArray1D<I, DIMENSION>()
  {
    assert(DIMENSION == 1);
    FixedArray1D<I, DIMENSION>::operator [] (0) = i0;
  }

  template <class I, unsigned int DIMENSION>
  MultiIndex<I, DIMENSION>::MultiIndex(const I& i0, const I& i1)
    : FixedArray1D<I, DIMENSION>()
  {
    assert(DIMENSION == 2);
    FixedArray1D<I, DIMENSION>::operator [] (0) = i0;
    FixedArray1D<I, DIMENSION>::operator [] (1) = i1;
  }

  template <class I, unsigned int DIMENSION>
  bool MultiIndex<I, DIMENSION>::operator == (const MultiIndex& lambda) const
  {
    for (unsigned int i(0); i < DIMENSION; i++)
      if (FixedArray1D<I, DIMENSION>::operator [] (i) != lambda[i]) return false;
    return true;
  }

  template <class I, unsigned int DIMENSION>
  inline
  bool MultiIndex<I, DIMENSION>::operator < (const MultiIndex& lambda) const
  {
         
      unsigned int r1(FixedArray1D<I,DIMENSION>::operator [] (0)),r2(multi_degree(lambda));
      for (unsigned int i(1); i < DIMENSION; i++)
      {
          r1 += FixedArray1D<I,DIMENSION>::operator [] (i);
      }
      return (r1<r2) || ((r1==r2) && (std::lexicographical_compare(FixedArray1D<I, DIMENSION>::begin(),
                                                                   FixedArray1D<I, DIMENSION>::end(),
                                                                   lambda.begin(), lambda.end()) ));
  }
  
  template <class I, unsigned int DIMENSION>
  inline
  bool MultiIndex<I, DIMENSION>::lex (const MultiIndex<I, DIMENSION>& lambda) const
  {
      return std::lexicographical_compare(FixedArray1D<I, DIMENSION>::begin(),
                                          FixedArray1D<I, DIMENSION>::end(),
                                          lambda.begin(), lambda.end());
  }
  
  template<class I, unsigned int DIMENSION>
  MultiIndex<I, DIMENSION>&
  MultiIndex<I, DIMENSION>::operator ++ ()
  {
  	bool is_zero = true;
  	for (int i(DIMENSION-1); i >= 0; i--)
  	{
  		I value (FixedArray1D<I, DIMENSION>::operator [] (i));
  		if (value != 0)
  		{
  			is_zero = false;
  			if (i != 0)
  			{
  				FixedArray1D<I,DIMENSION>::operator [] (i-1) = FixedArray1D<I,DIMENSION>::operator [] (i-1) +1;
  				FixedArray1D<I,DIMENSION>::operator [] (i) = 0;
  				FixedArray1D<I,DIMENSION>::operator [] (DIMENSION-1) = value -1;
  			} else
  			{
  				FixedArray1D<I,DIMENSION>::operator [] (0) = 0;
  				FixedArray1D<I,DIMENSION>::operator [] (DIMENSION-1) = value +1;
  			}
  			break;
  		}
  	}
  	if (is_zero==true) FixedArray1D<I,DIMENSION>::operator [] (DIMENSION-1) = 1;
  	return *this;
  }

    template<class I, unsigned int DIMENSION>
    unsigned long int MultiIndex<I, DIMENSION>::number() const
    {
        unsigned long int num;
        unsigned int level,a,b;
        switch (DIMENSION)
        {
            case 1:
                num = (FixedArray1D<I, DIMENSION>::operator [] (0));
                break;
            case 2:
                level = multi_degree(*this);
                num = (level*(level+1)/2 + (FixedArray1D<I, DIMENSION>::operator [] (0)) );
                break;
            case 3:
                //lambda = (a,b,c)
                level = multi_degree(*this);
                a = FixedArray1D<I, DIMENSION>::operator [] (0);
                b = FixedArray1D<I, DIMENSION>::operator [] (1);
                num = level*(level+1)*(level+2)/6 + a*(2*level-a+3)/2 + b;
                break;
            default:
                // much slower, but works in every dimension
                MultiIndex temp;
                num = 0;
                while ((temp) < (*this))
                {
                    ++temp;
                    num++;
                }
                break;
        }
        return num;
    };
    
    
  template<class I, unsigned int DIMENSION>
  std::map<MultiIndex<I,DIMENSION>,int>
  indexmapping(const MultiIndex<I, DIMENSION>& lambda)
  {

      typedef std::map<MultiIndex<I, DIMENSION>,int> map_type;
      typedef typename map_type::iterator map_it;

      MultiIndex<I, DIMENSION> temp;
      int num(0);
      map_type r;
      //pair <map_it,bool> res;

      map_it it;
      it=r.begin();
      //typedef typename set_type::value_type value_type;
      //cout << "computing indexnumbers fr lambda = "<<lambda<<endl;
      it = r.insert(it, std::pair<MultiIndex<I,DIMENSION>,int>(temp,num));
      while (temp < lambda)
      {
          //cout << "inserting temp,num=" << temp << " , " << num << endl;
          it=r.insert(it, std::pair<MultiIndex<I,DIMENSION>,int>(temp,num));
          ++temp;
          num++;
      }
      return r;
  }

  template<class I, unsigned int DIMENSION>
  std::set<MultiIndex<I, DIMENSION> >
  cuboid_indices(const MultiIndex<I, DIMENSION>& alpha,
		 const MultiIndex<I, DIMENSION>& beta)
  {
    typedef std::set<MultiIndex<I, DIMENSION> > set_type;
    set_type r;

    MultiIndex<I, DIMENSION> gamma(alpha);
    for (I index(alpha[0]); index <= beta[0]; index++)
      {
	gamma[0] = index;
	r.insert(gamma);
      }
    for (unsigned int i(1); i < DIMENSION; i++)
      {
	set_type sofar(r);
	for (typename set_type::const_iterator it(sofar.begin()); it != sofar.end(); ++it)
	  {
	    gamma = *it;
	    for (I index(alpha[i]); index <= beta[i]; index++)
	      {
		gamma[i] = index;
		r.insert(gamma);
	      }
	  }
      }

    return r;
  }

  template <class I, unsigned int DIMENSION>
  std::set<MultiIndex<I, DIMENSION> >
  degree_indices(const unsigned int k)
  {
    typedef std::set<MultiIndex<I, DIMENSION> > set_type;
    set_type r;

    MultiIndex<I, DIMENSION> alpha;
    r.insert(alpha);
    for (unsigned int i(1); i <= k; i++)
      {
	set_type sofar(r);
	r.clear();
	for (typename set_type::const_iterator it(sofar.begin()); it != sofar.end(); ++it)
	  {
	    alpha = *it;
	    for (unsigned int j(0); j < DIMENSION; j++)
	      {
		alpha[j] += 1;
		r.insert(alpha);
		alpha[j] -= 1;
	      }
	  }	
      }

    return r;
  }

  template <class I, unsigned int DIMENSION>
  unsigned int multi_degree(const MultiIndex<I, DIMENSION>& alpha)
  {
    unsigned int r(alpha[0]);
    for (unsigned int i(1); i < DIMENSION; i++)
      r += alpha[i];
    return r;
  }

  template<class I, unsigned int DIMENSION>
  unsigned int multi_factorial(const MultiIndex<I, DIMENSION>& alpha)
  {
    unsigned int r(factorial(alpha[0]));
    for (unsigned int i(1); i < DIMENSION; i++)
      r *= factorial(alpha[i]);
    return r;
  }

  template <class I, unsigned int DIMENSION>
  int multi_power(const MultiIndex<I, DIMENSION>& alpha,
		  const MultiIndex<I, DIMENSION>& beta)
  {
    int r(intpower(alpha[0], beta[0]));
    for (unsigned int i(1); i < DIMENSION; i++)
      r *= intpower(alpha[i], beta[i]);
    return r;
  }

  template <class I, unsigned int DIMENSION>
  int multi_binomial(const MultiIndex<I, DIMENSION>& beta,
		     const MultiIndex<I, DIMENSION>& alpha)
  {
    int r(binomial(beta[0], alpha[0]));
    for (unsigned int i(1); i < DIMENSION; i++)
      r *= binomial(beta[i], alpha[i]);
    return r;
  }

}
