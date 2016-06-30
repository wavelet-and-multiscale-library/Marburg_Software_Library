// implementation for ring_gramian.h

namespace WaveletTL
{
  template <int d, int dt, int s0, int s1>
  RingGramian<d,dt,s0,s1>::RingGramian(const RingBasis<d,dt,s0,s1>& basis,
				       const InfiniteVector<double,Index>& y)
    : basis_(basis), y_(y)
  {
  }
  
  template <int d, int dt, int s0, int s1>
  inline
  double
  RingGramian<d,dt,s0,s1>::a(const Index& lambda,
			     const Index& mu) const
  {
    Index0 lambda0(lambda.j(), lambda.e()[0], lambda.k()[0]),
      mu0(mu.j(), mu.e()[0], mu.k()[0]);
    Index1 lambda1(lambda.j(), lambda.e()[1], lambda.k()[1]),
      mu1(mu.j(), mu.e()[1], mu.k()[1]);
    
    double i0 = 0, i1 = 0;

    // lookup 1D integrals in the cache

    typename One_D_IntegralCache0::iterator col_lb0(one_d_integrals0.lower_bound(lambda0));
    typename One_D_IntegralCache0::iterator col_it0(col_lb0);
    if (col_lb0 == one_d_integrals0.end() ||
	one_d_integrals0.key_comp()(lambda0,col_lb0->first))
      {
	// insert a new column
	typedef typename One_D_IntegralCache0::value_type value_type;
	col_it0 = one_d_integrals0.insert(col_lb0, value_type(lambda0, Column1D_0()));
      }
    
    Column1D_0& col0(col_it0->second);
    
    typename Column1D_0::iterator lb0(col0.lower_bound(mu0));
    typename Column1D_0::iterator it0(lb0);
    if (lb0 == col0.end() ||
	col0.key_comp()(mu0, lb0->first))
      {
	// cache miss
	i0 = basis_.basis0().integrate(lambda0, mu0);
	typedef typename Column1D_0::value_type value_type;
	it0 = col0.insert(lb0, value_type(mu0, i0));
      }
    else
      {
	// cache hit
	i0 = it0->second;
      }

    if (fabs(i0) > 1e-15)
      {
	typename One_D_IntegralCache1::iterator col_lb1(one_d_integrals1.lower_bound(lambda1));
	typename One_D_IntegralCache1::iterator col_it1(col_lb1);
	if (col_lb1 == one_d_integrals1.end() ||
	    one_d_integrals1.key_comp()(lambda1,col_lb1->first))
	  {
	    // insert a new column
	    typedef typename One_D_IntegralCache1::value_type value_type;
	    col_it1 = one_d_integrals1.insert(col_lb1, value_type(lambda1, Column1D_1()));
	  }
    
	Column1D_1& col1(col_it1->second);
    
	typename Column1D_1::iterator lb1(col1.lower_bound(mu1));
	typename Column1D_1::iterator it1(lb1);
	if (lb1 == col1.end() ||
	    col1.key_comp()(mu1, lb1->first))
	  {
	    // cache miss
	    i1 = basis_.basis1().integrate(lambda1, mu1);
	    typedef typename Column1D_1::value_type value_type;
	    it1 = col1.insert(lb1, value_type(mu1, i1));
	  }
	else
	  {
	    // cache hit
	    i1 = it1->second;
	  }
	
	return i0 * i1; // == basis_.integrate(lambda, mu);
      }
    
    return 0; // for the compiler
  }

}
