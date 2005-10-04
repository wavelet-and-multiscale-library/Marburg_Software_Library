// implementation for cached_problem.h

namespace WaveletTL
{
  template <class PROBLEM>
  CachedProblem<PROBLEM>::CachedProblem(const PROBLEM& P)
    : PROBLEM(P)
  {
  }

  template <class PROBLEM>
  double
  CachedProblem<PROBLEM>::a(const typename WaveletBasis::Index& lambda,
			    const typename WaveletBasis::Index& nu) const
  {
    // first check whether the lambda-th column already exists in the cache
    typename ColumnCache::iterator col_lb(entries_cache.lower_bound(lambda));
    typename ColumnCache::iterator col_it(col_lb);
    if (col_lb == entries_cache.end() ||
	entries_cache.key_comp()(lambda, col_lb->first))
      {
	// insert a new cache column
	typedef typename ColumnCache::value_type value_type;
	col_it = entries_cache.insert(col_lb, value_type(lambda, Column()));
      }
    
    // for simplicity, extract the column
    Column& col(col_it->second);

    double r = 0;

    // check whether the entry has already been calculated
    typename Column::iterator lb(col.lower_bound(nu));
    typename Column::iterator it(lb);
    if (lb == col.end() ||
	col.key_comp()(nu, lb->first))
      {
	// compute the entry ...
	r = PROBLEM::a(lambda, nu);
	// ... and insert it in the cache
	typedef typename Column::value_type value_type;
	it = col.insert(lb, value_type(nu, r));
      }
    else
      r = it->second;
    
    return r;
  }
}
