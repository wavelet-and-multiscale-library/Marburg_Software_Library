// implementation for compression.h

#include <map>
#include <list>

namespace WaveletTL
{
  template <class PROBLEM>
  void
  add_compressed_column(const PROBLEM& P,
			const double factor,
			const typename PROBLEM::WaveletBasis::Index& lambda,
			const int J,
			InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& w,
			const int jmax,
			const CompressionStrategy strategy)
  {
    typedef typename PROBLEM::WaveletBasis WaveletBasis;
    typedef typename WaveletBasis::Index Index;
    typedef typename WaveletBasis::Support Support;
    typedef std::list<Index> IntersectingList;
    
    if (P.local_operator())
      {
	// differential operators
	
	if (strategy == CDD1) {
	  // [CDD1] strategy:
	  // active row indices nu have to fulfill ||nu|-|lambda|| <= J/d and
	  // the supports of psi_lambda and psi_nu have to intersect
	  
	  const int maxlevel = std::min(lambda.j()+J/P.space_dimension(), jmax);
	  IntersectingList nus;
	  for (int level = std::max(P.basis().j0()-1, lambda.j()-J/P.space_dimension());
	       level <= maxlevel; level++)
	    {
	      // compute all wavelets on level j, such that supp(psi_lambda) and supp(psi_nu) intersect
	      intersecting_wavelets(P.basis(), lambda,
				    std::max(level, P.basis().j0()),
				    level == (P.basis().j0()-1),
				    nus);
	      
	      // traverse the matrix block and update the result
	      const double d1 = P.D(lambda);
	      for (typename IntersectingList::const_iterator it(nus.begin()), itend(nus.end());
		   it != itend; ++it) {
		const double entry = P.a(*it, lambda) / (d1*P.D(*it));
		w.set_coefficient(*it,
				  w.get_coefficient(*it)
				  + entry * factor);
	      }
	    }
	}

	if (strategy == St04a) {
	  // [St04a] strategy:
	  // active row indices nu have to fulfill ||nu|-|lambda|| <= J/d and
	  // the support of psi_lambda has to intersect the singular support of psi_nu
	  
 	  const int maxlevel = std::min(lambda.j()+J/P.space_dimension(), jmax);
 	  IntersectingList nus;
 	  for (int level = std::max(P.basis().j0()-1, lambda.j()-J/P.space_dimension());
 	       level <= maxlevel; level++)
 	    {
 	      // compute all wavelets on level j, such that supp(psi_lambda) and singsupp(psi_nu) intersect
	      // (and vice versa)
 	      relevant_wavelets(P.basis(), lambda,
				std::max(level, P.basis().j0()),
				level == (P.basis().j0()-1),
				nus);
	      
 	      // traverse the matrix block and update the result
 	      const double d1 = P.D(lambda);
 	      for (typename IntersectingList::const_iterator it(nus.begin()), itend(nus.end());
 		   it != itend; ++it) {
 		const double entry = P.a(*it, lambda) / (d1*P.D(*it));
 		w.set_coefficient(*it,
 				  w.get_coefficient(*it)
 				  + entry * factor);
 	      }
 	    }
	}
      }
    else 
      {
	// integral operators: branch is not implemented so far
      }
  }
}
