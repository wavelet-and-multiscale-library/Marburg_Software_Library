// implementation for compression.h

#include <map>
#include <list>

namespace WaveletTL
{
  template <class PROBLEM>
  void
  add_compressed_column(const PROBLEM& P,
			const double factor,
			const typename PROBLEM::Index& lambda,
			const int J,
			InfiniteVector<double, typename PROBLEM::Index>& w,
			const int jmax,
			const CompressionStrategy strategy)
  {
    //typedef typename PROBLEM::WaveletBasis WaveletBasis;
    typedef typename PROBLEM::Index Index;
    //typedef typename WaveletBasis::Support Support;
    typedef std::list<Index> IntersectingList;
    
    if (P.local_operator())
      {
	// differential operators
	
	if (strategy == CDD1) {
	  // [CDD1] strategy:
	  // active row indices nu have to fulfill ||nu|-|lambda|| <= J/d and
	  // the supports of psi_lambda and psi_nu have to intersect
	  
	  const int maxlevel = std::min(lambda.j()+(J/P.space_dimension()), jmax);
	  for (int level = std::max(P.basis().j0()-1, lambda.j()-(J/P.space_dimension()));
	       level <= maxlevel; level++)
	    {
	      P.add_level(lambda,w,level,factor,J,strategy);
	    }
	}

	if (strategy == St04a) {
	  // [St04a] strategy:
	  // active row indices nu have to fulfill the following conditions:
	  //   ( ||nu|-|lambda|| <= k(j,d), where k(j,d) = j/(d-1) for d>1 and
	  //     j<=k(j,1)<=2^j and k(j,1)>j*min(t+mT,sigma)/(gamma-t) )
          //     and
	  //   ( ||nu|-|lambda|| <= j/d or supp(psi_lambda) intersects singsupp(psi_nu) (for |lambda|>|nu|) )
	  // Here gamma is the Sobolev regularity of the primal basis, t the order of the operator and
	  // sigma some exponent such that L,L':H^{t+sigma}\to H^{-t+sigma} are bounded.
	  // For the moment, we neglect sigma here.
	  
	  const double kjd = (P.space_dimension() == 1
			      ? std::min((double)J, ceil(J*(P.operator_order()+P.basis().primal_vanishing_moments()) / 
							 ((double) P.basis().primal_regularity()-P.operator_order())))
			      : J / (P.space_dimension()-1.));
	  //cout << P.operator_order() << endl;
	  //cout << "kjd = " << kjd << endl;
 	  const int maxlevel = std::min((int)floor(lambda.j()+kjd), jmax);
 	  IntersectingList nus;
 	  for (int level = std::max(P.basis().j0()-1, (int)ceil(lambda.j()-kjd));
 	       level <= maxlevel; level++)
 	    {
	      P.add_level(lambda,w,level,factor,J,strategy);
	    }
	}
      }
    else 
      {
	// integral operators: branch is not implemented so far
      }
  }
}
