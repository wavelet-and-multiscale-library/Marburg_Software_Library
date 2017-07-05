
#include "ldomain_frame_index.h"

// implementation for ldomain_frame_support.h

namespace WaveletTL
{
  template <class IFRAME>
  inline
  void support(const LDomainFrame<IFRAME>& frame,
	       const typename LDomainFrame<IFRAME>::Index& lambda,
	       typename LDomainFrame<IFRAME>::Support& supp)
  {
    frame.support(lambda, supp);
  }

  template <class IFRAME>
  bool intersect_supports(const LDomainFrame<IFRAME>& frame,
			  const typename LDomainFrame<IFRAME>::Index& lambda,
			  const typename LDomainFrame<IFRAME>::Support& supp2,
			  typename LDomainFrame<IFRAME>::Support& supp)
  {
    typedef typename LDomainFrame<IFRAME>::Support Support;

    Support supp1;
    frame.support(lambda, supp1);
    
    // quickly return false if the supports do not intersect at all
    if (supp1.xmin[0] == -1 && supp2.xmin[1] == -1 && supp2.xmin[2] == -1) return false;
    if (supp2.xmin[0] == -1 && supp1.xmin[1] == -1 && supp2.xmin[2] == -1) return false;
    if (supp2.xmin[0] == -1 && supp2.xmin[1] == -1 && supp1.xmin[2] == -1) return false;
    if (supp1.xmin[0] == -1 && supp1.xmin[1] == -1 && supp2.xmin[2] == -1) return false;
    if (supp1.xmin[0] == -1 && supp2.xmin[1] == -1 && supp1.xmin[2] == -1) return false;
    if (supp2.xmin[0] == -1 && supp1.xmin[1] == -1 && supp1.xmin[2] == -1) return false;
    
    
    bool r = false;

#if 1
    supp.j[0] = supp1.j[0];
    supp.j[1] = supp1.j[1];
    int xdiff1 = 0;
    int xdiff2 = supp.j[0]-supp2.j[0];
    int ydiff1 = 0;
    int ydiff2 = supp.j[1]-supp2.j[1];
    if (supp2.j[0] > supp1.j[0]) {
      supp.j[0] = supp2.j[0];
      xdiff1 = supp.j[0]-supp1.j[0];
      xdiff2 = 0;
    }
    if (supp2.j[1] > supp1.j[1]) {
      supp.j[1] = supp2.j[1];
      ydiff1 = supp.j[1]-supp1.j[1];
      ydiff2 = 0;
    }
#else
    // old version, slightly slower
    supp.j[0] = std::max(supp1.j[0], supp2.j[0]);
    supp.j[1] = std::max(supp1.j[1], supp2.j[1]);
    const int xdiff1 = supp.j[0]-supp1.j[0];
    const int xdiff2 = supp.j[0]-supp2.j[0];
    const int ydiff1 = supp.j[1]-supp1.j[1];
    const int ydiff2 = supp.j[1]-supp2.j[1];
#endif
    
    for (int patch = 0; patch <= 2; patch++) {
      if (supp1.xmin[patch] != -1 && supp2.xmin[patch] != -1) {
	// intersection of two nontrivial sets on patch p
	
	const int xmin1 = supp1.xmin[patch] << xdiff1;
	const int xmax1 = supp1.xmax[patch] << xdiff1;
	const int xmin2 = supp2.xmin[patch] << xdiff2;
	const int xmax2 = supp2.xmax[patch] << xdiff2;
	
	if (xmin1 < xmax2 && xmax1 > xmin2) {
	  // nontrivial intersection in x direction

	  const int ymin1 = supp1.ymin[patch] << ydiff1;
	  const int ymax1 = supp1.ymax[patch] << ydiff1;
	  const int ymin2 = supp2.ymin[patch] << ydiff2;
	  const int ymax2 = supp2.ymax[patch] << ydiff2;
	  
	  if (ymin1 < ymax2 && ymax1 > ymin2) {
	    // nontrivial intersection in y direction

	    supp.xmin[patch] = std::max(xmin1, xmin2);
	    supp.xmax[patch] = std::min(xmax1, xmax2);
	    supp.ymin[patch] = std::max(ymin1, ymin2);
	    supp.ymax[patch] = std::min(ymax1, ymax2);
	    
	    r = true;
	  } else {
	    supp.xmin[patch] = -1;
	  }
	} else {
	  supp.xmin[patch] = -1;
	}
      } else {
	supp.xmin[patch] = -1;
      }
    }
    
    return r;
  }

  template <class IFRAME>
  inline
  bool intersect_supports(const LDomainFrame<IFRAME>& frame,
			  const typename LDomainFrame<IFRAME>::Index& lambda1,
			  const typename LDomainFrame<IFRAME>::Index& lambda2,
			  typename LDomainFrame<IFRAME>::Support& supp)
  {
    typedef typename LDomainFrame<IFRAME>::Support Support;

    Support supp2;
    support(frame, lambda2, supp2);

    return intersect_supports(frame, lambda1, supp2, supp);
  }
  
  template <class IFRAME>
  void intersecting_quarklets(const LDomainFrame<IFRAME>& frame,
			     const typename LDomainFrame<IFRAME>::Index& lambda,
			     const typename LDomainFrameIndex<IFRAME>::level_type& j,
			     std::list<typename LDomainFrame<IFRAME>::Index>& intersecting,
                             const typename LDomainFrameIndex<IFRAME>::polynomial_type& p)
  {
    typedef typename LDomainFrame<IFRAME>::Index Index;
    typedef typename LDomainFrame<IFRAME>::Support Support;

    intersecting.clear();

    Support supp, supp_lambda;
    support(frame, lambda, supp_lambda);
    
    // a brute force solution
    if (j==frame.j0()) {
      Index last_qua(last_quarklet<IFRAME>(&frame, j, p));
      for (Index mu = first_generator<IFRAME>(&frame, j, p);; ++mu) {
	if (intersect_supports(frame, mu, supp_lambda, supp))
	  intersecting.push_back(mu);
	if (mu == last_qua) break;
      }
    } else {
      Index last_qua(last_quarklet<IFRAME>(&frame, j, p));
//      cout << "First quarklet: " << first_quarklet<IFRAME>(&frame, j, p) << endl;
//      cout << "Last quarklet: " << last_qua << endl;
      for (Index mu = first_quarklet<IFRAME>(&frame, j, p);; ++mu) {
//          cout << "Mu: " << mu << ", " << mu.number() << endl;
	if (intersect_supports(frame, mu, supp_lambda, supp))
	  intersecting.push_back(mu);
	if (mu == last_qua) break;
      }
    }
  }
  
  template <class IFRAME>
  bool intersect_singular_support(const LDomainFrame<IFRAME>& frame,
				  const typename LDomainFrame<IFRAME>::Index& lambda,
				  const typename LDomainFrame<IFRAME>::Index& mu)
  {
    // we cheat a bit here: we return true if already the supports intersect (overestimate)
    typedef typename LDomainFrame<IFRAME>::Support Support;
    Support supp;
    return intersect_supports(frame, lambda, mu, supp);
  }

}
