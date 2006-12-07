// implementation for map_tools.h

namespace MathTL
{
  template <class K, class C>
  void
  add_maps(const std::map<K,C>& m1, const std::map<K,C>& m2,
	   std::map<K,C>& result,
	   const double factor1 = 1.0, const double factor2 = 1.0)
  {
    // cf. set_union() from stl_algo.h:

    result.clear();

    typename std::map<K,C>::const_iterator it1(m1.begin()), it1end(m1.end());
    typename std::map<K,C>::const_iterator it2(m2.begin()), it2end(m2.end());
    typename std::map<K,C>::iterator hint(result.begin()), hint2(result.begin());
    
    while (it1 != it1end && it2 != it2end)
      {
	if (it1->first < it2->first)
	  {
	    hint2 = result.insert(hint, *it1);
	    hint2->second *= factor1;
	    hint = hint2;
	    ++it1;
	  }
	else
	  {
	    if (it2->first < it1->first)
	      {
		hint2 = result.insert(hint, *it2);
		hint2->second *= factor2;
		hint = hint2;
		++it2;
	      }
	    else
	      {
		const C value(factor1*it1->second + factor2*it2->second);
 		if (value != C(0)) {
		  hint2 = result.insert(hint, std::pair<K,C>(it1->first, value));
 		  hint = hint2;
 		}
		++it1;
		++it2;
	      }
	  }
      }

    while (it1 != it1end)
      {
 	hint2 = result.insert(hint, *it1);
	hint2->second *= factor1;
 	hint = hint2;
 	++it1;
      }

    while (it2 != it2end)
      {
 	hint2 = result.insert(hint, *it2);
	hint2->second *= factor2;
 	hint = hint2;
 	++it2;
      }
  }
}
