#include <iostream>

#include <interval/spline_basis.h>
#include <interval/spline_support.h>

using namespace std;
using namespace WaveletTL;

int main()
{
  cout << "Testing support calculation of SplineBasis generators and wavelets..." << endl;

#if 0
  typedef SplineBasis<2,2,P_construction> Basis;
  Basis basis("",0,0,0,0); // PBasis, no b.c.'s
#else
  typedef SplineBasis<3,3,DS_construction> Basis;
  typedef Basis::Index Index;
  Basis basis("bio5",0,0,0,0); // DSBasis, no b.c.'s
//   Basis basis("bio5-energy",0,0,0,0); // DSBasis, no b.c.'s
#endif

  typedef Basis::Index Index;
  typedef Basis::Support Support;


#if 1
  for (int level = basis.j0(); level <= basis.j0()+1; level++) {
    cout << "- computing the supports of all generators and wavelets on level j=" << level << ":" << endl;
    
    Index lambda(basis.first_generator(level));
    for (;; ++lambda) {
      int k1, k2;
      support(basis, lambda, k1, k2);
      cout << "  lambda=" << lambda << ", supp(psi_lambda)=2^{-"
	   << lambda.j()+lambda.e()
	   << "}["
	   << k1
	   << ","
	   << k2
	   << "]"
	   << endl;
      
      if (lambda == basis.last_wavelet(level)) break;
    }
  }
#endif
 
#if 0
  cout << "- calculating some support intersections:" << endl;
  for (Index lambda = first_generator(&basis, basis.j0());; ++lambda)
    {
      Support supp;
      support(basis, lambda, supp.k1, supp.k2);
      cout << "psi_lambda, lambda=" << lambda << " has support 2^{-"
	   << lambda.j()+lambda.e()
	   << "}["
	   << supp.k1
	   << ","
	   << supp.k2
	   << "]"
	   << endl;

      cout << "support intersection with first generator on level j0: ";
      bool inter = intersect_supports(basis, lambda, first_generator(&basis, basis.j0()), supp);
      if (inter)
	cout << "2^{-" << supp.j << "}[" << supp.k1 << "," << supp.k2 << "]" << endl;
      else
	cout << "none" << endl;
      
      if (lambda == last_wavelet(&basis, basis.j0())) break;
    }
#endif

#if 1
  cout << "- compute all intersecting wavelets:" << endl;
  for (Index lambda = first_generator(&basis, basis.j0());; ++lambda)
    {
      cout << "  * for lambda=" << lambda << ":" << endl;
      typedef std::list<std::pair<Index, Basis::Support> > SupportList;
      SupportList nus;
      intersecting_wavelets(basis, lambda, basis.j0(), true, nus);
      for (SupportList::const_iterator it(nus.begin()); it != nus.end(); ++it) {
	cout << "    nu=" << it->first 
	     << " with support intersection "
	     << "2^{-" << it->second.j << "}[" << it->second.k1 << "," << it->second.k2 << "]" << endl;
      }
      for (int level = basis.j0(); level <= basis.j0()+2; level++) {
	intersecting_wavelets(basis, lambda, level, false, nus);
	for (SupportList::const_iterator it(nus.begin()); it != nus.end(); ++it) {
	  cout << "    nu=" << it->first 
	       << " with support intersection "
	       << "2^{-" << it->second.j << "}[" << it->second.k1 << "," << it->second.k2 << "]" << endl;
	}
      }
      
      if (lambda == last_wavelet(&basis, basis.j0())) break;
    }
#endif  

  return 0;
}
