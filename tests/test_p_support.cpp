#include <iostream>

#include <interval/p_basis.h>
#include <interval/p_support.h>

using namespace std;
using namespace WaveletTL;

int main()
{
  cout << "Testing support calculation of [P] generators and wavelets..." << endl;

  const int d = 3;
  const int dT = 3;

  typedef PBasis<d,dT> Basis;
  typedef Basis::Index Index;
  typedef Basis::Support Support;

//   Basis basis(0, 0); // no b.c.'s
//   Basis basis(1, 0); // complementary b.c. at x=0
//   Basis basis(0, 1); // complementary b.c. at x=1
  Basis basis(1, 1); // complementary b.c.'s

#if 1
  for (int level = basis.j0(); level <= basis.j0()+1; level++) {
    cout << "- computing the supports of all generators and wavelets on level j=" << level << ":" << endl;
    
    Index lambda(first_generator(&basis, level));
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
      
      if (lambda == last_wavelet(&basis, level)) break;
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

#if 0
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

#if 0
  cout << "- checking intersection of singular supports:" << endl;
  for (Index lambda = first_generator(&basis, basis.j0()+2);; ++lambda)
    {
      Support supp;
      support(basis, lambda, supp.k1, supp.k2);
      cout << "psi_lambda, lambda=" << lambda << " has the support 2^{-"
	   << lambda.j()+lambda.e()
	   << "}["
	   << supp.k1
	   << ","
	   << supp.k2
	   << "]"
	   << endl;
      
      Support supp_0;
      support(basis, first_generator(&basis, basis.j0()), supp_0.k1, supp_0.k2);
      cout << "* first generator on level j0 has the support 2^{-"
	   << basis.j0()
	   << "}["
	   << supp_0.k1
	   << ","
	   << supp_0.k2
	   << "]"
	   << endl;

      cout << "* support intersection with first generator on level j0         : ";
      bool inter = intersect_supports(basis, lambda, first_generator(&basis, basis.j0()), supp);
      if (inter)
	cout << "2^{-" << supp.j << "}[" << supp.k1 << "," << supp.k2 << "]" << endl;
      else
	cout << "none" << endl;

      cout << "* singular support intersection with first generator on level j0: ";
      inter = intersect_singular_support(basis, lambda, first_generator(&basis, basis.j0()), supp.j, supp.k1, supp.k2);
      if (inter)
	cout << "2^{-" << supp.j << "}[" << supp.k1 << "," << supp.k2 << "]" << endl;
      else
	cout << "none" << endl;

      if (lambda == last_wavelet(&basis, basis.j0()+2)) break;
    }
#endif

  return 0;
}
