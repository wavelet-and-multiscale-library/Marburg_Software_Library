#include <iostream>

#include <Ldomain/ldomain_jl_basis.h>
#include <Ldomain/ldomain_jl_support.h>

using namespace std;
using namespace WaveletTL;

int main()
{
  cout << "Testing support calculations for LDomainJLBasis..." << endl;
  
  typedef LDomainJLBasis Basis;
  typedef Basis::Index Index;
  typedef LDomainJLBasis::Support Support;

  Basis basis;

#if 0
  for (int level = basis.j0(); level <= basis.j0()+1; level++) {
    for (Index lambda(first_generator(level)); lambda != last_wavelet(level); ++lambda) {
      Support supp;
      support(basis, lambda, supp);
      cout << "  supp(psi_{" << lambda << "})="
	   << "2^{-" << supp.j << "}"
	   << "[" << supp.xmin << "," << supp.xmax << "]"
	   << "x[" << supp.ymin << "," << supp.ymax << "]" << endl;
    }
  }
#endif

#if 1
  Index lambda = first_generator(basis.j0());
  Support supp_lambda;
  support(basis, lambda, supp_lambda);
  cout << "  supp(psi_{" << lambda << "})="
       << "2^{-" << supp_lambda.j << "}"
       << "[" << supp_lambda.xmin << "," << supp_lambda.xmax << "]"
       << "x[" << supp_lambda.ymin << "," << supp_lambda.ymax << "]" << endl;
  for (Index mu = first_generator(basis.j0()); mu != last_wavelet(basis.j0()+1); ++mu) {
    Support supp_intersect;
//     bool inter = intersect_supports(basis, mu, supp_lambda, supp_intersect);
    bool inter = intersect_supports(basis, lambda, mu, supp_intersect);
    cout << "- support intersection of psi_lambda and "
	 << "psi_{" << mu << "}: ";
    if (inter) {
      cout << "2^{-" << supp_intersect.j << "}"
	   << "[" << supp_intersect.xmin << "," << supp_intersect.xmax << "]"
	   << "x[" << supp_intersect.ymin << "," << supp_intersect.ymax << "]" << endl;
    } else {
      cout << "none!" << endl;
    }
  }
#endif

#if 1
  cout << "* for lambda=" << lambda << ", compute all intersecting generators on the coarsest level:" << endl;
  typedef std::list<Index> SupportList;
  SupportList nus;
  intersecting_wavelets(basis, lambda, basis.j0(), true, nus);
  for (SupportList::const_iterator it(nus.begin()); it != nus.end(); ++it) {
    cout << "    nu=" << *it << endl;
  }
#endif

  return 0;
}
