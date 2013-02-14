#include <iostream>

#include <interval/ds_basis.h>
#include <interval/ds_support.h>

using namespace std;
using namespace WaveletTL;

int main()
{
  cout << "Testing support calculation of DS generators and wavelets..." << endl;

  const int d = 2;
  const int dT = 4;

  typedef DSBasis<d,dT> Basis;
  typedef Basis::Index Index;
  typedef Basis::Support Support;
  
  Basis basis;
  Index lambda(first_generator(&basis, basis.j0()));
//   ++lambda;
//   for (int i = 1; i <= 4; i++, ++lambda);
//   lambda = last_generator(&basis, basis.j0());
//   lambda = first_wavelet(&basis, basis.j0()+1);
//   for (int i = 1; i <= 4; i++, ++lambda);
//   lambda = last_wavelet(&basis, basis.j0()+1);

  cout << "get_intersecting_wavelets_on_level does NOT work!!" << endl;
  cout << "The following code was moved from test_tbasis_support and may not compile correctly"
  abort();
#if 0
    /*
     * get_intersecting_wavelets_on_level does not work for DSBasis!
     * The following code is copied to the dsbasis_testfile
     */
    int tempA, tempB;
    //Index temp_mu(4,1,6,TBasisArray[0]->bases()[1]);
    Index first(basis.get_wavelet(0));
    //typedef Basis1D::Support Support1D;
    //Support1D supp1d;
    basis.support(first, tempA, tempB);
    cout << "wavelet = " << first << "; 1d support = (" << tempA << ", " << tempB << ")" << endl;
            
            
    Index temp_mu(basis.get_wavelet(20));
    cout << "temp_mu = " << temp_mu << endl;
    
    get_intersecting_wavelets_on_level(basis,
                    temp_mu,
            4,true,tempA,tempB);
    cout << "temp_mu -> get_intersecting_wavelets (4, true) (min, max) = (" << tempA << ", " << tempB << ")" << endl;
    get_intersecting_wavelets_on_level(basis,
                    temp_mu,
            4,false,tempA,tempB);
    cout << "temp_mu -> get_intersecting_wavelets (4, true) (min, max) = (" << tempA << ", " << tempB << ")" << endl;
     
    cout << "evaluate wavelet mu = " << temp_mu << endl;
    for (unsigned int i=0; i< 100; ++i)
    {
        cout << "f(" << (0.01*i) << ") = " << basis.evaluate(0, temp_mu, 0.01*i) << endl;
    }

    cout << "min = " << tempA << "; max = " << tempB << endl;
    for (unsigned int i = 0; i< 40; ++i)
    {
        Index temp_index(basis.get_wavelet(i));
        basis.support(temp_index, tempA, tempB);
        cout << "N = " << i << "; lam = " << temp_index << "; k1 = " << tempA << "; k2 = " << tempB << endl;
    }
    abort();
#endif
    
#if 0
  cout << "- point values at dyadic points of the " << lambda << " generator/wavelet:" << endl;
  evaluate(basis, lambda, true, 6).matlab_output(cout);
#endif

#if 0
  int k1, k2;
  support(basis, lambda, k1, k2);
  cout << "- support of psi_lambda is 2^{-"
       << lambda.j()+lambda.e()
       << "}["
       << k1
       << ","
       << k2
       << "]"
       << endl;
#endif
  
#if 1
  cout << "- calculating some support intersections:" << endl;
  for (lambda = first_generator(&basis, basis.j0());; ++lambda)
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
  for (lambda = first_generator(&basis, basis.j0());; ++lambda)
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
      for (int level = basis.j0(); level <= basis.j0()+1; level++) {
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

  cout << "- checking intersection of singular supports:" << endl;
  for (lambda = first_generator(&basis, basis.j0()+2);; ++lambda)
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

  return 0;
}
