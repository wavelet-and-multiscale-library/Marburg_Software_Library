#include <iostream>

#include <algebra/infinite_vector.h>
#include <algebra/sparse_matrix.h>
#include <utils/array1d.h>
#include <utils/fixed_array1d.h>

#include <interval/ds_basis.h>
#include <interval/p_basis.h>
#include <interval/jl_basis.h>
#include <interval/jl_support.h>
#include <interval/jl_evaluate.h>

#include <interval/i_index.h>
#include <cube/cube_basis.h>

using namespace std;
using namespace WaveletTL;

int main()
{
  cout << "Testing support calculations for wavelet bases on the cube..." << endl;

//   typedef DSBasis<3,3> Basis1D;
  typedef PBasis<2,2> Basis1D;
//   typedef JLBasis Basis1D;

  const unsigned int DIM = 2;
  typedef CubeBasis<Basis1D,DIM> Basis;
  typedef Basis::Index Index;

#if 1
  Basis basis;
#else
  FixedArray1D<bool,4> bc;
  bc[0] = bc[1] = true;
  bc[2] = bc[3] = true;
  Basis basis(bc);
#endif

#if 0
  cout << "- testing calculation of supports:" << endl;
  Basis::Support supp;
  for (Index lambda(first_generator<Basis1D,2,Basis>(&basis, basis.j0()));; ++lambda) {
//     support<Basis1D,2,Basis>(basis, lambda, supp);
    support(basis, lambda, supp);
    cout << lambda << " has support 2^{-" << supp.j << "}"
	 << "[" << supp.a[0] << "," << supp.b[0]
	 << "]x[" << supp.a[1] << "," << supp.b[1] << "]"
	 << endl;
    if (lambda == last_wavelet<Basis1D,2,Basis>(&basis, basis.j0()+1)) break;
  }  
#endif

#if 1
  cout << "- compute all intersecting wavelets:" << endl;
  for (Index lambda = first_generator<Basis1D,2,Basis>(&basis, basis.j0());; ++lambda)
    {
      cout << "  * for lambda=" << lambda << ":" << endl;
      typedef std::list<Index> SupportList;
      SupportList nus;
      intersecting_wavelets<Basis1D,2>(basis, lambda, basis.j0(), true, nus);
      for (SupportList::const_iterator it(nus.begin()); it != nus.end(); ++it) {
	Basis::Support supp;
	intersect_supports(basis, lambda, *it, supp);
	cout << "    nu=" << *it 
	     << " with support intersection "
	     << "2^{-" << supp.j << "}"
	     << "[" << supp.a[0] << "," << supp.b[0]
	     << "]x[" << supp.a[1] << "," << supp.b[1] << "]"
	     << endl;
      }
      for (int level = basis.j0(); level <= basis.j0()+2; level++) {
	intersecting_wavelets<Basis1D,2>(basis, lambda, level, false, nus);
	for (SupportList::const_iterator it(nus.begin()); it != nus.end(); ++it) {
	  Basis::Support supp;
	  intersect_supports(basis, lambda, *it, supp);
	  cout << "    nu=" << *it 
	       << " with support intersection "
	       << "2^{-" << supp.j << "}"
	       << "[" << supp.a[0] << "," << supp.b[0]
	       << "]x[" << supp.a[1] << "," << supp.b[1] << "]"
	       << endl;
	}
      }
      
      if (lambda == last_wavelet<Basis1D,2,Basis>(&basis, basis.j0())) break;
    }
#endif  
}
