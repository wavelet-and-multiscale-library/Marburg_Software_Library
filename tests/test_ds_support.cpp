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
  
  Basis basis;
  Index lambda(first_generator(&basis, basis.j0()));
//   ++lambda;
//   for (int i = 1; i <= 4; i++, ++lambda);
//   lambda = last_generator(&basis, basis.j0());
//   lambda = first_wavelet(&basis, basis.j0()+1);
//   for (int i = 1; i <= 4; i++, ++lambda);
//   lambda = last_wavelet(&basis, basis.j0()+1);

  cout << "- point values at dyadic points of the " << lambda << " generator/wavelet:" << endl;
  evaluate(basis, lambda, true, 6).matlab_output(cout);

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
  
//   cout << "- calculating some support intersections:" << endl;
//   for (lambda = basis.firstGenerator(basis.j0());; ++lambda)
//     {
//       int j, k1, k2;
//       support(basis, lambda, k1, k2);
//       cout << "psi_lambda, lambda=" << lambda << " has support 2^{-"
// 	   << lambda.j()+lambda.e()
// 	   << "}["
// 	   << k1
// 	   << ","
// 	   << k2
// 	   << "]"
// 	   << endl;

//       cout << "support intersection with first generator on level j0: ";
//       bool inter = intersect_supports(basis, lambda, basis.firstGenerator(basis.j0()), j, k1, k2);
//       if (inter)
// 	cout << "2^{-" << j << "}[" << k1 << "," << k2 << "]" << endl;
//       else
// 	cout << "none" << endl;
      
//       if (lambda == basis.lastWavelet(basis.j0())) break;
//     }

  return 0;
}
