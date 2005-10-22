#include <iostream>

#include <algebra/vector.h>
#include <algebra/infinite_vector.h>
#include <utils/function.h>
#include <utils/fixed_array1d.h>
#include <numerics/bvp.h>

#include <interval/ds_basis.h>
#include <cube/cube_basis.h>
#include <galerkin/cube_equation.h>

using namespace MathTL;
using namespace WaveletTL;

int main()
{
  cout << "Testing wavelet-Galerkin solution of an elliptic equation on the cube ..." << endl;

  ConstantFunction<2> rhs(Vector<double>(1, "1.0"));
  PoissonBVP<2> poisson(&rhs);

  const int d  = 2;
  const int dT = 2; // be sure to use a continuous dual here, otherwise the RHS test will fail
  typedef DSBasis<d,dT> Basis1D;
  typedef CubeBasis<Basis1D,2> CBasis;
  typedef CBasis::Index Index;

  FixedArray1D<bool,4> bc;
  bc[0] = bc[1] = bc[2] = bc[3] = true;

  CubeEquation<Basis1D,2,CBasis> eq(poisson, bc);

  InfiniteVector<double, Index> coeffs;

#if 0
  coeffs[first_generator<Basis1D,2,CBasis>(&eq.basis(), eq.basis().j0())] = 1.0;
  coeffs[last_generator<Basis1D,2,CBasis>(&eq.basis(), eq.basis().j0())] = 2.0;
  coeffs[first_wavelet<Basis1D,2,CBasis>(&eq.basis(), eq.basis().j0())] = 3.0;
  coeffs[last_wavelet<Basis1D,2,CBasis>(&eq.basis(), eq.basis().j0())] = 4.0;
  coeffs[first_wavelet<Basis1D,2,CBasis>(&eq.basis(), eq.basis().j0()+1)] = 5.0;
  cout << "- a coefficient set:" << endl
       << coeffs << endl;
  eq.rescale(coeffs, -1);
  cout << "- after rescaling with D^{-1}:" << endl
       << coeffs << endl;
#endif

//   eq.RHS(1e-2, coeffs);

  return 0;
}
