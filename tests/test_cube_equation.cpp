#include <iostream>

#include <algebra/vector.h>
#include <utils/function.h>
#include <numerics/bvp.h>
#include <galerkin/cube_equation.h>

using namespace MathTL;
using namespace WaveletTL;

int main()
{
  cout << "Testing wavelet-Galerkin solution of an elliptic equation on the cube ..." << endl;

  ConstantFunction<2> rhs(Vector<double>(1, "1.0"));
  PoissonBVP<2> poisson(&rhs);



  return 0;
}
