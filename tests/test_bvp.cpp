#include <iostream>
#include <numerics/bvp.h>
#include <geometry/chart.h>
#include <geometry/atlas.h>
#include <algebra/vector.h>

using std::cout;
using std::endl;
using namespace MathTL;

int main()
{
  cout << "Testing MathTL::EllipticBVP ..." << endl;

//   AffineLinearMapping<2> I;
//   Array1D<Chart<2>*> kappas(1);
//   kappas[0] = &I;
//   SymmetricMatrix<bool> S(1);
//   S(0, 0) = true;
//   Atlas<2> A(kappas, S);
//   FixedArray1D<int,4> bc;
//   bc[0] = bc[1] = bc[2] = bc[3] = 1; // homogeneous Dirichlet conditions everywhere
//   Array1D<FixedArray1D<int,4> > bchelp(1);
//   bchelp[0] = bc;

//   /*
//     the Poisson equation on any domain
//       -Delta u = 1 on Omega
//   */
//   ConstantFunction<2> one(Vector<double>(1, "1.0"));
//   ZeroFunction<2> zero;
//   EllipticBVP<2> bvp(&one, &zero, &one);

//   cout << bvp.a(Point<2>(1.0, 2.9)) << endl;
//   cout << bvp.q(Point<2>(1.0, 2.9)) << endl;
//   cout << bvp.f(Point<2>(1.0, 2.9)) << endl;

//   /*
//     the same equation
//    */
//   PoissonBVP<2> poisson(&one);
//   cout << poisson.a(Point<2>(1.0, 2.9)) << endl;
//   cout << poisson.q(Point<2>(1.0, 2.9)) << endl;
//   cout << poisson.f(Point<2>(1.0, 2.9)) << endl;

  return 0;
}
