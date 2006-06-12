#include <iostream>
#include <fstream>
#include <set>

#include <algebra/vector.h>
#include <algebra/infinite_vector.h>
#include <utils/function.h>
#include <utils/fixed_array1d.h>
#include <numerics/bvp.h>
#include <numerics/corner_singularity.h>
#include <geometry/sampled_mapping.h>

#include <interval/ds_basis.h>
#include <interval/p_basis.h>

#define _WAVELETTL_LDOMAINBASIS_VERBOSITY 1
#include <Ldomain/ldomain_basis.h>
#include <galerkin/ldomain_equation.h>

using namespace std;
using namespace MathTL;
using namespace WaveletTL;

int main()
{
  cout << "Testing wavelet-Galerkin solution of an elliptic equation on the L-shaped domain ..." << endl;

  const int d  = 2;
  const int dT = 2;
//   typedef DSBasis<d,dT> Basis1D;
  typedef PBasis<d,dT> Basis1D;
  typedef LDomainBasis<Basis1D> LBasis;
  typedef LBasis::Index Index;

  CornerSingularity    u_sing(Point<2>(0,0), 0.5, 1.5);
  CornerSingularityRHS f_sing(Point<2>(0,0), 0.5, 1.5);
  PoissonBVP<2> poisson(&f_sing);

  LDomainEquation<Basis1D> eq(&poisson);

  InfiniteVector<double, Index> coeffs;
  eq.RHS(1e-8, coeffs);

  return 0;
}
