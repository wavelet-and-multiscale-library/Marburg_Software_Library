#include <iostream>
#include <fstream>
#include <sstream>
#include <algebra/vector.h>
#include <utils/function.h>
#include <interval/spline_basis.h>
#include <galerkin/full_gramian.h>
#include <Rd/cdf_utils.h>
#include <numerics/iteratsolv.h>
#include <numerics/quadrature.h>
#include <numerics/schoenberg_splines.h>

using namespace std;
using namespace MathTL;
using namespace WaveletTL;

int main()
{
  cout << "Testing FullGramian ..." << endl;

  const unsigned int d = 2;
  const unsigned int dT = 2;

  SplineBasis<d,dT> basis("P","",1,1,0,0); // PBasis, complementary b.c.'s
  FullGramian<d,dT> G(basis);

  cout << "* Gramian matrix on coarsest level j0=" << basis.j0() << ":" << endl
       << G;
  
  G.set_level(basis.j0()+1);
  cout << "* Gramian matrix on next level j0+1=" << basis.j0()+1 << ":" << endl
       << G;
  

  return 0;
}
