#include <iostream>
#include <fstream>
#include <interval/spline_basis.h>
#include <galerkin/full_laplacian.h>

using namespace std;
using namespace WaveletTL;

int main()
{
  cout << "Testing FullLaplacian ..." << endl;

  SplineBasis<2,2> basis("P","",1,1,0,0); // PBasis, complementary b.c.'s
  FullLaplacian<2,2> delta(basis);

  cout << "* stiffness matrix on coarsest level j0=" << basis.j0() << ":" << endl
       << delta;
  
  delta.set_level(basis.j0()+1);
  cout << "* stiffness matrix on next level j0+1=" << basis.j0()+1 << ":" << endl
       << delta;
  
  

  return 0;
}
