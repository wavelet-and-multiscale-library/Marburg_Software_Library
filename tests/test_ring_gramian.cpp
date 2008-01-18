#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <utils/function.h>
#include <ring/ring_basis.h>
#include <galerkin/ring_gramian.h>
#include "ring_functions.h"

using namespace std;
using namespace MathTL;
using namespace WaveletTL;


int main()
{
  cout << "Testing wavelet bases on the ring-shaped domain ..." << endl;

  const int d  = 2;
  const int dt = 2;
  const int s0 = 0;
  const int s1 = 0;

  const double r0 = 0.5;
  const double r1 = 2.0;

  typedef RingBasis<d,dt,s0,s1> Basis;
  Basis basis(r0, r1);

  typedef Basis::Index Index;

  Function<2>* u = 0;
  const int testcase = 2;
  switch(testcase) {
  case 2:
    u = new RingFunction2(r0, r1);
    break;
  default: // also for testcase 1
    u = new RingFunction1(r0, r1);
  }

  RingGramian<d,dt,s0,s1> G(basis, InfiniteVector<double,Index>());


  if (u) delete u;

  return 0;
}
