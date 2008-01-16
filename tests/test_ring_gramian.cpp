#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <utils/function.h>
#include <ring/ring_basis.h>
#include <galerkin/ring_gramian.h>

using namespace std;
using namespace MathTL;
using namespace WaveletTL;


int main()
{
  cout << "Testing wavelet bases on the ring-shaped domain ..." << endl;

  const int d  = 2;
  const int dt = 2;

  typedef RingBasis<d,dt,1,1> Basis;
  Basis basis;

  typedef Basis::Index Index;

  
  return 0;
}
