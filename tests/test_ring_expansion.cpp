#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algebra/infinite_vector.h>
#include <interval/periodic.h>
#include <utils/function.h>
#include <ring/ring_basis.h>

using namespace std;
using namespace MathTL;
using namespace WaveletTL;

// the constant function f(x)=23
class RingFunction1
  : public Function<2>
{
public:
  inline double value(const Point<2>& p, const unsigned int component = 0) const {
    return 23;
  }
  void vector_value(const Point<2> &p, Vector<double>& values) const {
    values.resize(1, false);
    values[0] = value(p);
  }
};

int main()
{
  cout << "Testing wavelet bases on the ring-shaped domain ..." << endl;

  const int d  = 2;
  const int dt = 2;

  typedef RingBasis<d,dt,1,0> Basis;
  Basis basis;

  typedef Basis::Index Index;

  Function<2>* u = 0;
  const int testcase = 1;
  switch(testcase) {
  default:
    u = new RingFunction1();
  }
  
  InfiniteVector<double,Index> coeffs;
  
  const int j0 = basis.j0();
  const int jmax = j0+1;
  
  basis.expand(u, true, j0, coeffs);
  cout << "- integrals of the test function u against all wavelets up to level jmax="
       << j0 << ":" << endl
       << coeffs << endl;
  
  if (u) delete u;
  
  return 0;
}
