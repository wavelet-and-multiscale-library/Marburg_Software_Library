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

// the radial hat function
class RingFunction2
  : public Function<2>
{
public:
  inline double value(const Point<2>& p, const unsigned int component = 0) const {
    const double r = sqrt(p[0]*p[0]+p[1]*p[1]);
    return max(1-2*fabs(1.5-r),0.);
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

  typedef RingBasis<d,dt,1,1> Basis;
  Basis basis;

  typedef Basis::Index Index;

  Function<2>* u = 0;
  const int testcase = 2;
  switch(testcase) {
  case 2:
    u = new RingFunction2();
    break;
  default: // also for testcase 1
    u = new RingFunction1();
  }

  RingGramian<d,dt,1,1> G(basis, InfiniteVector<double,Index>());


  if (u) delete u;

  return 0;
}
