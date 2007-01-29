#include <iostream>
#include "algebra/fixed_vector.h"
#include "algebra/vector_norms.h"

using std::cout;
using std::endl;
using namespace MathTL;

int main()
{
  cout << "Testing the fixed vector class ..." << endl;

  const unsigned int dim(5);
  FixedVector<double,dim> v;
  cout << "- a zero fixed vector of dimension " << v.size() << ":" << endl
       << v << endl;

  cout << "- writing access on v:" << endl;
  v(1) = 1; v[3] = 42;
  cout << v << endl;

  cout << "- copy constructor w(v):" << endl;
  FixedVector<double,dim> w(v);
  cout << w << endl;

  cout << "- are the two vectors equal?" << endl;
  if (w == v)
    cout << "  ... yes!" << endl;
  else
    cout << "  ... no!" << endl;

  cout << "- in place summation v+=w:" << endl;
  v += w;
  cout << v << endl;

  cout << "- in place subtraction w-=v:" << endl;
  w -= v;
  cout << w << endl;

  cout << "- in place multiplication v*=2:" << endl;
  v *= 2;
  cout << v << endl;
  
  cout << "- in place division v/=3:" << endl;
  v /= 3;
  cout << v << endl;

  cout << "- external arithmetic functionality:" << endl;
  FixedVector<double,dim> a, b;
  a(1) = 23; a(2) = 10; b(1) = -1.5; b(2) = 3; b(4) = 8;
  cout << "  a=" << a << ", b=" << b << endl;
  cout << "  check lexicographical order a < b: ";
  if (a < b)
    cout << "true" << endl;
  else
    cout << "false" << endl;

  return 0;
}
