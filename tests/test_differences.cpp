#include <iostream>
#include <cmath>
#include <algebra/polynomial.h>
#include <algebra/infinite_vector.h>
#include <numerics/differences.h>

using std::cout;
using std::endl;

using namespace MathTL;

int main()
{
  cout << "Testing finite difference operators ..." << endl;
  
  Polynomial<double> p;
  p.set_coefficient(2, 1.0); // x^2

  InfiniteVector<double, int> values;
  for (int i = -4; i <= 4; i++)
    values[i] = p.value(i);

  cout << "- a sampled polynomial on a 1D grid:" << endl;
  cout << values;

  cout << "- result of forward_difference<0>:" << endl;
  cout << forward_difference<0>(values);

  cout << "- result of forward_difference<1>:" << endl;
  cout << forward_difference<1>(values);

  cout << "- result of forward_difference<2>:" << endl;
  cout << forward_difference<2>(values);

  cout << "- result of forward_difference<3>:" << endl;
  cout << forward_difference<3>(values);

  cout << "- result of backward_difference<0>:" << endl;
  cout << backward_difference<0>(values);

  cout << "- result of backward_difference<1>:" << endl;
  cout << backward_difference<1>(values);

  cout << "- result of backward_difference<2>:" << endl;
  cout << backward_difference<2>(values);

  return 0;
}
