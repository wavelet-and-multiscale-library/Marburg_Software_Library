#include <iostream>
#include <cmath>
#include <algebra/polynomial.h>
#include <algebra/infinite_vector.h>
#include <numerics/differences.h>
#include <numerics/multi_differences.h>
#include <utils/multiindex.h>

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

  Polynomial<double> q;
  q.set_coefficient(3, -1.0); // x^3
 
  InfiniteVector<double, MultiIndex<int, 2> > values2;
  for (int x = -2; x <= 3; x++)
    for (int y = -1; y <= 2; y++)
      values2.set_coefficient(MultiIndex<int, 2>(x, y), p.value(x)*q.value(y));

  cout << "- a sampled multivariate polynomial:" << endl;
  cout << values2;

  cout << "- result of multivariate_forward_difference<2, 0>:" << endl;
  cout << multivariate_forward_difference<2, 0>(values2);

  cout << "- result of multivariate_forward_difference<2, 1>:" << endl;
  cout << multivariate_forward_difference<2, 1>(values2);

  cout << "- result of multivariate_forward_difference<2, 0>:" << endl;
  cout << multivariate_backward_difference<2, 0>(values2);

  cout << "- result of multivariate_forward_difference<2, 1>:" << endl;
  cout << multivariate_backward_difference<2, 1>(values2);

  return 0;
}
