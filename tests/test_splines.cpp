#include <iostream>
#include <numerics/splines.h>

using std::cout;
using std::endl;
using namespace MathTL;

int main()
{
  cout << "Testing MathTL::Spline ..." << endl;

  Spline<2> f;
  f.dump();
}
