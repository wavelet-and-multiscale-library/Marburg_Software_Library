#include <iostream>
#include <geometry/point.h>

using std::cout;
using std::endl;
using namespace MathTL;

int main()
{
  cout << "Testing MathTL::Point ..." << endl;
  cout << "- a Point<3>: " << Point<3>() << endl;
  cout << "- a nontrivial Point<2>: " << Point<2>(1.0, 2.0) << endl;
}
