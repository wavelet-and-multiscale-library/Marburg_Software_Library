#include <iostream>
#include <utils/fixed_array1d.h>
#include <io/vector_io.h>

using std::cout;
using std::endl;
using namespace MathTL;

class TestClass
{
 public:
  TestClass()
  {
    cout << "TestClass::TestClass() called" << endl;
  }

  ~TestClass()
  {
    cout << "TestClass::~TestClass() called" << endl;
  }
};

int main()
{
  cout << "Testing MathTL::FixedArray1D ..." << endl;
  FixedArray1D<double,2> a;
  a[0] = 2; a[1] = -1;
  cout << "- a FixedArray1D<double,2>: " << a << endl;
  FixedArray1D<TestClass,3> b;

  return 0;
}
