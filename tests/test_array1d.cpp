#include <iostream>
#include <utils/array1d.h>

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
  cout << "Testing MathTL::Array1D ..." << endl;
  Array1D<double> a(2);
  Array1D<TestClass> b(3);
}
