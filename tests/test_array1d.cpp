#include <iostream>
#include <utils/array1d.h>
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
  cout << "Testing MathTL::Array1D ..." << endl;
  Array1D<double> a(2);
  cout << "- an Array1D<double>(2): " << a << endl;
  Array1D<TestClass> b(3);

  return 0;
}
