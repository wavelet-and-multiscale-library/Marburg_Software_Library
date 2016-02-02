#include <iostream>
#include <algebra/tensor.h>

using std::cout;
using std::endl;
using namespace MathTL;

int main()
{
  cout << "Testing MathTL::Tensor ..." << endl;
  cout << "- a Tensor<1,1>: " << Tensor<1,1>() << endl;
  cout << "- a Tensor<1,2>: " << Tensor<1,2>() << endl;
  cout << "- a Tensor<2,1>: " << Tensor<2,1>() << endl;
  cout << "- a Tensor<2,2>: " << Tensor<2,2>() << endl;
  cout << "- memory consumption of a Tensor<1,1>: " << Tensor<1,1>().memory_consumption() << endl;
  cout << "- memory consumption of a Tensor<1,2>: " << Tensor<1,2>().memory_consumption() << endl;
  cout << "- memory consumption of a Tensor<2,1>: " << Tensor<2,1>().memory_consumption() << endl;
  cout << "- memory consumption of a Tensor<2,2>: " << Tensor<2,2>().memory_consumption() << endl;
}
