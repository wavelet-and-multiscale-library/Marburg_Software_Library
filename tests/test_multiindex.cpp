#include <iostream>
#include <utils/array1d.h>
#include <utils/multiindex.h>

using std::cout;
using std::endl;
using namespace MathTL;

int main()
{
  cout << "Testing MathTL::MultiIndex ..." << endl;
  
  MultiIndex<int,2> lambda;
  cout << lambda << endl;

  return 0;
}
