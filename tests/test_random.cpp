#include <iostream>
#include <map>
#include <utils/random.h>

using std::cout;
using std::endl;
using std::map;

using namespace MathTL;

int main()
{
  cout << "Testing Random class ..." << endl;
  
  for (int i(0); i < 10; i++)
    cout << random_double() << endl;
  
  for (int i(0); i < 10; i++)
    cout << random_integer(0,10) << endl;

  return 0;
}
