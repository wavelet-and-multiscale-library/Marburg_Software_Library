#include <cmath>
#include <iostream>
#include <utils/tiny_tools.h>

using namespace std;

int main()
{
  cout << "Testing some of the tiny tools..." << endl;

  for (int j = 0; j <= 10; j++) {
    double x1 = sqrt((double) (1<<j));
    double x2 = twotothejhalf(j);
    cout << "* 2^{" << j << "/2}= " << x2 << ", error: " << fabs(x2-x1) << endl;
  }

  cout << endl
       << "- checking positive dyadic modulo routine:" << endl;
  
  const int expo = 3;
  for (int j = -10; j <= 10; j++) {
    cout << j << " modulo 2^" << expo << "="
	 << dyadic_modulo(j,expo)
	 << endl;
  }

  return 0;
}
