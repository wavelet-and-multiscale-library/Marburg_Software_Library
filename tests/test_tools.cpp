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

  return 0;
}
