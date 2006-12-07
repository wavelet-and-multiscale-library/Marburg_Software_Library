#include <iostream>
#include <map>
#include <utils/map_tools.h>

using namespace std;
using namespace MathTL;

int main()
{
  cout << "Testing map_tools..." << endl;

  map<int,double> m1, m2, m3;
  m1[0] = 1;
  m1[2] = -3.14;
  m1[3] = -1;
  m2[1] = 2;
  m2[2] = 1;
  m2[10] = 1;
 
  cout << "* two maps, m1=" << endl;
  for (map<int,double>::const_iterator it(m1.begin()); it != m1.end(); ++it) {
    cout << it->first << ": " << it->second << endl;
  }
  cout << "  m2=" << endl;
  for (map<int,double>::const_iterator it(m2.begin()); it != m2.end(); ++it) {
    cout << it->first << ": " << it->second << endl;
  }

  add_maps(m1, m2, m3, 2.0, 3.0);

  cout << "  m3=2*m1+3*m2=" << endl;
  for (map<int,double>::const_iterator it(m3.begin()); it != m3.end(); ++it) {
    cout << it->first << ": " << it->second << endl;
  }
  
  return 0;
}
