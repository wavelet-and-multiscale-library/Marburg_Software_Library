#include <cstdlib>
#include <iostream>
#include <algebra/infinite_vector.h>

using std::cout;
using std::endl;
using namespace MathTL;

int main()
{
  cout << "Testing the InfiniteVector class ..." << endl;

  InfiniteVector<float,long int> s;
  cout << "- a zero vector:" << endl
       << s << endl;

  cout << "- writing access on s:" << endl;
  s[1] = 2;
//  cout << "  (size after writing the first element: " << s.size() << ")" << endl;
  s[3] = 42;
//   cout << "  (size after writing the second element: " << s.size() << ")" << endl;
  cout << s;

  cout << "- copy constructor t(s):" << endl;
  InfiniteVector<float,long int> t(s);
  cout << t;

  cout << "- are the two vectors equal?" << endl;
  if (t == s)
    cout << "  ... yes!" << endl;
  else
    cout << "  ... no!" << endl;

  cout << "- are the two vectors inequal?" << endl;
  if (t != s)
    cout << "  ... yes!" << endl;
  else
    cout << "  ... no!" << endl;

  cout << "- in place summation s+=t:" << endl;
  s += t;
  cout << s;

  cout << "- in place subtraction t-=s:" << endl;
  t -= s;
  cout << t;

  cout << "- in place multiplication s*=2:" << endl;
  s *= 2;
  cout << s;
  
  cout << "- in place division s/=3:" << endl;
  s /= 3;
  cout << s;

  cout << "- ell_p norms of s:" << endl;
  cout << "  ||x||_2 = " << l2_norm(s)
       << ", ||x||_1 = " << l1_norm(s)
       << ", ||x||_infinity = " << linfty_norm(s) << endl;

  cout << "- external arithmetic functionality:" << endl;
  InfiniteVector<float,long int> sa, sb;
  sa[1] = 23; sa[2] = sa[3] = 10; sb[1] = -1.5; sb[3] = 3; sb[4] = 8;
  cout << "  a=" << endl << sa
       << "  b=" << endl << sb;
  swap(sa,sb);
  cout << "  after swapping, a=" << endl << sa << "  b=" << endl << sb;
  cout << "  a+b=" << endl << sa+sb
       << "  a-b=" << endl << sa-sb;
  cout << "  a*b=" << sa*sb << endl;
  cout << "  mean value of a: " << mean_value(sa) << endl;

//   cout << "- preparing a large random vector for the NCOARSE routine with size ";
//   InfiniteVector<float,int> v,w;
//   for (int i=0; i < 10000; i++) v(i) = (double)rand()/(double)RAND_MAX;
//   cout << v.size() << endl;
//   double eps = 0.1;
//   cout << "- NCOARSE(" << eps << ",w) yields w with ";
//   v.NCOARSE(eps,w);
//   cout << w.size() << " entries and ||v-w||_2=" << norm(v-w)
//        << " (||v||=" << norm(v) << ")" << endl;
//   eps = 1.0;
//   cout << "- NCOARSE(" << eps << ",w) yields w with ";
//   v.NCOARSE(eps,w);
//   cout << w.size() << " entries and ||v-w||_2=" << norm(v-w)
//        << " (||v||=" << norm(v) << ")" << endl;
//   eps = 10.0;
//   cout << "- NCOARSE(" << eps << ",w) yields w with ";
//   v.NCOARSE(eps,w);
//   cout << w.size() << " entries and ||v-w||_2=" << norm(v-w)
//        << " (||v||=" << norm(v) << ")" << endl;
  
  
//   cout << "- some weak \ell_\tau norms of v:" << endl;
//   for (double tau(1.8); tau >= 0.2; tau -= 0.2)
//     {
//       cout << "  tau=" << tau << ", ||v||_{\\ell^w_\\tau}=" << v.weak_norm(tau) << endl;
//     }

  return 0;
}
