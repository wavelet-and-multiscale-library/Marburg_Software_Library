#include <cstdlib>
#include <cmath>
#include <set>
#include <iostream>
#include <algebra/infinite_vector.h>

using std::cout;
using std::endl;
using namespace MathTL;

class Squares
  : public InfiniteDiagonalMatrix<float>
{
public:
  double diag(const int& i) const
  {
    return pow(i,2);
  }
};

int main()
{
  cout << "Testing the InfiniteVector class ..." << endl;

  InfiniteVector<float,long int> s;
  cout << "- a zero vector:" << endl
       << s << endl;

  cout << "- writing access on s:" << endl;
  s[1] = 2;
  cout << "  (size after writing the first element: " << s.size() << ")" << endl;
  s[3] = 42;
  cout << "  (size after writing the second element: " << s.size() << ")" << endl;
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

  cout << "- preparing a large random vector for the NCOARSE routine with size ";
  InfiniteVector<float> v, w;
  for (unsigned int i=0; i < 1000; i++)
    {
      v[i] = (double)rand()/(double)RAND_MAX;
    }
  cout << v.size()
       << " and ||v||_2=" << l2_norm(v) << endl;

  double eps = 0.1;
  cout << "- COARSE(" << eps << ",w) yields w with ";
  v.COARSE(eps,w);
  cout << w.size() << " entries and ||v-w||_2=" << l2_norm(v-w) << endl;
  eps = 1.0;
  cout << "- COARSE(" << eps << ",w) yields w with ";
  v.COARSE(eps,w);
  cout << w.size() << " entries and ||v-w||_2=" << l2_norm(v-w) << endl;
  eps = 10.0;
  cout << "- COARSE(" << eps << ",w) yields w with ";
  v.COARSE(eps,w);
  cout << w.size() << " entries and ||v-w||_2=" << l2_norm(v-w) << endl;
  
  cout << "- some weak \\ell_\\tau norms of v:" << endl;
  for (double tau(1.8); tau >= 0.2; tau -= 0.2)
    {
      cout << "  tau=" << tau << ", ||v||_{\\ell^w_\\tau}=" << v.weak_norm(tau) << endl;
    }

  v.clear();
  v[0] = 1e-10;
  v[1] = 0.5;
  v[2] = -1e-5;
  cout << "- another vector v:" << endl << v;
  v.compress(1e-2);
  cout << "- compressing with eta=1e-2:" << endl << v;

  v.clear();
  v[0] = 123;
  v[2] = 345;
  v[4] = -678;
  cout << "- another vector v:" << endl << v;
  std::set<int> supp;
  v.support(supp);
  cout << "- v has the support" << endl;
  for (std::set<int>::const_iterator it = supp.begin(); it != supp.end(); ++it)
    cout << *it << endl;
  std::set<int> Lambda;
  Lambda.insert(2);
  Lambda.insert(0);
  Lambda.insert(-1);
  v.clip(Lambda);
  cout << "- v clipped to an index set:" << endl << v;

  v.add_coefficient(0, 1.5);
  cout << "- added something to the first coefficient of v:" << endl << v;

  v.add_coefficient(2, -345);
  cout << "- added something to the second coefficient of v:" << endl << v;

  v.clear();
  v[0] = 1;
  v[1] = 2;
  v[2] = 3;
  v[3] = 4;
  w.clear();
  w[0] = 1;
  w[1] = -1;
  w[2] = -1;
  w[3] = 1;
  const double atol = 1;
  const double rtol = 1;
  cout << "- vectors v=" << endl << v << "  and w=" << endl << w;
  cout << "  weighted root mean square norm of v ("
       << "atol=" << atol << ", rtol=" << rtol << "): "
       << v.wrmsqr_norm(atol, rtol, w, w) << endl;

  Squares S;
  w.scale(&S);
  cout << "w weighted with an instance of Squares: " << endl << w;
  
  return 0;
}
