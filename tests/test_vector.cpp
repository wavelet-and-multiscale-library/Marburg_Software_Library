#include <iostream>
#include "algebra/vector.h"

using std::cout;
using std::endl;
using namespace MathTL;

int main()
{
  cout << "Testing the vector class ..." << endl;

  int dim(5);
  Vector<double> v(dim+1);
  cout << "- a zero vector of dimension " << v.size() << ":" << endl
       << v << endl;

//   v.resize(dim);
//   cout << "- resized to dimension " << v.dim() << ":" << endl
//        << v << endl;

//   cout << "- writing access on v:" << endl;
//   v(1) = 1; v(3) = 42;
//   cout << v << endl;

//   cout << "- copy constructor w(v):" << endl;
//   RawVector<double> w(v);
//   cout << w << endl;

//   cout << "- reference counter of w:" << endl
//        << w.refc() << endl;

//   cout << "- are the two vectors equal?" << endl;
//   if (w == v)
//     cout << "  ... yes!" << endl;
//   else
//     cout << "  ... no!" << endl;

//   cout << "- in place summation v+=w:" << endl;
//   v += w;
//   cout << v << endl;

//   cout << "- in place subtraction w-=v:" << endl;
//   w -= v;
//   cout << w << endl;

//   cout << "- in place multiplication v*=2:" << endl;
//   v *= 2;
//   cout << v << endl;
  
//   cout << "- in place division v/=3:" << endl;
//   v /= 3;
//   cout << v << endl;

//   cout << "- ell_p norms of v:" << endl;
//   cout << "  ||x||_2 = " << norm(v)
//        << ", ||x||_1 = " << norm(v,1)
//        << ", ||x||_infinity = " << maxnorm(v) << endl;

//   cout << "- external arithmetic functionality:" << endl;
//   RawVector<double> a(dim), b(dim);
//   a(1) = 23; a(2) = 10; b(1) = -1.5; b(2) = 3; b(4) = 8;
//   cout << "  a=" << a << ", b=" << b << endl;
//   cout << "  a+b=" << a+b << ", a-b=" << a-b << endl;
//   cout << "  a*b=" << a*b << endl;

//   cout << endl;
//   cout << "- testing the sparse (map) representation of vectors:" << endl;
//   RawVector<double,SparseArray1D<double> > s(dim);
//   cout << "- a zero vector of dimension " << s.dim() << " with size " << s.size() << ":" << endl
//        << s << endl;

//   cout << "- writing access on s:" << endl;
//   s(1) = 1;
//   cout << "  (size after writing the first element: " << s.size() << ")" << endl;
//   s(3) = 42;
//   cout << "  (size after writing the second element: " << s.size() << ")" << endl;
//   cout << s << endl;

//   cout << "- copy constructor t(s):" << endl;
//   RawVector<double,SparseArray1D<double> > t(s);
//   cout << t << endl;

//   cout << "- reference counter of t:" << endl
//        << t.refc() << endl;

//   cout << "- are the two vectors equal?" << endl;
//   if (t == s)
//     cout << "  ... yes!" << endl;
//   else
//     cout << "  ... no!" << endl;

//   cout << "- in place summation s+=t:" << endl;
//   s += t;
//   cout << s << endl;

//   cout << "- in place subtraction t-=s:" << endl;
//   t -= s;
//   cout << t << endl;

//   cout << "- in place multiplication s*=2:" << endl;
//   s *= 2;
//   cout << s << endl;
  
//   cout << "- in place division s/=3:" << endl;
//   s /= 3;
//   cout << s << endl;

//   cout << "- ell_p norms of s:" << endl;
//   cout << "  ||x||_2 = " << norm(s)
//        << ", ||x||_1 = " << norm(s,1)
//        << ", ||x||_infinity = " << maxnorm(s) << endl;

//   cout << "- external arithmetic functionality:" << endl;
//   RawVector<double,SparseArray1D<double> > sa(dim), sb(dim);
//   sa(1) = 23; sa(2) = 10; sb(1) = -1.5; sb(2) = 3; sb(4) = 8;
//   cout << "  a=" << sa << ", b=" << sb << endl;
//   cout << "  a+b=" << sa+sb << ", a-b=" << sa-sb << endl;
//   cout << "  a*b=" << sa*sb << endl;

  return 0;
}
