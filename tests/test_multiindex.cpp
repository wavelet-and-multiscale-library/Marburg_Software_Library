#include <iostream>
#include <set>
#include <map>
#include <utils/array1d.h>
#include <utils/multiindex.h>
#include <algebra/infinite_vector.h>

using namespace std;
using namespace MathTL;

int main()
{
  cout << "Testing MathTL::MultiIndex ..." << endl;
  
  typedef MultiIndex<int,3> MIndex;
  MIndex lambda;
  cout << "- an empty3-index: " << lambda << endl;
  MIndex lambda2;
  lambda2[0] = 3;
  lambda2[1] = 1;
  lambda2[2] = 4;
  cout << "- a nontrivial 3-index lambda2: " << lambda2 << endl;

  lambda[0] = 3;
  lambda[1] = 1;
  cout << "  * changed lambda into " << lambda << endl;
  cout << "  * check equality: " << (lambda==lambda2) << endl;
  cout << "  * check non-equality: " << (lambda!=lambda2) << endl;
  cout << "  * check lex. order operator <: "
       << (lambda2 < lambda) << ", " << (lambda < lambda2) << endl;

  cout << "  * making a set out of the two multiindices:" << endl;
  set<MIndex> M;
  M.insert(lambda);
  M.insert(lambda2);
  for (set<MIndex>::const_iterator it(M.begin()); it != M.end(); ++it)
    cout << *it << endl;

  cout << "  * making a map:" << endl;
  map<MIndex,double> c;
  MIndex lambda3;
  lambda3[0] = -2;
  c.insert(std::make_pair<MIndex,double>(lambda, 42.0));
  c.insert(std::make_pair<MIndex,double>(lambda2, 23.0));
  c.insert(std::make_pair<MIndex,double>(lambda3, 3.14));
  for (map<MIndex,double>::const_iterator it(c.begin()); it != c.end(); ++it)
    cout << it->first << ", " << it->second << endl;

  cout << "  * making an InfiniteVector:" << endl;
  InfiniteVector<double, MIndex> v;
  v[lambda] = 42.0;
  v[lambda2] = 23.0;
  v[lambda3] = 3.14;
  cout << v;

  cout << "  * making a map based on a 2-index:" << endl;
  typedef MultiIndex<int,2> MIndex2;
  MIndex2 n;
  n[0] = 2; n[1] = 4;
  map<MIndex2,double> d;
  d[n] = 1.0;
  n[0] = -6; n[1] = 6;
  d[n] = -1.0;
  for (map<MIndex2,double>::const_iterator it(d.begin()); it != d.end(); ++it)
    cout << it->first << ", " << it->second << endl;

  return 0;
}
