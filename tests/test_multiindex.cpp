#include <iostream>
#include <set>
#include <map>
#include <utils/array1d.h>
#include <utils/multiindex.h>

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
  cout << "  * making a map out of the two multiindices:" << endl;
  map<MIndex,double> c;
  c.insert(std::make_pair<MIndex,double>(lambda,42.0));
  c.insert(std::make_pair<MIndex,double>(lambda2,23.0));
  for (map<MIndex,double>::const_iterator it(c.begin()); it != c.end(); ++it)
    cout << it->first << ", " << it->second << endl;

  return 0;
}
