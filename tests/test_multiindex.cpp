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
  
  typedef MultiIndex<int,3> MIndex3;
  MIndex3 lambda;
  cout << "- an empty 3-index: " << lambda << endl;
  MIndex3 lambda2;
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
  set<MIndex3> M;
  M.insert(lambda);
  M.insert(lambda2);
  for (set<MIndex3>::const_iterator it(M.begin()); it != M.end(); ++it)
    cout << *it << endl;

  cout << "  * making a map:" << endl;
  map<MIndex3,double> c;
  MIndex3 lambda3;
  lambda3[0] = -2;
  c.insert(std::make_pair<MIndex3,double>(lambda, 42.0));
  c.insert(std::make_pair<MIndex3,double>(lambda2, 23.0));
  c.insert(std::make_pair<MIndex3,double>(lambda3, 3.14));
  for (map<MIndex3,double>::const_iterator it(c.begin()); it != c.end(); ++it)
    cout << it->first << ", " << it->second << endl;

  cout << "  * making an InfiniteVector:" << endl;
  InfiniteVector<double, MIndex3> v;
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

  typedef MultiIndex<int,1> MIndex1;
  MIndex1 alpha1, beta1;
  alpha1[0] = 0;
  beta1[0] = 2;
  cout << "- creating the cuboid between " << alpha1 << " and " << beta1 << ":" << endl;
  set<MIndex1> cuboid1(cuboid_indices<int, 1>(alpha1, beta1));
  for (set<MIndex1>::const_iterator it(cuboid1.begin()); it != cuboid1.end(); ++it)
    cout << *it << endl;

  MIndex2 alpha, beta;
  alpha[0] = alpha[1] = 0;
  beta[0] = 2; beta[1] = 1;
  cout << "- creating the cuboid between " << alpha << " and " << beta << ":" << endl;
  set<MIndex2> cuboid(cuboid_indices<int, 2>(alpha, beta));
  for (set<MIndex2>::const_iterator it(cuboid.begin()); it != cuboid.end(); ++it)
    cout << *it << endl;

  MIndex3 alpha3, beta3;
  alpha3[0] = alpha3[1] = alpha3[2] = 0;
  beta3[0] = 2; beta3[1] = 1; beta3[2] = 2;
  cout << "- creating the cuboid between " << alpha3 << " and " << beta3 << ":" << endl;
  set<MIndex3> cuboid3(cuboid_indices<int, 3>(alpha3, beta3));
  for (set<MIndex3>::const_iterator it(cuboid3.begin()); it != cuboid3.end(); ++it)
    cout << *it << endl;

  typedef MultiIndex<unsigned int, 2> MI;
  cout << "- all 2-indices with degree 4:" << endl;
  set<MI> degree4(degree_indices<2>(4));
  for (set<MI>::const_iterator it(degree4.begin()); it != degree4.end(); ++it)
    cout << *it << endl;

  MI gamma, delta;
  gamma[0] = 3; gamma[1] = 2;
  delta[0] = 1; delta[1] = 4;
  cout << "- the 2-index gamma=" << gamma << " has degree "
       << multi_degree(gamma) << " and faculty " << multi_faculty(gamma) << endl;
  cout << "- another 2-index delta=" << delta << " yields gamma^delta="
       << multi_power(gamma, delta) << endl;

  return 0;
}
