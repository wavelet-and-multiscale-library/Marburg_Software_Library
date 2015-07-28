#include <iostream>
#include <set>
#include <map>
#include <utils/array1d.h>
#include <utils/multiindex.h>
#include <algebra/infinite_vector.h>
#include <algebra/vector.h>

using namespace std;
using namespace MathTL;


int main()
{
  cout << "Testing MathTL::MultiIndex ..." << endl;
  
  typedef MultiIndex<int,1> MIndex1;
  typedef MultiIndex<int,2> MIndex2;
  typedef MultiIndex<int,3> MIndex3;
  
  MIndex3 lambda;
  cout << "- an empty 3-index: " << lambda << endl;
  MIndex3 lambda2;
  lambda2[0] = 3;
  lambda2[1] = 1;
  lambda2[2] = 4;
  cout << "- a nontrivial 3-index lambda2: " << lambda2 << endl;

  lambda[0] = 0;
  lambda[1] = 0;
  lambda[2] = 9;
  cout << "  * changed lambda into " << lambda << endl;
  cout << "  * check equality h == h2 : " << (lambda==lambda2) << endl;
  cout << "  * check non-equality h!=h2 : " << (lambda!=lambda2) << endl;
  cout << "  * check lex. ordering. h2.lex(h) : "
       << (lambda2.lex(lambda)) << ", h.lex(h2) : " << (lambda.lex(lambda2)) << endl;
  cout << "  * check ordering operator <. h2 < h : "
       << (lambda2 < lambda) << ", h < h2 : " << (lambda < lambda2) << endl;
  cout << "  * check increment:" << endl;
  cout << "lambda = " << lambda << "." << endl;
  cout << "++" << lambda; 
  cout << " = " << (++ lambda) << "." << endl;
  cout << "lambda = " << lambda << "." << endl;
  cout << lambda << "++ = " ;
  cout << (lambda++) << "." << endl;
  cout << "lambda = " << lambda << "." << endl;
  
  cout << "  * check number() :" << endl;
  cout << "  * 1D : " << endl;
  MIndex1 a;
  a[0]=2;
  for (unsigned int i=0; i<= 20; ++i)
  {
      cout << "i = " << i << " " << a << ".num = " << a.number() << endl;
      ++a;
  }
  
  MIndex2 aa;
  aa[0]=4; aa[1]=0;
  for (unsigned int i=0; i<= 20; ++i)
  {
      cout << "i = " << i << " " << aa << ".num = " << aa.number() << endl;
      ++aa;
  }
  
  MIndex3 aaa;
  aaa[0]=0; aaa[1]=3; aaa[2]=0;
  for (unsigned int i=0; i<= 20; ++i)
  {
      cout << "i = " << i << " " << aaa << ".num = " << aaa.number() << endl;
      ++aaa;
  }
  
#if 0
  cout << "  * check number() :" << endl;
  cout << "  * 1D : " << endl;
  MIndex1 a;
  a[0]=2;
  for (unsigned int i=0; i<= 20; ++i)
  {
      cout << "i = " << i << " " << a << ".num // .num2 = " << a.number() << " // " << a.number2() << endl;
      ++a;
  }
  
  MIndex2 aa;
  aa[0]=4; aa[1]=0;
  for (unsigned int i=0; i<= 20; ++i)
  {
      cout << "i = " << i << " " << aa << ".num // .num2 = " << aa.number() << " // " << aa.number2() << endl;
      ++aa;
  }
  
  MIndex3 aaa;
  aaa[0]=0; aaa[1]=3; aaa[2]=0;
  for (unsigned int i=0; i<= 20; ++i)
  {
      cout << "i = " << i << " " << aaa << ".num // .num2 = " << aaa.number() << " // " << aaa.number2() << endl;
      ++aaa;
  }
  
  clock_t tstart, tend;
  double time_small1(0), time_big1(0), time_small2(0), time_big2(0);
  int number_of_runs(30), small(500), big(8000), collect;
  
  // dim =1
  collect = 0;
  for (unsigned int run=0; run < number_of_runs; run++)
  {
      // test number and number2 on index 0 ... small
      tstart = clock();
      a[0]=0;
      for (unsigned int i=0; i < small; i++)
      {
          collect += a.number();
          a++;
      }
      tend = clock();
      time_small1 += (double)(tend-tstart)/CLOCKS_PER_SEC;
      
      tstart = clock();
      a[0]=0;
      for (unsigned int i=0; i < big; i++)
      {
          collect += a.number();
          a++;
      }
      tend = clock();
      time_big1 += (double)(tend-tstart)/CLOCKS_PER_SEC;
      
      tstart = clock();
      a[0]=0;
      for (unsigned int i=0; i < small; i++)
      {
          collect += a.number2();
          a++;
      }
      tend = clock();
      time_small2 += (double)(tend-tstart)/CLOCKS_PER_SEC;
      
      tstart = clock();
      a[0]=0;
      for (unsigned int i=0; i < big; i++)
      {
          collect += a.number2();
          a++;
      }
      tend = clock();
      time_big2 += (double)(tend-tstart)/CLOCKS_PER_SEC;
  }
  cout << "Speed up?? Dim = 1" << endl;
  cout <<" time_small1 // time_small2 = " << time_small1  << " // " << time_small2 << " seconds" << endl;
  cout <<" time_big1 // time_big2 = " << time_big1  << " // " << time_big2 << " seconds" << endl;
  
  // dim =2
  time_small1 = 0; time_small2 = 0;  time_big1 = 0; time_big2 = 0;
  collect = 0;
  for (unsigned int run=0; run < number_of_runs; run++)
  {
      // test number and number2 on index 0 ... small
      tstart = clock();
      aa[0]=0; aa[1]=0;
      for (unsigned int i=0; i < small; i++)
      {
          collect += aa.number();
          aa++;
      }
      tend = clock();
      time_small1 += (double)(tend-tstart)/CLOCKS_PER_SEC;
      
      tstart = clock();
      aa[0]=0; aa[1]=0;
      for (unsigned int i=0; i < big; i++)
      {
          collect += aa.number();
          aa++;
      }
      tend = clock();
      time_big1 += (double)(tend-tstart)/CLOCKS_PER_SEC;
      
      tstart = clock();
      aa[0]=0; aa[1]=0;
      for (unsigned int i=0; i < small; i++)
      {
          collect += aa.number2();
          aa++;
      }
      tend = clock();
      time_small2 += (double)(tend-tstart)/CLOCKS_PER_SEC;
      
      tstart = clock();
      aa[0]=0; aa[1]=0;
      for (unsigned int i=0; i < big; i++)
      {
          collect += aa.number2();
          aa++;
      }
      tend = clock();
      time_big2 += (double)(tend-tstart)/CLOCKS_PER_SEC;
  }
  cout << "Speed up?? Dim = 2" << endl;
  cout <<" time_small1 // time_small2 = " << time_small1  << " // " << time_small2 << " seconds" << endl;
  cout <<" time_big1 // time_big2 = " << time_big1  << " // " << time_big2 << " seconds" << endl;
  
  // dim =2
  time_small1 = 0; time_small2 = 0;  time_big1 = 0; time_big2 = 0;
  collect = 0;
  for (unsigned int run=0; run < number_of_runs; run++)
  {
      // test number and number2 on index 0 ... small
      tstart = clock();
      aaa[0]=0; aaa[1]=0; aaa[2]=0;
      for (unsigned int i=0; i < small; i++)
      {
          collect += aaa.number();
          aaa++;
      }
      tend = clock();
      time_small1 += (double)(tend-tstart)/CLOCKS_PER_SEC;
      
      tstart = clock();
      aaa[0]=0; aaa[1]=0; aaa[2]=0;
      for (unsigned int i=0; i < big; i++)
      {
          collect += aaa.number();
          aaa++;
      }
      tend = clock();
      time_big1 += (double)(tend-tstart)/CLOCKS_PER_SEC;
      
      tstart = clock();
      aaa[0]=0; aaa[1]=0; aaa[2]=0;
      for (unsigned int i=0; i < small; i++)
      {
          collect += aaa.number2();
          aaa++;
      }
      tend = clock();
      time_small2 += (double)(tend-tstart)/CLOCKS_PER_SEC;
      
      tstart = clock();
      aaa[0]=0; aaa[1]=0; aaa[2]=0;
      for (unsigned int i=0; i < big; i++)
      {
          collect += aaa.number2();
          aaa++;
      }
      tend = clock();
      time_big2 += (double)(tend-tstart)/CLOCKS_PER_SEC;
  }
  cout << "Speed up?? Dim = 3" << endl;
  cout <<" time_small1 // time_small2 = " << time_small1  << " // " << time_small2 << " seconds" << endl;
  cout <<" time_big1 // time_big2 = " << time_big1  << " // " << time_big2 << " seconds" << endl;
    
#endif
  lambda[0]=0; lambda[1]=1; lambda[2]=2;
    cout << "  * check indexmapping :" << endl;
  map<MIndex3, int> ind_map(indexmapping(lambda));
  cout << "lambda = " << lambda << endl
          << "indexmapping(" << lambda << ") = " << endl;
  for (map<MIndex3,int>::const_iterator it(ind_map.begin()); it!=ind_map.end(); ++it)
  {
      MIndex3 temp_ind (it->first);
      cout << it->first << ", " << it->second << "; .number() = " << temp_ind.number() << endl;
      if (it==ind_map.end()) break;
  }
    
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

  cout << "observe the ordering: " << endl;
  cout << lambda << " < " << lambda2 << " = " << (lambda < lambda2) << endl;
  cout << lambda << " < " << lambda3 << " = " << (lambda < lambda3) << endl;
  cout << lambda2 << " < " << lambda3 << " = " << (lambda2 < lambda3) << endl;
  
  cout << "  * making an InfiniteVector:" << endl;
  InfiniteVector<double, MIndex3> v;
  v[lambda] = 42.0;
  v[lambda2] = 23.0;
  v[lambda3] = 3.14;
  cout << v;

  cout << "  * making a map based on a 2-index:" << endl;
  
  MIndex2 n;
  map<MIndex2,double> d;
   n[0] = 6; n[1] = 6;
  d[n] = 2.0;
  n[0] = 2; n[1] = 4;
  d[n] = 1.0;
 
  for (map<MIndex2,double>::const_iterator it(d.begin()); it != d.end(); ++it)
    cout << it->first << ", " << it->second << endl;

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
  MI gamma, delta;
  gamma[0] = 3; gamma[1] = 2;
  delta[0] = 1; delta[1] = 4;
  cout << "- the 2-index gamma=" << gamma << " has degree "
       << multi_degree(gamma) << " and factorial " << multi_factorial(gamma) << endl;
  cout << "- another 2-index delta=" << delta << " yields gamma^delta="
       << multi_power(gamma, delta) << endl;
  gamma[0] = 4; gamma[1] = 5; delta[0]=1; delta[1] = 2;
    cout << "- multi-binomial ("<< gamma << " " << delta << ") = " 
       << multi_binomial(gamma, delta) << endl;

  cout << "- all 2-indices with degree 4:" << endl;
  set<MI> degree4(degree_indices<unsigned int, 2>(4));
  for (set<MI>::const_iterator it(degree4.begin()); it != degree4.end(); ++it)
    cout << *it << endl;
  return 0;
}
