#include <iostream>
#include <fstream>
#include <set>
#include <map>
#include <Rd/haar_mask.h>
#include <Rd/cdf_mask.h>
#include <Rd/r_basis.h>
#include <Rd/r_index.h>
#include <Rd/r_mw_index.h>
#include <Rd/cdf_basis.h>
#include <Rd/cdf_utils.h>

using namespace std;
using namespace WaveletTL;

int main()
{
  cout << "Testing wavelet bases over R^d..." << endl;

  RIndex lambda;
  cout << "- testing wavelet index class RIndex:" << endl;
  cout << "  * default index: " << lambda << endl;
  lambda = RIndex(3,1,4);
  cout << "  * some index: " << lambda << endl;
  RIndex lambda2(lambda);
  cout << "  * testing copy constructor: " << lambda2 << endl;
  cout << "  * check equality: " << (lambda==lambda2) << endl;
  cout << "  * check non-equality: " << (lambda!=lambda2) << endl;
  ++lambda;
  cout << "  * check preincrement operator ++: " << lambda << endl;
  cout << "  * check lex. order operator <: "
       << (lambda2 < lambda) << ", " << (lambda < lambda2) << endl;
  cout << "  * making a set out of the two wavelet indices:" << endl;
  set<RIndex> M;
  M.insert(lambda);
  M.insert(lambda2);
  for (set<RIndex>::const_iterator it(M.begin()); it != M.end(); ++it)
    cout << *it << endl;
  cout << "  * making a map out of the two wavelet indices and some coeffs:" << endl;
  map<RIndex,double> c;
  c.insert(std::make_pair<RIndex,double>(lambda,42.0));
  c.insert(std::make_pair<RIndex,double>(lambda2,23.0));
  for (map<RIndex,double>::const_iterator it(c.begin()); it != c.end(); ++it)
    cout << it->first << ", " << it->second << endl;

  cout << "- constructing a CDF basis with d=";
  const int d = 2, dt = 2;
  cout << d << ", dt=" << dt << endl;
  RBasis<CDFMask_primal<d>, CDFMask_dual<d, dt> > basis;
  cout << "  + primal mask:" << endl;
  cout << basis.a() << endl;
  cout << "  + dual mask:" << endl;
  cout << basis.aT() << endl;
  cout << "  + primal wavelet mask:" << endl;
  for (MultivariateLaurentPolynomial<double, 1>::const_iterator it(basis.b().begin());
       it != basis.b().end(); ++it)
    cout << "k=" << it.index()[0] << ": " << *it << endl;
  cout << "  + dual wavelet mask:" << endl;
  for (MultivariateLaurentPolynomial<double, 1>::const_iterator it(basis.bT().begin());
       it != basis.bT().end(); ++it)
    cout << "k=" << it.index()[0] << ": " << *it << endl;

  InfiniteVector<double, RIndex> coeff;
  coeff[RIndex(2,0,0)] = 1.0;
  coeff[RIndex(2,0,1)] = -3.14;
  coeff[RIndex(2,0,3)] = 2.0;
  cout << "  * a small but nontrivial coefficient set:" << endl << coeff;
  cout << "  * result of DECOMPOSE:" << endl;
  InfiniteVector<double,RIndex> wcoeff;
  basis.decompose(coeff,0,wcoeff);
  cout << wcoeff;
  cout << "  * RECONSTRUCT that:" << endl;
  InfiniteVector<double,RIndex> rcoeff;
  basis.reconstruct(wcoeff,2,rcoeff);
  cout << rcoeff;

  cout << "  * again a small but nontrivial coefficient set:" << endl << coeff;
  cout << "  * result of DECOMPOSEt:" << endl;
  wcoeff.clear();
  basis.decompose_t(coeff,0,wcoeff);
  cout << wcoeff;
  cout << "  * RECONSTRUCTt that:" << endl;
  rcoeff.clear();
  basis.reconstruct_t(wcoeff,2,rcoeff);
  cout << rcoeff;

#if 0
  cout << "- evaluating some primal and dual CDF functions:" << endl;
  CDFBasis<2, 2> basis2;
  lambda = RIndex(1,0,2);
  basis2.evaluate(0, lambda, true, -2, 2, 2).matlab_output(cout);
  basis2.evaluate(0, RIndex(1,1,0), true, -2, 2, 2).matlab_output(cout);
  
  InfiniteVector<double,RIndex> evalcoeffs;
  evalcoeffs[RIndex(0,0,0)] = 1.0;
  evalcoeffs[RIndex(1,0,2)] = 1.0;
  evalcoeffs[RIndex(0,1,0)] = 1.0;
  basis2.evaluate(0, evalcoeffs, true, -2, 2, 2).matlab_output(cout);
#endif

#if 0
  cout << "- evaluating primal CDF functions: write file 'primalwavelet.m'" << endl;
  std::ofstream fs("primalwavelet.m");
  lambda = RIndex(2,1,-2);
  int k1, k2;
  support<2,2>(lambda, true, k1, k2);
  basis2.evaluate(0, lambda, true,
		  (int)floor(ldexp(1.0, -lambda.j())*k1), (int)ceil(ldexp(1.0, -lambda.j())*k2),
		  10).matlab_output(fs);
  fs.close();

  cout << "- evaluating dual CDF functions: write file 'dualwavelet.m'" << endl;
  CDFBasis<1, 3> basis3;

  fs.open("dualwavelet.m");
  lambda = RIndex(2,1,-2);
  support<1,3>(lambda, false, k1, k2);
  basis3.evaluate(0, lambda, false,
		  (int)floor(ldexp(1.0, -lambda.j())*k1), (int)ceil(ldexp(1.0, -lambda.j())*k2),
		  10).matlab_output(fs);
#endif

  cout << "- testing multiwavelet index class RMWIndex:" << endl;
  RMWIndex rmwindex;
  cout << "  * default index: " << rmwindex << endl;
  rmwindex = RMWIndex(42, 1, 1, -23);
  cout << "  * a crude index: " << rmwindex << endl;

  return 0;
}
