#include <iostream>
#include <fstream>
#include <set>
#include <map>
#include <Rd/r_q_index.h>
#include <Rd/r_index.h>
#include <Rd/r_frame.h>
#include <Rd/cdf_mask.h>
#include <Rd/quarklet_frame.h>
#include <numerics/cardinal_splines.h>
#include <Rd/cdf_basis.h>
#include <Rd/quarklet_mask.h>


using namespace std;
using namespace WaveletTL;


int main(){
  

  ofstream ausgabe("fktplot.m");
  
  //Ausgabe auf Matlab-file
  cout << "- evaluating some primal and dual CDF functions:" << endl;
  CDFBasis<2, 2> basis2;
  QuarkletFrame<3, 3> frame4;
   frame4.evaluate(0, RQIndex(1,0,1,0), true, -4, 4, 5).matlab_output(ausgabe);
   ausgabe << "figure;\nplot(x,y);"
              << "title('Superplot');" << endl;
  ausgabe.close();
  
  CDFBasis<2, 2>  basis;
  cout << "CDF_evaluate-Wert: " << basis.evaluate(0,RIndex(0,1,0),0.4) << endl;
  
  QuarkletFrame<2, 2> frame3;
  cout << "QFrame-Wert: " << frame3.evaluate(0,RQIndex(0,0,1,0),0.4) << endl;
  cout << "Splinewert: " << EvaluateCardinalBSpline_td<2>(0,0,0.0) << endl;
  
    

  InfiniteVector<double, RQIndex> coeff;
  coeff[RQIndex(0,2,0,0)] = 1.0;
  coeff[RQIndex(1,2,0,1)] = -3.14;
  coeff[RQIndex(2,2,0,3)] = 2.0;
  cout << "  * a small but nontrivial coefficient set:" << endl << coeff;
  cout << "  * result of DECOMPOSE:" << endl;
  InfiniteVector<double,RQIndex> wcoeff;
  frame3.decompose(coeff,0,wcoeff);
  cout << wcoeff;
  cout << "  * RECONSTRUCT that:" << endl;
  InfiniteVector<double,RQIndex> rcoeff;
  frame3.reconstruct(wcoeff,2,rcoeff);
  cout << rcoeff;

  cout << "  * again a small but nontrivial coefficient set:" << endl << coeff;
  cout << "  * result of DECOMPOSEt:" << endl;
  wcoeff.clear();
  frame3.decompose_t(coeff,0,wcoeff);
  cout << wcoeff;
  cout << "  * RECONSTRUCTt that:" << endl;
  rcoeff.clear();
  frame3.reconstruct_t(wcoeff,2,rcoeff);
  cout << rcoeff;


  
  
  
  InfiniteVector<double,RIndex> evalcoeffs;
  evalcoeffs[RIndex(0,0,0)] = 1.0;
  evalcoeffs[RIndex(1,0,2)] = 1.0;
  evalcoeffs[RIndex(0,1,0)] = 1.0;
  basis2.evaluate(0, evalcoeffs, true, -2, 2, 2).matlab_output(cout);

  cout << "Testing wavelet bases over R^d..." << endl;

  RQIndex lambda;
  cout << "- testing wavelet index class RIndex:" << endl;
  cout << "  * default index: " << lambda << endl;
  lambda = RQIndex(3,1,4,0);
  cout << "  * some index: " << lambda << endl;
  RQIndex lambda2(lambda);
  cout << "  * testing copy constructor: " << lambda2 << endl;
  cout << "  * check equality: " << (lambda==lambda2) << endl;
  cout << "  * check non-equality: " << (lambda!=lambda2) << endl;
  ++lambda;
  cout << "  * check preincrement operator ++: " << lambda << endl;
  cout << "  * check lex. order operator <: "
       << (lambda2 < lambda) << ", " << (lambda < lambda2) << endl;
  cout << "  * making a set out of the two wavelet indices:" << endl;
  set<RQIndex> M;
  M.insert(lambda);
  M.insert(lambda2);
  for (set<RQIndex>::const_iterator it(M.begin()); it != M.end(); ++it)
    cout << *it << endl;
  cout << "  * making a map out of the two wavelet indices and some coeffs:" << endl;
  map<RQIndex,double> c;
  c.insert(std::make_pair<RQIndex,double>(lambda,42.0));
  c.insert(std::make_pair<RQIndex,double>(lambda2,23.0));
  for (map<RQIndex,double>::const_iterator it(c.begin()); it != c.end(); ++it)
  cout << it->first << ", " << it->second << endl;
  
  

  QuarkletFrame<2, 2> frame2;
  
  const RQIndex lambda3(0,0,1,0);
  int k1, k2;
  int j = lambda3.j()+lambda3.e();
  frame2.support(lambda3,k1,k2);
  cout << "Träger: " << "[" <<pow(2,-j)*k1 << 
  "," << pow(2,-j)*k2 << "]" << endl;

 

}
