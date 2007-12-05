#include <iostream>
#include <Rd/cdf_mask.h>
#include <Rd/dm_mask.h>
#include <Rd/refinable.h>

using namespace std;
using namespace WaveletTL;

int main()
{
  cout << "Testing the Dahmen/Micchelli mask ..." << endl;

  typedef CDFMask_primal<1> Haar;

  cout << "- integral of two shifted Haar functions is the hat function:" << endl;
  MultivariateRefinableFunction<DMMask1<Haar, Haar>, 1> dm;
  cout << dm << endl;
  
  cout << "- evaluate this at the integers:" << endl;
  cout << dm.evaluate();

  cout << "- integral of two Hat functions against each other:" << endl;
  MultivariateRefinableFunction<DMMask1<CDFMask_primal<2>, CDFMask_primal<2> >, 1> dm2;
  cout << dm2 << endl;
  
  cout << "- evaluate this at the integers:" << endl;
  cout << dm2.evaluate();

  cout << "- DM multivariate mask for integrals of three Haar fcts.:" << endl;
  MultivariateRefinableFunction<DMMask2<Haar, Haar, Haar>, 2> dm_multi;
  cout << dm_multi << endl;
  
  cout << "- evaluate this at the integers:" << endl;
  cout << dm_multi.evaluate();

  cout << "- DM multivariate mask for integrals of a Haar fct. and two translated hat fcts.:" << endl;
  MultivariateRefinableFunction<DMMask2<Haar, CDFMask_primal<2>, CDFMask_primal<2> >, 2> dm_multi2;
  cout << dm_multi2 << endl;
 
  cout << "- evaluate this at the integers:" << endl;
  cout << dm_multi2.evaluate();

  cout << "- DM multivariate mask for integrals of a Haar fct. and two translated (2,2) CDF fcts.:" << endl;
  MultivariateRefinableFunction<DMMask2<Haar, CDFMask_primal<2>, CDFMask_dual<2,2> >, 2> dm_multi3;
  cout << dm_multi3 << endl;
 
  cout << "- evaluate this at the integers:" << endl;
  cout << dm_multi3.evaluate();

  cout << "- DM multivariate mask for integrals of a Haar fct. and two translated (3,5) CDF fcts.:" << endl;
  MultivariateRefinableFunction<DMMask2<Haar, CDFMask_primal<3>, CDFMask_dual<3,5> >, 2> dm_multi4;
  cout << dm_multi4 << endl;
 
  cout << "- evaluate this at the integers:" << endl;
  cout << dm_multi4.evaluate();

  return 0;
}
