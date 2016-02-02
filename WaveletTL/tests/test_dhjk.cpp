#include <iostream>
#include <Rd/dhjk_mask.h>
#include <Rd/refinable.h>

using namespace std;
using namespace WaveletTL;

int main()
{
  cout << "Testing the Dahmen/Han/Jia/Kunoth mask ..." << endl;

  cout << "- primal multigenerator:" << endl;
  RefinableFunction<DHJKMask_primal> primal;
  cout << primal << endl;

  cout << "- dual multigenerator:" << endl;
  RefinableFunction<DHJKMask_dual> dual;
  cout << dual << endl;

//  cout << "point evaluation: \\tilde{\\phi}_1(0.5) = ";
//  cout << dual.multi_evaluate(0) << endl;
//  cout << "plotting the generators" << endl;
  
  return 0;
}
