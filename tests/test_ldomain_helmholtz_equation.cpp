#include <iostream>

#include <utils/function.h>
#include <Ldomain/ldomain_jl_basis.h>
#include <galerkin/ldomain_jl_helm_eq.h>

#include "ldomain_solutions.h"

using namespace std;
using namespace MathTL;
using namespace WaveletTL;

int main()
{
  cout << "Testing LDomainJLHelmholtzEquation ..." << endl;

  typedef LDomainJLBasis Basis;
  typedef Basis::Index Index;

  Basis basis;
  
  const int solution = 1;
  Function<2> *uexact = 0, *f = 0;
  switch(solution) {
  case 1:
    uexact = new EigenSolution();
    f = new EigenRHS();
    break;
  default:
    break;
  }

  if (uexact) delete uexact;
  if (f) delete f;

  return 0;
}
