#include <iostream>

#include <utils/function.h>
#include <Ldomain/ldomain_jl_basis.h>
#include <Ldomain/ldomain_jl_expansion.h>
#include <galerkin/ldomain_jl_gramian.h>

#include "ldomain_solutions.h"

using namespace std;
using namespace MathTL;
using namespace WaveletTL;

int main()
{
  cout << "Testing LDomainJLGramian ..." << endl;

  typedef LDomainJLBasis Basis;
  typedef Basis::Index Index;

  Basis basis;
  
  const int solution = 1;
  Function<2> *uexact = 0, *f = 0;
  switch(solution) {
  case 1:
    uexact = new PolySolution();
    f = new PolyRHS();
    break;
  case 2:
    uexact = new EigenSolution();
    f = new EigenRHS();
    break;
  default:
    break;
  }

#if 0
  // temporary hack, choose f=1 just to see the inner products
  delete f;
  f = new ConstantFunction<2>(Vector<double>(1, "1"));
#endif

  const int jmax = 1;
  
  typedef LDomainJLGramian Problem;
  Problem problem(basis, InfiniteVector<double,Index>());

  InfiniteVector<double,Index> fcoeffs;
  expand(f, basis, true, jmax, fcoeffs);
  problem.set_rhs(fcoeffs);

  cout << "- integrals of f against the primal wavelets:" << endl
       << fcoeffs << endl;


  if (uexact) delete uexact;
  if (f) delete f;

  return 0;
}
