#define _WAVELETTL_GALERKINUTILS_VERBOSITY 0
#define _WAVELETTL_CDD1_VERBOSITY 1

#include <iostream>
#include <map>
#include <time.h>

#include <algebra/infinite_vector.h>

#include <interval/i_index.h>
#include <interval/ds_basis.h>
#include <interval/ds_expansion.h>
#include <interval/p_basis.h>
#include <interval/p_expansion.h>

#include <galerkin/gramian.h>
#include <galerkin/cached_problem.h>

#include <adaptive/cdd1.h>

using namespace std;
using namespace WaveletTL;

using MathTL::CG;

/*
  test function
  N=1: f(x)=1
  N=2: f(x)=exp(x)
  N=3: f(x)=chi_{[0,0.5)}(x)
*/
template <unsigned int N>
class TestFunction
  : public Function<1>
{
  public:
  inline double value(const Point<1>& p,
		      const unsigned int component = 0) const
  {
    switch (N) {
    case 1:
      return 1.0;
      break;
    case 2:
      return exp(p[0]);
      break;
    case 3:
      return (p[0] < 0.5 ? 1.0 : 0.0);
      break;
    default:
      return 0.0;
      break;
    }
  }
  
  void vector_value(const Point<1> &p,
		    Vector<double>& values) const
  {
    values.resize(1, false);
    values[0] = value(p);
  }
};

int main()
{
  cout << "Testing adaptive wavelet-Galerkin methods for the identity operator with CDD1_SOLVE ..." << endl;

#if 1
  const int d  = 2;
  const int dT = 2;
//   typedef DSBasis<d,dT> Basis;
  typedef PBasis<d,dT> Basis;
#else
  typedef JLBasis Basis;
#endif
  typedef Basis::Index Index;

  Basis basis(0, 0); // no b.c.'s
//   Basis basis(1, 0); // complementary b.c. at x=0
//   Basis basis(0, 1); // complementary b.c. at x=1
//   Basis basis(1, 1); // complementary b.c.'s

  TestFunction<2> f;

  const int j0   = basis.j0();
  const int jmax = 8;

  IntervalGramian<Basis> problem(basis, InfiniteVector<double,Index>());

  InfiniteVector<double,Index> fcoeffs;
  expand(&f, basis, true, jmax, fcoeffs); // expand in the dual (!) basis
//   cout << "fcoeffs" << endl << fcoeffs;

  problem.set_rhs(fcoeffs);
//   CachedProblem<SturmEquation<Basis> > cproblem(&problem);

  // initialization with some precomputed DSBasis eigenvalue bounds:
//   CachedProblem<IntervalGramian<Basis> > cproblem(&problem, 1.4986, 47.7824); // d=2, dT=2 (no precond.)
//   CachedProblem<IntervalGramian<Basis> > cproblem(&problem, 1.28638, 705.413); // d=3, dT=3 (no precond.)

  // initialization with some precomputed PBasis eigenvalue bounds:
  CachedProblem<IntervalGramian<Basis> > cproblem(&problem, 0.928217, 24.7998); // d=2, dT=2 (no precond.)
//   CachedProblem<IntervalGramian<Basis> > cproblem(&problem, 0.978324, 19.8057); // d=2, dT=4 (no precond.)

//   double normA = problem.norm_A();
//   double normAinv = problem.norm_Ainv();

//   cout << "* estimate for normA: " << normA << endl;
//   cout << "* estimate for normAinv: " << normAinv << endl;

  InfiniteVector<double, Index> u_epsilon;
  CDD1_SOLVE(cproblem, 1e-4, u_epsilon, 8);
//   CDD1_SOLVE(cproblem, 1e-4, u_epsilon, 12);
//   CDD1_SOLVE(cproblem, 1e-7, u_epsilon, 12);
//   CDD1_SOLVE(cproblem, 1e-5, u_epsilon, 20);
//   CDD1_SOLVE(cproblem, 1e-10, u_epsilon, 20);
//   CDD1_SOLVE(cproblem, 1e-4, u_epsilon, 20);
//   CDD1_SOLVE(cproblem, 1e-2, u_epsilon, 10);
//   CDD1_SOLVE(cproblem, 1e-4, u_epsilon, 10);
//   CDD1_SOLVE(cproblem, 1e-4, u_epsilon);
  
  return 0;
}
