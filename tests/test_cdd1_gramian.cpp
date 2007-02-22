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
#include <interval/spline_basis.h>
#include <interval/spline_expansion.h>

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
  N=4: f(x)=x*(1-x)
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
    case 4:
      return p[0]*(1-p[0]);
      break;
    case 5:
      return p[0]*p[0]*(1-p[0]);
      break;
    default:
      break;
    }
    return 0.0;
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

  // for DSBasis+PBasis+SplineBasis
  const int d  = 3;
  const int dT = 3;
//   typedef DSBasis<d,dT> Basis;
//   typedef PBasis<d,dT> Basis;
  typedef SplineBasis<d,dT,P_construction> Basis;
  typedef Basis::Index Index;

#if 0
  // for DSBasis+PBasis
  Basis basis(0, 0); // no b.c.'s
//   Basis basis(1, 0); // complementary b.c. at x=0
//   Basis basis(0, 1); // complementary b.c. at x=1
//   Basis basis(1, 1); // complementary b.c.'s
#else
  // for SplineBasis
  Basis basis("",1,1,0,0); // PBasis, complementary b.c.'s
#endif

  TestFunction<5> f;

  const int jmax = 12;
  basis.set_jmax(12);

  IntervalGramian<Basis> problem(basis, InfiniteVector<double,Index>());

  InfiniteVector<double,Index> fcoeffs;
  expand(&f, basis, true, jmax, fcoeffs); // expand in the dual (!) basis
//   cout << "fcoeffs" << endl << fcoeffs;

  problem.set_rhs(fcoeffs);
  CachedProblem<IntervalGramian<Basis> > cproblem(&problem);

  // initialization with some precomputed DSBasis eigenvalue bounds:
//   CachedProblem<IntervalGramian<Basis> > cproblem(&problem, 1.4986, 47.7824); // d=2, dT=2 (no precond.)
//   CachedProblem<IntervalGramian<Basis> > cproblem(&problem, 1.28638, 705.413); // d=3, dT=3 (no precond.)

  // initialization with some precomputed PBasis eigenvalue bounds:
//   CachedProblem<IntervalGramian<Basis> > cproblem(&problem, 0.928217, 24.7998); // d=2, dT=2 (no precond.)
//   CachedProblem<IntervalGramian<Basis> > cproblem(&problem, 0.978324, 19.8057); // d=2, dT=4 (no precond.)

//   double normA = problem.norm_A();
//   double normAinv = problem.norm_Ainv();

//   cout << "* estimate for normA: " << normA << endl;
//   cout << "* estimate for normAinv: " << normAinv << endl;

  InfiniteVector<double, Index> u_epsilon;
//   CDD1_SOLVE(problem, 1e-4, u_epsilon, jmax);
//   CDD1_SOLVE(cproblem, 1e-4, u_epsilon, jmax);
//   CDD1_SOLVE(cproblem, 1e-4, u_epsilon, 12);
//   CDD1_SOLVE(problem, 1e-7, u_epsilon, jmax);
  CDD1_SOLVE(cproblem, 1e-4, u_epsilon, jmax);
//   CDD1_SOLVE(cproblem, 1e-10, u_epsilon, 20);
//   CDD1_SOLVE(cproblem, 1e-4, u_epsilon, 20);
//   CDD1_SOLVE(cproblem, 1e-2, u_epsilon, 10);
//   CDD1_SOLVE(cproblem, 1e-4, u_epsilon, 10);
//   CDD1_SOLVE(cproblem, 1e-4, u_epsilon);

  
  // insert coefficients into a dense vector
  Vector<double> wcoeffs(problem.basis().Deltasize(jmax+1));
  for (InfiniteVector<double,Index>::const_iterator it(u_epsilon.begin()),
	 itend(u_epsilon.end()); it != itend; ++it) {
    // determine number of the wavelet
    typedef Vector<double>::size_type size_type;
    size_type number = 0;
    if (it.index().e() == 0) {
      number = it.index().k()-problem.basis().DeltaLmin();
    } else {
      number = problem.basis().Deltasize(it.index().j())+it.index().k()-problem.basis().Nablamin();
    }
    wcoeffs[number] = *it;
  }
  
  // switch to generator representation
  Vector<double> gcoeffs(wcoeffs.size(), false);
  if (jmax+1 == problem.basis().j0())
    gcoeffs = wcoeffs;
  else
    problem.basis().apply_Tj(jmax, wcoeffs, gcoeffs);
  
  // evaluate Galerkin solution
  const unsigned int N = 100;
  const double h = 1./N;
  Vector<double> ulambda_values(N+1);
  for (unsigned int i = 0; i <= N; i++) {
    const double x = i*h;
    SchoenbergIntervalBSpline_td<d> sbs(jmax+1,0);
    for (unsigned int k = 0; k < gcoeffs.size(); k++) {
      sbs.set_k(problem.basis().DeltaLmin()+k);
      ulambda_values[i] += gcoeffs[k] * sbs.value(Point<1>(x));
    }
  }

  // evaluate exact solution
  Vector<double> uexact_values(N+1);
  for (unsigned int i = 0; i <= N; i++) {
    const double x = i*h;
    uexact_values[i] = f.value(Point<1>(x));
  }
  
  const double Linfty_error = linfty_norm(ulambda_values-uexact_values);
  cout << "  L_infinity error on a subgrid: " << Linfty_error << endl;


  return 0;
}
