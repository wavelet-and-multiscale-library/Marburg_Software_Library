#include <iostream>
#include <map>
#include <time.h>

#include <algebra/symmetric_matrix.h>
#include <algebra/sparse_matrix.h>
#include <numerics/bvp.h>
#include <numerics/iteratsolv.h>
#include <numerics/corner_singularity.h>

#include <interval/i_index.h>
#include <interval/ds_basis.h>
#include <interval/p_basis.h>

#define _WAVELETTL_LDOMAINBASIS_VERBOSITY 1
#define _WAVELETTL_GALERKINUTILS_VERBOSITY 1

#include <Ldomain/ldomain_basis.h>
#include <galerkin/cached_problem.h>
#include <galerkin/ldomain_equation.h>

using namespace std;
using namespace MathTL;
using namespace WaveletTL;

/*
  A test problem for the Poisson equation on the L--shaped domain with homogeneous Dirichlet b.c.'s:
    -Delta u(x,y) = 2*pi^2*sin(pi*x)*sin(pi*y)
  with exact solution
    u(x,y) = sin(pi*x)*sin(pi*y)
*/
class myRHS
  : public Function<2,double>
{
public:
  virtual ~myRHS() {};
  double value(const Point<2>& p, const unsigned int component = 0) const {
    return 2*M_PI*M_PI*sin(M_PI*p[0])*sin(M_PI*p[1]);
  }
  void vector_value(const Point<2>& p, Vector<double>& values) const {
    values[0] = value(p);
  }
};

class mySolution
  : public Function<2,double>
{
public:
  virtual ~mySolution() {};
  double value(const Point<2>& p, const unsigned int component = 0) const {
    return sin(M_PI*p[0])*sin(M_PI*p[1]);
  }
  void vector_value(const Point<2>& p, Vector<double>& values) const {
    values[0] = value(p);
  }
};

int main()
{
  cout << "Testing adaptive wavelet-Galerkin solution of a Poisson problem on the L-shaped domain with CDD1_SOLVE ..." << endl;

  const int d  = 3;
  const int dT = 3;
  typedef DSBasis<d,dT> Basis1D;
//   typedef PBasis<d,dT> Basis1D;

  typedef LDomainBasis<Basis1D> LBasis;
  typedef LBasis::Index Index;

//   CornerSingularity    u_sing(Point<2>(0,0), 0.5, 1.5);
//   CornerSingularityRHS f_sing(Point<2>(0,0), 0.5, 1.5);

  myRHS rhs;
  PoissonBVP<2> poisson(&rhs);

  typedef LDomainEquation<Basis1D> Problem;
  Problem problem(&poisson);
  //   CachedProblem<Problem> cproblem(&problem);

  // initialization with some precomputed DSBasis eigenvalue bounds: ### update the values!!!
  //   CachedProblem<Problem> cproblem(&problem, 19.97  ,    6.86044); // d=2, dT=2
  //   CachedProblem<Problem> cproblem(&problem, 29.8173,   25.6677 ); // d=2, dT=4
  //   CachedProblem<Problem> cproblem(&problem, 8.51622, 10000); //6311.51   ); // d=3, dT=3 

  // initialization with some precomputed PBasis eigenvalue bounds: ### update the values!!!
  //   CachedProblem<Problem> cproblem(&problem, 10.6941 ,   3.4127); // d=2, dT=2 (2^j-precond.)
  //   CachedProblem<Problem> cproblem(&problem,  2.77329,  11.1314); // d=2, dT=2 (diag. precond.)
  //   CachedProblem<Problem> cproblem(&problem, 37.9188 ,  14.6577); // d=2, dT=4 (2^j-precond.)
  //   CachedProblem<Problem> cproblem(&problem,  4.45301,  213.333); // d=2, dT=4 (diag. precond.)
  //   CachedProblem<Problem> cproblem(&problem,  7.15276, 9044.08 ); // d=2, dT=6 (diag. precond.)
  //   CachedProblem<Problem> cproblem(&problem,  2.35701,  80.8879); // d=3, dT=3 (2^j-precond.)
//   CachedProblem<Problem> cproblem(&problem,  4.91237,  23.5086); // d=3, dT=3 (diag. precond.)
  //   CachedProblem<Problem> cproblem(&problem,  2.4999 ,  67.5863); // d=3, dT=5 (2^j-precond., not exact)
  //   CachedProblem<Problem> cproblem(&problem,  5.49044, 124.85  ); // d=3, dT=5 (diag. precond.)

  double normA = problem.norm_A();
  double normAinv = problem.norm_Ainv();

  cout << "* estimate for normA: " << normA << endl;
  cout << "* estimate for normAinv: " << normAinv << endl;

//   InfiniteVector<double, Index> u_epsilon;

//   CDD1_SOLVE(cproblem, 1e-2, u_epsilon, 5);
// //   CDD1_SOLVE(cproblem, 1e-4, u_epsilon, 7);
// //   CDD1_SOLVE(cproblem, 1e-4, u_epsilon, 10);
// //   CDD1_SOLVE(cproblem, 1e-4, u_epsilon, 6, CDD1);
// //   CDD1_SOLVE(cproblem, 1e-4, u_epsilon);
  
  return 0;
}
