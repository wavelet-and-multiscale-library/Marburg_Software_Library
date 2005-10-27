#include <iostream>
#include <map>
#include <time.h>

#include <algebra/symmetric_matrix.h>
#include <algebra/sparse_matrix.h>
#include <numerics/bvp.h>
#include <numerics/iteratsolv.h>

#include <interval/i_index.h>
#include <interval/ds_basis.h>
#include <cube/cube_basis.h>
#include <galerkin/cached_problem.h>
#include <galerkin/cube_equation.h>

#include <adaptive/cdd1.h>

using namespace std;
using namespace WaveletTL;
using namespace MathTL;

/*
  Some test problems for the Poisson equation on the cube with homogeneous Dirichlet b.c.'s:
  1: u(x,y) = x(1-x)y(1-y), -Delta u(x,y) = 2(x(1-x)+y(1-y))
*/
template <unsigned int N>
class TestRHS
  : public Function<2,double>
{
public:
  virtual ~TestRHS() {};
  double value(const Point<2>& p, const unsigned int component = 0) const {
    switch(N) {
    case 1:
    default:
      return 2*(p[0]*(1-p[0])+p[1]*(1-p[1]));
      break;
    }
  }
  void vector_value(const Point<2>& p, Vector<double>& values) const {
    values[0] = value(p);
  }
};

int main()
{
  cout << "Testing adaptive wavelet-Galerkin solution of a Poisson problem on the cube with CDD1_SOLVE ..." << endl;

  const int d  = 2;
  const int dT = 2;
  typedef DSBasis<d,dT> Basis1D;
  typedef CubeBasis<Basis1D,2> Basis;
  typedef Basis::Index Index;

  TestRHS<1> rhs;
  PoissonBVP<2> poisson(&rhs);

  FixedArray1D<bool,4> bc;
  bc[0] = bc[1] = bc[2] = bc[3] = true;
  typedef CubeEquation<Basis1D,2,Basis> Problem;
  Problem problem(&poisson, bc);
  CachedProblem<Problem> cproblem(problem);

  InfiniteVector<double, Index> u_epsilon;
  CDD1_SOLVE(cproblem, 1e-4, u_epsilon, 8);
//   CDD1_SOLVE(cproblem, 1e-4, u_epsilon, 10);
//   CDD1_SOLVE(cproblem, 1e-4, u_epsilon, 6, CDD1);
//   CDD1_SOLVE(cproblem, 1e-4, u_epsilon);
  
  return 0;
}
