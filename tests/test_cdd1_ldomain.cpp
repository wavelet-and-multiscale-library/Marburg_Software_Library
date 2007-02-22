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
#include <interval/spline_basis.h>

#define _WAVELETTL_LDOMAINBASIS_VERBOSITY 0
#define _WAVELETTL_GALERKINUTILS_VERBOSITY 0

#include <Ldomain/ldomain_basis.h>
#include <Ldomain/ldomain_expansion.h>
#include <galerkin/cached_problem.h>
#include <galerkin/ldomain_equation.h>
#include <galerkin/ldomain_helmholtz_equation.h>

#define _WAVELETTL_CDD1_VERBOSITY 1
#include <adaptive/cdd1.h>

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

class mySum
  : public Function<2,double>
{
public:
  virtual ~mySum() {};
  double value(const Point<2>& p, const unsigned int component = 0) const {
    return (1+2*M_PI*M_PI)*sin(M_PI*p[0])*sin(M_PI*p[1]);
  }
  void vector_value(const Point<2>& p, Vector<double>& values) const {
    values[0] = value(p);
  }
};

int main()
{
  cout << "Testing adaptive wavelet-Galerkin solution of a Poisson problem on the L-shaped domain with CDD1_SOLVE ..." << endl;

  const int d  = 2;
  const int dT = 2;
  
  typedef SplineBasis<d,dT,DS_construction> Basis1D;
  Basis1D basis1d("bio5-energy",0,0,0,0);
  
  typedef LDomainBasis<Basis1D> Basis;
  typedef Basis::Index Index;
  Basis basis(basis1d);


  const int jmax = 5;
  basis.set_jmax(jmax);

//   CornerSingularity    u_sing(Point<2>(0,0), 0.5, 1.5);
//   CornerSingularityRHS f_sing(Point<2>(0,0), 0.5, 1.5);

//   myRHS rhs; // alpha=0
  mySum rhs; // alpha=1

  InfiniteVector<double,Index> rhs_coeffs;
  cout << "expanding right-hand side" << endl;
  expand(&rhs, basis, true, jmax, rhs_coeffs);
  cout << "done expanding right-hand side" << endl;

  typedef LDomainHelmholtzEquation<d,dT> Problem;
  Problem problem(basis,
		  "DS_B_2_2_5_G",
		  "DS_B_2_2_5_A",
		  jmax,
		  1.0,
		  rhs_coeffs);

//   double normA = problem.norm_A();
//   double normAinv = problem.norm_Ainv();

//   cout << "* estimate for normA: " << normA << endl;
//   cout << "* estimate for normAinv: " << normAinv << endl;

  InfiniteVector<double, Index> u_epsilon;

  CDD1_SOLVE(problem, 1e-2, u_epsilon, u_epsilon, 1.0, 1.0, jmax);
//   CDD1_SOLVE(problem, 1e-2, u_epsilon, u_epsilon, 4);

//   CDD1_SOLVE(cproblem, 1e-2, u_epsilon, 5);
// //   CDD1_SOLVE(cproblem, 1e-4, u_epsilon, 7);
// //   CDD1_SOLVE(cproblem, 1e-4, u_epsilon, 10);
// //   CDD1_SOLVE(cproblem, 1e-4, u_epsilon, 6, CDD1);
// //   CDD1_SOLVE(cproblem, 1e-4, u_epsilon);
  
  return 0;
}
