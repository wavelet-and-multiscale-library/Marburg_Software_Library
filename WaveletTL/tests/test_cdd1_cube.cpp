#include <iostream>
#include <map>
#include <time.h>

#include <algebra/symmetric_matrix.h>
#include <algebra/sparse_matrix.h>
#include <numerics/bvp.h>
#include <numerics/iteratsolv.h>

#include <interval/i_index.h>
#include <interval/ds_basis.h>
#include <interval/p_basis.h>
#include <interval/jl_basis.h>
#include <interval/jl_support.h>
#include <interval/jl_evaluate.h>

#include <cube/cube_indexplot.h>

#define _WAVELETTL_GALERKINUTILS_VERBOSITY 0

#include <cube/cube_basis.h>
#include <galerkin/cached_problem.h>
#include <galerkin/cube_equation.h>

#define _WAVELETTL_CDD1_VERBOSITY 1
#include <adaptive/cdd1.h>

using namespace std;
using namespace WaveletTL;
using namespace MathTL;

/*
  Some test problems for the Poisson equation on the cube with homogeneous Dirichlet b.c.'s:
  1: u(x,y) = x(1-x)y(1-y), -Delta u(x,y) = 2(x(1-x)+y(1-y))
  2: u(x,y) = exp(-50*((x-0.5)^2+(y-0.5)^2)), -Delta u(x,y)= (200-(100x-50)^2-(100y-50)^2)*u(x,y)
  3: u(x,y) = x(1-x)^2y^2(1-y), -Delta u(x,y)= 4*(1-x)*y^2*(1-y)-2*x*y^2*(1-y)-2*x*(1-x)^2*(1-y)+4*x*(1-x)^2*y
*/
template <unsigned int N>
class TestRHS
  : public Function<2,double>
{
public:
  virtual ~TestRHS() {};
  double value(const Point<2>& p, const unsigned int component = 0) const {
    switch(N) {
    case 2:
      return
	(200.-(100.*p[0]-50.)*(100.*p[0]-50.)-(100.*p[1]-50.)*(100.*p[1]-50.))
	* exp(-50.*((p[0]-0.5)*(p[0]-0.5)+(p[1]-0.5)*(p[1]-0.5)));
      break;
    case 3:
      return
	4*(1-p[0])*p[1]*p[1]*(1-p[1])
	- 2*p[0]*p[1]*p[1]*(1-p[1])
	- 2*p[0]*(1-p[0])*(1-p[0])*(1-p[1])
	+ 4*p[0]*(1-p[0])*(1-p[0])*p[1];
      break;
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

#if 1
  const int d  = 3;
  const int dT = 3;
  const int jmax = 3;
//   typedef DSBasis<d,dT> Basis1D;
  typedef PBasis<d,dT> Basis1D;
#else
  typedef JLBasis Basis1D; // does not work at the moment
#endif
  typedef CubeBasis<Basis1D,2> Basis;
  typedef Basis::Index Index;

  TestRHS<1> rhs;
  PoissonBVP<2> poisson(&rhs);

  FixedArray1D<bool,4> bc;
  bc[0] = bc[1] = bc[2] = bc[3] = true;
  typedef CubeEquation<Basis1D,2,Basis> Problem;
  Problem problem(&poisson, bc, jmax);
//   CachedProblem<Problem> cproblem(&problem);

  // initialization with some precomputed DSBasis eigenvalue bounds:
//   CachedProblem<Problem> cproblem(&problem, 19.97  ,    6.86044); // d=2, dT=2
//   CachedProblem<Problem> cproblem(&problem, 29.8173,   25.6677 ); // d=2, dT=4
//   CachedProblem<Problem> cproblem(&problem, 8.51622, 10000); //6311.51   ); // d=3, dT=3 

  // initialization with some precomputed PBasis eigenvalue bounds:
//   CachedProblem<Problem> cproblem(&problem, 10.6941 ,   3.4127); // d=2, dT=2 (2^j-precond.)
//   CachedProblem<Problem> cproblem(&problem,  2.77329,  11.1314); // d=2, dT=2 (diag. precond.)
//   CachedProblem<Problem> cproblem(&problem, 37.9188 ,  14.6577); // d=2, dT=4 (2^j-precond.)
//   CachedProblem<Problem> cproblem(&problem,  4.45301,  213.333); // d=2, dT=4 (diag. precond.)
//   CachedProblem<Problem> cproblem(&problem,  7.15276, 9044.08 ); // d=2, dT=6 (diag. precond.)
   CachedProblem<Problem> cproblem(&problem,  2.35701,  80.8879); // d=3, dT=3 (2^j-precond.)
//   CachedProblem<Problem> cproblem(&problem,  4.91237,  23.5086); // d=3, dT=3 (diag. precond.)
//   CachedProblem<Problem> cproblem(&problem,  2.4999 ,  67.5863); // d=3, dT=5 (2^j-precond., not exact)
//   CachedProblem<Problem> cproblem(&problem,  5.49044, 124.85  ); // d=3, dT=5 (diag. precond.)

  // initialization with some precomputed JLBasis eigenvalue bounds:
//   CachedProblem<Problem> cproblem(&problem, 0.3, 400); // (2^j-precond.)
//   CachedProblem<Problem> cproblem(&problem, 0.328494, 385.552); // (2^j-precond.)
//  CachedProblem<Problem> cproblem(&problem, 2.65769, 6.16906); // (diag. precond.)

//   double normA = problem.norm_A();
//   double normAinv = problem.norm_Ainv();

//   cout << "* estimate for normA: " << normA << endl;
//   cout << "* estimate for normAinv: " << normAinv << endl;

   InfiniteVector<double, Index> u_epsilon, f, w, r;


//   Index ind;


//   for (int i = 0; i < 100; i++) {
//     APPLY(cproblem, u_epsilon, 0.01, w, 3, CDD1);
//     cproblem.RHS(0.01, f);
//     r = f-w;
//     //r.COARSE(0.0001,w);
//     r = w;
//     r *= 0.2;
//     u_epsilon = u_epsilon + r;
//     cout << l2_norm(r) << endl;
//   }




   CDD1_SOLVE(cproblem, 1e-2, u_epsilon, jmax);
//   CDD1_SOLVE(cproblem, 1e-4, u_epsilon, jmax, CDD1);
//   CDD1_SOLVE(cproblem, 1e-4, u_epsilon, 7);
//   CDD1_SOLVE(cproblem, 1e-4, u_epsilon, 10);
//   CDD1_SOLVE(cproblem, 1e-4, u_epsilon, 6, CDD1);
//   CDD1_SOLVE(cproblem, 1e-4, u_epsilon);
   
   std::ofstream plotstream;
   plotstream.open("cube_indexplot_cdd1.m");
   plot_indices(&cproblem.basis(), u_epsilon, jmax, plotstream, "jet", true, true);

   plotstream.close();

   u_epsilon.scale(&cproblem, -1);
   cout << "saving computed solution to u_adaptCDD1_cube.m" << endl;
   SampledMapping<2> s(evaluate(cproblem.basis(), u_epsilon, true, d+2+2+1));
   std::ofstream u_stream("u_adaptCDD1_cube.m");
   s.matlab_output(u_stream);
   u_stream.close();
  
  return 0;
}
