#include <cmath>
#include <iostream>
#include <fstream>
#include <set>

#include <algebra/vector.h>
#include <algebra/infinite_vector.h>
#include <utils/function.h>
#include <utils/fixed_array1d.h>
#include <numerics/bvp.h>
#include <numerics/corner_singularity.h>
#include <geometry/sampled_mapping.h>

#include <interval/ds_basis.h>
#include <interval/p_basis.h>

// #define _WAVELETTL_LDOMAINBASIS_VERBOSITY 1
#include <Ldomain/ldomain_basis.h>
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
  cout << "Testing wavelet-Galerkin solution of an elliptic equation on the L-shaped domain ..." << endl;

  const int d  = 2;
  const int dT = 2;
//   typedef DSBasis<d,dT> Basis1D;
  typedef PBasis<d,dT> Basis1D;
  typedef LDomainBasis<Basis1D> LBasis;
  typedef LBasis::Index Index;

  CornerSingularity    u_sing(Point<2>(0,0), 0.5, 1.5);
  CornerSingularityRHS f_sing(Point<2>(0,0), 0.5, 1.5);

#if 0
  PoissonBVP<2> poisson(&f_sing);
#else
  myRHS rhs;
  PoissonBVP<2> poisson(&rhs);
#endif

  LDomainEquation<Basis1D> eq(&poisson);

  InfiniteVector<double, Index> coeffs;
  eq.RHS(1e-8, coeffs);

  set<Index> Lambda;
  for (Index lambda = eq.basis().first_generator(eq.basis().j0());; ++lambda) {
    Lambda.insert(lambda);
    if (lambda == eq.basis().last_wavelet(eq.basis().j0()+1)) break;
  }

//   cout << "- set up stiffness matrix with respect to the index set Lambda=" << endl;
//   for (set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it)
//     cout << *it << endl;

  cout << "- set up (preconditioned) stiffness matrix..." << endl;
  clock_t tstart, tend;
  double time;
  tstart = clock();

  SparseMatrix<double> A;
  setup_stiffness_matrix(eq, Lambda, A);
//   std::ofstream ofs("stiff_out.m");
//   print_matrix(A,ofs);
//   ofs.close();

  tend = clock();
  time = (double)(tend-tstart)/CLOCKS_PER_SEC;
  cout << "  ... done, time needed: " << time << " seconds" << endl;
//   cout << "- (preconditioned) stiffness matrix A=" << endl << A << endl;

//   cout << "- set up right-hand side..." << endl;
//   tstart = clock();
//   Vector<double> b;
//   setup_righthand_side(eq, Lambda, b);
//   tend = clock();
//   time = (double)(tend-tstart)/CLOCKS_PER_SEC;
//   cout << "  ... done, time needed: " << time << " seconds" << endl;
// //   cout << "- right hand side: " << b << endl;

//   Vector<double> x(Lambda.size()), err(Lambda.size()); x = 0;
//   unsigned int iterations;
//   CG(A, b, x, 1e-8, 100, iterations);
  
// //   cout << "- solution coefficients: " << x;
//   cout << " with residual (infinity) norm ";
//   A.apply(x, err);
//   err -= b;
//   cout << linfty_norm(err) << endl;


  return 0;
}
