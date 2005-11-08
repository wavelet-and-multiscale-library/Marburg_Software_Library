#include <iostream>
#include <fstream>
#include <set>

#include <algebra/vector.h>
#include <algebra/infinite_vector.h>
#include <utils/function.h>
#include <utils/fixed_array1d.h>
#include <numerics/bvp.h>
#include <geometry/sampled_mapping.h>

#include <interval/ds_basis.h>
#include <cube/cube_basis.h>
#include <galerkin/cube_equation.h>

using namespace std;
using namespace MathTL;
using namespace WaveletTL;

/*
  A test problem for the Poisson equation on the cube with homogeneous Dirichlet b.c.'s:
    -Delta u(x,y) = 2(x(1-x)+y(1-y))
  with exact solution
    u(x,y) = x(1-x)y(1-y)
*/
class myRHS
  : public Function<2,double>
{
public:
  virtual ~myRHS() {};
  double value(const Point<2>& p, const unsigned int component = 0) const {
    return 2*(p[0]*(1-p[0])+p[1]*(1-p[1]));
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
    return p[0]*(1-p[0])*p[1]*(1-p[1]);
  }
  void vector_value(const Point<2>& p, Vector<double>& values) const {
    values[0] = value(p);
  }
};

int main()
{
  cout << "Testing wavelet-Galerkin solution of an elliptic equation on the cube ..." << endl;

  ConstantFunction<2> constant_rhs(Vector<double>(1, "1.0"));
  PoissonBVP<2> poisson(&constant_rhs);

  const int d  = 2;
  const int dT = 2; // be sure to use a continuous dual here, otherwise the RHS test will fail
  typedef DSBasis<d,dT> Basis1D;
  typedef CubeBasis<Basis1D,2> CBasis;
  typedef CBasis::Index Index;

  FixedArray1D<bool,4> bc;
  bc[0] = bc[1] = bc[2] = bc[3] = true;

  CubeEquation<Basis1D,2,CBasis> eq(&poisson, bc);

  InfiniteVector<double, Index> coeffs;

#if 0
  coeffs[first_generator<Basis1D,2,CBasis>(&eq.basis(), eq.basis().j0())] = 1.0;
  coeffs[last_generator<Basis1D,2,CBasis>(&eq.basis(), eq.basis().j0())] = 2.0;
  coeffs[first_wavelet<Basis1D,2,CBasis>(&eq.basis(), eq.basis().j0())] = 3.0;
  coeffs[last_wavelet<Basis1D,2,CBasis>(&eq.basis(), eq.basis().j0())] = 4.0;
  coeffs[first_wavelet<Basis1D,2,CBasis>(&eq.basis(), eq.basis().j0()+1)] = 5.0;
  cout << "- a coefficient set:" << endl
       << coeffs << endl;
  eq.rescale(coeffs, -1);
  cout << "- after rescaling with D^{-1}:" << endl
       << coeffs << endl;
#endif
  
  eq.RHS(1e-8, coeffs);

#if 0
//   cout << "- approximate coefficient set of the right-hand side:" << endl
//        << coeffs << endl;
  cout << "- check expansion of the right-hand side in the dual basis:" << endl;
  eq.rescale(coeffs, 1);
  SampledMapping<2> S(evaluate<Basis1D,2>(eq.basis(), coeffs, false, 6));
//   std::ofstream rhs_stream("constant_rhs.m");
//   S.matlab_output(rhs_stream);
//   rhs_stream.close();
//   cout << "  ... done, see file 'constant_rhs.m'" << endl;
  S.add(-1.0, SampledMapping<2>(S, constant_rhs));
  cout << "  ... done, pointwise error: " << row_sum_norm(S.values()) << endl;
  eq.rescale(coeffs, -1);
#endif  

  set<Index> Lambda;
  for (Index lambda = first_generator<Basis1D,2,CBasis>(&eq.basis(), eq.basis().j0());; ++lambda) {
    Lambda.insert(lambda);
    if (lambda == last_wavelet<Basis1D,2,CBasis>(&eq.basis(), eq.basis().j0()+1)) break;
  }

//   cout << "- set up stiffness matrix with respect to the index set Lambda=" << endl;
//   for (set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it)
//     cout << *it << endl;

  // choose another rhs
  myRHS rhs;
  poisson.set_f(&rhs);
  eq.set_bvp(&poisson);
  eq.RHS(1e-4, coeffs);
  
#if 1
  cout << "- set up (preconditioned) stiffness matrix..." << endl;
  clock_t tstart, tend;
  double time;
  tstart = clock();
  SparseMatrix<double> A;
// <<<<<<< test_cube_equation.cpp
//   setup_stiffness_matrix(eq_test, Lambda, A);
//   std::ofstream ofs("stiff_out.m");
//   print_matrix(A,ofs);
//   ofs.close();
// =======
//   setup_stiffness_matrix(eq, Lambda, A);
// >>>>>>> 1.11
  setup_stiffness_matrix(eq, Lambda, A);
  tend = clock();
  time = (double)(tend-tstart)/CLOCKS_PER_SEC;
  cout << "  ... done, time needed: " << time << " seconds" << endl;
//   cout << "- (preconditioned) stiffness matrix A=" << endl << A << endl;

  cout << "- set up right-hand side..." << endl;
  tstart = clock();
  Vector<double> b;
  setup_righthand_side(eq, Lambda, b);
  tend = clock();
  time = (double)(tend-tstart)/CLOCKS_PER_SEC;
  cout << "  ... done, time needed: " << time << " seconds" << endl;
//   cout << "- right hand side: " << b << endl;

  Vector<double> x(Lambda.size()), err(Lambda.size()); x = 0;
  unsigned int iterations;
  CG(A, b, x, 1e-8, 100, iterations);
  
//   cout << "- solution coefficients: " << x;
  cout << " with residual (infinity) norm ";
  A.apply(x, err);
  err -= b;
  cout << linfty_norm(err) << endl;
  
  cout << "- point values of the solution:" << endl;
  InfiniteVector<double,Index> u;
  unsigned int i = 0;
  for (set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it, ++i)
    u.set_coefficient(*it, x[i]);
  
  eq.rescale(u, -1);
  SampledMapping<2> s(evaluate(eq.basis(), u, true, 6));
//   std::ofstream u_Lambda_stream("u_lambda.m");
//   s.matlab_output(u_Lambda_stream);
//   u_Lambda_stream.close();
//   cout << "  ... done, see file 'u_lambda.m'" << endl;
  mySolution u_Lambda;
  s.add(-1.0, SampledMapping<2>(s, u_Lambda));
  cout << "  ... done, pointwise error: " << row_sum_norm(s.values()) << endl;
#endif

  return 0;
}
