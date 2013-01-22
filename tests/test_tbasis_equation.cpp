#include <iostream>
#include <fstream>
#include <set>

#include <algebra/vector.h>
#include <algebra/infinite_vector.h>
#include <utils/function.h>
#include <utils/fixed_array1d.h>
#include <numerics/bvp.h>
#include <numerics/bezier.h>
#include <geometry/sampled_mapping.h>

//#include <interval/ds_basis.h>
#include <interval/p_basis.h>
//#include <interval/jl_basis.h>
//#include <interval/jl_support.h>
//#include <interval/jl_evaluate.h>

#define _WAVELETTL_GALERKINUTILS_VERBOSITY 1

#include <cube/tbasis.h>
#include <galerkin/tbasis_equation.h>
#include <cube/tbasis_evaluate.h>


using namespace std;
using namespace MathTL;
using namespace WaveletTL;

/*
  Some test problem for the Poisson equation on the cube with homogeneous Dirichlet b.c.'s:
  1: -Delta u(x,y) = 2(x(1-x)+y(1-y)), u(x,y) = x(1-x)y(1-y)
  2: -Delta u(x,y) = 2*pi^2*sin(pi*x)*sin(pi*y), u(x,y) = sin(pi*x)*sin(pi*y)
  3: u cubic Hermite interpolant
*/
template <unsigned int N>
class myRHS
  : public Function<2,double>
{
public:
  virtual ~myRHS() {};
  double value(const Point<2>& p, const unsigned int component = 0) const {
    switch(N) {
    case 1:
      return 2*(p[0]*(1-p[0])+p[1]*(1-p[1]));
      break;
    case 2:
      return 2*M_PI*M_PI*sin(M_PI*p[0])*sin(M_PI*p[1]);
      break;
    case 3:
      return
	-((p[0]<=0.5 ? 24-96*p[0] : 96*p[0]-72)
	  * (p[1]<=0.5 ? 4*p[1]*p[1]*(3-4*p[1]) : (2-2*p[1])*(2-2*p[1])*(4*p[1]-1))
	  +(p[0]<=0.5 ? 4*p[0]*p[0]*(3-4*p[0]) : (2-2*p[0])*(2-2*p[0])*(4*p[0]-1))
	  * (p[1]<=0.5 ? 24-96*p[1] : 96*p[1]-72));
      break;
    default:
      return 1;
      break;
    }
  }
  void vector_value(const Point<2>& p, Vector<double>& values) const {
    values[0] = value(p);
  }
};

template <unsigned int N>
class mySolution
  : public Function<2,double>
{
public:
  virtual ~mySolution() {};
  double value(const Point<2>& p, const unsigned int component = 0) const {
    switch(N) {
    case 1:
      return p[0]*(1-p[0])*p[1]*(1-p[1]);
      break;
    case 2:
      return sin(M_PI*p[0])*sin(M_PI*p[1]);
      break;
    case 3:
      return
	(p[0]<=0.5 ? 4*p[0]*p[0]*(3-4*p[0]) : (2-2*p[0])*(2-2*p[0])*(4*p[0]-1))
	* (p[1]<=0.5 ? 4*p[1]*p[1]*(3-4*p[1]) : (2-2*p[1])*(2-2*p[1])*(4*p[1]-1));
      break;
    default:
      return 0;
      break;
    }
  }
  void vector_value(const Point<2>& p, Vector<double>& values) const {
    values[0] = value(p);
  }
};

int main()
{
  cout << "Testing wavelet-Galerkin solution of an elliptic equation on the cube ..." << endl;



#if 1
  const int d  = 4;
  const int dT = 4; // be sure to use a continuous dual here, otherwise the RHS test will fail
  const int dim = 2;
  const int radius = 0; // range for the Level in the indexset of the stiffness matrix

  ConstantFunction<dim> constant_rhs(Vector<double>(1, "1.0"));
  PoissonBVP<dim> poisson(&constant_rhs);
  //typedef DSBasis<d,dT> Basis1D;
  typedef PBasis<d,dT> Basis1D;
#else
  typedef JLBasis Basis1D; // does not work at the moment
#endif
  typedef TensorBasis<Basis1D,dim> Basis;
  
  typedef Basis::Index Index;

  FixedArray1D<bool,(2*dim)> bc;
  if (dim==1)
  {
      bc[0] = bc[1] = true;
  }
  else
  {
      bc[0] = bc[1] = bc[2] = bc[3] = true;
  }

  TensorEquation<Basis1D,dim,Basis> eq(&poisson, bc);
  
  InfiniteVector<double, Index> coeffs;

#if 0
  coeffs[first_generator<Basis1D,2,Basis>(&eq.basis())] = 1.0;
  coeffs[last_generator<Basis1D,2,Basis>(&eq.basis())] = 2.0;
  coeffs[first_wavelet<Basis1D,2,Basis>(&eq.basis(), eq.basis().j0())] = 3.0;
  Index lambda = last_wavelet<Basis1D,2,Basis>(&eq.basis(), eq.basis().j0());
  coeffs[lambda] = 4.0;
  ++lambda;
  coeffs[lambda] = 5.0;
  cout << "- a coefficient set:" << endl
       << coeffs << endl;
  coeffs.scale(&eq, -1);
  cout << "- after rescaling with D^{-1}:" << endl
       << coeffs << endl;
#endif

  eq.RHS(1e-8, coeffs);

#if 0
  
//   cout << "- approximate coefficient set of the right-hand side:" << endl
//        << coeffs << endl;
  cout << "- check expansion of the right-hand side in the dual basis:" << endl;
  coeffs.scale(&eq, 1);
  SampledMapping<2> S(evaluate<Basis1D,2>(eq.basis(), coeffs, false, 6));
//   std::ofstream rhs_stream("constant_rhs.m");
//   S.matlab_output(rhs_stream);
//   rhs_stream.close();
//   cout << "  ... done, see file 'constant_rhs.m'" << endl;
  S.add(-1.0, SampledMapping<2>(S, constant_rhs));
  cout << "  ... done, pointwise error: " << row_sum_norm(S.values()) << endl;
  coeffs.scale(&eq, -1);
#endif

  set<Index> Lambda;
  for (Index lambda = first_generator<Basis1D,dim,Basis>(&eq.basis()), itend = last_wavelet<Basis1D,dim,Basis>(&eq.basis(), multi_degree(eq.basis().j0())+radius);; ++lambda)
  //for (Index lambda = first_generator<Basis1D,2,Basis>(&eq.basis()), itend = last_wavelet<Basis1D,2,Basis>(&eq.basis(), eq.basis().j0());; ++lambda)
  {
      // Test fuer Richtigkeit der Reihenfolge die durch ++ gegeben wird (unzureichend!)
      //Index tempindex(lambda);
      //++tempindex;
      //cout << "lambda  = "<<lambda << " < " << tempindex << " = "<<(lambda<tempindex);
      //cout << "  #" << lambda.number() <<", #" << tempindex.number()<<endl;
    Lambda.insert(lambda);
    if (lambda == itend) break;
  }

//   cout << "- set up stiffness matrix with respect to the index set Lambda=" << endl;
//   for (set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it)
//     cout << *it << endl;

  // choose another rhs
#if 1
  const unsigned int N = 2;
  //myRHS<N> rhs;
  //poisson.set_f(&rhs);
  eq.set_bvp(&poisson);
  eq.RHS(1e-8, coeffs);


  cout << "- set up (preconditioned) stiffness matrix..." << endl;
  clock_t tstart, tend;
  double time;
  tstart = clock();
  SparseMatrix<double> A;
  setup_stiffness_matrix(eq, Lambda, A);
  A.compress(1e-15);

//   cout << "  output in file \"stiff_out.m\"..." << endl;
//   std::ofstream ofs("stiff_out.m");
//   ofs << "A=";
//   print_matrix(A, ofs);
//   ofs << ";" << endl;
//   ofs.close();
//   cout << "  ...done!" << endl;

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
#endif
  // test normA and normAinv (from cached_tproblem)
      cout << "computing norm_Ainv() ..." << endl;
      //Index tempindex(problem->basis().last_wavelet(multi_degree(problem->basis().j0())+5));
      //cout << "number of last wavelet is "<<tempindex.number()<<endl;

      set<Index> Lambda2;
      cout << "hier sollte alles gleich aussehen: "<<endl;

      cout << first_generator<Basis1D,dim,Basis>(&eq.basis()) << " =? " ;
      cout <<( eq.basis().first_generator() )<<endl;
      cout << last_wavelet<Basis1D,dim,Basis>(&eq.basis(), multi_degree(eq.basis().j0())+radius) << " =? ";
      cout << (eq.basis().last_wavelet(multi_degree(eq.basis().j0())+radius)) << endl;
/*
      for (Index lambda ( eq.basis().first_generator() ), itend(eq.basis().last_wavelet(multi_degree(eq.basis().j0())+radius));; ++lambda) {
	Lambda2.insert(lambda);
	if (lambda == itend) break;
      }
      SparseMatrix<double> A_Lambda2;
      setup_stiffness_matrix(eq, Lambda2, A_Lambda2);
 * */
      double help, normA;
      unsigned int iterations2;
      LanczosIteration(A, 1e-6, help, normA, 200, iterations2);
      double normAinv ( 1./help);
      cout << "normA = "<<normA<<endl;
      cout << "normAinv = "<<normAinv<<endl;
      cout << "kond = "<<(normA*normAinv)<<endl;
      //LanczosIteration(A_Lambda2, 1e-6, help, normA, 200, iterations2);
      //normAinv = ( 1./help);
      //cout << "normA = "<<normA<<endl;
      //cout << "normAinv = "<<normAinv<<endl;


  cout << "- point values of the solution:" << endl;
  InfiniteVector<double,Index> u;
  unsigned int i = 0;
  for (set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it, ++i)
    u.set_coefficient(*it, x[i]);

  
  u.scale(&eq, -1);
  // evaluate(eq.basis(), u, true, 6);
  SampledMapping<2> s(evaluate(eq.basis(), u, true, 6));
  std::ofstream u_Lambda_stream("u_lambda_t.m");
  s.matlab_output(u_Lambda_stream);
  u_Lambda_stream.close();
  cout << "  ... done, see file 'u_lambda_t.m'" << endl;

  mySolution<N> u_Lambda;
  s.add(-1.0, SampledMapping<2>(s, u_Lambda));
  std::ofstream u_Lambda_error_stream("u_lambda_t_error.m");
  s.matlab_output(u_Lambda_error_stream);
  u_Lambda_error_stream.close();
  cout << "  ... done, see file 'u_lambda_t_error.m', pointwise maximal error: " << maximum_norm(s.values()) << endl;

   SampledMapping<2> sexact(s, u_Lambda);
   std::ofstream u_exact_stream("u_exact_t.m");
   sexact.matlab_output(u_exact_stream);
   u_exact_stream.close();
   cout << " ... done writing the exact solution to file 'u_exact_t.m'"<<endl;

//   SampledMapping<2> srhs(s, rhs);
//   std::ofstream rhs_stream("rhs.m");
//   srhs.matlab_output(rhs_stream);
//   rhs_stream.close();
   /*
 */


  return 0;
}

