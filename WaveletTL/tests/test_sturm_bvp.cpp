#include <iostream>
#define _WAVELETTL_GALERKINUTILS_VERBOSITY 1
#include <fstream>
#include <map>
#include <time.h>

#include <algebra/symmetric_matrix.h>
#include <algebra/sparse_matrix.h>
#include <numerics/sturm_bvp.h>
#include <numerics/iteratsolv.h>

#include <interval/i_index.h>
#include <interval/i_indexplot.h>
#include <interval/i_q_index.h>
#include <interval/i_q_indexplot.h>
#include <interval/ds_basis.h>
#include <interval/p_basis.h>
#include <interval/pq_frame.h>
#include <interval/jl_basis.h>
#include <interval/jl_support.h>
#include <interval/jl_evaluate.h>
#include <interval/spline_basis.h>
#undef BASIS
#define FRAME
#undef DELTADIS

#include <galerkin/sturm_equation.h>

using namespace std;
using namespace WaveletTL;

using MathTL::SimpleSturmBVP;
using MathTL::CG;

/*
  different test problems with homogeneous Dirichlet b.c.'s
  1: y(t)=x*(1-x), -y''(t)=2
 */
template <unsigned int N>
class TestProblem
  : public SimpleSturmBVP
{
public:
  double p(const double t) const {
    switch(N) {
    case 1:
      return 1;
      break;
    case 2:
      return 0;
      break;
    case 3:
      return 1;
      break;
    case 4:
      return 1;
      break;
    case 5:
       return 1;
       break;
    default:
      return 0;
      break;
    }
  }
  double p_prime(const double t) const {
    switch(N) {
    case 1:
      return 0;
      break;
    case 2:
      return 0;
      break;
    case 3:
      return 0;
      break;
    case 4:
      return 0;
      break;
    case 5:
      return 0;
      break;
    default:
      return 0;
      break;
    }
  }
  double q(const double t) const {
    switch(N) {
    case 1:
      return 0;
      break;
    case 2:
      return 1;
      break;
    case 3:
      return 0;
      break;
    case 4:
      return 0;
      break;
    case 5:
      return 0;
      break;  
    default:
      return 0;
      break;
    }
  }
  double g(const double t) const {
    switch(N) {
    case 1:
      return 2;
      break;
    case 2:
      return t*(1-t);
      break;
    case 3:
      //return exp(-100*(t-0.5)*(t-0.5));
      return (200-(200*t-100)*(200*t-100))*exp(-100*(t-0.5)*(t-0.5));
      break;
    case 4:
      return ( -100*exp(5*t)*(1-(exp(5*t)-1)/(exp(5.)-1)) /
             (exp(5.)-1)+200*exp(10*t)/((exp(5.)-1) *
             (exp(5.)-1))+100*(exp(5*t)-1)*exp(5*t) /
             ((exp(5.)-1)*(exp(5.)-1)) );
      break; 
    case 5:
      return -9*M_PI*M_PI*sin(3*M_PI*t)-4;  
      break;
    default:
      return 0;
      break;
    }
  }
  bool bc_left() const { return true; }
  bool bc_right() const { return true; }
};


int main()
{
  cout << "Testing wavelet-Galerkin solution of a Sturm b.v.p. ..." << endl;

  TestProblem<1> T;

  const int d  = 3;
  const int dT = 3; // be sure to use a continuous dual here, otherwise the RHS test will fail
  
  const int jmax = 4;
#ifdef FRAME
  const int pmax = 1;
#endif

  // typedef DSBasis<d,dT> Basis; //Basis basis(true, true);
  // typedef PBasis<d,dT> Basis;
  // typedef JLBasis Basis;
  // typedef SplineBasis<d,dT,P_construction,1,1,0,0,SplineBasisData_j0<d,dT,P_construction,1,1,0,0>::j0> Basis;
  
#ifdef BASIS
  typedef PBasis<d,dT> Basis;

  Basis basis(1,1);
  typedef Basis::Index Index;
  const char* basis_type = "Primbs basis";
  basis.set_jmax(jmax);
#endif
  
#ifdef FRAME
  typedef PQFrame<d,dT> Basis;

  Basis basis(1,1, false);
  typedef Basis::Index Index;
  const char* basis_type = "Primbs quarklet frame";
  basis.set_jpmax(jmax,pmax);
  
#endif
  
  SturmEquation<Basis> eq(T, basis);
  const int j0=eq.basis().j0();
  InfiniteVector<double, Index> coeffs;

#if 0
  coeffs.set_coefficient(eq.basis().first_generator(eq.basis().j0()), 1.0);
  coeffs.set_coefficient(eq.basis().last_generator(eq.basis().j0()), 2.0);
  coeffs.set_coefficient(eq.basis().first_wavelet(eq.basis().j0()), 3.0);
  coeffs.set_coefficient(eq.basis().last_wavelet(eq.basis().j0()), 4.0);
  coeffs.set_coefficient(eq.basis().first_wavelet(eq.basis().j0()+1), 5.0);
  cout << "- a coefficient set:" << endl
       << coeffs << endl;
  coeffs.scale(&eq, -1);
  cout << "- after rescaling with D^{-1}:" << endl
       << coeffs << endl;
#endif

  eq.RHS(1e-2, coeffs);

#if 1
  cout << "- approximate coefficient set of the right-hand side:" << endl
       << coeffs << endl;
//   cout << "- check expansion of the right-hand side in the dual basis:" << endl;
//   coeffs.scale(&eq, 1);
//   evaluate(eq.basis(), coeffs, false, 8).matlab_output(cout);
//   coeffs.scale(&eq, -1);
#endif  

  //Implementation of the index set
  set<Index> Lambda;
#ifdef FRAME
    int p = 0;
      
      for (Index lambda = eq.basis().first_generator(j0,0);;) {
	Lambda.insert(lambda);
	if (lambda == eq.basis().last_wavelet(jmax,pmax)) break;
        //if (i==7) break;
        if (lambda == eq.basis().last_wavelet(jmax,p)){
            ++p;
            lambda = eq.basis().first_generator(j0,p);
        }
        else
            ++lambda;
      }
#else
  
//  eq.basis().setWavelets();
  for (Index lambda = eq.basis().first_generator(j0);; ++lambda) {
    Lambda.insert(lambda);
    if (lambda == eq.basis().last_wavelet(jmax)) break;    
  }
#endif
  
  
  

  
     cout << "- set up stiffness matrix with respect to the index set Lambda=" << endl;
     for (set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it)
       cout << *it << endl;

#if 1
  cout << "- set up (preconditioned) stiffness matrix (j0=" << eq.basis().j0() << ",jmax=" << jmax << ")..." << endl;
  clock_t tstart, tend;
  double time;
  tstart = clock(); 
  SparseMatrix<double> A;
  setup_stiffness_matrix(eq, Lambda, A);
  tend = clock();
  time = (double)(tend-tstart)/CLOCKS_PER_SEC;
  cout << "  ... done, time needed: " << time << " seconds" << endl;
//   cout << "- (preconditioned) stiffness matrix A=" << endl << A << endl;

//   cout << "- writing A to the file stiffmat.m ..." << endl;
//   std::ofstream Astream("stiffmat.m");
//   Astream << "A=";
//   print_matrix(A, Astream);
//   Astream << ";" << endl;
//   cout << "  ... done!" << endl;

  cout << "- set up right-hand side..." << endl;
  tstart = clock();
  Vector<double> b;
  setup_righthand_side(eq, Lambda, b);
  tend = clock();
  time = (double)(tend-tstart)/CLOCKS_PER_SEC;
  cout << "  ... done, time needed: " << time << " seconds" << endl;

  cout << "- right hand side: " << b << endl << endl;

  Vector<double> x(Lambda.size()), err(Lambda.size()), err2(Lambda.size()); x = 0;
  unsigned int iterations;
  
  CG(A, b, x, 1e-8, 200, iterations);
  
  
  A.apply(x, err);
  //cout << "Ax " << err << endl << endl;
  
    
  err -= b;
  //cout << "residual " << err << endl;
  cout << " residual (infinity) norm " << linfty_norm(err) << endl;
  
  
//  cout << "- point values of the solution:" << endl;
  InfiniteVector<double,Index> u,v;
  unsigned int i = 0;
  for (set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it, ++i)
    u.set_coefficient(*it, x[i]);
  
  
  u.COARSE(1e-6,v);
  
  
  
  cout << "solution: " << endl << u << endl;
  cout << "coarsed solution: " << endl << v << endl;

  
  const char* filenameSolution1 = "sturm_bvp_solution.m";
  
  
  /* plot solution coefficients */
  
#ifdef FRAME
  string filenameCoefficients[4] = {"sturm_bvp_solution_coefficients_p_0.m", "sturm_bvp_solution_coefficients_p_1.m",
  "sturm_bvp_solution_coefficients_p_2.m", "sturm_bvp_solution_coefficients_p_3.m"};
  
//  string filenameCoefficients[1] = {"sturm_bvp_solution_coefficients_p_0.m"};
  
  for(int p=0;p<=pmax;p++){
  const char* cstr = filenameCoefficients[p].c_str();
  cout << filenameCoefficients[p] << endl;
  std::ofstream coeff_stream1 (cstr);
  coeff_stream1 << "figure;" << endl;
  plot_indices(&basis, u, jmax, coeff_stream1, p, "jet", false, true, -8);
  coeff_stream1 << "title('coefficients on the level p=" << p <<" of the test problem ("
                  << basis_type << " basis)');" << endl;
  coeff_stream1.close();    
  }
#else
    const char* filenameCoefficients1 = "sturm_bvp_solution_coefficients.m";
    std::ofstream coeff_stream1;
    coeff_stream1.open(filenameCoefficients1);
    coeff_stream1 << "figure;" << endl;
    plot_indices(&basis, v, jmax, coeff_stream1, "jet", false, true, -8);
    coeff_stream1 << "title('Sturm_bvp: solution coefficients of the test problem ("
                  << basis_type << ")');" << endl;
    coeff_stream1.close();
#endif
    
  
 /* plot solution*/ 
    
  v.scale(&eq, -1);
  SampledMapping<1> s(evaluate(eq.basis(), v, true, 7));
  std::ofstream u_stream1(filenameSolution1);
  s.matlab_output(u_stream1);
  u_stream1 << "figure;\nplot(x,y);"
            << "title('Sturm bvp: solution to test problem ("
            << basis_type << ")');" << endl;
  u_stream1.close();

#endif

#if 0
  cout << "- estimate for ||D^{-1}AD^{-1}||: " << eq.norm_A() << endl;
  cout << "- estimate for ||(D^{-1}AD^{-1})^{-1}||: " << eq.norm_Ainv() << endl;
  // von Ulli importiert:
      double help, normA;
      unsigned int iterations2;
      LanczosIteration(A, 1e-6, help, normA, 200, iterations2);
      double normAinv ( 1./help);
      cout << "normA = "<<normA<<endl;
      cout << "normAinv = "<<normAinv<<endl;
      cout << "kond = "<<(normA*normAinv)<<endl;
#endif

#if 0
  cout << "- checking add_column and the compression strategy:" << endl;
  InfiniteVector<double,Index> v, w;
  for (set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it) {
    v.clear();
    w.clear();
    eq.add_column(1.0, *it, 99, w);
    cout << "lambda=" << *it << endl;
//     cout << "* result of add_column:" << endl << w << endl;
    const int jmax = 8;
    cout << "* checking correctness up to level " << jmax << "..." << endl;
    for (Index nu = first_generator(&eq.basis(), eq.basis().j0());; ++nu) {
      v.set_coefficient(nu, eq.a(nu, *it)/(eq.D(nu)*eq.D(*it)));
      if (nu == last_wavelet(&eq.basis(), jmax)) break;
    }
//     cout << "  ... done, error:" << endl;
//     InfiniteVector<double,Index> error = v - w;
//     error.compress();
//     cout << error << endl;

    cout << "  ... done, absolute error " << linfty_norm(v-w) << endl;

    break;
  }

  cout << "- check add_column for a second time to test the cache:" << endl;
  for (set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it) {
    v.clear();
    w.clear();
    eq.add_column(1.0, *it, 99, w);
    cout << "lambda=" << *it << endl;
//     cout << "* result of add_column:" << endl << w << endl;
    const int jmax = 8;
    cout << "* checking correctness up to level " << jmax << "..." << endl;
    for (Index nu = first_generator(&eq.basis(), eq.basis().j0());; ++nu) {
      v.set_coefficient(nu, eq.a(nu, *it)/(eq.D(nu)*eq.D(*it)));
      if (nu == last_wavelet(&eq.basis(), jmax)) break;
    }
    cout << "  ... done, absolute error " << linfty_norm(v-w) << endl;

    break;
  }
#endif

//   InfiniteVector<double,Index> coeffsp;
//   eq.left_preconditioner()->apply_preconditioner(coeffs,coeffsp);

  return 0;
}
