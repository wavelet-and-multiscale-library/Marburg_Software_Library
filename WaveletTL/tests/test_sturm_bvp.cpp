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
#undef BASIS
#define FRAME
#undef DELTADIS

#include <galerkin/sturm_equation.h>
#include <galerkin/galerkin_utils.h>

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
      return 1;
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

  TestProblem<2> T;

  const int d  = 2;
  const int dT = 2; // be sure to use a continuous dual here, otherwise the RHS test will fail
  
  const int jmax = 3;
#ifdef FRAME
  const int pmax = 2;
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

  Basis basis(0,0, false);
  typedef Basis::Index Index;
  const char* basis_type = "Primbs quarklet frame";
  basis.set_jpmax(jmax,pmax);
  
#endif
  
  SturmEquation<Basis> eq(T, basis);
  const int j0=eq.basis().j0();
  InfiniteVector<double, Index> coeffs;

  
#if 1
 //Testing some entrys of the stiffness matrix
  Index mu(1,3,0,2, &basis);
  Index nu (1,3,0,2, &basis);
  double entry=eq.a(mu,nu);
  cout << "Entry in stiffness matrix A for indices " << mu << "and " << nu << ": " << entry << endl;
#endif

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
  
  
  

#if 1  
     cout << "- set up stiffness matrix with respect to the index set Lambda=" << endl;
     for (set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it)
       cout << *it << endl;


  cout << "- set up (preconditioned) stiffness matrix (j0=" << eq.basis().j0() << ",jmax=" << jmax << ")..." << endl;
  clock_t tstart, tend;
  double time;
  tstart = clock(); 
  SparseMatrix<double> A;
  setup_stiffness_matrix(eq, Lambda, A);
  tend = clock();
  time = (double)(tend-tstart)/CLOCKS_PER_SEC;
  cout << "  ... done, time needed: " << time << " seconds" << endl;
   //cout << "- (preconditioned) stiffness matrix A=" << endl << A << endl;

   cout << "- writing A to the file stiffmat.m ..." << endl;
   std::ofstream Astream("stiffmat.m");
   Astream << "A=";
   print_matrix(A, Astream);
   Astream << ";" << endl;
   Astream.close();
   cout << "  ... done!" << endl;

  cout << "- set up right-hand side..." << endl;
  tstart = clock();
  Vector<double> b;
  setup_righthand_side(eq, Lambda, b);
  tend = clock();
  time = (double)(tend-tstart)/CLOCKS_PER_SEC;
  cout << "  ... done, time needed: " << time << " seconds" << endl;
  
  cout << "- writing b to the file rhs.m ..." << endl;
   std::ofstream bstream("rhs.m");
   bstream << "b=";
   print_vector(b, bstream);
   bstream << ";" << endl;
   bstream.close();
   cout << "  ... done!" << endl;

  
  b.compress(1e-14);
  cout << "- right hand side: " << b << endl << endl;

  Vector<double> x(Lambda.size()), err(Lambda.size()), err2(Lambda.size()); x = 0;
  unsigned int iterations;
  
  
#ifdef FRAME
//  Matrix<double> evecs;
//  Vector<double> evals;
//  SymmEigenvalues(A, evals, evecs); // for safety, we compute ALL eigenvalues
////   double lambdamax = PowerIteration(A, xk, 1e-6, 1000, iterations);
////   double lambdamin = 1./InversePowerIteration(A, xk, 1e-6, 1000, iterations);
//  double lambdamax = evals[evals.size()-1];
//  int j=0;
//  while(evals[j]<1e-1) j++;
//  double lambdamin = evals[j];
//  cout << "  lambdamax=" << lambdamax
//       << ", lambdamin=" << lambdamin << endl;
//  double omegastar = 2./(lambdamin+lambdamax);
//  cout << "omegastar: " << omegastar << endl;
//  Richardson(A, b, x, omegastar, 1e-8, 10000, iterations);
//  GaussSeidel(A, b, x, 1e-8, 10000, iterations);
  CG(A, b, x, 1e-8, 200, iterations);
#else
  CG(A, b, x, 1e-8, 200, iterations);
#endif
  
  
 
  
//  Vector<double> x2(Lambda.size()," -90.69813778142046 277.6049681626802 137.763924751933 3740.815469402151 -553.5220285733511 -2571.078754156906 -880.6581000170487 -90.82424114257628 -29.46472526365101 -32.48153042189944 340.9701025624454 -556.4234231003472 -313.6400421185935 67.8221987536722 -32.59435947723999 -29.43250221197568 -2708.989044118656 -6729.763648330345 -35809.96876058168 -1453.428184644519 14946.34582699998 6561.842913774125 -908.4681731417184 2406.360859896055 464.522704015511 1277.593500413668 -5807.643139152126 -12422.6705934792 55752.00389593931 55790.39412933643 -12372.14758315131 -5783.017058896096 -554.5402161403698 -4307.107220873145 -4307.086264666619 -553.7583062933044"); 
 Vector<double> x2(Lambda.size()," 0.4995196796211758 0.4995930793194577 0.4993781557876391 0.4995930793195578 0.499519679621447 0.0001961663434695258 0.0004941575858113331 0.0004941575858145899 0.0001961663432735496 7.059617586933921e-05 4.386870953756039e-05 6.423344769786667e-05 7.433745784499037e-05 7.433745776920209e-05 6.42334477370391e-05 4.386870945398857e-05 7.059617587131275e-05 -0.00296136681093653 -3.799663503842528e-12 0.002961366805766139 -0.001841220606581801 0.001841220605553969 -1.056830778084215e-05 -0.0001740815697862949 -0.0001232956657777421 0.0001232956657720723 0.0001740815697222662 1.056830781413941e-05 0.00203775954392417 0.03041258804115368 0.00203775955151338 -0.0003098495771240994 -0.0003098495767476009 -3.977715111333296e-05 2.724649018440909e-05 -1.726262057569815e-05 -1.726262080100297e-05 2.724649042325838e-05 -3.977715130052654e-05");
//  readVectorFromFile(x2,"myfile.m");
  //cout << "x2: " << x2 << endl;
  
  A.apply(x2,err2);
  
  cout << "Ax2" << err2 << endl << endl;
  
  //x.compress(1e-4);
  
  A.apply(x, err);
  cout << "Ax " << err << endl << endl;
  
    
  err -= b;
  err2-=b;
  cout << "residual " << err << endl;
  cout << "residual Matlab " << err2 << endl;
  cout << " residual (infinity) norm " << linfty_norm(err) << endl;
  cout << " residual (infinity) norm Matlab " << linfty_norm(err2) << endl;
  
//  cout << "- point values of the solution:" << endl;
  InfiniteVector<double,Index> u,v;
  unsigned int i = 0;
  for (set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it, ++i)
//    u.set_coefficient(*it, x[i]);
    u.set_coefficient(*it, x2[i]);
  
  u.COARSE(1e-6,v);
  //v.set_coefficient(Index(0,2,0,0,&basis),1);
  
  
  
  cout << "solution: " << endl << u << endl;
  cout << "coarsed solution: " << endl << v << endl;
  cout << "Number of iterations: " << iterations << endl;

  
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
    
  u.scale(&eq, -1);
  SampledMapping<1> s(evaluate(eq.basis(), u, true, 7));
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
