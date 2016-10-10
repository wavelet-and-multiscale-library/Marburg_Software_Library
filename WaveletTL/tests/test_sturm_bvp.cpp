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

#undef ADAPTIVE
#define NONADAPTIVE

#include <galerkin/sturm_equation.h>
#include <galerkin/galerkin_utils.h>
#ifdef BASIS
#include <galerkin/cached_problem.h>
#endif
#ifdef FRAME
#include <galerkin/cached_quarklet_problem.h>
#endif
#include <galerkin/TestProblem.h>
#include <galerkin/TestFunctions.h>
#include <adaptive/cdd2.h>

using namespace std;
using namespace WaveletTL;

using MathTL::SimpleSturmBVP;
using MathTL::CG;


//Auslagern in galerkin/TestProblem.h
/*
  different test problems with homogeneous Dirichlet b.c.'s
  1: y(t)=x*(1-x), -y''(t)=2
  2: y(t) = exp(-100*(x-0.5)^2), -y''(t) = (200-(200x-100)^2)*exp(-100*(x-0.5)^2)
  3: 1D example from [BBCCDDU]
  4: y(t) = -sin(3*M_PI*t), -y''(t) = -9*M_PI*M_PI*sin(3*M_PI*t)
  5: y(t) = t*t*(1-t) 
  6: "Hat function as gramian problem (for test purposes)";  
  7: "Another hat function with bigger support as gramian problem (for test purposes)";
  8: "Quark function";*/

int main()
{
  cout << "Testing wavelet-Galerkin solution of a Sturm b.v.p. ..." << endl;

  const unsigned int testcase=1;
  TestProblem<testcase> T;
  Function<1>* uexact = 0;
    switch(testcase) {
        case 1:
            uexact = new Function2();
            break;
        case 2:
            uexact = new Function2b();
            break;
        case 3:
            uexact = new FunctionBBCCDDU();
            break;
        case 4:
            uexact = new Function4();
            break;
        case 5:
            uexact = new Function2a();
            break;
        case 6:
            uexact = new scaledHat2a();
            break;
        case 7:
            uexact = new scaledHat2b();
            break;
        case 8:
            uexact = new scaledQuark();
            break;
        default:
            break;
    }
    
    // evaluate exact solution
    const unsigned int N = 100;
    const double h = 1./N;
    Vector<double> uexact_values(N+1);
    for (unsigned int i = 0; i <= N; i++) {
      const double point = i*h;
      uexact_values[i] = uexact->value(Point<1>(point));
    }
  
  

  const int d  = 2;
  const int dT = 2; // be sure to use a continuous dual here, otherwise the RHS test will fail
  
  const int jmax = 6;
  const int pmax = 2;
  
  
#ifdef BASIS
  typedef PBasis<d,dT> Basis;
  Basis basis(1,1);
  typedef Basis::Index Index;
  const char* basis_type = "Primbs basis";
  basis.set_jmax(jmax);
#endif
  
#ifdef FRAME
  typedef PQFrame<d,dT> Basis;
  Basis basis(true, true, false);
  typedef Basis::Index Index;
  const char* basis_type = "Primbs quarklet frame";
  basis.set_jpmax(jmax,pmax);
  
#endif
  
  SturmEquation<Basis> eq(T, basis);
  const int j0=eq.basis().j0();
  

#ifdef NONADAPTIVE
  //nonadaptive setting
  
#if 0
 //Testing some entrys of the stiffness matrix
  Index mu(1,3,0,2, &basis);
  Index nu (1,3,0,2, &basis);
  double entry=eq.a(mu,nu);
  cout << "Entry in stiffness matrix A for indices " << mu << "and " << nu << ": " << entry << endl;
#endif



  

#if 0
  InfiniteVector<double, Index> coeffs;
  eq.RHS(1e-2, coeffs);
  cout << "- approximate coefficient set of the right-hand side:" << endl
       << coeffs << endl;

#endif  

  //Implementation of the index set
  set<Index> Lambda;
#ifdef FRAME
    
#if 1
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
#elif 0  
      for (Index lambda = eq.basis().first_generator(j0,0);;) {
         Lambda.insert(lambda);
         if (lambda == eq.basis().last_generator(j0,0)) break; 
         else
             ++lambda;
      }
#else
    int p = 0;
    
    for (Index lambda = eq.basis().first_generator(j0,0);;) {
	Lambda.insert(lambda);
	if (lambda == eq.basis().last_generator(j0,pmax)) break;
        if (lambda == eq.basis().last_generator(j0,p)){
            ++p;
            lambda = eq.basis().first_generator(j0,p);
        }
        else
            ++lambda;
      }
#endif
#else
  

  for (Index lambda = eq.basis().first_generator(j0);; ++lambda) {
    Lambda.insert(lambda);
    if (lambda == eq.basis().last_wavelet(jmax)) break;    
  }
#endif
  
  
  


    
    
    
  //gives the index set
//     cout << "- set up stiffness matrix with respect to the index set Lambda=" << endl;
//     for (set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it)
//       cout << *it << endl;


  cout << "- set up (preconditioned) stiffness matrix (j0=" << eq.basis().j0() << ",jmax=" << jmax << ",pmax=" << pmax << ")..." << endl;
  SparseMatrix<double> A;
  setup_stiffness_matrix(eq, Lambda, A);
  
  
  


   cout << "- writing A to the file stiffmat.m ..." << endl;
   std::ofstream Astream("stiffmat.m");
   Astream << "A=";
   print_matrix(A, Astream);
   Astream << ";" << endl;
   Astream.close();
   cout << "  ... done!" << endl;

  cout << "- set up right-hand side..." << endl;
  Vector<double> b;
  setup_righthand_side(eq, Lambda, b);
  
  
  cout << "- writing b to the file rhs.m ..." << endl;
   std::ofstream bstream("rhs.m");
   bstream << "b=";
   print_vector(b, bstream);
   bstream << ";" << endl;
   bstream.close();
   cout << "  ... done!" << endl;

  
  b.compress(1e-14);
  

  Vector<double> x(Lambda.size()), x2(Lambda.size()), err(Lambda.size()), err2(Lambda.size()); x = 0, x2 = 0;
  unsigned int iterations;
  
  
CG(A, b, x, 1e-8, 5000, iterations);

A.apply(x, err);
  //cout << "Ax " << err << endl << endl;
  
    
  err -= b;
  //cout << "residual " << err << endl;
  cout << " residual (infinity) norm " << linfty_norm(err) << endl;
  
#if 1
  
  // evaluate approximate solution on a grid
    Vector<double> u_values(N+1);
    for (unsigned int i = 0; i <= N; i++) {
      const double point = i*h;
      int id = 0;
      for (set<Index>::const_iterator it(Lambda.begin()); it != Lambda.end(); ++it, ++id) 
	u_values[i] += x[id] * basis.evaluate(0, *it, point)*1./eq.D(*it);
    }
     

    
    
    
    
    
     

    // compute some error-norms
    const double Linfty_error = linfty_norm(u_values-uexact_values);
    cout << "  L_infinity error on a subgrid: " << Linfty_error << endl;
    
    
    const double L2_error = sqrt(l2_norm_sqr(u_values-uexact_values)*h);
    cout << "  L_2 error on a subgrid: " << L2_error << endl;
    
#endif
  
  
//  cout << "- point values of the solution:" << endl;
  InfiniteVector<double,Index> u,v,w,testv;
  unsigned int i = 0;
  for (set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it, ++i)
    u.set_coefficient(*it, x[i]);
    
  
  u.COARSE(1e-6,v);
  //v.set_coefficient(Index(0,2,0,0,&basis),1);
  
  //optional outputs
//  Matrix<double> evecs;
//  Vector<double> evals;
//  SymmEigenvalues(A,evals,evecs);
//  cout<< "Eigenwerte A: " << evals << endl;
//  cout << "- (preconditioned) stiffness matrix A=" << endl << A << endl;
//  cout << "- right hand side: " << b << endl << endl;
//  cout << "  point values of Galerkin solution: " << u_values << endl;
//  cout << "  point values of exact solution: " << uexact_values << endl;
//  cout << " pointwise error (exact solution): " << u_values-uexact_values << endl;
//  cout << "solution: " << endl << u << endl;
//  cout << "coarsed solution: " << endl << v << endl;
  cout << "Number of iterations: " << iterations << endl;

  
  //plots
  
  
  
  /* plot solution coefficients */
  // output for Frame has to be renewed
  
#ifdef FRAME
#if 0
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
#endif
#else
    const char* filenameCoefficients1 = "sturm_bvp_solution_coefficients.m";
    std::ofstream coeff_stream1;
    coeff_stream1.open(filenameCoefficients1);
    coeff_stream1 << "figure;" << endl;
    plot_indices(&basis, v, jmax, coeff_stream1, "jet", false, true, -8);
    coeff_stream1 << "title('Sturm bvp: solution coefficients of the test problem ("
                  << basis_type << ")');" << endl;
    coeff_stream1.close();
#endif
    
  
 /* plot solution*/ 
  const char* filenameSolution1 = "sturm_bvp_solution.m";  
  u.scale(&eq, -1);
  //u.set_coefficient(Index(1,3,0,1,&basis), u[Index(1,3,0,1,&basis)]+1);
//  w.set_coefficient(Index(1,3,1,1,&basis), 1);
//  basis.reconstruct(w,4,testv);
//  cout << "reconstruction: " << testv << endl;
//  cout << "scaled solution: " << endl << u << endl;
  SampledMapping<1> s(evaluate(eq.basis(), u, true, 7));
//  SampledMapping<1> s(evaluate(eq.basis(), Index(1,3,0,1,&basis), true, 7));
  std::ofstream u_stream1(filenameSolution1);
  s.matlab_output(u_stream1);
  u_stream1 << "figure;\nplot(x,y);"
            << "title('Sturm bvp: solution to test problem ("
            << basis_type << "), " << "pmax= " << pmax << ", " << "d= " << d << "');" << endl;
  u_stream1.close();
  
#endif
  
#ifdef ADAPTIVE
  
  //adaptive setting (CDD2)
  
#ifdef FRAME  
  CachedQuarkletProblem<SturmEquation<Basis> > ceq(&eq);
#endif
#ifdef BASIS
  CachedProblem<SturmEquation<Basis> > ceq(&eq);
#endif
  InfiniteVector<double, Index> F_eta;
  ceq.RHS(1e-6, F_eta);
  const double nu = ceq.norm_Ainv() * l2_norm(F_eta);
  double epsilon = 1e-6;
  InfiniteVector<double, Index> u_epsilon;
#ifdef FRAME
  CDD2_SOLVE(ceq, nu, epsilon, u_epsilon, jmax, DKR, pmax, 2, 2);
#endif
#ifdef BASIS
  CDD2_SOLVE(ceq, nu, epsilon, u_epsilon, jmax);
#endif



  
  // evaluate approximate solution on a grid
    Vector<double> u_values_ad(N+1);
    for (unsigned int i = 0; i <= N; i++) {
      const double point = i*h;
      for (InfiniteVector<double, Index>::const_iterator it(u_epsilon.begin()); it != u_epsilon.end(); ++it) 
	u_values_ad[i] += *it * basis.evaluate(0, it.index(), point)*1./ceq.D(it.index());
    }
    
    // compute some error-norms
    const double Linfty_error_ad = linfty_norm(u_values_ad-uexact_values);
    cout << "  L_infinity error on a subgrid (adaptive): " << Linfty_error_ad << endl;
    
    
    const double L2_error_ad = sqrt(l2_norm_sqr(u_values_ad-uexact_values)*h);
    cout << "  L_2 error on a subgrid (adaptive): " << L2_error_ad << endl;
    
    
/* plot solution coefficients */
  // output for Frame has to be renewed
  
#ifdef FRAME
#if 0
  string filenameCoefficients2[4] = {"sturm_bvp_solution_coefficients_p_0_ad.m", "sturm_bvp_solution_coefficients_p_1_ad.m",
  "sturm_bvp_solution_coefficients_p_2_ad.m", "sturm_bvp_solution_coefficients_p_3_ad.m"};
  
//  string filenameCoefficients2[1] = {"sturm_bvp_solution_coefficients_p_0_ad.m"};
  
  for(int p=0;p<=pmax;p++){
  const char* cstr = filenameCoefficients2[p].c_str();
  cout << filenameCoefficients2[p] << endl;
  std::ofstream coeff_stream2 (cstr);
  coeff_stream2 << "figure;" << endl;
  plot_indices(&basis, u_epsilon, jmax, coeff_stream2, p, "jet", false, true, -8);
  coeff_stream2 << "title('adaptive coefficients on the level p=" << p <<" of the test problem ("
                  << basis_type << " basis)');" << endl;
  coeff_stream2.close();   
  }
#endif
#else
    const char* filenameCoefficients2 = "sturm_bvp_solution_coefficients_adaptive.m";
    std::ofstream coeff_stream2;
    coeff_stream2.open(filenameCoefficients2);
    coeff_stream2 << "figure;" << endl;
    plot_indices(&basis, u_epsilon, jmax, coeff_stream2, "jet", false, true, -8);
    coeff_stream2 << "title('Sturm bvp: adaptive solution coefficients of the test problem ("
                  << basis_type << ")');" << endl;
    coeff_stream2.close();
#endif    
  
//plot of the adaptive solution
  const char* filenameSolution2 = "sturm_bvp_solution_adaptive.m";  
  u_epsilon.scale(&ceq, -1);
  SampledMapping<1> s2(evaluate(ceq.basis(), u_epsilon, true, 7));
  std::ofstream u_stream2(filenameSolution2);
  s2.matlab_output(u_stream2);
  u_stream2 << "figure;\nplot(x,y);"
            << "title('Sturm bvp: adaptive solution to test problem ("
            << basis_type << "), " << "pmax= " << pmax << ", " << "d= " << d << "');" << endl;
  
  u_stream2.close();

#endif
  return 0;
}
