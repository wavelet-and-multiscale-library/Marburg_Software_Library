#include <iostream>
#define _WAVELETTL_GALERKINUTILS_VERBOSITY 1

#undef DYADIC
#undef TRIVIAL
#define ENERGY

#define JACO
#undef COGRAD
#undef OTHER

#define JMAX 6
#define PMAX 1
#define ONE_D
#define _DIM 1

#undef BASIS
#define FRAME
#undef DELTADIS

#define NONADAPTIVE


#define PRIMALORDER 2
#define DUALORDER   2
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


#include <galerkin/sturm_equation.h>
#include <galerkin/galerkin_utils.h>
#ifdef BASIS
#include <interval/p_expansion.h>
#include <galerkin/cached_problem.h>
#endif
#ifdef FRAME
#include <galerkin/cached_quarklet_problem.h>
#endif
#include <galerkin/TestProblem.h>
#include <galerkin/TestFunctions.h>



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

    const unsigned int testcase=8;
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
          case 9:
              //dummy. needs to be replacedwith the correct function. Not yet implemented
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
  
  

  const int d  = PRIMALORDER;
  const int dT = DUALORDER; // be sure to use a continuous dual here, otherwise the RHS test will fail
  
  const int jmax = JMAX;
  const int pmax = PMAX;

  
  
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
  cout << "setup equation.." << endl;
  SturmEquation<Basis> eq(T, basis);
  cout << "end setup equation!" << endl;
  
 

#ifdef NONADAPTIVE
  //nonadaptive setting
  
//  const int j0=eq.basis().j0();
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

    

  
  for (int i=0; i<basis.degrees_of_freedom();i++) {
//  for (int i=0; i<928;i++) {    
      if(i>=0){
        Lambda.insert(*(basis.get_quarklet(i)));
        cout << *(basis.get_quarklet(i)) << endl;
      }
  }


  
  
  
  
  
  
  
  


    
    
    
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
  
  
  
//  abort();

  
  
  
  

  Vector<double> x(Lambda.size()), x2(Lambda.size()), err(Lambda.size()), err2(Lambda.size()); x = 0, x2 = 0;
//  Vector<double> y(Lambda.size()); y = 0;
//  y[15]=1;
  
#ifdef JACO
  unsigned int iterations;
//GaussSeidel(A, b, x, 1e-8, 5000, iterations);
  double omega=0.5;
Richardson(A, b, x, omega, 1e-12, 1e+6, iterations);
cout << "Number of iterations: " << iterations << endl;
#endif
#ifdef COGRAD  
unsigned int iterations;
CG(A, b, x, 1e-8, 5000, iterations);
cout << "Number of iterations: " << iterations << endl;
#endif
#ifdef OTHER
double lmax = 1;
double alpha_n = 2. / lmax - 0.001;
  
  Vector<double> resid(x.size());
  Vector<double> help(x.size());
  for (int i = 0; i < 10000;i++) {
    A.apply(x,help);
    resid = b - help;
    cout << "i = " << i << " " << sqrt(resid*resid) << endl;
    A.apply(resid,help);
    alpha_n = (resid * resid) * (1.0 / (resid * help));
    resid *= alpha_n;
    x = x + resid;
    if(sqrt(resid*resid)<1e-8)
        break;
  }
#endif

  

//A.apply(x, err);
  A.apply(x,err);
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
    
//    std::ofstream u_stream1(filenameSolution1);
//  s.matlab_output(u_stream1);
//  u_stream1 << "figure;\nplot(x,y);"
//            << "title('Sturm bvp: solution to test problem ("
//            << basis_type << "), " << "pmax= " << pmax << ", " << "jmax= " << jmax << ", " << "d= " << d << "');" << endl;
//  u_stream1.close();
    
     

    // compute some error-norms
    const double Linfty_error = linfty_norm(u_values-uexact_values);
    cout << "  L_infinity error on a subgrid: " << Linfty_error << endl;
    
    
    const double L2_error = sqrt(l2_norm_sqr(u_values-uexact_values)*h);
    cout << "  L_2 error on a subgrid: " << L2_error << endl;
    
    cout << " \ell_1 Norm Solution: " << l1_norm(x) << endl;
    
#endif
  
  
//  cout << "- point values of the solution:" << endl;
  InfiniteVector<double,Index> u,v,w,testv;
  unsigned int i = 0;
  for (set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it, ++i)
    u.set_coefficient(*it, x[i]);
//  u.set_coefficient(*(basis.get_quarklet(15)),1);
    
  
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
  

  
  //plots
  
  
  
  /* plot solution coefficients */
  // output for Frame has to be renewed
  
#ifdef FRAME

  for(int p=0;p<=pmax;p++){
    char filenameCoefficients2[128];
    sprintf(filenameCoefficients2, "%s%d%s", "sturm_bvp_solution_coefficients_p_" , p , ".m");
    std::ofstream coeff_stream2 (filenameCoefficients2);
  coeff_stream2 << "figure;" << endl;
  plot_indices_iq(&basis, u, jmax, coeff_stream2, p, "jet", false, true, -8);
  coeff_stream2 << "title('coefficients on the level p=" << p <<" of the test problem ("
                  << basis_type << " basis)');" << endl;
  coeff_stream2.close();   
  }
  
  cout << v << endl;
  
  
  

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
            << basis_type << "), " << "pmax= " << pmax << ", " << "jmax= " << jmax << ", " << "d= " << d << "');" << endl;
  u_stream1.close();
  
#endif
  

  return 0;
}
