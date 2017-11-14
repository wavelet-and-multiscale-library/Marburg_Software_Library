#include <iostream>
#define _WAVELETTL_GALERKINUTILS_VERBOSITY 1

#define PARALLEL 0

#undef DYADIC
#define ENERGY
#undef TRIVIAL

#undef SD
#define RICHARDSON
#undef CONGRAD

#undef SHRINKAGE


#define JMAX 15
#define PMAX 0
#define ONE_D

#define PRIMALORDER 3
#define DUALORDER   3
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
#define DELTADIS

#define ADAPTIVE
#undef NONADAPTIVE

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
#include <adaptive/cdd2.h>
#include <adaptive/steepest_descent_ks.h>
#include <adaptive/duv.h>
#include <adaptive/apply.h>
#include <adaptive/compression.h>


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

    const unsigned int testcase=9;
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
          case 10:
              uexact = new SolutionToDeltadis();
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
  Basis basis(true, true, true);
  typedef Basis::Index Index;
  const char* basis_type = "Primbs quarklet frame";
  basis.set_jpmax(jmax,pmax);
  
#endif
  cout << "setup equation.." << endl;
  SturmEquation<Basis> eq(T, basis);
  cout << "end setup equation!" << endl;
  
#if 0
  for(int p=0; p<100; p++){
      Index lambda4(p,basis.j0(),1,3,&basis);
      Index lambda5(0,4,1,5,&basis);
      cout <<eq.a(lambda4,lambda5) <<endl;;
  }
  
  /* plot solution*/ 
  InfiniteVector<double,Index> indplot;
  const char* filenameInd = "indexplot.m";  
////  u.scale(&eq, -1);
//  //u.set_coefficient(Index(1,3,0,1,&basis), u[Index(1,3,0,1,&basis)]+1);
  indplot.set_coefficient(Index(0,6,1,10,&basis), 1);
////  basis.reconstruct(w,4,testv);
////  cout << "reconstruction: " << testv << endl;
////  cout << "scaled solution: " << endl << u << endl;
  SampledMapping<1> indsamp(evaluate(eq.basis(), indplot, true, 7));
//  SampledMapping<1> s(evaluate(eq.basis(), Index(1,3,0,1,&basis), true, 7));
  std::ofstream ind_stream1(filenameInd);
  indsamp.matlab_output(ind_stream1);
  ind_stream1 << "figure;\nplot(x,y);"
            << "title('Indexplot');" << endl;
  ind_stream1.close();
  
  abort();
  
  int mysize=PMAX;
  Vector<double> Squarenorm;
  Squarenorm.resize(mysize);
  for(int p=0; p<mysize; p++){
      Index lambdag(p,basis.j0(),0,3,&basis);
      Squarenorm[p]=eq.a(lambdag,lambdag);
  }
  
  cout << "- writing Squarenorm to the file squarenorm.m ..." << endl;
   std::ofstream sqstream("squarenorm.m");
   sqstream << "sqnorm=";
   print_vector(Squarenorm, sqstream);
   sqstream << ";" << endl;
   sqstream.close();
   cout << "  ... done!" << endl;
#endif
  
//  abort();
 

#ifdef NONADAPTIVE
  //nonadaptive setting
  
  const int j0=eq.basis().j0();
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
  for (int i=0; i<basis.degrees_of_freedom();i++) {
//  for (int i=0; i<928;i++) {    

//    if(((*(frame.get_quarklet(i))).k()[1]>1 && (*(frame.get_quarklet(i))).k()[1]<6) || (*(frame.get_quarklet(i))).p()[1]==0){
//    if(i<7 || (i>=15 && i<22)){
        Lambda.insert(*(basis.get_quarklet(i)));
//        cout << *(basis.get_quarklet(i)) << endl;
//    }
//    }
//    }
//    cout << discrete_poisson.a(*(frame.get_quarklet(i)),*(frame.get_quarklet(i)))<<endl;


  }
//  abort();
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
//    if (lambda == eq.basis().last_generator(j0)) break;
  }
#endif
  
  
  


    
    
    
  //gives the index set
//     cout << "- set up stiffness matrix with respect to the index set Lambda=" << endl;
//     for (set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it)
//       cout << *it << endl;


  cout << "- set up (preconditioned) stiffness matrix (j0=" << eq.basis().j0() << ",jmax=" << jmax << ",pmax=" << pmax << ")..." << endl;
  SparseMatrix<double> A;
  setup_stiffness_matrix(eq, Lambda, A);
//  cout << "A size: " << A.size() << endl;
  Vector<double> A_diagonal;
  A_diagonal.resize(Lambda.size());
  for(int i=0; i<Lambda.size();i++){
//      cout << i << endl;
      A_diagonal[i]=A.get_entry(i,i);
  }
  
  cout << "- writing diagonal of A to the file A_diagonal.m ..." << endl;
   std::ofstream diagstream("A_diagonal.m");
   diagstream << "A_diag=";
   print_vector(A_diagonal, diagstream);
   diagstream << ";" << endl;
   diagstream.close();
   cout << "  ... done!" << endl;
  


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
   
//   abort();

  
  b.compress(1e-14);
  

  Vector<double> x(Lambda.size()), x2(Lambda.size()), err(Lambda.size()), err2(Lambda.size()); x = 0, x2 = 0;
  
  unsigned int iterations;
  
//  cout << b << endl;
//CG(A, b, x, 1e-8, 10000, iterations);
//Landweber(A,b,x, 0.5, 1e-8, 10000, iterations);
  double mu=0.1;
  int potenz = 1;
#ifdef SHRINKAGE
  
  while (mu > 1e-10){
      
//    x=b;
      x=0;
      
#ifdef RICHARDSON
    Richardson_sparse(A,b,x, 0.5, mu, 1e-8, 1000, iterations);    
//    Richardson_sparse(A,b,x, 0.3, 0, 1e-8, 1000, iterations);
#endif
#ifdef CONGRAD
  //  CG(A, b, x, 1e-8, 10000, iterations);
#endif
  //  Richardson(A,b,x, 0.1, 1e-8, 5000, iterations);
  //  Landweber
    
  
#else
    bool myboolean = true;
    while(myboolean){ 
#ifdef RICHARDSON
        Richardson(A,b,x, 0.05, 1e-8, 20000, iterations);
#endif
#ifdef CONGRAD
       CG(A, b, x, 1e-8, 10000, iterations); 
#endif
        myboolean=false;
#endif
//cout << "Solution: " << endl << x << endl;
int nontrivial(0);
for(int i=0; i<x.size();i++){
    if(x[i]!=0)
        nontrivial++;
}
cout << "Number of nontrivial entries: " << nontrivial << endl;

cout << "- writing x to the file solution.m ..." << endl;
std::ofstream xstream("solution.m");
xstream << "x=";
print_vector(x, xstream);
xstream << ";" << endl;
xstream.close();
cout << "  ... done!" << endl;


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
     

    
    
    
    
    
#ifndef DELTADIS

    // compute some error-norms
    const double Linfty_error = linfty_norm(u_values-uexact_values);
    cout << "  L_infinity error on a subgrid: " << Linfty_error << endl;
    
    
    const double L2_error = sqrt(l2_norm_sqr(u_values-uexact_values)*h);
    cout << "  L_2 error on a subgrid: " << L2_error << endl;
#endif
    
#endif
  
  
//  cout << "- point values of the solution:" << endl;
  InfiniteVector<double,Index> u,v,w,testv;
  unsigned int i = 0;
  for (set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it, ++i)
    u.set_coefficient(*it, x[i]);
  
//  cout << "u vor shrinkage: " << u << endl;
//  u.shrinkage(0.1);
//  cout << "u nach shrinkage: " << u << endl;
    
  
  
//  u.COARSE(1e-6,v);
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
#if 1
  
  
  for(int p=0;p<=pmax;p++){
    char filenameCoefficients[128];
    sprintf(filenameCoefficients, "%s%d%s%d%s", "sturm_bvp_solution_coefficients_p_" , p , "_mu_", potenz , ".m");
    std::ofstream coeff_stream (filenameCoefficients);
  coeff_stream << "figure;" << endl;
  plot_indices(&basis, u, jmax, coeff_stream, p, "jet", true, true, -8);
  coeff_stream << "title('nonadaptive coefficients on the level p=" << p <<" of the test problem ("
                  << basis_type << ")');" << endl;
  coeff_stream.close();   
  }
#endif
#else
    const char* filenameCoefficients1 = "sturm_bvp_solution_coefficients.m";
    std::ofstream coeff_stream1;
    coeff_stream1.open(filenameCoefficients1);
    coeff_stream1 << "figure;" << endl;
    plot_indices(&basis, u, jmax, coeff_stream1, "jet", true, true, -8);
    coeff_stream1 << "title('Sturm bvp: solution coefficients of the test problem ("
                  << basis_type << ")');" << endl;
    coeff_stream1.close();
#endif
    
  
 /* plot solution*/ 
  char filenameSolution1[128];
  sprintf(filenameSolution1, "%s%d%s", "sturm_bvp_mu_", potenz , ".m");
//  const char* filenameSolution1 = "sturm_bvp_solution.m";  
  u.scale(&eq, -1);
  //u.set_coefficient(Index(1,3,0,1,&basis), u[Index(1,3,0,1,&basis)]+1);
//  w.set_coefficient(Index(0,3,1,3,&basis), 1);
//  basis.reconstruct(w,4,testv);
//  cout << "reconstruction: " << testv << endl;
//  cout << "scaled solution: " << endl << u << endl;
  SampledMapping<1> s(evaluate(eq.basis(), u, true, 7));
//  SampledMapping<1> s(evaluate(eq.basis(), Index(1,3,0,1,&basis), true, 7));
  std::ofstream u_stream1(filenameSolution1);
  s.matlab_output(u_stream1);
  u_stream1 << "figure;\nplot(x,y);"
            << "title('Sturm bvp: solution to test problem ("
            << basis_type << "), " << "pmax= " << pmax << ", jmax= " << jmax << ", d= " << d << "');" << endl;
  u_stream1.close();
  
  mu*=0.1;
  potenz++;
  }
  
#endif
  
#ifdef ADAPTIVE
  
  //adaptive setting (CDD2)
  
#ifdef FRAME  
  cout << "setup cached equation.." << endl;
//  CachedQuarkletProblem<SturmEquation<Basis> > ceq(&eq, 3.7, 5);
  CachedQuarkletProblem<SturmEquation<Basis> > ceq(&eq, 1., 1.);
  cout << "end setup cached equation.." << endl;
#endif

#ifdef BASIS
//  CachedProblem<SturmEquation<Basis> > ceq(&eq, 1.0, 1.0);
  CachedProblem<SturmEquation<Basis> > ceq(&eq);
  
  //calculate the exact expansion and testing APPLY
  InfiniteVector<double, Index> exact_solution, exact_solution_coarsed, apply_result, approx_apply_result, right_side;
//  expand(uexact, basis, 0, jmax, exact_solution);
//  exact_solution.scale(&ceq, 1);
//  APPLY(ceq, exact_solution, 1e-9, apply_result, jmax, CDD1);
//  ceq.RHS(1e-10,right_side);
//  cout.precision(17);
//  cout << "exact_solution: " << endl << exact_solution << endl;
//  cout << "A*(exact_solution): " << endl << apply_result << endl;
//  cout << "RHS: " << endl << right_side << endl;
  
#endif
  InfiniteVector<double, Index> F_eta;
  ceq.RHS(1e-6, F_eta);
  double epsilon = 1e-6;
  InfiniteVector<double,Index> u_epsilon;
  InfiniteVector<double,int> u_epsilon_int;
  clock_t tic = clock();

#ifdef FRAME
//  const double norminv = ceq.norm_Ainv();  
//  const double nu = norminv*l2_norm(F_eta); 
//  CDD2_QUARKLET_SOLVE(ceq, nu, epsilon, u_epsilon, jmax, DKR, pmax, 2, 2);
//  DUV_QUARKLET_SOLVE_SD(ceq, nu, epsilon, u_epsilon, CDD1, pmax, jmax, 2, 2);
//  steepest_descent_ks_QUARKLET_SOLVE(ceq, epsilon, u_epsilon, DKR, 2, 2);
#ifdef SD
  steepest_descent_ks_QUARKLET_SOLVE(ceq, epsilon, u_epsilon_int, DKR, 2, 2);
#endif
//  CDD2_QUARKLET_SOLVE(ceq, nu, epsilon, u_epsilon, jmax, DKR, pmax, 2, 2);
  
  
//  richardson_QUARKLET_SOLVE(ceq,epsilon,u_epsilon,DKR, 2, 2);
#ifdef RICHARDSON
#ifdef SHRINKAGE
  const double shrink = 0.01;
  richardson_QUARKLET_SOLVE(ceq,epsilon,u_epsilon_int,DKR, 1, 1, shrink);  
#else
  richardson_QUARKLET_SOLVE(ceq,epsilon,u_epsilon_int,DKR, 1, 1);
  //  CDD2_QUARKLET_SOLVE(ceq, nu, epsilon, u_epsilon, jmax, DKR, pmax, 2, 2);
#endif

#endif
for (typename InfiniteVector<double,int>::const_iterator it(u_epsilon_int.begin()),
 	   itend(u_epsilon_int.end()); it != itend; ++it){
        u_epsilon.set_coefficient(*(basis.get_quarklet(it.index())), *it);
    }
  
#endif
#ifdef BASIS
  CDD2_SOLVE(ceq, nu, epsilon, u_epsilon, jmax);
//  DUV_SOLVE_SD(ceq, nu, epsilon, u_epsilon, jmax);
//  Array1D<InfiniteVector<double, Index> > approximations(1);

//  steepest_descent_ks_SOLVE(ceq, epsilon, u_epsilon);
//  cout << "Solution:" << endl << u_epsilon << endl;
  
  
#endif
  


clock_t toc = clock();
double time = (double)(toc-tic);
cout << "\nTime taken: " << (time/CLOCKS_PER_SEC) << " s";



  
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
    
#ifdef BASIS
    cout << "||u-u_epsilon||_{l_2}: " << l2_norm(u_epsilon-exact_solution) << endl;
    cout << "||u-u_epsilon||_{l_infty}: " << linfty_norm(u_epsilon-exact_solution) << endl;
    APPLY(ceq, exact_solution, 1e-9, apply_result, jmax, CDD1);
    APPLY(ceq, u_epsilon, 1e-9, approx_apply_result, jmax, CDD1);
    exact_solution.COARSE(1e-6,exact_solution_coarsed);
    cout << "||Au-Au_epsilon||_{l_2}: " << l2_norm(apply_result-approx_apply_result) << endl;
    cout << "||Au-Au_epsilon||_{l_infty}: " << linfty_norm(apply_result-approx_apply_result) << endl;
//    cout << "Expansion of the exact solution" << endl << exact_solution_coarsed << endl;
//    cout << "Expansion of the approximate solution" << endl << u_epsilon << endl;
#endif
    
    
    
    
    
/* plot solution coefficients */
  // output for Frame has to be renewed
  
#ifdef FRAME

//  string filenameCoefficients2[4] = {"sturm_bvp_solution_coefficients_p_0_ad.m", "sturm_bvp_solution_coefficients_p_1_ad.m",
//  "sturm_bvp_solution_coefficients_p_2_ad.m", "sturm_bvp_solution_coefficients_p_3_ad.m"};
//  
////  string filenameCoefficients2[1] = {"sturm_bvp_solution_coefficients_p_0_ad.m"};
//  
//  for(int p=0;p<=pmax;p++){
//  const char* cstr = filenameCoefficients2[p].c_str();
//  cout << filenameCoefficients2[p] << endl;
//  std::ofstream coeff_stream2 (cstr);
//  coeff_stream2 << "figure;" << endl;
//  plot_indices(&basis, u_epsilon, jmax, coeff_stream2, p, "jet", true, true, -6);
//  coeff_stream2 << "title('adaptive coefficients on the level p=" << p <<" of the test problem ("
//                  << basis_type << ")');" << endl;
//  coeff_stream2.close();   
//  }
  
  for(int p=0;p<=pmax;p++){
    char filenameCoefficients2[128];
    sprintf(filenameCoefficients2, "%s%d%s", "sturm_bvp_solution_coefficients_p_" , p , "_ad.m");
    std::ofstream coeff_stream2 (filenameCoefficients2);
  coeff_stream2 << "figure;" << endl;
  plot_indices(&basis, u_epsilon, jmax, coeff_stream2, p, "jet", true, true, -8);
  coeff_stream2 << "title('adaptive coefficients on the level p=" << p <<" of the test problem ("
                  << basis_type << ")');" << endl;
  coeff_stream2.close();   
  }
  
  

#else
    const char* filenameCoefficients2 = "./sd_results33_basis/sturm_bvp_solution_coefficients_adaptive.m";
    std::ofstream coeff_stream2;
    coeff_stream2.open(filenameCoefficients2);
    coeff_stream2 << "figure;" << endl;
    plot_indices(&basis, u_epsilon, jmax, coeff_stream2, "jet", true, true, -6);
    coeff_stream2 << "title('Sturm bvp: adaptive solution coefficients of the test problem ("
                  << basis_type << ")');" << endl;
    coeff_stream2.close();
#endif    
 
//plot of the adaptive solution
  const char* filenameSolution2 = "sturm_bvp_solution_adaptive.m";  
  u_epsilon.scale(&ceq, -1);
#ifdef FRAME
  SampledMapping<1> s2(evaluate(ceq.frame(), u_epsilon, true, 7));
#else
  SampledMapping<1> s2(evaluate(ceq.basis(), u_epsilon, true, 7));
#endif
  std::ofstream u_stream2(filenameSolution2);
  s2.matlab_output(u_stream2);
  u_stream2 << "figure;\nplot(x,y);"
            << "title('Sturm bvp: adaptive solution to test problem ("
            << basis_type << "), " << "pmax= " << pmax << ", jmax= " << jmax << ", d= " << d << ", dT= " << dT << "');" << endl;
  
  u_stream2.close();

#endif


  return 0;
}
