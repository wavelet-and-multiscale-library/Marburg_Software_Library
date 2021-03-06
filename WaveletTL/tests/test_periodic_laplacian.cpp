/* 
 * File:   test_periodic_basis.cpp
 * Author: keding
 *
 * Created on July 6, 2016, 2:33 PM
 */
#define PARALLEL 0
#define PRIMALORDER 3
#define DUALORDER 3
#define JMAX 8
#define PMAX 0
#define _DIM 1

#include <cstdlib>
#define _WAVELETTL_GALERKINUTILS_VERBOSITY 0
#include <iostream>
#include <fstream>
#include <utils/array1d.h>
#include <utils/function.h>
#include <Rd/r_basis.h>
#include <Rd/r_frame.h>
#include <Rd/quarklet_frame.h>
#include <Rd/cdf_basis.h>
#include <interval/periodic.h>
#include <interval/periodic_frame.h>
#include <interval/i_indexplot.h>
#include <numerics/periodicgr.h>
#undef PERIODIC_CDFBASIS
#ifndef PERIODIC_CDFBASIS
#define PERIODIC_FRAME
#endif

#define ADAPTIVE
#undef NONADAPTIVE
#undef  DELTADIS

#ifdef PERIODIC_CDFBASIS
#include <galerkin/cached_problem.h>
#include <galerkin/periodic_laplacian.h>
#endif
#include <galerkin/galerkin_utils.h>
#include <galerkin/Periodic_TestProblem.h>

#ifdef PERIODIC_FRAME
#include <galerkin/cached_quarklet_problem.h>
#include <galerkin/periodic_frame_laplacian.h>
#endif
#include <numerics/eigenvalues.h>
#include <adaptive/compression.h>
#include <adaptive/cdd2.h>
#include <adaptive/apply.h>
#include <galerkin/TestFunctions.h>
#include <interval/i_q_indexplot.h>


using namespace std;
using namespace WaveletTL;

/*testcase -u''=f, u(0)=u(1), \int_0^1 u=0
    different laplacian test problems with periodic b.c.'s -u''=f, u(0)=u(1), \int_0^1 u=0
    1: u(t)=sin(2*M_PI*t)
    2: u(t)=cos(2*M_PI*t)
    3: u(t)= f(t)=9*M_PI*M_PI*sin(3*M_PI*t)-6*M_PI 
    4: u(t)= f(t)=t-0.5 */

int main(int argc, char** argv) {
    
    const int d  = PRIMALORDER;
    const int dT = DUALORDER;
    const int jmax = JMAX;
    const bool normalization = 1;
    const int pmax = PMAX;

    
    
    
     
    const unsigned int testcase=4;
    PeriodicTestProblem2<testcase> tper;
    Function<1>* uexact = 0;
    switch(testcase) {
        case 1:
            uexact = new Function8();
            break;
        case 2:
            uexact = new Function3();
            break;
        case 3:
            uexact = new Function7();
            break;
        case 4:
            uexact = new Function6();
            break;
        default:
            break;
    }
    
    
    
#ifdef PERIODIC_CDFBASIS
    typedef CDFBasis<d,dT> RBasis;
    typedef PeriodicBasis<RBasis> Basis;
    typedef Basis::Index Index;
    Basis basis;
    basis.set_jmax(jmax);
    typedef PeriodicIntervalLaplacian<RBasis> Problem;
    Problem L(tper, basis);
    typedef CachedProblem<Problem> CPROBLEM;
    CPROBLEM cachedL(&L);
#endif
    
#ifdef PERIODIC_FRAME
    
    typedef QuarkletFrame<d,dT> RBasis;
    typedef PeriodicFrame<RBasis> Basis;
    typedef Basis::Index Index;
    Basis basis;
    basis.set_jpmax(jmax, pmax);
    typedef PeriodicFrameIntervalLaplacian<RBasis> Problem;
    Problem L(tper, basis);
    typedef CachedQuarkletProblem<Problem> CPROBLEM;
    CPROBLEM cachedL(&L);
#endif
  
    
      
    // evaluate exact solution
    const unsigned int N = 100;
    const double h = 1./N;
    Vector<double> uexact_values(N+1);
    for (unsigned int i = 0; i <= N; i++) {
      const double point = i*h;
      uexact_values[i] = uexact->value(Point<1>(point));
      
    }
    
   
#if 0
    //Plot of wavelet with coefficient mu
    Index mu(3,0,7);
    SampledMapping<1> sm1(basis.evaluate(mu, 8, normalization));
    std::ofstream u_stream1("plot_periodicwavelet.m");
    sm1.matlab_output(u_stream1);
    u_stream1 << "figure;\nplot(x,y);"
            << "title('periodic wavelet/generator with index: " << mu << "');" << endl;
    u_stream1.close();
    cout << "Integralwert: " << 1-basis.rintegrate(mu)*basis.rintegrate(mu) << endl;
    
#endif
#ifdef NONADAPTIVE
   //nonadaptive setting

    const int j0= L.basis().j0();
    set<Index> Lambda;
    Vector<double> x;
    
  for (int j = j0; j <= jmax; j++) {
        Lambda.clear();
        
    cout << "  j=" << j << ":" << endl;
    //Implementation of the index set
    
#ifdef PERIODIC_FRAME
    int p=0;
    for (Index lambda = L.basis().first_generator(j0,0);;) {
	Lambda.insert(lambda);
	if (lambda == L.basis().last_wavelet(j,pmax)) break;
        //if (i==7) break;
        if (lambda == L.basis().last_wavelet(j,p)){
            ++p;
            lambda = L.basis().first_generator(j0,p);
        }
        else
            ++lambda;
      }
#else
    
    for (Index lambda = L.basis().first_generator(j0);; ++lambda) {
        Lambda.insert(lambda);
        if (lambda == L.basis().last_wavelet(j)) break;    
    }
#endif
    SparseMatrix<double> A;
    Matrix<double> evecs;
    setup_stiffness_matrix(L, Lambda, A);
//    cout << "A: " << endl << A << endl;
    Vector<double> evals;
//    SymmEigenvalues(A,evals,evecs);
//    cout<< "Eigenwerte A: " << evals << endl;
    
    
    Vector<double> b;
    setup_righthand_side(L, Lambda, b);
//    cout << "- right hand side: " << b << endl << endl;
    
    
    x.resize(Lambda.size()); x = 0;
    Vector<double> residual(Lambda.size()); 
    unsigned int iterations;
    
    CG(A, b, x, 1e-8, 5000, iterations);
    cout << "  Galerkin system solved with residual (infinity) norm ";
    A.apply(x, residual);
    residual -= b;
    cout << linfty_norm(residual)
	 << " (" << iterations << " iterations needed)" << endl;
    
    // evaluate approximate solution on a grid
    
    Vector<double> u_values(N+1);
    for (unsigned int i = 0; i <= N; i++) {
      const double point = i*h;
      int id = 0;
      for (set<Index>::const_iterator it(Lambda.begin()); it != Lambda.end(); ++it, ++id) 
	u_values[i] += x[id] * basis.evaluate(0, *it, point, normalization)*1./L.D(*it);
    }
//     cout << "  point values of Galerkin solution: " << u_values << endl;

    
//     cout << "  point values of exact solution: " << uexact_values << endl;

    // compute some errors
    const double Linfty_error = linfty_norm(u_values-uexact_values);
    cout << "  L_infinity error on a subgrid: " << Linfty_error << endl;
    
    
    const double L2_error = sqrt(l2_norm_sqr(u_values-uexact_values)*h);
    cout << "  L_2 error on a subgrid: " << L2_error << endl << endl;
  }
    
    //Solution plot
    InfiniteVector<double,Index> u;
        unsigned int i = 0;
    for (set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it, ++i)
        u.set_coefficient(*it, x[i]);
//   cout << "Solution: " << endl << u << endl;    
        
    
    u.scale(&cachedL,-1);
    SampledMapping<1> sm2(basis.evaluate(u, 8, normalization));
    std::ofstream u_stream2("plot_periodicsolution.m");
    sm2.matlab_output(u_stream2);
    u_stream2 << "figure;\nplot(x,y);"
            << "title('solution of the laplacian, pmax= " << pmax << ", jmax=" << jmax << ", d=" << d << ", dT=" << dT << "')" << endl;
    u_stream2.close();
#endif
    
#ifdef ADAPTIVE
    
    //adaptive setting
    
    InfiniteVector<double, Index> F_eta;
    cachedL.RHS(1e-6, F_eta);
    const double nu = cachedL.norm_Ainv() * l2_norm(F_eta);
    InfiniteVector<double, Index> u_epsilon;
    const double epsilon= 1e-6;
#ifdef PERIODIC_CDFBASIS 
    CDD2_SOLVE(cachedL, nu,  epsilon,  u_epsilon, jmax);
#else
    CDD2_QUARKLET_SOLVE(cachedL, nu,  epsilon,  u_epsilon, jmax, DKR, pmax, 2, 2);
#endif

    #if 1
  
  
//  string filenameCoefficients2[1] = {"sturm_bvp_solution_coefficients_p_0_ad.m"};
  
  for(int p=0;p<=pmax;p++){
    char filenameCoefficients2[128];
    sprintf(filenameCoefficients2, "%s%d%s", "Matlab_outputs/laplacian_solution_coefficients_p_" , p , "_ad.m");
    std::ofstream coeff_stream2 (filenameCoefficients2);
  coeff_stream2 << "figure;" << endl;
  plot_indices_iq2(&basis, u_epsilon, jmax, coeff_stream2, p, "jet", false, true, -8);
  coeff_stream2 << "title('adaptive coefficients on the level p=" << p <<" of the test problem ("
                  <<  "periodic frame)');" << endl;
  coeff_stream2.close();   
  }
#endif
    
    
     u_epsilon.scale(&cachedL,-1);
    SampledMapping<1> sm3(basis.evaluate(u_epsilon, 8, normalization));
    std::ofstream u_stream3("Matlab_outputs/plot_adaptiveperiodicsolution.m");
    sm3.matlab_output(u_stream3);
    u_stream3 << "figure;\nplot(x,y);"
            << "title('adaptive solution of the laplacian, pmax= " << pmax << ", jmax=" << jmax << ", d=" << d << ", dT=" << dT << "')" << endl;
    u_stream3.close();
#endif

    
 
    

    
}

