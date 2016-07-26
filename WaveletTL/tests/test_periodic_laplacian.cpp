/* 
 * File:   test_periodic_basis.cpp
 * Author: keding
 *
 * Created on July 6, 2016, 2:33 PM
 */

#include <cstdlib>
#define _WAVELETTL_GALERKINUTILS_VERBOSITY 0
#include <iostream>
#include <fstream>
#include <utils/array1d.h>
#include <utils/function.h>
#include <Rd/r_basis.h>
#include <Rd/cdf_basis.h>
#include <interval/periodic.h>
#include <interval/i_indexplot.h>
#include <numerics/periodicgr.h>
#define PERIODIC_CDFBASIS

#include <galerkin/periodic_laplacian.h>
#include <galerkin/galerkin_utils.h>
#include <galerkin/Periodic_TestProblem.h>
#include <galerkin/cached_problem.h>
#include "TestFunctions.h"


using namespace std;
using namespace WaveletTL;

int main(int argc, char** argv) {
    
    const int d  = 3;
    const int dT = 3;
    const int jmax = 10;
    const bool normalization = 1;//choose 1 for Laplacian
    
    const unsigned int testcase=3;
    PeriodicTestProblem2<testcase> tper;
    Function<1>* uexact = 0;
    switch(testcase) {
        case 1:
            uexact = new Function2();
            break;
        case 2:
            uexact = new Function3();
//            break;
        case 3:
            uexact = new Function4();
//            break;
//        case 4:
//            uexact = new Function4();
//            break;
//        case 5:
//            uexact = new Hat();    
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
    const int j0= L.basis().j0();
    
#if 1
    //Plot of wavelet with coefficient mu
    Index mu(3,0,1);
    SampledMapping<1> sm1(basis.evaluate(mu, 8, normalization));
    std::ofstream u_stream1("plot_periodicwavelet.m");
    sm1.matlab_output(u_stream1);
    u_stream1 << "figure;\nplot(x,y);"
            << "title('periodic wavelet/generator with index: " << mu << "');" << endl;
    u_stream1.close();
#endif
    
    set<Index> Lambda;
    Vector<double> x;
    
  for (int j = j0; j <= jmax; j++) {
        Lambda.clear();
        
    cout << "  j=" << j << ":" << endl;
    //Implementation of the index set
    
    for (Index lambda = L.basis().first_generator(j0);; ++lambda) {
        Lambda.insert(lambda);
        if (lambda == L.basis().last_wavelet(j)) break;    
    }
    SparseMatrix<double> A;
    cout << "- set up stiffness matrix..." << endl;
    setup_stiffness_matrix(L, Lambda, A);
//    cout << "A: " << endl << A << endl;
    
    
    cout << "- set up right-hand side..." << endl;
    Vector<double> b;
    setup_righthand_side(L, Lambda, b);
//    cout << "- right hand side: " << b << endl << endl;
    
    
    x.resize(Lambda.size()); x = 0;
    Vector<double> residual(Lambda.size()); 
    unsigned int iterations;
    
    CG(A, b, x, 1e-15, 250, iterations);
    cout << "  Galerkin system solved with residual (infinity) norm ";
    A.apply(x, residual);
    residual -= b;
    cout << linfty_norm(residual)
	 << " (" << iterations << " iterations needed)" << endl;
    
    // evaluate approximate solution on a grid
    const unsigned int N = 100;
    const double h = 1./N;
    Vector<double> u_values(N+1);
    for (unsigned int i = 0; i <= N; i++) {
      const double point = i*h;
      int id = 0;
      for (set<Index>::const_iterator it(Lambda.begin()); it != Lambda.end(); ++it, ++id) 
	u_values[i] += x[id] * basis.evaluate(0, *it, point, normalization)*1./L.D(*it);
    }
//     cout << "  point values of Galerkin solution: " << u_values << endl;

    // evaluate exact solution
    Vector<double> uexact_values(N+1);
    for (unsigned int i = 0; i <= N; i++) {
      const double point = i*h;
      uexact_values[i] = uexact->value(Point<1>(point));
    }
//     cout << "  point values of exact solution: " << uexact_values << endl;

    // compute some errors
    const double Linfty_error = linfty_norm(u_values-uexact_values);
    cout << "  L_infinity error on a subgrid: " << Linfty_error << endl;
    
    
    const double L2_error = sqrt(l2_norm_sqr(u_values-uexact_values)*h);
    cout << "  L_2 error on a subgrid: " << L2_error << endl;
  }
    
    //Solution plot
    InfiniteVector<double,Index> u;
        unsigned int i = 0;
    for (set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it, ++i)
        u.set_coefficient(*it, x[i]);
//    cout << "Solution: " << endl << u << endl;    
        
    
    u.scale(&cachedL,-1);
    SampledMapping<1> sm2(basis.evaluate(u, 8, normalization));
    std::ofstream u_stream2("plot_periodicsolution.m");
    sm2.matlab_output(u_stream2);
    u_stream2 << "figure;\nplot(x,y);"
            << "title('solution of the laplacian')" << endl;
    u_stream2.close();
    
    
    
    
  
    
#endif
    

    
}

