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

#include <galerkin/periodic_gramian.h>
#include <galerkin/galerkin_utils.h>
#include <galerkin/Periodic_TestProblem.h>
#include <galerkin/TestFunctions.h>


using namespace std;
using namespace WaveletTL;

/*
  different test problems with periodic b.c.'s y(0)=y(1)
  1: y(t)=1;
  2: y(t) = t*(1-t)
  3: y(t)= cos(2*M_PI*t);
  4: y(t) = -sin(3*M_PI*t)
  5: y(t) = "hat function"
  6: y(t) = t; */

int main(int argc, char** argv) {
    
    const int d  = 2;
    const int dT = 2;
    const int jmax = 10;
    const bool normalization = 0;
    
    const unsigned int testcase=5;
    PeriodicTestProblem<testcase> tper;
    Function<1>* uexact = 0;
    switch(testcase) {
        case 1:
            uexact = new Function1();
            break;
        case 2:
            uexact = new Function2();
            break;
        case 3:
            uexact = new Function3();
            break;
        case 4:
            uexact = new Function4();
            break;
        case 5:
            uexact = new Hat(); 
            break;
        /*case 6:
            uexact = new Hat();
            break;*/
        default:
            break;
    }
    
    
    
#ifdef PERIODIC_CDFBASIS
    typedef CDFBasis<d,dT> RBasis;
    typedef PeriodicBasis<RBasis> Basis;
    typedef Basis::Index Index;
    Basis basis;
    basis.set_jmax(jmax);
    
    typedef PeriodicIntervalGramian<RBasis> Problem;
    Problem G(tper, basis);
    const int j0= G.basis().j0();
    
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
    
    for (Index lambda = G.basis().first_generator(j0);; ++lambda) {
        Lambda.insert(lambda);
        if (lambda == G.basis().last_wavelet(j)) break;    
    }
    SparseMatrix<double> A;
    Matrix<double> evecs;
    cout << "- set up stiffness matrix..." << endl;
    setup_stiffness_matrix(G, Lambda, A);
//    cout << "A: " << endl << A << endl;
    Vector<double> evals;
//    SymmEigenvalues(A,evals,evecs);
//    cout<< "Eigenwerte A: " << evals << endl;
    
    
    cout << "- set up right-hand side..." << endl;
    Vector<double> b;
    setup_righthand_side(G, Lambda, b);
//    cout << "- right hand side: " << b << endl << endl;
    
    
    x.resize(Lambda.size()); x = 0;
    Vector<double> residual(Lambda.size()); 
    unsigned int iterations;
    
    CG(A, b, x, 1e-15, 5000, iterations);
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
	u_values[i] += x[id] * basis.evaluate(0, *it, point, normalization);
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
        
    SampledMapping<1> sm2(basis.evaluate(u, 8, 0));
    std::ofstream u_stream2("plot_periodicsolution.m");
    sm2.matlab_output(u_stream2);
    u_stream2 << "figure;\nplot(x,y);"
            << "title('solution of the gramian')" << endl;
    u_stream2.close();
    
    
    
    
  
    
#endif
    

    
}

