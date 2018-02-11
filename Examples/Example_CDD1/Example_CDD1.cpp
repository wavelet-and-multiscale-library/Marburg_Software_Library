/**
 * Example_CDD1.cpp (version 1.0)
 *
 * Intended to be an implementation of an 
 * adaptive wavelet algorithm based on [CDD1] 
 * for Sturm boundary value problems.
 *
 * This software is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,
 * either expressed or implied.
 *
 * Literatur:
 *
 * [CDD1] A. Cohen, W. Dahmen, R. DeVore.
 * Adaptive wavelet methods for elliptic operator equations: Convergence rates.
 * Math. Comput. 70, No.233, 27-75 (2001).
 *
 * Contact: AG Numerik, Philipps University Marburg
 * http://www.mathematik.uni-marburg.de/~numerik/
 */

#define PARALLEL 0
#define _DIM 1

#include <iostream>

#include <interval/i_index.h>
#include <interval/i_indexplot.h>
#include <interval/spline_basis.h>
#include <interval/ds_basis.h>
#include <interval/p_basis.h>

#include <algebra/symmetric_matrix.h>
#include <algebra/sparse_matrix.h>
#include <numerics/iteratsolv.h>
#include <numerics/sturm_bvp.h>
#include <galerkin/sturm_equation.h>
#include <galerkin/cached_problem.h>
#include <./TestProblem.h>
#define _WAVELETTL_GALERKINUTILS_VERBOSITY 0
#define _WAVELETTL_CDD1_VERBOSITY 1
#include <adaptive/cdd1.h>

using namespace std;
using namespace WaveletTL;
using MathTL::SimpleSturmBVP;



/*
 * Select exactly one wavelet basis:
 * _SPLINE_BASIS
 * _DS_BASIS
 * _PRIMBS_BASIS
 */
#undef _DS_BASIS
#define _PRIMBS_BASIS



/**
 * Example CDD1: main()
 *
 *
 */
int main()
{
    cout << "Example [CDD1]" << endl
         << "\nComputation of adaptive wavelet-Galerkin solution of"
         << " Sturm boundary value problems based on [CDD1] ..." << endl;

    /*
     * Setup wavelet basis
     */
    const int d  = 3;    /* ... of wavelet basis */
    const int dT = 3;    /* vanishing moments of wavelet basis */

#ifdef _DS_BASIS
    cout << "\nSetup DS wavelet basis (d = "<< d <<", dT = "<< dT <<") ..." << endl;
    typedef DSBasis<d,dT> Basis;
    const char* basis_type = "DS";
#endif
#ifdef _PRIMBS_BASIS
    cout << "\nSetup Primbs wavelet basis (d = "<< d <<", dT = "<< dT <<") ..." << endl;
    typedef PBasis<d,dT> Basis;
    const char* basis_type = "Primbs";
#endif

    Basis basis(true, true);     /* create instance of wavelet basis */
    typedef Basis::Index Index;  /* (just an abbreviation) */

    /* set maximal available level within wavelet basis */
    const int jmax = 8;
    basis.set_jmax(jmax);

    cout << "Wavelet basis setup done (maximal level = "<< jmax <<")."<< endl;



    /*
     * First test problem
     */
    double epsilon1 = 1e-5;  /* set tolerance of the algorithm */
    const char* filenameConvergence1 = "Example_CDD1_Problem_1_convergence_logs.m";
    const char* filenameCoefficients1 = "Example_CDD1_Problem_1_solution_coefficients.m";
    const char* filenameSolution1 = "Example_CDD1_Problem_1_solution.m";

    TestProblem<1> testProblem1;

    cout << "\nSetup of first test problem:" << endl
         << testProblem1.toString() << endl;

    SturmEquation<Basis> problem1(testProblem1, basis);
    CachedProblem<SturmEquation<Basis> > cproblem1(&problem1);

    cout << "Done with setup of first test problem." << endl
         << "\nBegin to solve first test problem (up to tolerance " << epsilon1 << ") ..." << endl;

    MathTL::ConvergenceLogger logger1;
    InfiniteVector<double, Index> solution1_epsilon;
    CDD1_SOLVE(cproblem1, epsilon1, solution1_epsilon, logger1, jmax);

    cout << "Done solving first test problem." << endl
         << "\nWrite results of first test problem to disk ..." << endl;

    /* writing convergence logs to disk */
    std::ofstream conv_stream1(filenameConvergence1);
    conv_stream1 << "figure;" << endl;
    logger1.writeConvergencePlots(conv_stream1);
    conv_stream1.close();

    /* plot solution coefficients */
    std::ofstream coeff_stream1;
    coeff_stream1.open(filenameCoefficients1);
    coeff_stream1 << "figure;" << endl;
    plot_indices(&basis, solution1_epsilon, jmax, coeff_stream1, "jet", false, true, -8);
    coeff_stream1 << "title('Example [CDD1]: solution coefficients of the first test problem ("
                  << basis_type << " basis)');" << endl;
    coeff_stream1.close();

    /* plot solution graph */
    solution1_epsilon.scale(&cproblem1, -1); /* scaling because ... */
    SampledMapping<1> sm1(evaluate(cproblem1.basis(), solution1_epsilon, true, jmax+1)); //Increase last parameter, if "Assertion `resolution >= 0' failed."
    std::ofstream u_stream1(filenameSolution1);
    sm1.matlab_output(u_stream1);
    u_stream1 << "figure;\nplot(x,y);"
              << "title('Example [CDD1]: solution to first test problem ("
              << basis_type << " basis)');" << endl;
    u_stream1.close();

    cout << "Output of first test problem written to files: " << endl
         << filenameCoefficients1 << endl
         << filenameSolution1 << endl
         << filenameConvergence1 << endl;



    /*
     * Second test problem
     */
    double epsilon2 = 1e-5;  /* set tolerance of the algorithm */
    const char* filenameConvergence2 = "Example_CDD1_Problem_2_convergence_logs.m";
    const char* filenameCoefficients2 = "Example_CDD1_Problem_2_solution_coefficients.m";
    const char* filenameSolution2 = "Example_CDD1_Problem_2_solution.m";

    TestProblem<2> testProblem2;

    cout << "\nSetup second test problem:" << endl
         << testProblem2.toString() << endl;

    SturmEquation<Basis> problem2(testProblem2, basis);
    CachedProblem<SturmEquation<Basis> > cproblem2(&problem2);

    cout << "Done with setup of second test problem." << endl
         << "\nBegin to solve second test problem (up to tolerance " << epsilon2 << ") ..." << endl;

    MathTL::ConvergenceLogger logger2;
    InfiniteVector<double, Index> solution2_epsilon;
    CDD1_SOLVE(cproblem2, epsilon2, solution2_epsilon, logger2, jmax);

    cout << "Done solving second test problem." << endl
         << "\nWrite results of second test problem to disk ..." << endl;

    /* writing convergence logs to disk */
    std::ofstream conv_stream2(filenameConvergence2);
    conv_stream2 << "figure;" << endl;
    logger2.writeConvergencePlots(conv_stream2);
    conv_stream2.close();

    /* plot solution coefficients */
    std::ofstream coeff_stream2;
    coeff_stream2.open(filenameCoefficients2);
    coeff_stream2 << "figure;" << endl;
    plot_indices(&basis, solution2_epsilon, jmax, coeff_stream2, "jet", false, true, -8);
    coeff_stream2 << "title('Example [CDD1]: solution coefficients of the second test problem ("
                  << basis_type << " basis)');" << endl;
    coeff_stream2.close();

    /* plot solution graph */
    solution2_epsilon.scale(&cproblem2, -1); /* scaling because ... */
    SampledMapping<1> sm2(evaluate(cproblem2.basis(), solution2_epsilon, true, jmax+1)); //Increase last parameter, if "Assertion `resolution >= 0' failed."
    std::ofstream u_stream2(filenameSolution2);
    sm2.matlab_output(u_stream2);
    u_stream2 << "figure;\nplot(x,y);"
              << "title('Example [CDD1]: solution to second test problem ("
                  << basis_type << " basis)');" << endl;
    u_stream2.close();

    cout << "Output of second test problem written to files: " << endl
         << filenameCoefficients2 << endl
         << filenameSolution2 << endl
         << filenameConvergence2 << endl;



    /*
     * Third test problem
     */
    double epsilon3 = 1e-5;  /* set tolerance of the algorithm */
    const char* filenameConvergence3 = "Example_CDD1_Problem_3_convergence_logs.m";
    const char* filenameCoefficients3 = "Example_CDD1_Problem_3_solution_coefficients.m";
    const char* filenameSolution3 = "Example_CDD1_Problem_3_solution.m";

    TestProblem<3> testProblem3;

    cout << "\nSetup third test problem:" << endl
         << testProblem3.toString() << endl;

    SturmEquation<Basis> problem3(testProblem3, basis);
    CachedProblem<SturmEquation<Basis> > cproblem3(&problem3);

    cout << "Done with setup of third test problem." << endl
         << "\nBegin to solve third test problem (up to tolerance " << epsilon3 << ") ..." << endl;

    MathTL::ConvergenceLogger logger3;
    InfiniteVector<double, Index> solution3_epsilon;
    CDD1_SOLVE(cproblem3, epsilon3, solution3_epsilon, logger3, jmax);

    cout << "Done solving third test problem." << endl
         << "\nWrite results of third test problem to disk ..." << endl;

    /* writing convergence logs to disk */
    std::ofstream conv_stream3(filenameConvergence3);
    conv_stream3 << "figure;" << endl;
    logger3.writeConvergencePlots(conv_stream3);
    conv_stream3.close();

    /* plot solution coefficients */
    std::ofstream coeff_stream3;
    coeff_stream3.open(filenameCoefficients3);
    coeff_stream3 << "figure;" << endl;
    plot_indices(&basis, solution3_epsilon, jmax, coeff_stream3, "jet", false, true, -8);
    coeff_stream3 << "title('Example [CDD1]: solution coefficients of the third test problem ("
                  << basis_type << " basis)');" << endl;
    coeff_stream3.close();

    /* plot solution graph */
    solution3_epsilon.scale(&cproblem3, -1); /* scaling because ... */
    SampledMapping<1> sm3(evaluate(cproblem3.basis(), solution3_epsilon, true, jmax+1)); //Increase last parameter, if "Assertion `resolution >= 0' failed."
    std::ofstream u_stream3(filenameSolution3);
    sm3.matlab_output(u_stream3);
    u_stream3 << "figure;\nplot(x,y);"
              << "title('Example [CDD1]: solution to third test problem ("
              << basis_type << " basis)');" << endl;
    u_stream3.close();

    cout << "Output of third test problem written to files: " << endl
         << filenameCoefficients3 << endl
         << filenameSolution3 << endl
         << filenameConvergence3 << endl;


    /**
     * Contact: AG Numerik, Philipps University Marburg
     * http://www.mathematik.uni-marburg.de/~numerik/
     */
    cout << "\nContact: AG Numerik, Philipps University Marburg" << endl
         << "http://www.mathematik.uni-marburg.de/~numerik/" << endl
         << "\nEnd of Example [CDD1]." << endl;

    return EXIT_SUCCESS;
}
