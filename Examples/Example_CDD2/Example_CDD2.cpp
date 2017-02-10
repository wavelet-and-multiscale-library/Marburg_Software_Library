/**
 * Example_CDD2.cpp (version 0.9)
 *
 * Intended to be an implementation of an 
 * adaptive wavelet algorithm based on [CDD2] 
 * for Sturm boundary value problems.
 *
 * This software is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,
 * either expressed or implied.
 *
 * Literatur:
 *
 * [CDD2] A. Cohen, W. Dahmen, R. DeVore.
 * Adaptive wavelet methods. II: Beyond the elliptic case.
 * Found. Comput. Math. 2, No. 3, 203-245 (2002).
 *
 * Contact: AG Numerik, Philipps University Marburg
 * http://www.mathematik.uni-marburg.de/~numerik/
 */

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
#include <adaptive/cdd2.h>

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
 * Example CDD2: main()
 *
 *
 */
int main()
{
    cout << "Example [CDD2]" << endl
         << "\nComputation of adaptive wavelet-Galerkin solution of"
         << " Sturm boundary value problems based on [CDD2] ..." << endl;

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
    const int jmax = 12;
    basis.set_jmax(jmax);

    cout << "Wavelet basis setup done (maximal level = "<< jmax <<")."<< endl;



   
 /*
     * First test problem
     */
    double epsilon1 = 1e-6;  /* set tolerance of the algorithm */
    const char* filenameCoefficients1 = "Example_CDD2_Problem_1_solution_coefficients.m";
    const char* filenameSolution1 = "Example_CDD2_Problem_1_solution.m";

    TestProblem<1> testProblem1;

    cout << "\nSetup of first test problem:" << endl
         << testProblem1.toString() << endl;

    SturmEquation<Basis> problem1(testProblem1, basis);
    CachedProblem<SturmEquation<Basis> > cproblem1(&problem1);

    cout << "Done with setup of first test problem." << endl
         << "\nBegin to solve first test problem (up to tolerance " << epsilon1 << ") ..." << endl;

    InfiniteVector<double, Index> F1_eta;
    cproblem1.RHS(1e-6, F1_eta);
    const double nu1 = cproblem1.norm_Ainv() * l2_norm(F1_eta);
    cout << "nu = " << nu1 << endl;
    InfiniteVector<double, Index> solution1_epsilon;
    CDD2_SOLVE(cproblem1, nu1, epsilon1, solution1_epsilon, jmax);

    cout << "Done solving first test problem." << endl
         << "\nWrite results of first test problem to disk ..." << endl;

    /* plot solution coefficients */
    std::ofstream coeff_stream1;
    coeff_stream1.open(filenameCoefficients1);
    coeff_stream1 << "figure;" << endl;
    plot_indices(&basis, solution1_epsilon, jmax, coeff_stream1, "jet", false, true, -8);
    coeff_stream1 << "title('Example [CDD2]: solution coefficients of the first test problem ("
                  << basis_type << " basis)');" << endl;
    coeff_stream1.close();

    /* plot solution graph */
    solution1_epsilon.scale(&cproblem1, -1); /* scaling because ... */
    SampledMapping<1> sm1(evaluate(cproblem1.basis(), solution1_epsilon, true, 2*jmax)); // Increase last parameter, if "Assertion `resolution >= 0' failed."
    std::ofstream u_stream1(filenameSolution1);
    sm1.matlab_output(u_stream1);
    u_stream1 << "figure;\nplot(x,y);"
              << "title('Example [CDD2]: solution to first test problem ("
              << basis_type << " basis)');" << endl;
    u_stream1.close();

    cout << "Output of first test problem written to files: " << endl
         << filenameCoefficients1 << endl
         << filenameSolution1 << endl;



    /*
     * Second test problem
     */
    double epsilon2 = 1e-6;  /* set tolerance of the algorithm */
    const char* filenameCoefficients2 = "Example_CDD2_Problem_2_solution_coefficients.m";
    const char* filenameSolution2 = "Example_CDD2_Problem_2_solution.m";

    TestProblem<2> testProblem2;

    cout << "\nSetup second test problem:" << endl
         << testProblem2.toString() << endl;

    SturmEquation<Basis> problem2(testProblem2, basis);
    CachedProblem<SturmEquation<Basis> > cproblem2(&problem2);

    cout << "Done with setup of second test problem." << endl
         << "\nBegin to solve second test problem (up to tolerance " << epsilon2 << ") ..." << endl;

    InfiniteVector<double, Index> F2_eta;
    cproblem2.RHS(1e-6, F2_eta);
    const double nu2 = cproblem2.norm_Ainv() * l2_norm(F2_eta);
    cout << "nu = " << nu2 << endl;
    InfiniteVector<double, Index> solution2_epsilon;
    CDD2_SOLVE(cproblem2, nu2, epsilon2, solution2_epsilon, jmax);

    cout << "Done solving second test problem." << endl
         << "\nWrite results of second test problem to disk ..." << endl;

    /* plot solution coefficients */
    std::ofstream coeff_stream2;
    coeff_stream2.open(filenameCoefficients2);
    coeff_stream2 << "figure;" << endl;
    plot_indices(&basis, solution2_epsilon, jmax, coeff_stream2, "jet", false, true, -8);
    coeff_stream2 << "title('Example [CDD2]: solution coefficients of the second test problem ("
                  << basis_type << " basis)');" << endl;
    coeff_stream2.close();

    /* plot solution graph */
    solution2_epsilon.scale(&cproblem2, -1); /* scaling because ... */
    SampledMapping<1> sm2(evaluate(cproblem2.basis(), solution2_epsilon, true, 2*jmax)); //" Increase last parameter, if Assertion `resolution >= 0' failed."
    std::ofstream u_stream2(filenameSolution2);
    sm2.matlab_output(u_stream2);
    u_stream2 << "figure;\nplot(x,y);"
              << "title('Example [CDD2]: solution to second test problem ("
                  << basis_type << " basis)');" << endl;
    u_stream2.close();

    cout << "Output of second test problem written to files: " << endl
         << filenameCoefficients2 << endl
         << filenameSolution2 << endl;



    /*
     * Third test problem
     */
    double epsilon3 = 1e-6;  /* set tolerance of the algorithm */
    const char* filenameCoefficients3 = "Example_CDD2_Problem_3_solution_coefficients.m";
    const char* filenameSolution3 = "Example_CDD2_Problem_3_solution.m";

    TestProblem<3> testProblem3;

    cout << "\nSetup third test problem:" << endl
         << testProblem3.toString() << endl;

    SturmEquation<Basis> problem3(testProblem3, basis);
    CachedProblem<SturmEquation<Basis> > cproblem3(&problem3);

    cout << "Done with setup of third test problem." << endl
         << "\nBegin to solve third test problem (up to tolerance " << epsilon3 << ") ..." << endl;

    InfiniteVector<double, Index> F3_eta;
    cproblem3.RHS(1e-6, F3_eta);
    const double nu3 = cproblem3.norm_Ainv() * l2_norm(F3_eta);
    cout << "nu = " << nu3 << endl;
    InfiniteVector<double, Index> solution3_epsilon;
    CDD2_SOLVE(cproblem3, nu3, epsilon3, solution3_epsilon, jmax);

    cout << "Done solving third test problem." << endl
         << "\nWrite results of third test problem to disk ..." << endl;

    /* plot solution coefficients */
    std::ofstream coeff_stream3;
    coeff_stream3.open(filenameCoefficients3);
    coeff_stream3 << "figure;" << endl;
    plot_indices(&basis, solution3_epsilon, jmax, coeff_stream3, "jet", false, true, -8);
    coeff_stream3 << "title('Example [CDD2]: solution coefficients of the third test problem ("
                  << basis_type << " basis)');" << endl;
    coeff_stream3.close();

    /* plot solution graph */
    solution3_epsilon.scale(&cproblem3, -1); /* scaling because ... */
    SampledMapping<1> sm3(evaluate(cproblem3.basis(), solution3_epsilon, true, 2*jmax)); //" Increase last parameter, if Assertion `resolution >= 0' failed."
    std::ofstream u_stream3(filenameSolution3);
    sm3.matlab_output(u_stream3);
    u_stream3 << "figure;\nplot(x,y);"
              << "title('Example [CDD2]: solution to third test problem ("
              << basis_type << " basis)');" << endl;
    u_stream3.close();

    cout << "Output of third test problem written to files: " << endl
         << filenameCoefficients3 << endl
         << filenameSolution3 << endl;


    /**
     * Contact: AG Numerik, Philipps University Marburg
     * http://www.mathematik.uni-marburg.de/~numerik/
     */
    cout << "\nContact: AG Numerik, Philipps University Marburg" << endl
         << "http://www.mathematik.uni-marburg.de/~numerik/" << endl
         << "\nEnd of Example [CDD2]." << endl;

    return EXIT_SUCCESS;
}
