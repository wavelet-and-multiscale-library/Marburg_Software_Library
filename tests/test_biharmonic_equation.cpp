//#define _WAVELETTL_GALERKINUTILS_VERBOSITY 2

// resolution of the grid
#define RESOLUTION 7
// should the data be saved to files?
//#define SAVE_DATA

// choose which basis to use
#define BASIS_S 0
#define BASIS_P 1
#define BASIS_DS 2

#define BASIS BASIS_S


#define JMAX_END 15


#include <iostream>

#include <algebra/infinite_vector.h>
#include <algebra/vector_norms.h>
#include <algebra/matrix.h>
#include <algebra/sparse_matrix.h>
#include <algebra/matrix_norms.h>
#include <utils/random.h>
#include <ctime>

#if (BASIS == BASIS_S)
 #include <interval/s_basis.h>
 #include <interval/s_support.h>
 #include <interval/interval_evaluate.h>
 typedef WaveletTL::SBasis Basis;
 #define BASIS_NAME "s"
#elif (BASIS == BASIS_P)
 #include <interval/p_basis.h>
 #include <interval/p_support.h>
 #include <interval/p_evaluate.h>
 typedef WaveletTL::PBasis<4, 4> Basis;
 #define BASIS_NAME "p"
#elif (BASIS == BASIS_DS)
 #include <interval/ds_basis.h>
 #include <interval/ds_support.h>
 #include <interval/ds_evaluate.h>
 typedef WaveletTL::DSBasis<4, 4> Basis;
 #define BASIS_NAME "ds"
#endif

#include <numerics/bvp.h>
#include <galerkin/galerkin_utils.h>
#include <galerkin/biharmonic_equation.h>
#include <numerics/eigenvalues.h>
#include <numerics/iteratsolv.h>

using namespace std;
using namespace WaveletTL;
using namespace MathTL;


//! computes the square of x
double quad(double x)
{
  return x*x;
}


int main()
{
  cout << "Testing biharmonic equation ..." << endl;
  #if (BASIS == BASIS_S)
  Basis basis;
  #else
  Basis basis(2, 2);
  #endif
  Basis::Index lambda;
  set<Basis::Index> Lambda;
  
  for (int jmax = basis.j0(); jmax <= JMAX_END; jmax+=3) {
    cout << "jmax = " << jmax << endl;
    Lambda.clear();
    for (lambda = basis.first_generator(basis.j0()); lambda <= basis.last_wavelet(jmax); ++lambda)
      Lambda.insert(lambda);

    Vector<double> value(1);
    value[0] = 384;
    ConstantFunction<1> const_fkt(value);
//    BiharmonicBVP<1> biharmonic(&const_fkt);
    BiharmonicEquation1D<Basis> discrete_biharmonic(basis, &const_fkt);

    cout << "Setting up full stiffness matrix ..." << endl;
    SparseMatrix<double> stiff(Lambda.size());

    setup_stiffness_matrix(discrete_biharmonic, Lambda, stiff, true);
    #ifdef SAVE_DATA
    ostringstream filename;
    filename << BASIS_NAME << "_stiff_" << jmax;
    cout << "Saving stiffness matrix to " << filename << ".m" << endl;
    stiff.matlab_output(filename.str().c_str(), "A", 1);
    #endif
    double lambdamin, lambdamax;
    unsigned int iterations;
    LanczosIteration(stiff, 1e-7, lambdamin, lambdamax, 200, iterations);
    cout << "Extremal Eigenvalues of stiffness matrix: " << lambdamin << ", " << lambdamax << "; " << iterations << " iterations needed." << endl;
    cout << "Condition number of stiffness matrix: " << lambdamax/lambdamin << endl;

    cout << "Setting up full right-hand side ..." << endl;
    Vector<double> rh;
    setup_righthand_side(discrete_biharmonic, Lambda, rh);

    cout << "Computing solution ..." << endl;
    Vector<double> x(rh);
    CG(stiff, rh, x, 1e-15, 200, iterations);
    cout << "Solved linear equation system in " << iterations << " iterations. " << endl;
    InfiniteVector<double,Basis::Index> u;
    unsigned int i = 0;
    for (set<Basis::Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it, ++i)
      u.set_coefficient(*it, x[i]);
    u.scale(&discrete_biharmonic, -1); // undo preconditioning
    SampledMapping<1> s(evaluate(basis, u, true, RESOLUTION));

    #ifdef SAVE_DATA
    cout << "Writing solution to file biharmonic-solution.dat ..." << endl;
    std::ofstream fs("biharmonic-solution.dat");
    s.gnuplot_output(fs);
    fs.close();
    #endif

    cout << "Maximal difference to exact solution: ";
    Polynomial<double> solution(Vector<double>(5, "0 0 16 -32 16")); // 16 x^2 (1-x)^2 = 16 (x^2 - 2 x^3 + x^4)
    SampledMapping<1> difference(Grid<1>(s.points()), solution);
    difference.add(-1.0,s);
    Array1D<double> err(difference.values());
    cout << linfty_norm(err) << endl;
    cout << "L_2 error approximation: ";
    // use added trapezuidal rule
    double h = 1.0/(double)err.size();
    double err_L2 = quad(err[0])/2.0;
    for (i = 1; i < err.size()-1; i++) {
      err_L2 += quad(err[i]);
    }
    err_L2 += quad(err[err.size()-1])/2.0;
    err_L2 *= h;
    cout << err_L2 << endl;
  }

  return 0;
}
