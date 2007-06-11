//#define _WAVELETTL_GALERKINUTILS_VERBOSITY 2

// resolution of the grid
#define RESOLUTION 15
// should the data be saved to files?
//#define SAVE_DATA
#define SAVE_RESULTS

// choose which basis to use
#define BASIS_S 0
#define BASIS_P 1
#define BASIS_DS 2
#define BASIS_ADAPTED_S 3

#define BASIS BASIS_ADAPTED_S


#define JMAX_END 15
#define JMAX_STEP 1


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
#elif (BASIS == BASIS_ADAPTED_S)
 #define ADAPTED_BASIS
 #include <interval/s_basis.h>
 #include <interval/s_support.h>
 #include <interval/interval_evaluate.h>
 typedef WaveletTL::SBasis MultiBasis;
 #define BASIS_NAME "adapted_s"
#endif
#ifdef ADAPTED_BASIS
 #include <interval/adapted_basis.h>
 #include <interval/adapted_support.h>
 typedef WaveletTL::AdaptedBasis<MultiBasis> Basis;
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


//! computes the L_2 norm of a function given on a grid as an Array1D
double L2norm(const Array1D<double>& f)
{
  Array1D<double>::size_type i;
  double norm;

  // use added trapezuidal rule
  const double h = 1.0/(double)f.size();
  norm = quad(f[0])/2.0;
  for (i = 1; i < f.size()-1; i++) {
    norm += quad(f[i]);
  }
  norm += quad(f[f.size()-1])/2.0;
  norm *= h;

  return sqrt(norm);
}


//! differentiates a function known only on a grid, given as a vector
void differentiate(const Array1D<double>& f, Array1D<double>& derivative)
{
  assert(derivative.size() == f.size());
  const double h = 1.0/(double)f.size();
  Array1D<double>::size_type i;
  const Array1D<double>::size_type imax = f.size()-1;

  derivative[0] = (f[1]-f[0])/h;
  for (i = 1; i < imax; i++)
    derivative[i] = (f[i+1]-f[i-1])/(2*h);
  derivative[imax] = (f[imax]-f[imax-1])/h;

  return;
}



//! main routine
int main()
{
  cout << "Testing biharmonic equation ..." << endl;
  #if (BASIS == BASIS_S)
  Basis basis;
  #elif ((BASIS == BASIS_P) || (BASIS == BASIS_DS))
  Basis basis(2, 2);
  #elif defined(ADAPTED_BASIS)
  MultiBasis multi_basis;
  Basis basis(&multi_basis);
  #endif
  Basis::Index lambda;
  set<Basis::Index> Lambda;
  
  #ifdef SAVE_RESULTS
  std::ofstream fs;
  ostringstream filename;
  filename << "biharmonic_results_" << BASIS_NAME << ".dat";
  fs.open(filename.str().c_str()); // ("biharmonic-results.dat");
  #endif // SAVE_RESULTS
  for (int jmax = basis.j0(); jmax <= JMAX_END; jmax+=JMAX_STEP) {
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
    ostringstream matrix_filename;
    matrix_filename << BASIS_NAME << "_stiff_" << jmax;
    cout << "Saving stiffness matrix to " << matrix_filename << ".m" << endl;
    stiff.matlab_output(matrix_filename.str().c_str(), "A", 1);
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
    cout << "Evaluating solution on grid with resolution 2^{" << RESOLUTION << "} ... " << endl;
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
    double err_L2 = L2norm(err);
    cout << err_L2 << endl;
    cout << "H^2 error approximation: ";
    Array1D<double> d1(err.size());
    differentiate(err, d1);
    Array1D<double> d2(err.size());
    differentiate(d1,d2);
    double err_H2 = err_L2 + L2norm(d1) + L2norm(d2);
    cout << err_H2 << endl;
    #ifdef SAVE_RESULTS
    fs << jmax << "\t" << lambdamax/lambdamin << "\t" << linfty_norm(err) << "\t" << err_L2 << "\t" << err_H2 << endl;
    #endif // SAVE_RESULTS
  }
  #ifdef SAVE_RESULTS
  fs.close();
  #endif // SAVE_RESULTS

  return 0;
}
