#include <iostream>
#include <fstream>
#include <utils/array1d.h>
#include <Rd/r_basis.h>
#include <Rd/cdf_basis.h>
#include <interval/periodic.h>
#include <galerkin/periodic_gramian.h>
#include <galerkin/galerkin_utils.h>
#include <algebra/sparse_matrix.h>
#include <algebra/vector.h>

using namespace std;
using namespace WaveletTL;

//
// some 1-periodic test functions

// f(x)=1
class Function1 : public Function<1> {
public:
  inline double value(const Point<1>& p, const unsigned int component = 0) const {
    return 1.0;
  }
  
  void vector_value(const Point<1> &p, Vector<double>& values) const {
    values.resize(1, false);
    values[0] = value(p);
  }
};

// f(x)=x*(1-x)
class Function2 : public Function<1> {
public:
  inline double value(const Point<1>& p, const unsigned int component = 0) const {
    return p[0]*(1-p[0]);
  }
  
  void vector_value(const Point<1> &p, Vector<double>& values) const {
    values.resize(1, false);
    values[0] = value(p);
  }
};

template <unsigned int N>
class PeriodicTestProblem
  : public PeriodicLaplacianProblem
{
public:
  double g(const double t) const {
    switch(N) {
    case 1:
      //u(t) = 100*(t*t*t-1.5*t*t+0.5*t)
        //return 100*(t*t*t-1.5*t*t+0.5*t);
        return t*(1-t);
      //return 100*(-6*t+3);
      break;
    case 2:
      return sin(2*M_PI*t);
    case 3:
      return cos(2*M_PI*t);      
    default:
      return 0;
      break;
    }
  }
};

int main()
{
  cout << "Testing periodic Gramians ..." << endl;

  cout << "* a periodic CDF basis:" << endl;

  typedef CDFBasis<2,2> RBasis;
  typedef PeriodicBasis<RBasis> Basis;
  typedef Basis::Index Index;
  Basis basis;
  PeriodicTestProblem<1> tper;

  typedef PeriodicIntervalGramian<RBasis> Problem;
  Problem G(tper, basis);
  
  const unsigned int testcase=2;
  Function<1>* u = 0;

  switch(testcase) {
  case 1:
    u = new Function1();
    break;
  case 2:
    u = new Function2();
    break;
  default:
    break;
  }

  cout << "* compute approximate expansions of the test function for several levels..." << endl;
  const int jmin = basis.j0();
//   const int jmax = jmin;
  const int jmax = 4;
  Vector<double> js(jmax-jmin+1);
  Vector<double> Linfty_errors(jmax-jmin+1), L2_errors(jmax-jmin+1);

  Vector<double> rhs;
  InfiniteVector<double, Index> rhsinf;
  

  for (int j = jmin; j <= jmax; j++) {
    cout << "  j=" << j << ":" << endl;
    js[j-jmin] = j;

    // compute expansion coefficients of u in the dual basis
    basis.expand(u, true, j, rhs);
//    cout << "  inner products of u against all wavelets on level " << j << ":" << endl
//   	 << rhs << endl;
    
    basis.expand(u, false, j, rhsinf);
//    cout << "  inner products of u against all wavelets on level " << j << " (InfiniteVector-version):" << endl
//   	 << rhsinf << endl;
    
    // setup the set Lambda of active wavelets
    set<Index> Lambda;
    for (Index lambda = basis.first_generator(basis.j0());; ++lambda) {
      Lambda.insert(lambda);
      if (lambda == basis.last_wavelet(j)) break;
    }
//     cout << "  active wavelet set:" << endl;
//     for (set<Index>::const_iterator it(Lambda.begin()); it != Lambda.end(); ++it) {
//       cout << *it << endl;
//     }

    SparseMatrix<double> A;
    setup_stiffness_matrix(G, Lambda, A);
//     cout << "  stiffness matrix:" << endl
//  	 << A;
    
    // solve Gramian system
    Vector<double> uj(A.row_dimension()), residual(A.row_dimension(), false);
    InfiniteVector<double, Index> ujinf, residualinf;
    
    unsigned int iterations;
    CG(A, rhs, uj, 1e-15, 250, iterations);
    //CG(A, rhsinf, ujinf, 1e-15, 250, iterations);
//     cout << "  solution coefficients: " << uj << endl;
    cout << "  Galerkin system solved with residual (infinity) norm ";
    A.apply(uj, residual);
    residual -= rhs;
    cout << linfty_norm(residual)
	 << " (" << iterations << " iterations needed)" << endl;

    // evaluate approximate solution on a grid
    const unsigned int N = 100;
    const double h = 1./N;
    Vector<double> uj_values(N+1);
    for (unsigned int i = 0; i <= N; i++) {
      const double x = i*h;
      int id = 0;
      for (set<Index>::const_iterator it(Lambda.begin()); it != Lambda.end(); ++it, ++id) {
	uj_values[i] +=
	  uj[id] * basis.evaluate(0, *it, x);
      }
    }
//     cout << "  point values of Galerkin solution: " << uj_values << endl;

    // evaluate exact solution
    Vector<double> uexact_values(N+1);
    for (unsigned int i = 0; i <= N; i++) {
      const double x = i*h;
      uexact_values[i] = u->value(Point<1>(x));
    }
//     cout << "  point values of exact solution: " << uexact_values << endl;

    // compute some errors
    const double Linfty_error = linfty_norm(uj_values-uexact_values);
    cout << "  L_infinity error on a subgrid: " << Linfty_error << endl;
    Linfty_errors[j-jmin] = Linfty_error;
    
    const double L2_error = sqrt(l2_norm_sqr(uj_values-uexact_values)*h);
    cout << "  L_2 error on a subgrid: " << L2_error << endl;
    L2_errors[j-jmin] = L2_error;
    
//    cout << uj << endl;
//    cout << rhsinf << endl;
  }

  if (u) delete u;
  
  
  
  
  //solution1_epsilon.scale(&cproblem, -1);
  //cout << solution1_epsilon << endl;
  SampledMapping<1> sm1(basis.evaluate(rhsinf, 8));
  std::ofstream u_stream1("../../Desktop/plotthis.m");
  sm1.matlab_output(u_stream1);
  u_stream1 << "figure;\nplot(x,y);"
            << "title('mygraph');" << endl;
  u_stream1.close();

  return 0;
}