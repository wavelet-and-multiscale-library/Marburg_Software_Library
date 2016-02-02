#include <iostream>
#include <map>
#include <time.h>

#include <algebra/symmetric_matrix.h>
#include <algebra/sparse_matrix.h>
#include <numerics/sturm_bvp.h>
#include <numerics/iteratsolv.h>

#include <interval/i_index.h>
#include <interval/i_indexplot.h>
#include <interval/ds_basis.h>
#include <galerkin/sturm_equation.h>
#include <galerkin/cached_problem.h>

#include <adaptive/cdd2.h>

using namespace std;
using namespace WaveletTL;

using MathTL::SimpleSturmBVP;
using MathTL::CG;

/*
  different test problems with homogeneous Dirichlet b.c.'s
  1: y(t)=x*(1-x), -y''(t)=2
 */
template <unsigned int N>
class TestProblem
  : public SimpleSturmBVP
{
public:
  double p(const double t) const {
    switch(N) {
    case 1:
      return 1;
      break;
    default:
      return 0;
      break;
    }
  }
  double p_prime(const double t) const {
    switch(N) {
    case 1:
      return 0;
      break;
    default:
      return 0;
      break;
    }
  }
  double q(const double t) const {
    switch(N) {
    case 1:
      return 0;
      break;
    default:
      return 0;
      break;
    }
  }
  double g(const double t) const {
    switch(N) {
    case 1:
      return 2;
      break;
    default:
      return 0;
      break;
    }
  }
  bool bc_left() const { return true; }
  bool bc_right() const { return true; }
};


int main()
{
  cout << "Testing adaptive wavelet-Galerkin solution of a Sturm b.v.p. with CDD2_SOLVE ..." << endl;

  TestProblem<1> T;

  const int d  = 3;
  const int dT = 3;
  typedef DSBasis<d,dT> Basis;
  typedef Basis::Index Index;
  Basis basis(true, true);
  int jmax = 8;
  basis.set_jmax(jmax);

  SturmEquation<Basis> problem(T,basis);
  CachedProblem<SturmEquation<Basis> > cproblem(&problem);

  InfiniteVector<double, Index> F_eta;
  cproblem.RHS(1e-6, F_eta);
  const double nu = cproblem.norm_Ainv() * l2_norm(F_eta);
  cout << "nu = " << nu << endl;
  InfiniteVector<double, Index> u_epsilon;
  CDD2_SOLVE(cproblem, nu, 1, u_epsilon, jmax);
  
  
  /* plot solution coefficients */
    std::ofstream coeff_stream1("../../Arbeitsfläche/plotcoeff.m");
#if 1
    coeff_stream1 << "figure;" << endl;
    plot_indices(&basis, u_epsilon, jmax, coeff_stream1, "jet", false, true, -8);
    coeff_stream1 << "title('Example [CDD2]: solution coefficients of test problem');" << endl;
    coeff_stream1.close();

    /* plot solution graph */
    u_epsilon.scale(&cproblem, -1); /* scaling because ... */
    SampledMapping<1> sm1(evaluate(cproblem.basis(), u_epsilon, true, 2*jmax)); // Increase last parameter, if "Assertion `resolution >= 0' failed."
    std::ofstream u_stream1("../../Arbeitsfläche/plotgraph.m");
    sm1.matlab_output(u_stream1);
    u_stream1 << "figure;\nplot(x,y);"
              << "title('Example [CDD2]: solution to test problem');" << endl;
    u_stream1.close();
#endif
  
  return 0;
}




