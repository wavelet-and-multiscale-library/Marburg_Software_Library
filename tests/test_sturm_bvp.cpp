#include <iostream>
#include <map>

#include <algebra/symmetric_matrix.h>
#include <numerics/sturm_bvp.h>
#include <numerics/iteratsolv.h>

#include <interval/dku_index.h>
#include <interval/dku_basis.h>
#include <galerkin/sturm_equation.h>

using namespace std;
using namespace WaveletTL;

using MathTL::simpleSturmBVP;
using MathTL::CG;

/*
  different test problems with homogeneous Dirichlet b.c.'s
  1: y(t)=x*(1-x), -y''(t)=2
 */
template <unsigned int N>
class TestProblem
  : public simpleSturmBVP
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
  cout << "Testing wavelet-Galerkin solution of a Sturm b.v.p. ..." << endl;

  TestProblem<1> T;

  const int d  = 2;
  const int dT = 4;
  typedef DKUBasis<d,dT> Basis;
  typedef Basis::Index Index;

  SturmEquation<Basis> eq(T);

  InfiniteVector<double, Index> coeffs;

#if 0
  coeffs[eq.basis().firstGenerator(eq.basis().j0())] = 1.0;
  coeffs[eq.basis().lastGenerator(eq.basis().j0())] = 2.0;
  coeffs[eq.basis().firstWavelet(eq.basis().j0())] = 3.0;
  coeffs[eq.basis().lastWavelet(eq.basis().j0())] = 4.0;
  coeffs[eq.basis().firstWavelet(eq.basis().j0()+1)] = 5.0;
  cout << "- a coefficient set:" << endl
       << coeffs << endl;
  eq.rescale(coeffs, -1);
  cout << "- after rescaling with D^{-1}:" << endl
       << coeffs << endl;
#endif

  eq.RHS(coeffs, 1e-2);
//   cout << "- approximate coefficient set of the right-hand side:" << endl
//        << coeffs << endl;

#if 1
  cout << "- check expansion of the right-hand side in the dual basis:" << endl;
  eq.rescale(coeffs, 1);
  eq.basis().evaluate(coeffs, false, 7).matlab_output(cout);
  eq.rescale(coeffs, -1);
#endif  

  set<Index> Lambda;
  for (InfiniteVector<double,Index>::const_iterator it = coeffs.begin(); it != coeffs.end(); ++it)
    Lambda.insert(it.index());

//   cout << "- set up stiffness matrix with respect to the index set Lambda=" << endl;
//   for (set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it)
//     cout << *it << endl;

#if 0
  cout << "- set up (preconditioned) stiffness matrix..." << endl;
  SymmetricMatrix<double> A(Lambda.size());

  unsigned int i = 0;
  for (set<Index>::const_iterator it1 = Lambda.begin(); it1 != Lambda.end(); ++it1, ++i)
    {
      set<Index>::const_iterator it2 = Lambda.begin();
      for (; *it2 != *it1; ++it2);

      unsigned int j = i;
      for (; it2 != Lambda.end(); ++it2, ++j)
	{
	  A.set_entry(i, j, eq.a(*it2, *it1) / eq.D(*it1) / eq.D(*it2));
	}
    }
  cout << "  ... done!" << endl;
//   cout << "- (preconditioned) stiffness matrix A=" << endl << A << endl;

  Vector<double> b(coeffs.size());
  i = 0;
  for (InfiniteVector<double,Index>::const_iterator it = coeffs.begin(); it != coeffs.end(); ++it, ++i)
    b[i] = *it;

  cout << "- right hand side: " << b << endl;

  Vector<double> x(coeffs.size()), err(coeffs.size()); x = 0;
  unsigned int iterations;
  CG(A, b, x, 1e-8, 100, iterations);
  
  cout << "- solution coefficients: " << x;
  cout << " with residual (infinity) norm ";
  A.apply(x, err);
  err -= b;
  cout << linfty_norm(err) << endl;
  
  cout << "- point values of the solution:" << endl;
  InfiniteVector<double,Index> u;
  i = 0;
  for (InfiniteVector<double,Index>::const_iterator it = coeffs.begin(); it != coeffs.end(); ++it, ++i)
    u.set_coefficient(it.index(), x[i]);
  
  eq.rescale(u, -1);
  SampledMapping<1> s(eq.basis().evaluate(u, true, 7));
  s.matlab_output(cout);
#endif

  return 0;
}
