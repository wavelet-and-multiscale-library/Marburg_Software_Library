#include <iostream>
#include <map>

#include <algebra/symmetric_matrix.h>
#include <numerics/sturm_bvp.h>

#include <interval/dku_index.h>
#include <interval/dku_basis.h>
#include <galerkin/sturm_equation.h>

using namespace std;
using namespace WaveletTL;

using MathTL::simpleSturmBVP;

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
  SturmEquation<Basis> eq(T);

  typedef Basis::Index Index;
  set<Index> Lambda;
  for (Index lambda = eq.basis().firstGenerator(eq.basis().j0());;++lambda)
    {
      Lambda.insert(lambda);
      if (lambda == eq.basis().lastWavelet(eq.basis().j0()))
	break;
    }

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

  cout << "- stiffness matrix A=" << endl << A << endl;
  
  return 0;
}
