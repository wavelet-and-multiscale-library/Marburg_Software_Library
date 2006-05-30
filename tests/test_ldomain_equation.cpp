#include <iostream>
#include <fstream>
#include <set>

#include <algebra/vector.h>
#include <algebra/infinite_vector.h>
#include <utils/function.h>
#include <utils/fixed_array1d.h>
#include <numerics/bvp.h>
#include <geometry/sampled_mapping.h>

#include <interval/ds_basis.h>
#include <interval/p_basis.h>
#include <Ldomain/ldomain_basis.h>
#include <galerkin/ldomain_equation.h>

using namespace std;
using namespace MathTL;
using namespace WaveletTL;

/*
  A test problem for the Poisson equation on the cube with homogeneous Dirichlet b.c.'s:
    -Delta u(x,y) = 2(x(1-x)+y(1-y))
  with exact solution
    u(x,y) = x(1-x)y(1-y)
*/
class myRHS
  : public Function<2,double>
{
public:
  virtual ~myRHS() {};
  double value(const Point<2>& p, const unsigned int component = 0) const {
    return 2*(p[0]*(1-p[0])+p[1]*(1-p[1]));
  }
  void vector_value(const Point<2>& p, Vector<double>& values) const {
    values[0] = value(p);
  }
};

class mySolution
  : public Function<2,double>
{
public:
  virtual ~mySolution() {};
  double value(const Point<2>& p, const unsigned int component = 0) const {
    return p[0]*(1-p[0])*p[1]*(1-p[1]);
  }
  void vector_value(const Point<2>& p, Vector<double>& values) const {
    values[0] = value(p);
  }
};

int main()
{
  cout << "Testing wavelet-Galerkin solution of an elliptic equation on the L-shaped domain ..." << endl;

  const int d  = 2;
  const int dT = 2;
//   typedef DSBasis<d,dT> Basis1D;
  typedef PBasis<d,dT> Basis1D;
  typedef LDomainBasis<Basis1D> LBasis;
  typedef LBasis::Index Index;

  ConstantFunction<2> constant_rhs(Vector<double>(1, "1.0"));
  PoissonBVP<2> poisson(&constant_rhs);

  LDomainEquation<Basis1D> eq(&poisson);


  return 0;
}
