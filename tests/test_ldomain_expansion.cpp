#include <iostream>
#include <set>

#include <utils/function.h>
#include <utils/fixed_array1d.h>
#include <algebra/sparse_matrix.h>
#include <numerics/iteratsolv.h>
#include <numerics/corner_singularity.h>
#include <Ldomain/ldomain_basis.h>
#include <Ldomain/ldomain_evaluate.h>
#include <Ldomain/ldomain_expansion.h>

#include "ldomain_solutions.h"

using namespace std;
using namespace MathTL;
using namespace WaveletTL;

class Polyx : public Function<2,double> {
public:
  virtual ~Polyx() {};
  double value(const Point<2>& p, const unsigned int component = 0) const {
    return p[0];
  }
  void vector_value(const Point<2>& p, Vector<double>& values) const {
    values[0] = value(p);
  }
};

class Polyy : public Function<2,double> {
public:
  virtual ~Polyy() {};
  double value(const Point<2>& p, const unsigned int component = 0) const {
    return p[1];
  }
  void vector_value(const Point<2>& p, Vector<double>& values) const {
    values[0] = value(p);
  }
};

class Polyxy : public Function<2,double> {
public:
  virtual ~Polyxy() {};
  double value(const Point<2>& p, const unsigned int component = 0) const {
    return p[0]*p[1];
  }
  void vector_value(const Point<2>& p, Vector<double>& values) const {
    values[0] = value(p);
  }
};

int main()
{
  cout << "Testing expansion routines for LDomainBasis ..." << endl;

  const int d  = 2;
  const int dT = 2;

#if 0
  typedef DSBasis<d,dT,BernsteinSVD> Basis1D;
  Basis1D basis1d;
#else
  typedef SplineBasis<d,dT,DS_construction> Basis1D;
  Basis1D basis1d("bio5-energy",0,0,0,0);
#endif

  typedef LDomainBasis<Basis1D> Basis;
  Basis basis(basis1d);

  typedef Basis::Index Index;

  ConstantFunction<2> f(Vector<double>(1, "1"));
//   Polyx f;
//   Polyy f;
//   Polyxy f;
//   PolyRHS f;
//   CubicHermiteInterpolant2D_td f(1, 0, 0, -1, 1);
//   CubicHermiteInterpolant2D_td f(2, 0, 0, -1, 1);
//   CornerSingularity f(Point<2>(0.,0.), 0.5, 1.5);

  const int jmax = basis.j0();
  
  InfiniteVector<double,Index> fcoeffs;
  expand(&f, basis, true, jmax, fcoeffs);
  fcoeffs.compress(1e-12);

  cout << "- integrals of f against the primal wavelets:" << endl
       << fcoeffs << endl;


  return 0;
}
