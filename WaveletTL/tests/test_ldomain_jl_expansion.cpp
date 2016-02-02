#include <iostream>
#include <set>

#include <utils/function.h>
#include <utils/fixed_array1d.h>
#include <algebra/sparse_matrix.h>
#include <numerics/iteratsolv.h>
#include <numerics/corner_singularity.h>
#include <Ldomain/ldomain_jl_basis.h>
#include <Ldomain/ldomain_jl_evaluate.h>
#include <Ldomain/ldomain_jl_expansion.h>

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
  cout << "Testing expansion routines for LDomainJLBasis ..." << endl;

  typedef LDomainJLBasis Basis;
  typedef Basis::Index Index;

  Basis basis;
  
//   ConstantFunction<2> f(Vector<double>(1, "1"));
//   Polyx f;
//   Polyy f;
//   Polyxy f;
//   PolyRHS f;
//   CubicHermiteInterpolant2D_td f(1, 0, 0, -1, 1);
  CubicHermiteInterpolant2D_td f(2, 0, 0, -1, 1);
//   CornerSingularity f(Point<2>(0.,0.), 0.5, 1.5);

  const int jmax = 2;
  
  InfiniteVector<double,Index> fcoeffs;
  expand(&f, basis, true, jmax, fcoeffs);

  cout << "- integrals of f against the primal wavelets:" << endl
       << fcoeffs << endl;

  cout << "- checking the integrals:" << endl;
  for (Index lambda = basis.first_generator(basis.j0());; ++lambda) {
    const int N = 5000;
    const double h = 1./N;

    // sum up all the point values:
    double s = 0;

    FixedArray1D<Array1D<double>,2> values;
    FixedArray1D<Array1D<double>,2> knots;
    knots[0].resize(N+1);
    knots[1].resize(N+1);
    Point<2> x;

    // patch Omega_0 = [-1,0]x[0,1]
    if (lambda.k()[0] <= 0 && lambda.k()[1] >= 0) {
      for (int k0 = 0; k0 < N; k0++)
	knots[0][k0] = -1.0+k0*h;
      evaluate(0, lambda.j(), lambda.e()[0], lambda.c()[0], lambda.k()[0], knots[0], values[0]);
      for (int k1 = 0; k1 < N; k1++)
	knots[1][k1] = k1*h;
      evaluate(0, lambda.j(), lambda.e()[1], lambda.c()[1], lambda.k()[1], knots[1], values[1]);
      for (int k0 = 0; k0 < N; k0++) {
 	x[0] = knots[0][k0];
 	for (int k1 = 1; k1 < N; k1++) {
 	  x[1] = knots[1][k1];
 	  s += values[0][k0] * values[1][k1] * f.value(x);
 	}
      }
    }
    // patch Omega_1 = [-1,0]x[-1,0]
    if (lambda.k()[0] <= 0 && lambda.k()[1] <= 0) {
      for (int k0 = 0; k0 < N; k0++)
	knots[0][k0] = -1.0+k0*h;
      evaluate(0, lambda.j(), lambda.e()[0], lambda.c()[0], lambda.k()[0], knots[0], values[0]);
      for (int k1 = 0; k1 < N; k1++)
	knots[1][k1] = -1.0+k1*h;
      evaluate(0, lambda.j(), lambda.e()[1], lambda.c()[1], lambda.k()[1], knots[1], values[1]);
      for (int k0 = 0; k0 < N; k0++) {
 	x[0] = knots[0][k0];
 	for (int k1 = 1; k1 < N; k1++) {
 	  x[1] = knots[1][k1];
 	  s += values[0][k0] * values[1][k1] * f.value(x);
 	}
      }
    }
    // patch Omega_2 = [0,1]x[-1,0]
    if (lambda.k()[0] >= 0 && lambda.k()[1] <= 0) {
      for (int k0 = 0; k0 < N; k0++)
	knots[0][k0] = k0*h;
      evaluate(0, lambda.j(), lambda.e()[0], lambda.c()[0], lambda.k()[0], knots[0], values[0]);
      for (int k1 = 0; k1 < N; k1++)
	knots[1][k1] = -1.0+k1*h;
      evaluate(0, lambda.j(), lambda.e()[1], lambda.c()[1], lambda.k()[1], knots[1], values[1]);
      for (int k0 = 0; k0 < N; k0++) {
 	x[0] = knots[0][k0];
 	for (int k1 = 1; k1 < N; k1++) {
 	  x[1] = knots[1][k1];
 	  s += values[0][k0] * values[1][k1] * f.value(x);
 	}
      }
    }

    const double flambda = fcoeffs.get_coefficient(lambda);
    const double flambdaapprox = s*h*h;

    cout << "lambda=" << lambda << ": "
	 << "f(lambda): " << flambda
	 << ", approx.: " << flambdaapprox
	 << ", deviation: " << fabs(flambda-flambdaapprox) << endl;

    if (lambda == basis.last_wavelet(jmax)) break;
  }
  


  return 0;
}
