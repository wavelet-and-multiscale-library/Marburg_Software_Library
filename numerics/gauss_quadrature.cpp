// imlementation of gauss_quadrature.h inline functions

#include <cassert>
#include <cmath>
#include <iostream>
#include <numerics/gauss_data.h>
#include <algebra/vector.h>
#include <algebra/matrix.h>
#include <algebra/symmetric_matrix.h>
#include <numerics/eigenvalues.h>

namespace MathTL
{
  GaussLegendreRule::GaussLegendreRule(const unsigned int N)
  {
    assert(N <= 10);

    points_.resize(N);
    weights_.resize(N);
    for (unsigned int n(0); n < N; n++)
      {
	// the precomputed values live on [-1,1]
	points_[n]  = 0.5 + GaussPoints[N-1][n]/2.0;
	weights_[n] = GaussWeights[N-1][n];
      }
  }

  GaussRule::GaussRule(const OrthogonalPolynomial& poly,
		       const double a, const double b,
		       const unsigned int N)
  {
    points_.resize(N);
    weights_.resize(N);

    SymmetricMatrix<double> J(N);
    for (unsigned int n(0); n < N; n++)
      {
	J(n, n) = poly.a(n+1);
	if (n < N-1)
	  J(n, n+1) = sqrt(poly.b(n+2));
      }

    // the eigenvalues of J are the Gauss points,
    // the first component of the respective eigenvector is the
    // square root of the weight
    Matrix<double> evecs;
    Vector<double> evals;
    SymmEigenvalues(J,evals,evecs);

    // copy points and weights, rescale interval
    for (unsigned int n(0); n < N; n++)
      {
	points_[n] = (evals[n] - a)/(b-a);
	weights_[n] = evecs(0, n) * evecs(0, n);
      }
  }
}
