// imlementation of gauss_quadrature.h inline functions

#include <cassert>
#include <numerics/gauss_data.h>

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
}
