// implementation for regularity.h

#include <algebra/vector.h>
#include <algebra/matrix.h>
#include <numerics/eigenvalues.h>

using namespace MathTL;

namespace WaveletTL
{
  template <class MASK>
  AutocorrelationMask<MASK>::AutocorrelationMask()
  {
    MASK a;

    // b(k) = \sum_m a(k+m)*a(m)/2
    const int k1 = a.begin().index()[0];
    const int k2 = a.rbegin().index()[0];

    for (int k = k1-k2; k <= k2-k1; k++)
      for (int m = k1; m <= k2; m++)
	if (m+k >= k1 && m+k <= k2)
	  set_coefficient(MultiIndex<int,1>(k),
			  get_coefficient(MultiIndex<int,1>(k))
			  + a.get_coefficient(MultiIndex<int,1>(m))
			  * a.get_coefficient(MultiIndex<int,1>(m+k))/2.0);
  }
  
  template <class MASK>
  double
  Sobolev_regularity()
  {
    double r = 0;

    // setup the autocorrelation mask according to the given mask
    AutocorrelationMask<MASK> b;
//     cout << "* the autocorrelation mask of a=" << MASK() << " is " << b << endl;

    // setup A=(b_{2i-j})_{i,j\in\Omega}
    const int N = b.rbegin().index()[0]; // offset for A
    Matrix<double> A(2*N+1, 2*N+1);
    for (int i = -N; i <= N; i++)
      for (int j = -N; j <= N; j++)
	A(N+i, N+j) = b.get_coefficient(MultiIndex<int,1>(2*i-j));
    cout << "A=" << endl << A << endl;

    // TODO: solve nonsymmetric eigenvalue problem here!

    return r;
  }
}
