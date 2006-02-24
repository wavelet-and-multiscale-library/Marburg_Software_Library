// implementation for splines.h

#include <algebra/triangular_matrix.h>

namespace MathTL
{
  template <int d>
  void
  compute_Bspline_refinement_matrix(const KnotSequence* knots, Matrix<double>& M)
  {
    // generalization of Lemma 3.16 from [P]

    const int k0 = knots->k0();

    Matrix<double> B1(d-1-k0, -k0);
    for (int n = 0; n <= d-2-k0; n++)
      for (int k = 0; k < -k0; k++)
	B1.set_entry(n, k, evaluate_Bspline<d>(knots, k0+k, n/2.));

    UpperTriangularMatrix<double> B2(d-1-k0, d-1-k0);
    for (int k = 0; k <= d-2-k0; k++)
      for (int n = k; n <= d-2-k0; n++)
	B2.set_entry(k, n, evaluate_Bspline<d>(knots, k0+n, k));

    UpperTriangularMatrix<double> B2Inv;
    B2.inverse(B2Inv);

    M = Matrix<double>(B2Inv) * B1;
  }
}
