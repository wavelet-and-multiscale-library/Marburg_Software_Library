// implementation for boundary_gramian.h

#include <numerics/matrix_decomp.h>

namespace WaveletTL
{
  template <class MASK1, class MASK2>
  void
  compute_biorthogonal_boundary_gramian(const Matrix<double>& ML,
					const Matrix<double>& MLT,
					Matrix<double>& GammaL)
  {
    // The system for GammaL 
    //   GammaL = ThetaL * ThetaLTilde^T
    //          = 1/2 * M_L^T * (Gamma_L 0) * MTilde_L
    //                          (0       I)
    // can be written as a linear system of equations for col(GammaL), where
    //   col(A) = (a_1^T,...,a_n^T)^T, for A=(a_1,...,a_n).
    // More precisely, because of the identity col(ABC) = kron(C^T, A) * col(B),
    // and using corresponding partitions
    //   M_L = (A),  MTilde_L = (ATilde)
    //         (B)              (BTilde)
    // one can transform the equation for GammaL into
    //   col(GammaL) = 1/2 * col(A^T*GammaL*ATilde + B^T*BTilde)
    //               = 1/2 * (kron(ATilde^T, A^T)*col(GammaL) + col(B^T*BTilde))
    // i.e. the vector col(GammaL) solves the linear equation
    //   (I-1/2*kron(ATilde^T, A^T)) * col(GammaL) = 1/2*col(B^T*BTilde).

    const unsigned int n = ML.column_dimension(); // for readability
    Matrix<double> A(n, n), B(MLT.row_dimension()-n, n);
    ML.get_block(0, 0, n, n, A);
    ML.get_block(n, 0, ML.row_dimension()-n, n, B); // here we need the default param. of get_block()
    Matrix<double> AT(n, n), BT(MLT.row_dimension()-n, n);
    MLT.get_block(0, 0, n, n, AT);
    MLT.get_block(n, 0, MLT.row_dimension()-n, n, BT);
//     cout << "A=" << endl << A;
//     cout << "B=" << endl << B;
//     cout << "AT=" << endl << AT;
//     cout << "BT=" << endl << BT;

    KroneckerMatrix<double,Matrix<double>,Matrix<double> > K(transpose(AT), transpose(A));
    Matrix<double> H(K);
    H.scale(0.5);
    Matrix<double> I; I.diagonal(n*n, 1.0);
    Matrix<double> S = I - H;
//     cout << "S=" << endl << S;

    Matrix<double> R(transpose(B)*BT);
    R.scale(0.5);
    Vector<double> b;
    R.col(b);
//     cout << "b=" << b << endl;

    Vector<double> x;
    MathTL::QUDecomposition<double>(S).solve(b, x);
    GammaL.decol(x, n);
  }
}
