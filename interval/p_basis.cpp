// implementation for p_basis.h

#include <cassert>
#include <cmath>
#include <numerics/schoenberg_splines.h>
#include <algebra/triangular_matrix.h>
#include <utils/tiny_tools.h>
#include <interval/boundary_gramian.h>

namespace WaveletTL
{
  template <int d, int dT>
  PBasis<d,dT>::PBasis(const int s0, const int s1, const int sT0, const int sT1) {
    assert(std::max(s0,s1) < d && std::max(sT0,sT1) < dT);
        
    this->s0 = s0;
    this->s1 = s1;
    this->sT0 = sT0;
    this->sT1 = sT1;

    setup();
  }


  template <int d, int dT>
  void
  PBasis<d,dT>::setup() {
    // For simplicity, we do not implement a generic setup here
    // but only setup the parameters for several important special cases from [P].
    // Namely, we restrict the setting to the case d=2+s0, where no "additional"
    // interior dual generators have to be constructed. In a later version of
    // this class, we will fix this.

    
    // setup j0
    switch(d) {
    case 2:
      j0_ = 3;
      break;
    case 3:
      j0_ = 4;
      break;
    default:
      j0_ = 0;
      break;
    }

    // setup the refinement matrix block for all "real" primal boundary B-splines,
    // obeying the block structure (3.15)
    // (ignoring the current values of s0 and s1
    MathTL::SchoenbergKnotSequence<d> sknots;
    Matrix<double> ML_0;
    MathTL::compute_Bspline_refinement_matrix<d>(&sknots, ML_0);

    ML_.resize(3*dT+2*d-5, d+dT-2);
    ML_.set_block(0, 0, ML_0);
    for (int k = d; k <= d+dT-2; k++)
      for (int n = 2*k-d; n <= 2*k; n++)
	ML_.set_entry(n-1,k-1,cdf.a().get_coefficient(MultiIndex<int,1>(-(d/2)+n+d-2*k)));
    
    cout << "ML=" << endl << ML_;

    // setup the expansion coefficients of the unbiorthogonalized dual generators
    // w.r.t. the truncated CDF generators, see [P, Bem. 4.2]
    MathTL::LowerTriangularMatrix<double> D1(dT, dT), D2(dT, dT), D3(dT, dT);
    for (unsigned int r = 0; r < dT; r++)
      for (unsigned int i = 0; i <= r; i++)
 	D1.set_entry(r, i, binomial(r, i) * alpha(0, r-i)); // can be speeded up a bit...
    
    for (unsigned int i = 0; i < dT; i++)
      for (unsigned int n = 0; n <= i; n++)
	D2.set_entry(i, n, minus1power(n) * binomial(i, n) * intpower(-ell1T<d,dT>()-1, i-n));

    D3.set_entry(0, 0, 1.0); // the only nontrivial value for n=0
    for (unsigned int n = 1; n < dT; n++)
      for (unsigned int l = 1; l <= n; l++) {
	double help = 0;
	for (unsigned int k = 0; k < n; k++)
	  help += binomial(n, k) * D3.get_entry(k, l-1);
	D3.set_entry(n, l, help);
      }

    MathTL::LowerTriangularMatrix<double> DTilde = D1 * D2 * D3;
    MathTL::LowerTriangularMatrix<double> DTildeInv; DTilde.inverse(DTildeInv);
//     cout << "D1=" << endl << D1;
//     cout << "D2=" << endl << D2;
//     cout << "D3=" << endl << D3;
//     cout << "DTilde=" << endl << DTilde;

    MLT_.resize(3*dT+d-3, dT);
    MathTL::Matrix<double> muT(3*dT+d-3, dT), JdT(dT, dT);
    for (unsigned int n = 0; n < dT; n++)
      JdT.set_entry(n, dT-n-1, 1.0);
//     cout << "JdT=" << endl << JdT;
    
    Matrix<double> D0;
    D0.diagonal(3*dT+d-3, 1.0);
    D0.set_block(0, 0, JdT * transpose(Matrix<double>(DTilde)));
//     cout << "D0=" << endl << D0;

    for (unsigned int n = 0; n < dT; n++)
      muT.set_entry(n, n, ldexp(1.0, -n));
    for (unsigned int n = dT; n < 3*dT+d-3; n++)
      for (unsigned int k = 0; k < dT; k++)
	muT.set_entry(n, k, betaL(n+ell2<d>()-2, k)); // corrects a misprint in [P]
//     cout << "muT=" << endl << muT;
    
    MLT_ = D0 * muT * transpose(Matrix<double>(DTildeInv)) * JdT;
    MLT_.compress(1e-10);

    cout << "MLT=" << endl << MLT_;

    // for the biorthogonalization of the generators,
    // compute the gramian matrix of the primal and dual boundary generators
    compute_biorthogonal_boundary_gramian
      <CDFMask_primal<d>,CDFMask_dual<d,dT> >(ML_, MLT_, GammaL);
    cout << "GammaL=" << endl << GammaL;    
  }

  template <int d, int dT>
  const double
  PBasis<d,dT>::alpha(const int m, const unsigned int r) const {
    double result = 0;
    if (r == 0)
      return 1; // [DKU] (5.1.1)
    else {
      if (m == 0) {
	// [DKU] (5.1.3)
	for (int k = ell1<d>(); k <= ell2<d>(); k++) {
	  double help = 0;
	  for (unsigned int s = 0; s < r; s++)
	    help += binomial(r, s) * intpower(k, r-s) * alpha(0, s);
	  result += cdf.a().get_coefficient(MultiIndex<int,1>(k)) * help;
	}
	result /= ldexp(1.0, r+1) - 2.0;
      } else {
	// [DKU] (5.1.2)
	for (unsigned int i = 0; i <= r; i++)
	  result += binomial(r, i) * intpower(m, i) * alpha(0, r-i);
      }
    }
    return result;
  }
  
  template <int d, int dT>
  const double
  PBasis<d,dT>::betaL(const int m, const unsigned int r) const {
    // [DKU] (3.2.31)
    double result = 0;
    for (int q = (int)ceil((m-ell2T<d,dT>())/2.0); q < -ell1T<d,dT>(); q++)
      result += alpha(q, r) * cdf.aT().get_coefficient(MultiIndex<int,1>(m-2*q));
    return result;
  }

}
