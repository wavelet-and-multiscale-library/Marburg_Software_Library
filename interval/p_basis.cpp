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
  PBasis<d,dT>::PBasis(const int s0, const int s1) {
    assert(d >= 2 && d <= 3);
    assert(s0 >= d-2 && s1 >= d-2);

    this->s0 = s0;
    this->s1 = s1;

    setup();
  }


  template <int d, int dT>
  void
  PBasis<d,dT>::setup() {
    // For simplicity, we do not implement a fully generic setup here
    // but only setup the parameters for several important special cases from [P].
    // Namely, we restrict the setting to the case s0 >= d-2, where no "additional"
    // interior dual generators have to be constructed. In a later version of
    // this class, we will fix this.

    // setup j0 (TODO: use more generic formula?)
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
    // (ignoring the current values of s0 (and s1))
    MathTL::SchoenbergKnotSequence<d> sknots;
    Matrix<double> ML_0;
    MathTL::compute_Bspline_refinement_matrix<d>(&sknots, ML_0);
    
    // The total number of (left) boundary generators is always exactly dT
    // (to be able reproduce all polynomials by the dual generators).
    ML_.resize(3*dT+s0-1, dT);
    for (int column = s0; column < d-1; column++)
      for (int row = s0; row < 2*(d-1); row++)
	ML_.set_entry(row-s0, column-s0, ML_0.get_entry(row,column));
    for (int k = d; k <= dT+s0; k++)
      for (int n = 2*k-d; n <= 2*k; n++)
 	ML_.set_entry(n-1-s0,k-1-s0,cdf.a().get_coefficient(MultiIndex<int,1>(-(d/2)+n+d-2*k)));
    cout << "ML=" << endl << ML_;

    MR_.resize(3*dT+s1-1, dT);
    for (int column = s1; column < d-1; column++)
      for (int row = s1; row < 2*(d-1); row++)
	MR_.set_entry(row-s1, column-s1, ML_0.get_entry(row,column));
    for (int k = d; k <= dT+s1; k++)
      for (int n = 2*k-d; n <= 2*k; n++)
 	MR_.set_entry(n-1-s1,k-1-s1,cdf.a().get_coefficient(MultiIndex<int,1>(-(d/2)+n+d-2*k)));
    cout << "MR=" << endl << MR_;

    setup_Mj0(ML_, MR_, Mj0); // [DKU, (3.5.1)]

    // setup the expansion coefficients of the unbiorthogonalized dual generators
    // w.r.t. the truncated CDF generators, see [DKU, Lemma 3.1]
    Matrix<double> MLTp(3*dT+s0-1, dT);
    for (unsigned int r = 0; r < dT; r++) {
      MLTp.set_entry(r, r, ldexp(1.0, -r));
      for (int m = ellT_l(); m <= 2*ellT_l()+ell1T<d,dT>()-1; m++)
	MLTp.set_entry(dT+m-ellT_l(), r, alpha(m, r));
      for (int m = 2*ellT_l()+ell1T<d,dT>(); m <= 2*ellT_l()+ell2T<d,dT>()-2; m++)
	MLTp.set_entry(dT+m-ellT_l(), r, betaL(m, r));
    }
    cout << "MLTp=" << endl << MLTp;
        
    Matrix<double> MRTp(3*dT+s1-1, dT);
    for (unsigned int r = 0; r < dT; r++) {
      MRTp.set_entry(r, r, ldexp(1.0, -r));
      for (int m = ellT_r(); m <= 2*ellT_r()+ell1T<d,dT>()-1; m++)
	MRTp.set_entry(dT+m-ellT_r(), r, alpha(m, r));
      for (int m = 2*ellT_r()+ell1T<d,dT>(); m <= 2*ellT_r()+ell2T<d,dT>()-2; m++)
	MRTp.set_entry(dT+m-ellT_r(), r, betaR(m, r));
    }
    cout << "MRTp=" << endl << MRTp;
        
    // for the biorthogonalization of the generators,
    // compute the gramian matrix of the primal and dual boundary generators
    compute_biorthogonal_boundary_gramian
      <CDFMask_primal<d>,CDFMask_dual<d,dT> >(ML_, MLTp, GammaL);
    cout << "GammaL=" << endl << GammaL;

    compute_biorthogonal_boundary_gramian
      <CDFMask_primal<d>,CDFMask_dual<d,dT> >(MR_, MRTp, GammaR);
    cout << "GammaR=" << endl << GammaR;
   
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
    for (int q = (int)ceil((m-ell2T<d,dT>())/2.0); q < ellT_l(); q++)
      result += alpha(q, r) * cdf.aT().get_coefficient(MultiIndex<int,1>(m-2*q));
    return result;
  }

  template <int d, int dT>
  const double
  PBasis<d,dT>::betaR(const int m, const unsigned int r) const {
    // [DKU] (3.2.31)
    double result = 0;
    for (int q = (int)ceil((m-ell2T<d,dT>())/2.0); q < ellT_r(); q++)
      result += alpha(q, r) * cdf.aT().get_coefficient(MultiIndex<int,1>(m-2*q));
    return result;
  }

  template <int d, int dT>
  void
  PBasis<d,dT>::setup_Mj0(const Matrix<double>& ML, const Matrix<double>& MR, SparseMatrix<double>& Mj0) {
    // IGPMlib reference: I_Basis_Bspline_s::Mj0()
    // cf. [DKU section 3.5]

    const int nj  = Deltasize(j0());
    const int njp = Deltasize(j0()+1);
    Mj0.resize(njp, nj);

    for (unsigned int i = 0; i < ML.row_dimension(); i++)
      for (unsigned int k = 0; k < ML.column_dimension(); k++)
	Mj0.set_entry(i, k, ML.get_entry(i, k));
    
    for (unsigned int i = 0; i < MR.row_dimension(); i++)
      for (unsigned int k = 0; k < MR.column_dimension(); k++)
	Mj0.set_entry(njp-i-1, nj-k-1, MR.get_entry(i, k));
    
    int startrow = d+ell_l()+ell1<d>()-2*s0;
    for (int col = d-s0; col < nj-(d-s1); col++, startrow+=2) {
      int row = startrow;
      for (MultivariateLaurentPolynomial<double, 1>::const_iterator it(cdf.a().begin());
	   it != cdf.a().end(); ++it, row++)
	Mj0.set_entry(row, col, *it);
    }
    
    cout << "Mj0=" << endl << Mj0 << endl;
  }
}
