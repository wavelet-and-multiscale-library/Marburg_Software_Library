// implementation for ds_basis

#include <cassert>
#include <cmath>

#include <numerics/matrix_decomp.h>

#include <Rd/haar_mask.h>
#include <Rd/cdf_mask.h>
#include <Rd/dm_mask.h>
#include <Rd/refinable.h>
#include <Rd/r_index.h>

namespace WaveletTL
{
  template <int d, int dT>
  DSBasis<d,dT>::DSBasis(const int s0, const int s1, const int sT0, const int sT1,
			 DSBiorthogonalizationMethod bio) {
    assert(std::max(s0,s1) < d && std::max(sT0,sT1) < dT);
    
//     assert(((s0 == 0 && sT0 > 0) || (s0 > 0 && sT0 == 0))
// 	   && ((s1 == 0 && sT1 > 0) || (s1 > 0 && sT1 == 0)));
    
    this->s0 = s0;
    this->s1 = s1;
    this->sT0 = sT0;
    this->sT1 = sT1;
    this->bio = bio;

#if 0
    cout << "* some moments alpha(m,r):" << endl;
    for (unsigned int r = 1; r < d; r++)
      for (int m = -5; m <= 5; m++)
	cout << "m=" << m << ", r=" << r << ": " << alpha(m,r) << endl;
    cout << "* some moments alphaT(m,r):" << endl;
    for (unsigned int r = 1; r < dT; r++)
      for (int m = -5; m <= 5; m++)
	cout << "m=" << m << ", r=" << r << ": " << alphaT(m,r) << endl;
#endif

    setup_GammaLR();
    setup_CX_CXT();
    setup_CXA_CXAT();

    // IGPMlib reference: I_Basis_Bspline_s::Setup()

    setup_Cj();

    Matrix<double> ml = ML(); // (3.5.2)
    Matrix<double> mr = MR(); // (3.5.2)
    SparseMatrix<double> mj0;   setup_Mj0  (ml,   mr,   mj0);   // (3.5.1)

    Matrix<double> mltp = MLTp(); // (3.5.6)
    Matrix<double> mrtp = MRTp(); // (3.5.6)
    SparseMatrix<double> mj0tp; setup_Mj0Tp(mltp, mrtp, mj0tp); // (3.5.5)

    // biorthogonalize the generators
    Mj0  = transpose(inv_Cjp)  * mj0   * transpose(Cj); // (2.4.3)
    Mj0T = transpose(inv_CjpT) * mj0tp * transpose(CjT);
    
#if 1
    cout << "DSBasis(): check biorthogonality of Mj0, Mj0T:" << endl;
//     cout << "Mj0=" << endl << Mj0 << endl << "Mj0T=" << endl << Mj0T << endl;

    SparseMatrix<double> testbio0 = transpose(Mj0) * Mj0T;
//     cout << "Mj0^T*Mj0T=" << endl << testbio0 << endl;
    for (unsigned int i = 0; i < testbio0.row_dimension(); i++)
      testbio0.set_entry(i, i, testbio0.get_entry(i, i) - 1.0);
    cout << "* ||Mj0^T*Mj0T-I||_infty: " << row_sum_norm(testbio0) << endl;

    testbio0 = transpose(Mj0T) * Mj0;
//     cout << "Mj0T*Mj0^T=" << endl << testbio0 << endl;
    for (unsigned int i = 0; i < testbio0.row_dimension(); i++)
      testbio0.set_entry(i, i, testbio0.get_entry(i, i) - 1.0);
    cout << "* ||Mj0T^T*Mj0-I||_infty: " << row_sum_norm(testbio0) << endl;
#endif    

    // construction of the wavelet basis: initial stable completion, [DKU section 4.1]

    SparseMatrix<double> FF; F(FF);         // (4.1.14)
    SparseMatrix<double> PP; P(ml, mr, PP); // (4.1.22)
 
    SparseMatrix<double> A, H, Hinv;
    GSetup(A, H, Hinv); // (4.1.1), (4.1.13)

#if 1
    SparseMatrix<double> Aold(A); // for the checks below
#endif

    GElim (A, H, Hinv); // elimination (4.1.4)ff.
    SparseMatrix<double> BB; BT(A, BB); // (4.1.13)

#if 1
    cout << "DSBasis(): check properties (4.1.15):" << endl;
    SparseMatrix<double> test4115 = transpose(BB)*A;
    for (unsigned int i = 0; i < test4115.row_dimension(); i++)
      test4115.set_entry(i, i, test4115.get_entry(i, i) - 1.0);
    cout << "* ||Bj*Ajd-I||_infty: " << row_sum_norm(test4115) << endl;

    test4115 = transpose(FF)*FF;
    for (unsigned int i = 0; i < test4115.row_dimension(); i++)
      test4115.set_entry(i, i, test4115.get_entry(i, i) - 1.0);
    cout << "* ||Fj^T*Fj-I||_infty: " << row_sum_norm(test4115) << endl;    

    test4115 = transpose(BB)*FF;
    cout << "* ||Bj*Fj||_infty: " << row_sum_norm(test4115) << endl;    

    test4115 = transpose(FF)*A;
    cout << "* ||Fj^T*A||_infty: " << row_sum_norm(test4115) << endl;    
#endif

#if 1
    cout << "DSBasis(): check factorization of A:" << endl;
    SparseMatrix<double> testAfact = Aold - Hinv*A;
    cout << "* in infty-norm: " << row_sum_norm(testAfact) << endl;

    cout << "DSBasis(): check that H is inverse to Hinv:" << endl;
    SparseMatrix<double> testHinv = H*Hinv;
    for (unsigned int i = 0; i < testHinv.row_dimension(); i++)
      testHinv.set_entry(i, i, testHinv.get_entry(i, i) - 1.0);
    cout << "* in infty-norm: " << row_sum_norm(testHinv) << endl;
#endif

    SparseMatrix<double> mj1ih = PP * Hinv * FF; // (4.1.23)
    SparseMatrix<double> PPinv; InvertP(PP, PPinv);

#if 1
    cout << "DSBasis(): check that PPinv is inverse to PP:" << endl;
    SparseMatrix<double> testPinv = PP*PPinv;
    for (unsigned int i = 0; i < testPinv.row_dimension(); i++)
      testPinv.set_entry(i, i, testPinv.get_entry(i, i) - 1.0);
    cout << "* in infty-norm: " << row_sum_norm(testPinv) << endl;
#endif

    SparseMatrix<double> help = H * PPinv;
    SparseMatrix<double> gj0ih = transpose(BB) * help;
    SparseMatrix<double> gj1ih = transpose(FF) * help; // (4.1.24)

#if 1
    cout << "DSBasis(): check initial stable completion:" << endl;
    SparseMatrix<double> mj_initial(mj0.row_dimension(),
				    mj0.column_dimension() + mj1ih.column_dimension());
    for (unsigned int i = 0; i < mj0.row_dimension(); i++)
      for (unsigned int j = 0; j < mj0.column_dimension(); j++)	{
	const double help = mj0.get_entry(i, j);
	if (help != 0)
	  mj_initial.set_entry(i, j, help);
      }
    for (unsigned int i = 0; i < mj1ih.row_dimension(); i++)
      for (unsigned int j = 0; j < mj1ih.column_dimension(); j++) {
	const double help = mj1ih.get_entry(i, j);
	if (help != 0)
	  mj_initial.set_entry(i, j+mj0.column_dimension(), help);
      }
    
    SparseMatrix<double> gj_initial(gj0ih.row_dimension() + gj1ih.row_dimension(),
 				    gj0ih.column_dimension());
    for (unsigned int i = 0; i < gj0ih.row_dimension(); i++)
      for (unsigned int j = 0; j < gj0ih.column_dimension(); j++) {
	const double help = gj0ih.get_entry(i, j);
	if (help != 0)
	  gj_initial.set_entry(i, j, help);
      }
    for (unsigned int i = 0; i < gj1ih.row_dimension(); i++)
      for (unsigned int j = 0; j < gj1ih.column_dimension(); j++) {
	const double help = gj1ih.get_entry(i, j);
	if (help != 0)
	  gj_initial.set_entry(i+gj0ih.row_dimension(), j, help);
      }
    
    SparseMatrix<double> test_initial = mj_initial * gj_initial;
    for (unsigned int i = 0; i < test_initial.row_dimension(); i++)
      test_initial.set_entry(i, i, test_initial.get_entry(i, i) - 1.0);
    cout << "* ||Mj*Gj-I||_infty: " << row_sum_norm(test_initial) << endl;

    test_initial = gj_initial * mj_initial;
    for (unsigned int i = 0; i < test_initial.row_dimension(); i++)
      test_initial.set_entry(i, i, test_initial.get_entry(i, i) - 1.0);
    cout << "* ||Gj*Mj-I||_infty: " << row_sum_norm(test_initial) << endl;
#endif

    // construction of the wavelet basis: stable completion with basis transformations
    SparseMatrix<double> I; I.diagonal(Deltasize(j0()+1), 1.0);
    Mj1  = (I - (Mj0*transpose(Mj0T))) * (transpose(inv_Cjp) * mj1ih);
    Mj1T = Cjp * transpose(gj1ih);

    Mj0 .compress(1e-8);
    Mj1 .compress(1e-8);
    Mj0T.compress(1e-8);
    Mj1T.compress(1e-8);
    
#if 1
    cout << "DSBasis(): check new stable completion:" << endl;
    
    SparseMatrix<double> mj_new(Mj0.row_dimension(),
				Mj0.column_dimension() + Mj1.column_dimension());
    for (unsigned int i = 0; i < Mj0.row_dimension(); i++)
      for (unsigned int j = 0; j < Mj0.column_dimension(); j++) {
	const double help = Mj0.get_entry(i, j);
	if (help != 0)
	  mj_new.set_entry(i, j, help);
      }
    for (unsigned int i = 0; i < Mj1.row_dimension(); i++)
      for (unsigned int j = 0; j < Mj1.column_dimension(); j++) {
	const double help = Mj1.get_entry(i, j);
	if (help != 0)
	  mj_new.set_entry(i, j+Mj0.column_dimension(), help);
      }
    
    SparseMatrix<double> gj0_new = transpose(Mj0T); gj0_new.compress();
    SparseMatrix<double> gj1_new = transpose(Mj1T); gj1_new.compress();
    SparseMatrix<double> gj_new(gj0_new.row_dimension() + gj1_new.row_dimension(),
				gj0_new.column_dimension());
    for (unsigned int i = 0; i < gj0_new.row_dimension(); i++)
      for (unsigned int j = 0; j < gj0_new.column_dimension(); j++) {
	const double help = gj0_new.get_entry(i, j);
	if (help != 0)
	  gj_new.set_entry(i, j, help);
      }
    for (unsigned int i = 0; i < gj1_new.row_dimension(); i++)
      for (unsigned int j = 0; j < gj1_new.column_dimension(); j++) {
	const double help = gj1_new.get_entry(i, j);
	if (help != 0)
	  gj_new.set_entry(i+gj0_new.row_dimension(), j, help);
      }
    
    SparseMatrix<double> test_new = mj_new * gj_new;
    for (unsigned int i = 0; i < test_new.row_dimension(); i++)
      test_new.set_entry(i, i, test_new.get_entry(i, i) - 1.0);
//     cout << "Mj*Gj-I=" << endl << test_new << endl;
    cout << "* ||Mj*Gj-I||_infty: " << row_sum_norm(test_new) << endl;

    test_new = gj_new * mj_new;
    for (unsigned int i = 0; i < test_new.row_dimension(); i++)
      test_new.set_entry(i, i, test_new.get_entry(i, i) - 1.0);
//     cout << "Gj*Mj-I=" << endl << test_new << endl;
    cout << "* ||Gj*Mj-I||_infty: " << row_sum_norm(test_new) << endl;
        
    SparseMatrix<double> mjt_new(Mj0T.row_dimension(),
				 Mj0T.column_dimension() + Mj1T.column_dimension());
    for (unsigned int i = 0; i < Mj0T.row_dimension(); i++)
      for (unsigned int j = 0; j < Mj0T.column_dimension(); j++) {
	const double help = Mj0T.get_entry(i, j);
	if (help != 0)
	  mjt_new.set_entry(i, j, help);
      }
    for (unsigned int i = 0; i < Mj1T.row_dimension(); i++)
      for (unsigned int j = 0; j < Mj1T.column_dimension(); j++) {
	const double help = Mj1T.get_entry(i, j);
	if (help != 0)
	  mjt_new.set_entry(i, j+Mj0T.column_dimension(), help);
      }
    
    SparseMatrix<double> gjt0_new = transpose(Mj0); gjt0_new.compress();
    SparseMatrix<double> gjt1_new = transpose(Mj1); gjt1_new.compress();
    SparseMatrix<double> gjt_new(gjt0_new.row_dimension() + gjt1_new.row_dimension(),
				 gjt0_new.column_dimension());
    for (unsigned int i = 0; i < gjt0_new.row_dimension(); i++)
      for (unsigned int j = 0; j < gjt0_new.column_dimension(); j++) {
	const double help = gjt0_new.get_entry(i, j);
	if (help != 0)
	  gjt_new.set_entry(i, j, help);
      }
    for (unsigned int i = 0; i < gjt1_new.row_dimension(); i++)
      for (unsigned int j = 0; j < gjt1_new.column_dimension(); j++) {
	const double help = gjt1_new.get_entry(i, j);
	if (help != 0)
	  gjt_new.set_entry(i+gjt0_new.row_dimension(), j, help);
      }
    
    test_new = mjt_new * gjt_new;
    for (unsigned int i = 0; i < test_new.row_dimension(); i++)
      test_new.set_entry(i, i, test_new.get_entry(i, i) - 1.0);
//     cout << "MjT*GjT-I=" << endl << test_new << endl;
    cout << "* ||MjT*GjT-I||_infty: " << row_sum_norm(test_new) << endl;

    test_new = gjt_new * mjt_new;
    for (unsigned int i = 0; i < test_new.row_dimension(); i++)
      test_new.set_entry(i, i, test_new.get_entry(i, i) - 1.0);
//     cout << "GjT*MjT-I=" << endl << test_new << endl;
    cout << "* ||GjT*MjT-I||_infty: " << row_sum_norm(test_new) << endl;
#endif

    // construction of the wavelet basis: symmetrization for odd d and symmetric b.c.'s
    if (d%2 && (s0 == s1 && sT0 == sT1))
      {
	DS_symmetrization(Mj1, Mj1T);
#if 1
	{
	  cout << "DKUBasis(): check [DS] symmetrization:" << endl;
	  
	  SparseMatrix<double> mj_symm(Mj0.row_dimension(),
				       Mj0.column_dimension() + Mj1.column_dimension());
	  for (unsigned int i = 0; i < Mj0.row_dimension(); i++)
	    for (unsigned int j = 0; j < Mj0.column_dimension(); j++) {
	      const double help = Mj0.get_entry(i, j);
	      if (help != 0)
		mj_symm.set_entry(i, j, help);
	    }
	  for (unsigned int i = 0; i < Mj1.row_dimension(); i++)
	    for (unsigned int j = 0; j < Mj1.column_dimension(); j++) {
	      const double help = Mj1.get_entry(i, j);
	      if (help != 0)
		mj_symm.set_entry(i, j+Mj0.column_dimension(), help);
	    }
	
	  SparseMatrix<double> gj0_symm = transpose(Mj0T); gj0_symm.compress();
	  SparseMatrix<double> gj1_symm = transpose(Mj1T); gj1_symm.compress();
	  SparseMatrix<double> gj_symm(gj0_symm.row_dimension() + gj1_symm.row_dimension(),
				       gj0_symm.column_dimension());
	  for (unsigned int i = 0; i < gj0_symm.row_dimension(); i++)
	    for (unsigned int j = 0; j < gj0_symm.column_dimension(); j++) {
	      const double help = gj0_symm.get_entry(i, j);
	      if (help != 0)
		gj_symm.set_entry(i, j, help);
	    }
	  for (unsigned int i = 0; i < gj1_symm.row_dimension(); i++)
	    for (unsigned int j = 0; j < gj1_symm.column_dimension(); j++) {
	      const double help = gj1_symm.get_entry(i, j);
	      if (help != 0)
		gj_symm.set_entry(i+gj0_symm.row_dimension(), j, help);
	    }
	
	  SparseMatrix<double> test_symm = mj_symm * gj_symm;
	  for (unsigned int i = 0; i < test_symm.row_dimension(); i++)
	    test_symm.set_entry(i, i, test_symm.get_entry(i, i) - 1.0);
	  cout << "* ||Mj*Gj-I||_infty: " << row_sum_norm(test_symm) << endl;
	
	  test_symm = gj_symm * mj_symm;
	  for (unsigned int i = 0; i < test_symm.row_dimension(); i++)
	    test_symm.set_entry(i, i, test_symm.get_entry(i, i) - 1.0);
	  cout << "* ||Gj*Mj-I||_infty: " << row_sum_norm(test_symm) << endl;
        
	  SparseMatrix<double> mjt_symm(Mj0T.row_dimension(),
					Mj0T.column_dimension() + Mj1T.column_dimension());
	  for (unsigned int i = 0; i < Mj0T.row_dimension(); i++)
	    for (unsigned int j = 0; j < Mj0T.column_dimension(); j++) {
	      const double help = Mj0T.get_entry(i, j);
	      if (help != 0)
		mjt_symm.set_entry(i, j, help);
	    }
	  for (unsigned int i = 0; i < Mj1T.row_dimension(); i++)
	    for (unsigned int j = 0; j < Mj1T.column_dimension(); j++) {
	      const double help = Mj1T.get_entry(i, j);
	      if (help != 0)
		mjt_symm.set_entry(i, j+Mj0T.column_dimension(), help);
	    }
	
	  SparseMatrix<double> gjt0_symm = transpose(Mj0); gjt0_symm.compress();
	  SparseMatrix<double> gjt1_symm = transpose(Mj1); gjt1_symm.compress();
	  SparseMatrix<double> gjt_symm(gjt0_symm.row_dimension() + gjt1_symm.row_dimension(),
					gjt0_symm.column_dimension());
	  for (unsigned int i = 0; i < gjt0_symm.row_dimension(); i++)
	    for (unsigned int j = 0; j < gjt0_symm.column_dimension(); j++) {
	      const double help = gjt0_symm.get_entry(i, j);
	      if (help != 0)
		gjt_symm.set_entry(i, j, help);
	    }
	  for (unsigned int i = 0; i < gjt1_symm.row_dimension(); i++)
	    for (unsigned int j = 0; j < gjt1_symm.column_dimension(); j++) {
	      const double help = gjt1_symm.get_entry(i, j);
	      if (help != 0)
		gjt_symm.set_entry(i+gjt0_symm.row_dimension(), j, help);
	    }
	
	  test_symm = mjt_symm * gjt_symm;
	  for (unsigned int i = 0; i < test_symm.row_dimension(); i++)
	    test_symm.set_entry(i, i, test_symm.get_entry(i, i) - 1.0);
	  cout << "* ||MjT*GjT-I||_infty: " << row_sum_norm(test_symm) << endl;
	
	  test_symm = gjt_symm * mjt_symm;
	  for (unsigned int i = 0; i < test_symm.row_dimension(); i++)
	    test_symm.set_entry(i, i, test_symm.get_entry(i, i) - 1.0);
	  cout << "* ||GjT*MjT-I||_infty: " << row_sum_norm(test_symm) << endl;
	}
#endif
      }
    
    Mj0_t  = transpose(Mj0);
    Mj0T_t = transpose(Mj0T);
    Mj1_t  = transpose(Mj1);
    Mj1T_t = transpose(Mj1T);
  }
  
  template <int d, int dT>
  const double
  DSBasis<d,dT>::alpha(const int m, const unsigned int r) const {
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
	result /= (double)((1<<(r+1))-2);
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
  DSBasis<d,dT>::alphaT(const int m, const unsigned int r) const {
    double result = 0;
    if (r == 0)
      return 1; // [DKU] (5.1.1)
    else {
      if (m == 0) {
	// [DKU] (5.1.3)
	for (int k = ell1T<d,dT>(); k <= ell2T<d,dT>(); k++) {
	  double help = 0;
	  for (unsigned int s = 0; s < r; s++)
	    help += binomial(r, s) * intpower(k, r-s) * alphaT(0, s);
	  result += cdf.aT().get_coefficient(MultiIndex<int,1>(k)) * help;
	}
	result /= (double)((1<<(r+1))-2);
      } else {
	// [DKU] (5.1.2)
	for (unsigned int i = 0; i <= r; i++)
	  result += binomial(r, i) * intpower(m, i) * alphaT(0, r-i);
      }
    }
    return result;
  }

  template <int d, int dT>
  const double
  DSBasis<d,dT>::betaL(const int m, const unsigned int r) const {
    // [DKU] (3.2.31)
    double result = 0;
    for (int q = (int)ceil((m-ell2T<d,dT>())/2.0); q < ellT_l()-sT0; q++)
      result += alpha(q, r) * cdf.aT().get_coefficient(MultiIndex<int,1>(m-2*q));
    return result * M_SQRT1_2;
  }

  template <int d, int dT>
  const double
  DSBasis<d,dT>::betaLT(const int m, const unsigned int r) const {
    // [DKU] (3.2.31)
    double result = 0;
    for (int q = (int)ceil((m-ell2<d>())/2.0); q < ell_l()-s0; q++)
      result += alphaT(q, r) * cdf.a().get_coefficient(MultiIndex<int,1>(m-2*q));
    return result * M_SQRT1_2;
  }

  template <int d, int dT>
  const double
  DSBasis<d,dT>::betaR(const int m, const unsigned int r) const {
    // [DKU] (3.2.31)
    double result = 0;
    for (int q = (int)ceil((m-ell2T<d,dT>())/2.0); q < ellT_r()-sT1; q++)
      result += alpha(q, r) * cdf.aT().get_coefficient(MultiIndex<int,1>(m-2*q));
    return result * M_SQRT1_2;
  }

  template <int d, int dT>
  const double
  DSBasis<d,dT>::betaRT(const int m, const unsigned int r) const {
    // [DKU] (3.2.31)
    double result = 0;
    for (int q = (int)ceil((m-ell2<d>())/2.0); q < ell_r()-s1; q++)
      result += alphaT(q, r) * cdf.a().get_coefficient(MultiIndex<int,1>(m-2*q));
    return result * M_SQRT1_2;
  }

  template <int d, int dT>
  void
  DSBasis<d,dT>::setup_GammaLR() {
    // IGPMlib reference: I_Mask_Bspline::EvalGammaL(), ::EvalGammaR()

    const unsigned int GammaLsize = std::max(d-s0, dT-sT0);
    GammaL.resize(GammaLsize, GammaLsize);

    const unsigned int GammaRsize = std::max(d-s1, dT-sT1);
    GammaR.resize(GammaRsize, GammaRsize);

    // 1. compute the integrals
    //      z(s,t) = \int_0^1\phi(x-s)\tilde\phi(x-t)\,dx
    //    exactly with the [DM] trick
    MultivariateRefinableFunction<DMMask2<HaarMask, CDFMask_primal<d>, CDFMask_dual<d,dT> >,2> zmask;
    InfiniteVector<double, MultiIndex<int,2> > zvalues(zmask.evaluate());

//     cout << "z=" << endl << zvalues << endl;

    // 2. compute the integrals
    //      I(nu,mu) = \int_0^\infty\phi(x-\nu)\tilde\phi(x-\mu)\,dx
    //    exactly using the z(s,t) values

    const int I1Low  = -ell2<d>()+1;
    const int I1Up   = ell_l()-d+dT-1;
    const int I1Up_r = ell_r()-d+dT-1;
    const int I2Low  = -ell2T<d,dT>()+1;
    const int I2Up   = ellT_l()-1;
    const int I2Up_r = ellT_r()-1;

    Matrix<double> I(std::max(I1Up,I1Up_r)-I1Low+1,
		     std::max(I2Up,I2Up_r)-I2Low+1);
    for (int nu = I1Low; nu < -ell1<d>(); nu++)
      for (int mu = I2Low; mu < -ell1T<d,dT>(); mu++) {
  	double help(0);
  	const int diff(mu - nu);
   	const int sLow = std::max(-ell2<d>()+1,-ell2T<d,dT>()+1-diff);
  	for (int s = sLow; s <= nu; s++)
  	  help += zvalues.get_coefficient(MultiIndex<int,2>(s, s+diff));
  	I(-I1Low+nu, -I2Low+mu) = help; // [DKU] (5.1.7)
      }
    for (int nu = I1Low; nu <= std::max(I1Up,I1Up_r); nu++)
      for (int mu = I2Low; mu <= std::max(I2Up,I2Up_r); mu++) {
 	if ((nu >= -ell1<d>()) || ((nu <= ell_l()-1) && (mu >= -ell1T<d,dT>())))
 	  I(-I1Low+nu, -I2Low+mu) = (nu == mu ? 1 : 0); // [DKU] (5.1.6)
      }

//     cout << "I=" << endl << I << endl;

    // 3. finally, compute the Gramian GammaL
    for (int r = s0; r < d; r++) {
      // phiL's against phiTL's
      for (int k = sT0; k < dT; k++) {
     	double help = 0;
     	for (int nu = I1Low; nu < ell_l()-s0; nu++)
     	  for (int mu = I2Low; mu < ellT_l()-sT0; mu++)
     	    help += alphaT(nu, r) * alpha(mu, k) * I(-I1Low+nu, -I2Low+mu);
    	GammaL(r-s0, k-sT0) = help; // [DKU] (5.1.4)
      }
      // phiL's against phiT's
      for (int k = dT; k < d+sT0-s0; k++) {
	double help = 0;
	for (int nu = I1Low; nu < ell_l()-s0; nu++)
	  help += alphaT(nu, r) * I(-I1Low+nu, -I2Low+ellT_l()-dT-sT0+k);
	GammaL(r-s0, k-sT0) = help; // analogous to [DKU] (5.1.4) (!)
      }
    }
    // phi's against phiTL's
    for (int r = d; r < dT+s0-sT0; r++)
      for (int k = sT0; k < dT; k++) {
    	double help = 0;
    	for (int mu = I2Low; mu < ellT_l()-sT0; mu++)
    	  help += alpha(mu, k) * I(-I1Low+ell_l()-d-s0+r, -I2Low+mu);
    	GammaL(r-s0, k-sT0) = help; // [DKU] (5.1.5)
      }

//     cout << "GammaL=" << endl << GammaL << endl;

    // The same for GammaR:

    for (int r = s1; r < d; r++) {
      // phiR's against phiTR's
      for (int k = sT1; k < dT; k++) {
     	double help = 0;
     	for (int nu = I1Low; nu < ell_r()-s1; nu++)
     	  for (int mu = I2Low; mu < ellT_r()-sT1; mu++)
     	    help += alphaT(nu, r) * alpha(mu, k) * I(-I1Low+nu, -I2Low+mu);
    	GammaR(r-s1, k-sT1) = help; // [DKU] (5.1.4)
      }
      // phiR's against phiT's
      for (int k = dT; k < d+sT1-s1; k++) {
	double help = 0;
	for (int nu = I1Low; nu < ell_r()-s1; nu++)
	  help += alphaT(nu, r) * I(-I1Low+nu, -I2Low+ellT_r()-dT-sT1+k);
	GammaR(r-s1, k-sT1) = help; // analogous to [DKU] (5.1.4) (!)
      }
    }
    // phi's against phiTR's
    for (int r = d; r < dT+s1-sT1; r++)
      for (int k = sT1; k < dT; k++) {
    	double help = 0;
    	for (int mu = I2Low; mu < ellT_r()-sT1; mu++)
    	  help += alpha(mu, k) * I(-I1Low+ell_r()-d-s1+r, -I2Low+mu);
    	GammaR(r-s1, k-sT1) = help; // [DKU] (5.1.5)
      }

//     cout << "GammaR=" << endl << GammaR << endl;
  }

  template <int d, int dT>
  void
  DSBasis<d,dT>::setup_CX_CXT()
  {
    // IGPMlib reference: I_Mask_Bspline::EvalCL(), ::EvalCR()

    CL.resize(GammaL.row_dimension(), GammaL.column_dimension());
    CLT.resize(GammaL.row_dimension(), GammaL.column_dimension());
    CR.resize(GammaR.row_dimension(), GammaR.column_dimension());
    CRT.resize(GammaR.row_dimension(), GammaR.column_dimension());

    if (bio == none) {
      CL.diagonal(GammaL.row_dimension(), 1.0);
      Matrix<double> CLGammaLInv;
      QUDecomposition<double>(GammaL).inverse(CLGammaLInv);
      CLT = transpose(CLGammaLInv);

      CR.diagonal(GammaR.row_dimension(), 1.0);
      Matrix<double> CRGammaRInv;
      QUDecomposition<double>(GammaR).inverse(CRGammaRInv);
      CRT = transpose(CRGammaRInv);      
    }
    
    if (bio == SVD) {
      MathTL::SVD<double> svd(GammaL);
      Matrix<double> U, V;
      Vector<double> S;
      svd.getU(U);
      svd.getV(V);
      svd.getS(S);
      for (unsigned int i = 0; i < GammaL.row_dimension(); i++) {
 	S[i] = 1.0 / sqrt(S[i]);
 	for (unsigned int j = 0; j < GammaL.row_dimension(); j++)
 	  CL(i, j)  = S[i] * U(j, i);
      }
      CL.compress(1e-14);
      Matrix<double> CLGammaLInv;
      QUDecomposition<double>(CL * GammaL).inverse(CLGammaLInv);
      CLT = transpose(CLGammaLInv);
      CLT.compress(1e-14);
      
      MathTL::SVD<double> svd_r(GammaR);
      svd_r.getU(U);
      svd_r.getV(V);
      svd_r.getS(S);
      for (unsigned int i = 0; i < GammaR.row_dimension(); i++) {
 	S[i] = 1.0 / sqrt(S[i]);
 	for (unsigned int j = 0; j < GammaR.row_dimension(); j++)
 	  CR(i, j)  = S[i] * U(j, i);
      }
      CR.compress(1e-14);
      Matrix<double> CRGammaRInv;
      QUDecomposition<double>(CR * GammaR).inverse(CRGammaRInv);
      CRT = transpose(CRGammaRInv);      
      CRT.compress(1e-14);
    }
    
//     if (bio == Bernstein) {
//       double b(0);
//       if (d == 2)
//  	b = 0.7; // cf. [DKU]
//       else {
//  	b = 3.6;
//       }
      
//       // CL : transformation to Bernstein polynomials (triangular)
//       //   (Z_b^m)_{r,l} = (-1)^{r-1}\binom{m-1}{r}\binom{r}{l}b^{-r}, r>=l, 0 otherwise
//       for (int j(0); j < d; j++)
//  	for (int k(0); k <= j; k++)
//  	  CL_(k, j) = minus1power(j-k) * std::pow(b, -j) * binomial(d-1, j) * binomial(j, k);
//       for (int j(d); j < lupT-llowT+1; j++)
//  	CL_(j, j) = 1.0;
      
//       Matrix<double> CLGammaLInv;
//       QUDecomposition<double>(CL_ * GammaL_).inverse(CLGammaLInv);
//       CLT_ = transpose(CLGammaLInv);
      
//       CLT_.compress(1e-14);

//       for (int j(0); j < d; j++)
// 	for (int k(0); k <= j; k++)
// 	  CR_(k, j) = minus1power(j-k) * std::pow(b, -j) * binomial(d-1, j) * binomial(j, k);
//       for (int j(d); j < rlowTh-rupTh+1; j++)
// 	CR_(j, j) = 1.0;
      
//       Matrix<double> CRGammaRInv;
//       QUDecomposition<double>(CR_ * GammaR_).inverse(CRGammaRInv);
//       CRT_ = transpose(CRGammaRInv);
      
//       CRT_.compress(1e-14);
//     }
    
    if (bio == partialSVD) {
      MathTL::SVD<double> svd(GammaL);
      Matrix<double> U, V;
      Vector<double> S;
      svd.getU(U);
      svd.getV(V);
      svd.getS(S);
      Matrix<double> GammaLInv;
      QUDecomposition<double>(GammaL).inverse(GammaLInv);
      const double a = 1.0 / GammaLInv(0, 0);
      Matrix<double> R(GammaL.row_dimension(), GammaL.column_dimension());
      for (unsigned int i = 0; i < R.row_dimension(); i++)
 	S[i] = sqrt(S[i]);
      for (unsigned int j = 0; j < R.row_dimension(); j++)
 	R(0, j) = a * V(0, j) / S[j];
      for (unsigned int i = 1; i < R.row_dimension(); i++)
 	for (unsigned int j = 0; j < R.row_dimension(); j++)
 	  R(i, j) = U(j, i) * S[j];
      for (unsigned int i = 0; i < R.row_dimension(); i++)
 	for (unsigned int j = 0; j < R.row_dimension(); j++)
 	  U(i, j) /= S[i];
      CL = R*U;
      CL.compress(1e-14);
      Matrix<double> CLGammaLInv;
      QUDecomposition<double>(CL * GammaL).inverse(CLGammaLInv);
      CLT = transpose(CLGammaLInv);
      CLT.compress(1e-14);

      MathTL::SVD<double> svd_r(GammaR);
      svd_r.getU(U);
      svd_r.getV(V);
      svd_r.getS(S);
      Matrix<double> GammaRInv;
      QUDecomposition<double>(GammaR).inverse(GammaRInv);
      const double a_r = 1.0 / GammaRInv(0, 0);
      R.resize(GammaR.row_dimension(), GammaR.column_dimension());
      for (unsigned int i = 0; i < R.row_dimension(); i++)
	S[i] = sqrt(S[i]);
      for (unsigned int j = 0; j < R.row_dimension(); j++)
 	R(0, j) = a_r * V(0, j) / S[j];
      for (unsigned int i = 1; i < R.row_dimension(); i++)
 	for (unsigned int j = 0; j < R.row_dimension(); j++)
 	  R(i, j) = U(j, i) * S[j];
      for (unsigned int i = 0; i < R.row_dimension(); i++)
 	for (unsigned int j = 0; j < R.row_dimension(); j++)
 	  U(i, j) /= S[i];
      CR = R*U;
      CR.compress(1e-14);
      Matrix<double> CRGammaRInv;
      QUDecomposition<double>(CR * GammaR).inverse(CRGammaRInv);
      CRT = transpose(CRGammaRInv);
      CRT.compress(1e-14);
    }

//     if (bio_ == BernsteinSVD) {
//       double b(0);
//       if (d == 2)
// 	b = 0.7; // cf. [DKU]
//       else {
// 	b = 3.6;
//       }
      
//       // CL_ : transformation to Bernstein polynomials (triangular)
//       //   (Z_b^m)_{r,l} = (-1)^{r-1}\binom{m-1}{r}\binom{r}{l}b^{-r}, r>=l, 0 otherwise
//       for (int j(0); j < d; j++)
// 	for (int k(0); k <= j; k++)
// 	  CL_(k, j) = minus1power(j-k) * std::pow(b, -j) * binomial(d-1, j) * binomial(j, k);
//       for (int j(d); j < lupT-llowT+1; j++)
// 	CL_(j, j) = 1.0;
      
//       Matrix<double> GammaLNew(CL_ * GammaL_);
      
//       MathTL::SVD<double> svd(GammaLNew);
//       Matrix<double> U, V;
//       Vector<double> S;
//       svd.getU(U);
//       svd.getV(V);
//       svd.getS(S);
      
//       Matrix<double> GammaLNewInv;
//       QUDecomposition<double>(GammaLNew).inverse(GammaLNewInv);
//       const double a = 1.0 / GammaLNewInv(0, 0);
      
//       Matrix<double> R(lupT-llowT+1, lupT-llowT+1);
//       for (int i(0); i < lupT-llowT+1; i++)
// 	S[i] = sqrt(S[i]);
//       for (int j(0); j < lupT-llowT+1; j++)
// 	R(0, j) = a * V(0, j) / S[j];
//       for (int i(1); i < lupT-llowT+1; i++)
// 	for (int j(0); j < lupT-llowT+1; j++)
// 	  R(i, j) = U(j, i) * S[j];
      
//       for (int i(0); i < lupT-llowT+1; i++)
// 	for (int j(0); j < lupT-llowT+1; j++)
// 	  U(i, j) /= S[i];
//       CL_ = R*U*CL_;
      
//       CL_.compress(1e-14);

//       Matrix<double> CLGammaLInv;
//       QUDecomposition<double>(CL_ * GammaL_).inverse(CLGammaLInv);
//       CLT_ = transpose(CLGammaLInv);

//       CLT_.compress(1e-14);

//       for (int j(0); j < d; j++)
// 	for (int k(0); k <= j; k++)
// 	  CR_(k, j) = minus1power(j-k) * std::pow(b, -j) * binomial(d-1, j) * binomial(j, k);
//       for (int j(d); j < rlowTh-rupTh+1; j++)
// 	CR_(j, j) = 1.0;

//       Matrix<double> GammaRNew(CR_ * GammaR_);
      
//       MathTL::SVD<double> svd_r(GammaRNew);
//       svd_r.getU(U);
//       svd_r.getV(V);
//       svd_r.getS(S);

//       Matrix<double> GammaRNewInv;
//       QUDecomposition<double>(GammaRNew).inverse(GammaRNewInv);
//       const double a_r = 1.0 / GammaRNewInv(0, 0);
      
//       R.resize(rlowTh-rupTh+1,rlowTh-rupTh+1);
//       for (int i(0); i < rlowTh-rupTh+1; i++)
// 	S[i] = sqrt(S[i]);
//       for (int j(0); j < rlowTh-rupTh+1; j++)
// 	R(0, j) = a_r * V(0, j) / S[j];
//       for (int i(1); i < rlowTh-rupTh+1; i++)
// 	for (int j(0); j < rlowTh-rupTh+1; j++)
// 	  R(i, j) = U(j, i) * S[j];
      
//       for (int i(0); i < rlowTh-rupTh+1; i++)
// 	for (int j(0); j < rlowTh-rupTh+1; j++)
// 	  U(i, j) /= S[i];
//       CR_ = R*U*CR_;

//       CR_.compress(1e-14);

//       Matrix<double> CRGammaRInv;
//       QUDecomposition<double>(CR_ * GammaR_).inverse(CRGammaRInv);
//       CRT_ = transpose(CRGammaRInv);

//       CRT_.compress(1e-14);
//     }
    
#if 1
    // check biorthogonality of the matrix product CL * GammaL * (CLT)^T
//     cout << "GammaL=" << endl << GammaL << endl;
//     cout << "CL=" << endl << CL << endl;
//     cout << "CLT=" << endl << CLT << endl;
    Matrix<double> check(CL * GammaL * transpose(CLT));
    for (unsigned int i(0); i < check.row_dimension(); i++)
      check(i, i) -= 1;
    cout << "error for CLT: " << row_sum_norm(check) << endl;
#endif

    QUDecomposition<double>(CL).inverse(inv_CL);
    QUDecomposition<double>(CLT).inverse(inv_CLT);

#if 1
    // check biorthogonality of the matrix product CR * GammaR * (CRT)^T
//     cout << "GammaR=" << endl << GammaR << endl;
//     cout << "CR=" << endl << CR << endl;
//     cout << "CRT=" << endl << CRT << endl;
    Matrix<double> check2(CR * GammaR * transpose(CRT));
    for (unsigned int i(0); i < check2.row_dimension(); i++)
      check2(i, i) -= 1;
    cout << "error for CRT: " << row_sum_norm(check2) << endl;
#endif

    QUDecomposition<double>(CR).inverse(inv_CR);
    QUDecomposition<double>(CRT).inverse(inv_CRT);
  }

  template <int d, int dT>
  void DSBasis<d,dT>::setup_CXA_CXAT() {
    // IGPMlib reference: I_Mask_Bspline::EvalCL(), ::EvalCR()

    // setup CLA <-> AlphaT * (CL)^T
    CLA.resize(ell_l()+ell2<d>()-s0+std::max(0,dT-d+s0-sT0)-1, CL.row_dimension());
    for (int i = 1-ell2<d>(); i < ell_l()-s0; i++)
      for (unsigned int r = s0; r < d; r++)
	CLA(ell2<d>()-1+i, r-s0) = alphaT(i, r);
    for (unsigned int r = d; r < CL.row_dimension()+s0; r++)
      CLA(ell2<d>()-1+ell_l()-s0+r-d, r-s0) = 1.0;

//     cout << "CLA before biorthogonalization:" << endl << CLA << endl;
    CLA = CLA * transpose(CL);
    CLA.compress(1e-12);
//     cout << "CLA after biorthogonalization:" << endl << CLA << endl;

    // setup CLAT <-> Alpha * (CLT)^T
    CLAT.resize(ellT_l()+ell2T<d,dT>()-sT0+std::max(0,d-dT+sT0-s0)-1, CLT.row_dimension());
    for (int i = 1-ell2T<d,dT>(); i < ellT_l()-sT0; i++)
      for (unsigned int r = sT0; r < dT; r++)
	CLAT(ell2T<d,dT>()-1+i, r-sT0) = alpha(i, r);
    for (unsigned int r = dT; r < CLT.row_dimension()+sT0; r++)
      CLAT(ell2T<d,dT>()-1+ellT_l()-sT0+r-dT, r-sT0) = 1.0;

//     cout << "CLAT before biorthogonalization:" << endl << CLAT << endl;
    CLAT = CLAT * transpose(CLT);
    CLAT.compress(1e-12);
//     cout << "CLAT after biorthogonalization:" << endl << CLAT << endl;

    // the same for CRA, CRAT:
    CRA.resize(ell_r()+ell2<d>()-s1+std::max(0,dT-d+s1-sT1)-1, CR.row_dimension());
    for (int i = 1-ell2<d>(); i < ell_r()-s1; i++)
      for (unsigned int r = s1; r < d; r++)
	CRA(ell2<d>()-1+i, r-s1) = alphaT(i, r);
    for (unsigned int r = d; r < CR.row_dimension()+s1; r++)
      CRA(ell2<d>()-1+ell_r()-s1+r-d, r-s1) = 1.0;

//     cout << "CRA before biorthogonalization:" << endl << CRA << endl;
    CRA = CRA * transpose(CR);
    CRA.compress(1e-12);
//     cout << "CRA after biorthogonalization:" << endl << CRA << endl;

    CRAT.resize(ellT_r()+ell2T<d,dT>()-sT1+std::max(0,d-dT+sT1-s1)-1, CRT.row_dimension());
    for (int i = 1-ell2T<d,dT>(); i < ellT_r()-sT1; i++)
      for (unsigned int r = sT1; r < dT; r++)
	CRAT(ell2T<d,dT>()-1+i, r-sT1) = alpha(i, r);
    for (unsigned int r = dT; r < CRT.row_dimension()+sT1; r++)
      CRAT(ell2T<d,dT>()-1+ellT_r()-sT1+r-dT, r-sT1) = 1.0;

//     cout << "CRAT before biorthogonalization:" << endl << CRAT << endl;
    CRAT = CRAT * transpose(CRT);
    CRAT.compress(1e-12);
//     cout << "CRAT after biorthogonalization:" << endl << CRAT << endl;

#if 0
    cout << "CLA=" << endl << CLA << endl;
    cout << "CLAT=" << endl << CLAT << endl;
    cout << "CRA=" << endl << CRA << endl;
    cout << "CRAT=" << endl << CRAT << endl;
#endif
  }

  template <int d, int dT>
  void DSBasis<d,dT>::setup_Cj() {
    // IGPMlib reference: I_Basis_Bspline_s::setup_Cj(), ::put_Mat()

    // (5.2.5)
    Cj.diagonal(Deltasize(j0()), 1.0);
    Cj.set_block(0, 0, CL);
    Cj.set_block(Deltasize(j0())-CR.row_dimension(),
		 Deltasize(j0())-CR.column_dimension(),
		 CR, true);

    inv_Cj.diagonal(Deltasize(j0()), 1.0);
    inv_Cj.set_block(0, 0, inv_CL);
    inv_Cj.set_block(Deltasize(j0())-inv_CR.row_dimension(),
		     Deltasize(j0())-inv_CR.column_dimension(),
		     inv_CR, true);

    CjT.diagonal(Deltasize(j0()), 1.0);
    CjT.set_block(0, 0, CLT);
    CjT.set_block(Deltasize(j0())-CRT.row_dimension(),
		   Deltasize(j0())-CRT.column_dimension(),
		   CRT, true);

    inv_CjT.diagonal(Deltasize(j0()), 1.0);
    inv_CjT.set_block(0, 0, inv_CLT);
    inv_CjT.set_block(Deltasize(j0())-inv_CRT.row_dimension(),
		      Deltasize(j0())-inv_CRT.column_dimension(),
		      inv_CRT, true);

    Cjp.diagonal(Deltasize(j0()+1), 1.0);
    Cjp.set_block(0, 0, CL);
    Cjp.set_block(Deltasize(j0()+1)-CR.row_dimension(),
		  Deltasize(j0()+1)-CR.column_dimension(),
		  CR, true);

    inv_Cjp.diagonal(Deltasize(j0()+1), 1.0);
    inv_Cjp.set_block(0, 0, inv_CL);
    inv_Cjp.set_block(Deltasize(j0()+1)-inv_CR.row_dimension(),
		      Deltasize(j0()+1)-inv_CR.column_dimension(),
		      inv_CR, true);

    CjpT.diagonal(Deltasize(j0()+1), 1.0);
    CjpT.set_block(0, 0, CLT);
    CjpT.set_block(Deltasize(j0()+1)-CRT.row_dimension(),
		   Deltasize(j0()+1)-CRT.column_dimension(),
		   CRT, true);

    inv_CjpT.diagonal(Deltasize(j0()+1), 1.0);
    inv_CjpT.set_block(0, 0, inv_CLT);
    inv_CjpT.set_block(Deltasize(j0()+1)-inv_CRT.row_dimension(),
		       Deltasize(j0()+1)-inv_CRT.column_dimension(),
		       inv_CRT, true);

#if 1
    cout << "DSBasis: testing setup of Cj:" << endl;

    SparseMatrix<double> test1 = CjT * inv_CjT;
    for (unsigned int i = 0; i < test1.row_dimension(); i++)
      test1.set_entry(i, i, test1.get_entry(i, i) - 1.0);
    cout << "* ||CjT*inv_CjT-I||_infty: " << row_sum_norm(test1) << endl;

    SparseMatrix<double> test2 = Cj * inv_Cj;
    for (unsigned int i = 0; i < test2.row_dimension(); i++)
      test2.set_entry(i, i, test2.get_entry(i, i) - 1.0);
    cout << "* ||Cj*inv_Cj-I||_infty: " << row_sum_norm(test2) << endl;

    SparseMatrix<double> test3 = CjpT * inv_CjpT;
    for (unsigned int i = 0; i < test3.row_dimension(); i++)
      test3.set_entry(i, i, test3.get_entry(i, i) - 1.0);
    cout << "* ||CjpT*inv_CjpT-I||_infty: " << row_sum_norm(test3) << endl;

    SparseMatrix<double> test4 = Cjp * inv_Cjp;
    for (unsigned int i = 0; i < test4.row_dimension(); i++)
      test4.set_entry(i, i, test4.get_entry(i, i) - 1.0);
    cout << "* ||Cjp*inv_Cjp-I||_infty: " << row_sum_norm(test4) << endl;
#endif
  }

  template <int d, int dT>
  Matrix<double>
  DSBasis<d, dT>::ML() const {
    // IGPMlib reference: I_Basis_Bspline_s::ML()
    
    Matrix<double> ML(ell_l()+d-2*s0+ell2<d>()-1, d-s0);

    for (int row = s0; row < d; row++)
      ML(row-s0, row-s0) = 1.0 / sqrt(ldexp(1.0, 2*row+1));
    for (int m = ell_l()-s0; m <= 2*(ell_l()-s0)+ell1<d>()-1; m++)
      for (int k = s0; k < d; k++)
     	ML(-ell_l()+d+m, k-s0) = alphaT(m, k) / sqrt(ldexp(1.0, 2*k+1));
    for (int m = 2*(ell_l()-s0)+ell1<d>(); m <= 2*(ell_l()-s0)+ell2<d>()-2; m++)
      for (int k = s0; k < d; k++)
	ML(-ell_l()+d+m, k-s0) = betaLT(m, k);

//     cout << "ML=" << endl << ML << endl;
    
    return ML;
  }

  template <int d, int dT>
  Matrix<double>
  DSBasis<d, dT>::MR() const {
    // IGPMlib reference: I_Basis_Bspline_s::MR()

    Matrix<double> MR(ell_r()+d-2*s1+ell2<d>()-1, d-s1);

    for (int row = s1; row < d; row++)
      MR(row-s1, row-s1) = 1.0 / sqrt(ldexp(1.0, 2*row+1));
    for (int m = ell_r()-s1; m <= 2*(ell_r()-s1)+ell1<d>()-1; m++)
      for (int k = s1; k < d; k++)
     	MR(-ell_r()+d+m, k-s1) = alphaT(m, k) / sqrt(ldexp(1.0, 2*k+1));
    for (int m = 2*(ell_r()-s1)+ell1<d>(); m <= 2*(ell_r()-s1)+ell2<d>()-2; m++)
      for (int k = s1; k < d; k++)
	MR(-ell_r()+d+m, k-s1) = betaRT(m, k);

//     cout << "MR=" << endl << MR << endl;

    return MR;
  }

  template <int d, int dT>
  Matrix<double>
  DSBasis<d, dT>::MLTp() const {
    // IGPMlib reference: I_Basis_Bspline_s::MLts()

    Matrix<double> MLTp(ellT_l()+dT-2*sT0+ell2T<d,dT>()-1, dT-sT0);

    for (int row = sT0; row < dT; row++)
      MLTp(row-sT0, row-sT0) = 1.0 / sqrt(ldexp(1.0, 2*row+1));
    for (int m = ellT_l()-sT0; m <= 2*(ellT_l()-sT0)+ell1T<d,dT>()-1; m++)
      for (int k = sT0; k < dT; k++)
     	MLTp(-ellT_l()+dT+m, k-sT0) = alpha(m, k) / sqrt(ldexp(1.0, 2*k+1));
    for (int m = 2*(ellT_l()-sT0)+ell1T<d,dT>(); m <= 2*(ellT_l()-sT0)+ell2T<d,dT>()-2; m++)
      for (int k = sT0; k < dT; k++)
	MLTp(-ellT_l()+dT+m, k-sT0) = betaL(m, k);

//     cout << "MLTp=" << endl << MLTp << endl;

    return MLTp;
  }

  template <int d, int dT>
  Matrix<double>
  DSBasis<d, dT>::MRTp() const {
    // IGPMlib reference: I_Basis_Bspline_s::MRts()

    Matrix<double> MRTp(ellT_r()+dT-2*sT1+ell2T<d,dT>()-1, dT-sT1);

    for (int row = sT1; row < dT; row++)
      MRTp(row-sT1, row-sT1) = 1.0 / sqrt(ldexp(1.0, 2*row+1));
    for (int m = ellT_r()-sT1; m <= 2*(ellT_r()-sT1)+ell1T<d,dT>()-1; m++)
      for (int k = sT1; k < dT; k++)
     	MRTp(-ellT_r()+dT+m, k-sT1) = alpha(m, k) / sqrt(ldexp(1.0, 2*k+1));
    for (int m = 2*(ellT_r()-sT1)+ell1T<d,dT>(); m <= 2*(ellT_r()-sT1)+ell2T<d,dT>()-2; m++)
      for (int k = sT1; k < dT; k++)
	MRTp(-ellT_r()+dT+m, k-sT1) = betaR(m, k);

//     cout << "MRTp=" << endl << MRTp << endl;
 
    return MRTp;
  }

  template <int d, int dT>
  void
  DSBasis<d,dT>::setup_Mj0(const Matrix<double>& ML, const Matrix<double>& MR, SparseMatrix<double>& Mj0) {
    // IGPMlib reference: I_Basis_Bspline_s::Mj0()
    
    // TODO: enhance readability! (<-> [DKU section 3.5])

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
	Mj0.set_entry(row, col, M_SQRT1_2 * *it);
    }

//     cout << "Mj0=" << endl << Mj0 << endl;
  }
  
  template <int d, int dT>
  void
  DSBasis<d,dT>::setup_Mj0Tp(const Matrix<double>& MLTp, const Matrix<double>& MRTp, SparseMatrix<double>& Mj0Tp) {
    // IGPMlib reference: I_Basis_Bspline_s::Mj0ts()

    const int nj  = Deltasize(j0());
    const int njp = Deltasize(j0()+1);
    Mj0Tp.resize(njp, nj);
    
    for (unsigned int i = 0; i < MLTp.row_dimension(); i++)
      for (unsigned int k = 0; k < MLTp.column_dimension(); k++)
 	Mj0Tp.set_entry(i, k, MLTp.get_entry(i, k));
    
    for (unsigned int i = 0; i < MRTp.row_dimension(); i++)
      for (unsigned int k = 0; k < MRTp.column_dimension(); k++)
 	Mj0Tp.set_entry(njp-i-1, nj-k-1, MRTp.get_entry(i, k));
    
    int startrow = dT+ellT_l()+ell1T<d,dT>()-2*sT0;
    for (int col = dT-sT0; col < nj-(dT-sT1); col++, startrow+=2) {
      int row = startrow;
      for (MultivariateLaurentPolynomial<double, 1>::const_iterator it(cdf.aT().begin());
	   it != cdf.aT().end(); ++it, row++)
	Mj0Tp.set_entry(row, col, M_SQRT1_2 * *it);
    }
    
//     cout << "Mj0Tp=" << endl << Mj0Tp << endl;
  }

  template <int d, int dT>
  void
  DSBasis<d, dT>::F(SparseMatrix<double>& FF) {
    // IGPMlib reference: I_Basis_Bspline_s::F()
    
    const int FLow = ell_l()-s0+(d%2);       // start column index for F_j in (4.1.14)
    const int FUp  = FLow+(DeltaRmin(j0())-DeltaLmax())-1; // end column index for F_j in (4.1.14)
    
    // (4.1.14):
    
    FF.resize(Deltasize(j0()+1), 1<<j0());

    for (int r = 0; r < FLow; r++)
      FF.set_entry(r+d-s0, r, 1.0);
    
    int i = d+ell_l()+(d%2)-1-2*s0;
    for (int r = FLow; r <= FUp; r++) {
      FF.set_entry(i, r-1, 1.0);
      i += 2;
    } 
    
    i = Deltasize(j0()+1)-d-1+s1;
    for (int r = (1<<j0()); r >= FUp+1; r--) {
      FF.set_entry(i, r-1, 1.0);
      i--;
    }

//     cout << "F=" << endl << FF << endl;
  }

  template <int d, int dT>
  void
  DSBasis<d, dT>::P(const Matrix<double>& ML, const Matrix<double>& MR, SparseMatrix<double>& PP) {
    // IGPMlib reference: I_Basis_Bspline_s::P()
    
    // (4.1.22):

    PP.diagonal(Deltasize(j0()+1), 1.0);
    
    for (unsigned int i = 0; i < ML.row_dimension(); i++)
      for (unsigned int k = 0; k < ML.column_dimension(); k++)
	PP.set_entry(i, k, ML.get_entry(i, k));

    for (unsigned int i = 0; i < MR.row_dimension(); i++)
      for (unsigned int k = 0; k < MR.column_dimension(); k++)
	PP.set_entry(Deltasize(j0()+1)-i-1, Deltasize(j0()+1)-k-1, MR.get_entry(i, k));

//     cout << "P=" << endl << PP << endl;
  }

  template <int d, int dT>
  void
  DSBasis<d, dT>::GSetup(SparseMatrix<double>& A, SparseMatrix<double>& H, SparseMatrix<double>& Hinv) {
    // IGPMlib reference: I_Basis_Bspline_s::GSetup()

    // (4.1.13):
    
    const int nj  = Deltasize(j0());
    const int njp = Deltasize(j0()+1);
    A.resize(njp, nj);

    for (int row = 0; row < d-s0; row++)
      A.set_entry(row, row, 1.0);
    
    int startrow = d+ell_l()+ell1<d>()-2*s0;
    for (int col = d-s0; col < nj-(d-s1); col++, startrow+=2) {
      int row = startrow;
      for (MultivariateLaurentPolynomial<double, 1>::const_iterator it(cdf.a().begin());
 	   it != cdf.a().end(); ++it, row++) {
 	A.set_entry(row, col, M_SQRT1_2 * *it);
      }
    }
    
    for (int row = njp-1, col = nj-1; col > nj-1-(d-s1); row--, col--)
      A.set_entry(row, col, 1.0);

    // prepare H, Hinv for elimination process:
    H   .diagonal(njp, 1.0);
    Hinv.diagonal(njp, 1.0);

//     cout << "A=" << endl << A << endl;
  }

  template <int d, int dT>
  void
  DSBasis<d, dT>::GElim(SparseMatrix<double>& A, SparseMatrix<double>& H, SparseMatrix<double>& Hinv) {
    // IGPMlib reference: I_Basis_Bspline_s::gelim()
    
    // A_j=A_j^{(0)} in (4.1.1) is a q times p matrix
    const int firstcol = d-s0; // first column of A_j^{(d)} in Ahat_j^{(d)}
    const int lastcol  = (Deltasize(j0())-1)-(d-s1); // last column
    const int firstrow = d+ell_l()+ell1<d>()-2*s0; // first row of A_j^{(d)} in Ahat_j^{(d)}
    const int lastrow  = (Deltasize(j0()+1)-1)-(d+ell_r()-ell2<d>()+(d%2)-2*s1); // last row
    
    int p = lastcol-firstcol+1;
//     int q = lastrow-firstrow+1;

    SparseMatrix<double> help;
    
    // elimination (4.1.4)ff.:
    for (int i = 1; i <= d; i++) {
      help.diagonal(Deltasize(j0()+1), 1.0);

      const int elimrow = i%2 ? firstrow+(i-1)/2 : lastrow-(int)floor((i-1)/2.);

      const int HhatLow = i%2 ? elimrow : ell_l()+ell2<d>()+2-(d%2)-(i/2)-2*s0;
      const int HhatUp  = i%2 ? HhatLow + 2*p-1+(d+(d%2))/2 : elimrow;
      
      if (i%2) // i odd, elimination from above (4.1.4a)
	{
	  assert(fabs(A.get_entry(elimrow+1, firstcol)) >= 1e-10);
	  const double Uentry = -A.get_entry(elimrow, firstcol) / A.get_entry(elimrow+1, firstcol);
	  
	  // insert Uentry in Hhat
	  for (int k = HhatLow; k <= HhatUp; k += 2)
	    help.set_entry(k, k+1, Uentry);
	}
      else // i even, elimination from below (4.1.4b)
	{
	  assert(fabs(A.get_entry(elimrow-1, lastcol)) >= 1e-10);
	  const double Lentry = -A.get_entry(elimrow, lastcol) / A.get_entry(elimrow-1, lastcol);
	  
	  // insert Lentry in Hhat
	  for (int k = HhatLow; k <= HhatUp; k += 2)
	    help.set_entry(k+1, k, Lentry);
	}
      
      A = help * A;
      H = help * H;
      
      A.compress(1e-14);

      // invert help
      if (i%2) {
	for (int k = HhatLow; k <= HhatUp; k += 2)
	  help.set_entry(k, k+1, -help.get_entry(k, k+1));
      }	else {
	for (int k = HhatLow; k <= HhatUp; k += 2)
	  help.set_entry(k+1, k, -help.get_entry(k+1, k));
      }
      
      Hinv = Hinv * help;
    }
  }

  template <int d, int dT>
  void
  DSBasis<d, dT>::BT(const SparseMatrix<double>& A, SparseMatrix<double>& BB) {
    // IGPMlib reference: I_Basis_Bspline_s::Btr()
    
    const int p = (1<<j0()) - ell_l() - ell_r() - (d%2) + 1;
    const int llow = ell_l()-d;
    
    BB.resize(Deltasize(j0()+1), Deltasize(j0()));

    for (int r = 0; r < d-s0; r++)
      BB.set_entry(r, r, 1.0);

    const double help = 1./A.get_entry(d+ell_l()+ell1<d>()+ell2<d>(), d);

    for (int c = d-s0, r = d+ell_l()+ell1<d>()+ell2<d>(); c < d+p+s1; c++, r += 2)
      BB.set_entry(r-2*s0, c, help);

    for (int r = DeltaRmax(j0())-d+1+s1; r <= DeltaRmax(j0()); r++)
      BB.set_entry(-llow+r+DeltaRmax(j0()+1)-DeltaRmax(j0()), -llow+r, 1.0);
  }

  template <int d, int dT>
  void
  DSBasis<d, dT>::InvertP(const SparseMatrix<double>& PP, SparseMatrix<double>& PPinv) {
    // IGPMlib reference: I_Basis_Bspline_s::InverseP()
    
    PPinv.diagonal(PP.row_dimension(), 1.0);
    
    const int msize_l = d+ell_l()+ell2<d>()-1;
    const int msize_r = d+ell_r()+ell2<d>()-1;
    
    Matrix<double> ml;
    ml.diagonal(msize_l, 1.0);
    for (int i = 0; i < msize_l; i++)
      for (int k = 0; k <= d; k++)
	ml.set_entry(i, k, PP.get_entry(i, k));
    
    Matrix<double> mr;
    mr.diagonal(msize_r, 1.0);
    for (int i = 0; i < msize_r; i++)
      for (int k = 0; k <= d; k++)
	mr.set_entry(i, msize_r-d-1+k, PP.get_entry(PP.row_dimension()-msize_r+i, PP.column_dimension()-d-1+k));
    
    Matrix<double> mlinv, mrinv;
    QUDecomposition<double>(ml).inverse(mlinv);
    QUDecomposition<double>(mr).inverse(mrinv);
    
    for (int i = 0; i < msize_l; i++)
      for (int k = 0; k <= d; k++)
	PPinv.set_entry(i, k, mlinv.get_entry(i, k));
    for (int i = 0; i < msize_r; i++)
      for (int k = 0; k <= d; k++) {
	PPinv.set_entry(PP.row_dimension()-msize_r+i, PP.column_dimension()-d-1+k, mrinv.get_entry(i, msize_r-d-1+k));
      }
  }
  
  template <int d, int dT>
  void
  DSBasis<d, dT>::DS_symmetrization(SparseMatrix<double>& Mj1, SparseMatrix<double>& Mj1T) {
    // IGPMlib reference: I_Basis_Bspline::Modify()
    
    SparseMatrix<double> Hj1(Deltasize(j0()+1), 1<<j0()),
      Hj1T(Deltasize(j0()+1), 1<<j0());
    
    // copy left halves of Mj1, Mj1T, right halves are mirrored
    for (int i = 0; i < Deltasize(j0()+1); i++)
      for (int j = 0; j < 1<<(j0()-1); j++) {
	double help = Mj1.get_entry(i, j);
	if (help != 0) {
	  Hj1.set_entry(i, j, help);
	  Hj1.set_entry(Deltasize(j0()+1)-1-i, (1<<j0())-1-j, help);
	}
	
	help = Mj1T.get_entry(i, j);
	if (help != 0) {
	  Hj1T.set_entry(i, j, help);
	  Hj1T.set_entry(Deltasize(j0()+1)-1-i, (1<<j0())-1-j, help);
	}
      }
    
    SparseMatrix<double> Kj = transpose(Hj1T)*Hj1; Kj.compress(1e-12);
    Matrix<double> Kj_full; Kj.get_block(0, 0, Kj.row_dimension(), Kj.column_dimension(), Kj_full);
    Matrix<double> invKj_full;
    QUDecomposition<double>(Kj_full).inverse(invKj_full);
    SparseMatrix<double> invKj(Kj.row_dimension(), Kj.column_dimension());
    invKj.set_block(0, 0, invKj_full);
    
    Hj1 = Hj1 * invKj;
    Hj1.compress();
    
    Mj1  = Hj1;
    Mj1T = Hj1T;
  }

  template <int d, int dT>
  void
  DSBasis<d, dT>::decompose(const InfiniteVector<double, Index>& c,
			     const int jmin,
			     InfiniteVector<double, Index>& v) const {
    v.clear();
    for (typename InfiniteVector<double, Index>::const_iterator it(c.begin()), itend(c.end());
	 it != itend; ++it) {
      InfiniteVector<double, Index> help;
      decompose_1(it.index(), jmin, help);
      v += *it * help;
    }
  }
  
  template <int d, int dT>
  void
  DSBasis<d, dT>::decompose_t(const InfiniteVector<double, Index>& c,
			       const int jmin,
			       InfiniteVector<double, Index>& v) const {
    v.clear();
    for (typename InfiniteVector<double, Index>::const_iterator it(c.begin()), itend(c.end());
	 it != itend; ++it) {
      InfiniteVector<double, Index> help;
      decompose_t_1(it.index(), jmin, help);
      v += *it * help;
    }
  }
  
  template <int d, int dT>
  void
  DSBasis<d, dT>::reconstruct(const InfiniteVector<double, Index>& c,
			       const int j,
			       InfiniteVector<double, Index>& v) const {
    v.clear();
    for (typename InfiniteVector<double, Index>::const_iterator it(c.begin()), itend(c.end());
	 it != itend; ++it) {
      InfiniteVector<double, Index> help;
      reconstruct_1(it.index(), j, help);
      v += *it * help;
    }
  }

  template <int d, int dT>
  void
  DSBasis<d, dT>::reconstruct_t(const InfiniteVector<double, Index>& c,
				 const int j,
				 InfiniteVector<double, Index>& v) const {
    v.clear();
    for (typename InfiniteVector<double, Index>::const_iterator it(c.begin()), itend(c.end());
	 it != itend; ++it) {
      InfiniteVector<double, Index> help;
      reconstruct_t_1(it.index(), j, help);
      v += *it * help;
    }
  }

  template <int d, int dT>
  void
  DSBasis<d, dT>::decompose_1(const Index& lambda,
			       const int jmin,
			       InfiniteVector<double, Index>& c) const {
    assert(jmin >= j0());
    assert(lambda.j() >= jmin);
    
    c.clear();

    if (lambda.e() == 1) // wavelet
      c.set_coefficient(lambda, 1.0); // true wavelet coefficients don't have to be modified
    else // generator
      {
	if (lambda.j() == jmin) // generators on the coarsest level don't have to be modified
	  c.set_coefficient(lambda, 1.0);
	else // j > jmin
	  {
	    // For the multiscale decomposition of psi_lambda, we have to compute
	    // the corresponding column of the transformation matrix G_{j-1}=\tilde M_{j-1}^T,
	    // i.e. one row of G_{j-1}^T=(\tilde M_{j-1,0}, \tilde M_{j-1,1}).
	    
	    typedef Vector<double>::size_type size_type;

	    const size_type row = lambda.k() - DeltaLmin();

  	    // compute d_{j-1}
	    size_type row_j0 = row;
	    size_type offset = 0;

 	    if (lambda.j()-1 == j0()) {		
	      for (size_type k(0); k < Mj1T.entries_in_row(row_j0); k++)
		c.set_coefficient(Index(j0(), 1, Mj1T.get_nth_index(row_j0,k), this),
				  Mj1T.get_nth_entry(row_j0,k));
	    } else {
	      // Due to the [DS] symmetrization, we have to be a bit careful here.
	      const size_type third = Deltasize(j0()+1)*2/5;
	      if (row < third)
		for (size_type k(0); k < Mj1T.entries_in_row(row_j0); k++)
		  c.set_coefficient(Index(lambda.j()-1, 1, Mj1T.get_nth_index(row_j0,k), this),
				    Mj1T.get_nth_entry(row_j0,k));
	      else {
		const size_type bottom_third = Deltasize(lambda.j())-Deltasize(j0()+1)*2/5;
		if (row >= bottom_third) {
		  row_j0 = row+Deltasize(j0()+1)-Deltasize(lambda.j());
		  offset = (1<<(lambda.j()-1))-(1<<j0());
		  for (size_type k(0); k < Mj1T.entries_in_row(row_j0); k++)
		    c.set_coefficient(Index(lambda.j()-1, 1, Mj1T.get_nth_index(row_j0,k)+offset, this),
				      Mj1T.get_nth_entry(row_j0,k));
		} else {
		  // Left half of Mj1T:
		  
		  const int col_left = (1<<(j0()-1))-1;
		  
		  // The row ...
		  const int first_row_left = Mj1T_t.get_nth_index(col_left, 0);
		  // ... is the first nontrivial row in column (1<<(j0-1))-1 and the row ...
		  const int last_row_left  = Mj1T_t.get_nth_index(col_left, Mj1T_t.entries_in_row(col_left)-1);
		  // ... is the last nontrivial one,
		  // i.e. the row last_row_left begins at column (1<<(j0()-1))-1, so does the row last_row_left-1.
		  // So the row "row" starts at column ...
		  const int first_column_left = (1<<(j0()-1))-1+(int)floor(((int)row+1-last_row_left)/2.);
		  
		  for (int col = first_column_left, filter_row = last_row_left-abs(row-last_row_left)%2;
		       col < (1<<(lambda.j()-2)) && filter_row >= first_row_left; col++, filter_row -= 2)
		    c.set_coefficient(Index(lambda.j()-1, 1, col, this),
				      Mj1T_t.get_nth_entry(col_left, filter_row-first_row_left));
		  
		  // Analogous strategy for the right half:
		  
		  const int col_right = 1<<(j0()-1);
		  
		  const int offset_right = (Deltasize(lambda.j())-Deltasize(j0()+1))-(1<<(lambda.j()-1))+(1<<j0()); // row offset for the right half
		  const int first_row_right = Mj1T_t.get_nth_index(col_right, 0)+offset_right;
		  const int last_row_right  = Mj1T_t.get_nth_index(col_right, Mj1T_t.entries_in_row(col_right)-1)+offset_right;
		  
		  // The rows first_row_right and first_row_right+1 end at column 1<<(lambda.j()-2),
		  // so the row "row" ends at column ...
		  const int last_column_right = (1<<(lambda.j()-2))+(int)floor(((int)row-first_row_right)/2.);
		  
		  for (int col = last_column_right, filter_row = first_row_right-offset_right+abs(row-first_row_right)%2;
		       col >= 1<<(lambda.j()-2) && filter_row <= last_row_right-offset_right; col--, filter_row += 2)
		    c.set_coefficient(Index(lambda.j()-1, 1, col, this),
				      Mj1T_t.get_nth_entry(col_right, filter_row+offset_right-first_row_right));
		}
	      }	    
	    }
	    
   	    // compute c_{jmin} via recursion
	    row_j0 = row;
	    offset = 0;
	    if (lambda.j()-1 != j0()) {
	      const size_type rows_top = (int)ceil(Deltasize(j0()+1)/2.0);
	      if (row >= rows_top) {
		const size_type bottom = Deltasize(lambda.j())-Deltasize(j0()+1)/2;
		if (row >= bottom) {
		  row_j0 = row + rows_top - bottom;
		  offset = Deltasize(lambda.j()-1) - Deltasize(j0());
		} else {
		  row_j0 = rows_top-2+(row-rows_top)%2;
		  offset = (row-rows_top)/2+1;
		}
	      }
	    }
	    for (size_type k(0); k < Mj0T.entries_in_row(row_j0); k++) {
	      InfiniteVector<double, Index> dhelp;
	      decompose_1(Index(lambda.j()-1, 0, DeltaLmin()+Mj0T.get_nth_index(row_j0,k)+offset, this), jmin, dhelp);
	      c += Mj0T.get_nth_entry(row_j0,k) * dhelp;
	    }
	  }
      }
  }
  
  template <int d, int dT>
  void
  DSBasis<d, dT>::decompose_t_1(const Index& lambda,
				 const int jmin,
				 InfiniteVector<double, Index>& c) const {
    assert(jmin >= j0());
    assert(lambda.j() >= jmin);
    
    c.clear();
    
    if (lambda.e() == 1) // wavelet
      c.set_coefficient(lambda, 1.0); // true wavelet coefficients don't have to be modified
    else // generator
      {
	if (lambda.j() == jmin)
	  c.set_coefficient(lambda, 1.0);
	else // j > jmin
	  {
	    // For the multiscale decomposition of psi_lambda, we have to compute
	    // the corresponding column of the transformation matrix \tilde G_{j-1}=M_{j-1}^T,
	    // i.e. one row of \tilde G_{j-1}^T=(M_{j-1,0}, M_{j-1,1}).
	    
	    typedef Vector<double>::size_type size_type;
	    
	    const size_type row = lambda.k() - DeltaLmin();
	    
  	    // compute d_{j-1}
	    size_type row_j0 = row;
	    size_type offset = 0;
	    
 	    if (lambda.j()-1 == j0()) {		
	      for (size_type k(0); k < Mj1.entries_in_row(row_j0); k++)
		c.set_coefficient(Index(j0(), 1, Mj1.get_nth_index(row_j0,k), this),
				  Mj1.get_nth_entry(row_j0,k));
	    } else {
	      // Due to the [DS] symmetrization, we have to be a bit careful here.
	      const size_t third = Deltasize(j0()+1)*2/5;
	      if (row < third)
		for (size_type k(0); k < Mj1.entries_in_row(row_j0); k++)
		  c.set_coefficient(Index(lambda.j()-1, 1, Mj1.get_nth_index(row_j0,k), this),
				    Mj1.get_nth_entry(row_j0,k));
	      else {
		const size_t bottom_third = Deltasize(lambda.j())-Deltasize(j0()+1)*2/5;
		if (row >= bottom_third) {
		  row_j0 = row+Deltasize(j0()+1)-Deltasize(lambda.j());
		  offset = (1<<(lambda.j()-1))-(1<<j0());
		  for (size_type k(0); k < Mj1.entries_in_row(row_j0); k++)
		    c.set_coefficient(Index(lambda.j()-1, 1, Mj1.get_nth_index(row_j0,k)+offset, this),
				      Mj1.get_nth_entry(row_j0,k));
		} else {
		  // Left half of Mj1:
		  
		  const int col_left = (1<<(j0()-1))-1;
		  
		  // The row ...
		  const int first_row_left = Mj1_t.get_nth_index(col_left, 0);
		  // ... is the first nontrivial row in column (1<<(j0-1))-1 and the row ...
		  const int last_row_left  = Mj1_t.get_nth_index(col_left, Mj1_t.entries_in_row(col_left)-1);
		  // ... is the last nontrivial one,
		  // i.e. the row last_row_left begins at column (1<<(j0()-1))-1, so does the row last_row_left-1.
		  // So the row "row" starts at column ...
		  const int first_column_left = (1<<(j0()-1))-1+(int)floor(((int)row+1-last_row_left)/2.);
		  
		  for (int col = first_column_left, filter_row = last_row_left-abs(row-last_row_left)%2;
		       col < (1<<(lambda.j()-2)) && filter_row >= first_row_left; col++, filter_row -= 2)
		    c.set_coefficient(Index(lambda.j()-1, 1, col, this),
				      Mj1_t.get_nth_entry(col_left, filter_row-first_row_left));
		  
		  // Analogous strategy for the right half:
		  
		  const int col_right = 1<<(j0()-1);
		  
		  const int offset_right = (Deltasize(lambda.j())-Deltasize(j0()+1))-(1<<(lambda.j()-1))+(1<<j0()); // row offset for the right half
		  const int first_row_right = Mj1_t.get_nth_index(col_right, 0)+offset_right;
		  const int last_row_right  = Mj1_t.get_nth_index(col_right, Mj1_t.entries_in_row(col_right)-1)+offset_right;
		  
		  // The rows first_row_right and first_row_right+1 end at column 1<<(lambda.j()-2),
		  // so the row "row" ends at column ...
		  const int last_column_right = (1<<(lambda.j()-2))+(int)floor(((int)row-first_row_right)/2.);
		  
		  for (int col = last_column_right, filter_row = first_row_right-offset_right+abs(row-first_row_right)%2;
		       col >= 1<<(lambda.j()-2) && filter_row <= last_row_right-offset_right; col--, filter_row += 2)
		    c.set_coefficient(Index(lambda.j()-1, 1, col, this),
				      Mj1_t.get_nth_entry(col_right, filter_row+offset_right-first_row_right));
		}
	      }	    
	    }
	    
   	    // compute c_{jmin} via recursion
	    row_j0 = row;
	    offset = 0;
	    if (lambda.j()-1 != j0()) {
	      const size_type rows_top = (int)ceil(Deltasize(j0()+1)/2.0);
	      if (row >= rows_top) {
		const size_type bottom = Deltasize(lambda.j())-Deltasize(j0()+1)/2;
		if (row >= bottom) {
		  row_j0 = row + rows_top - bottom;
		  offset = Deltasize(lambda.j()-1) - Deltasize(j0());
		} else {
		  row_j0 = rows_top-2+(row-rows_top)%2;
		  offset = (row-rows_top)/2+1;
		}
	      }
	    }
	    for (size_type k(0); k < Mj0.entries_in_row(row_j0); k++) {
	      InfiniteVector<double, Index> dhelp;
	      decompose_t_1(Index(lambda.j()-1, 0, DeltaLmin()+Mj0.get_nth_index(row_j0,k)+offset, this), jmin, dhelp);
	      c += Mj0.get_nth_entry(row_j0,k) * dhelp;
	    }
	  }
      }
  }
  
  template <int d, int dT>
  void
  DSBasis<d, dT>::reconstruct_1(const Index& lambda,
				 const int j,
				 InfiniteVector<double, Index>& c) const {
    c.clear();
    
    if (lambda.j() >= j)
      c.set_coefficient(lambda, 1.0); 
    else {
      // For the reconstruction of psi_lambda, we have to compute
      // the corresponding column of the transformation matrix Mj=(Mj0, Mj1).
      
      // reconstruct by recursion
      typedef Vector<double>::size_type size_type;
      
      const size_type row = (lambda.e() == 0 ? lambda.k() - DeltaLmin() : lambda.k());
      
      const SparseMatrix<double>& M = (lambda.e() == 0 ? Mj0_t : Mj1_t);
      
      size_type row_j0 = row;
      size_type offset = 0;
      if (lambda.e() == 0) {
	if (lambda.j() > j0()) {
	  const size_type rows_top = (int)ceil(Deltasize(j0())/2.0);
	  if (row >= rows_top) {
	    const size_type bottom = Deltasize(lambda.j())-Deltasize(j0())/2;
	    if (row >= bottom) {
	      row_j0 = row+rows_top-bottom;
	      offset = Deltasize(lambda.j()+1)-Deltasize(j0()+1);
	    } else {
	      row_j0 = rows_top-1;
	      offset = 2*(row-rows_top)+2;
	    }
	  }
	}
      }
      else {
	if (lambda.j() > j0()) {
	  const size_type rows_top = 1<<(j0()-1);
	  if (row >= rows_top) {
	    const size_type bottom = (1<<lambda.j())-(1<<(j0()-1));
	    if (row >= bottom) {
	      row_j0 = row+rows_top-bottom;
	      offset = Deltasize(lambda.j()+1)-Deltasize(j0()+1);
	    } else {
	      if ((int)row < (1<<(lambda.j()-1))) {
		row_j0 = rows_top-1;
		offset = 2*(row-rows_top)+2;
	      } else {
		row_j0 = 1<<(j0()-1);
		offset = Deltasize(lambda.j()+1)-Deltasize(j0()+1)+2*((int)row-bottom);
	      }
	    }
	  }
	}
      }
      for (size_type k(0); k < M.entries_in_row(row_j0); k++) {
	InfiniteVector<double, Index> dhelp;
	reconstruct_1(Index(lambda.j()+1, 0, DeltaLmin()+M.get_nth_index(row_j0,k)+offset, this), j, dhelp);
	c += M.get_nth_entry(row_j0,k) * dhelp;
      }
    }
  }
  
  template <int d, int dT>
  void
  DSBasis<d, dT>::reconstruct_t_1(const Index& lambda,
				   const int j,
				   InfiniteVector<double, Index>& c) const {
    c.clear();
    
    if (lambda.j() >= j)
      c.set_coefficient(lambda, 1.0);
    else {
      // For the reconstruction of psi_lambda, we have to compute
      // the corresponding column of the transformation matrix \tilde Mj=(\tilde Mj0, \tilde Mj1).
      
      // reconstruct by recursion
      typedef Vector<double>::size_type size_type;
      
      const size_type row = (lambda.e() == 0 ? lambda.k() - DeltaLmin() : lambda.k());
      
      const SparseMatrix<double>& M = (lambda.e() == 0 ? Mj0T_t : Mj1T_t);
      
      size_type row_j0 = row;
      size_type offset = 0;
      if (lambda.e() == 0) {
	if (lambda.j() > j0()) {
	  const size_type rows_top = (int)ceil(Deltasize(j0())/2.0);
	  if (row >= rows_top) {
	    const size_type bottom = Deltasize(lambda.j())-Deltasize(j0())/2;
	    if (row >= bottom) {
	      row_j0 = row+rows_top-bottom;
	      offset = Deltasize(lambda.j()+1)-Deltasize(j0()+1);
	    } else {
	      row_j0 = rows_top-1;
	      offset = 2*(row-rows_top)+2;
	    }
	  }
	}
      } else {
	if (lambda.j() > j0()) {
	  const size_type rows_top = 1<<(j0()-1);
	  if (row >= rows_top) {
	    const size_type bottom = (1<<lambda.j())-(1<<(j0()-1));
	    if (row >= bottom) {
	      row_j0 = row+rows_top-bottom;
	      offset = Deltasize(lambda.j()+1)-Deltasize(j0()+1);
	    } else {
	      if ((int)row < (1<<(lambda.j()-1))) {
		row_j0 = rows_top-1;
		offset = 2*(row-rows_top)+2;
	      } else {
		row_j0 = 1<<(j0()-1);
		offset = Deltasize(lambda.j()+1)-Deltasize(j0()+1)+2*((int)row-bottom);
	      }
	    }
	  }
	}
      }
      for (size_type k(0); k < M.entries_in_row(row_j0); k++) {
	InfiniteVector<double, Index> dhelp;
	reconstruct_t_1(Index(lambda.j()+1, 0, DeltaLmin()+M.get_nth_index(row_j0,k)+offset, this), j, dhelp);
	c += M.get_nth_entry(row_j0,k) * dhelp;
      }
    }
  }

  template <int d, int dT>
  void
  DSBasis<d, dT>::assemble_Mj0(const int j, SparseMatrix<double>& mj0) const {
    if (j == j0())
      mj0 = Mj0;
    else {
      mj0.resize(Deltasize(j+1), Deltasize(j));
      
      const int rows       = Deltasize(j0()+1);
      const int cols_left  = (int)ceil(Deltasize(j0())/2.0);
      const int cols_right = Deltasize(j0()) - cols_left;
      
      // upper left block
      for (int row = 0; row < rows; row++)
	for (int col = 0; col < cols_left; col++)
	  mj0.set_entry(row, col, Mj0.get_entry(row, col));
      
      // lower right block
      for (int row = 0; row < rows; row++)
	for (int col = 0; col < cols_right; col++)
	  mj0.set_entry(Deltasize(j+1)-rows+row, Deltasize(j)-cols_right+col,
			Mj0.get_entry(row, col + cols_left));
      
      // central bands
      InfiniteVector<double, Vector<double>::size_type> v;
      Mj0_t.get_row(cols_left-1, v); // last column of left half
      for (typename InfiniteVector<double, Vector<double>::size_type>::const_iterator it(v.begin());
	   it != v.end(); ++it)
	for (int col = cols_left; col <= Deltasize(j) - cols_right - 1; col++)
	  mj0.set_entry(it.index()+2*(col-cols_left)+2, col, *it);
      
      mj0.compress();
    }
  }

  template <int d, int dT>
  void
  DSBasis<d, dT>::assemble_Mj0_t(const int j, SparseMatrix<double>& mj0_t) const {
    if (j == j0())
      mj0_t = Mj0_t;
    else {
      SparseMatrix<double> help;
      assemble_Mj0(j, help);
      mj0_t = transpose(help);
    }
  }

  template <int d, int dT>
  void
  DSBasis<d, dT>::assemble_Mj0T(const int j, SparseMatrix<double>& mj0T) const {
    if (j == j0())
      mj0T = Mj0T;
    else {
      mj0T.resize(Deltasize(j+1), Deltasize(j));
      
      const int rows       = Deltasize(j0()+1);
      const int cols_left  = (int)ceil(Deltasize(j0())/2.0);
      const int cols_right = Deltasize(j0()) - cols_left;
      
      // upper left block
      for (int row = 0; row < rows; row++)
	for (int col = 0; col < cols_left; col++)
	  mj0T.set_entry(row, col, Mj0T.get_entry(row, col));
      
      // lower right block
      for (int row = 0; row < rows; row++)
	for (int col = 0; col < cols_right; col++)
	  mj0T.set_entry(Deltasize(j+1)-rows+row, Deltasize(j)-cols_right+col,
			 Mj0T.get_entry(row, col + cols_left));
      
      // central bands
      InfiniteVector<double, Vector<double>::size_type> v;
      Mj0T_t.get_row(cols_left-1, v); // last column of left half
      for (typename InfiniteVector<double, Vector<double>::size_type>::const_iterator it(v.begin());
	   it != v.end(); ++it)
	for (int col = cols_left; col <= Deltasize(j) - cols_right - 1; col++)
	  mj0T.set_entry(it.index()+2*(col-cols_left)+2, col, *it);
      
      mj0T.compress();
    }
  }
  
  template <int d, int dT>
  void
  DSBasis<d, dT>::assemble_Mj0T_t(const int j, SparseMatrix<double>& mj0T_t) const {
    if (j == j0())
      mj0T_t = Mj0T_t;
    else {
      SparseMatrix<double> help;
      assemble_Mj0T(j, help);
      mj0T_t = transpose(help);
    }
  }

  template <int d, int dT>
  void
  DSBasis<d, dT>::assemble_Mj1(const int j, SparseMatrix<double>& mj1) const {
    if (j == j0())
      mj1 = Mj1;
    else {
      mj1.resize(Deltasize(j+1), 1<<j);
      
      const int rows       = Deltasize(j0()+1);
      const int cols_left  = 1<<(j0()-1);
      const int cols_right = cols_left;
      
      // upper left block
      for (int row = 0; row < rows; row++)
	for (int col = 0; col < cols_left; col++)
	  mj1.set_entry(row, col, Mj1.get_entry(row, col));
      
      // lower right block
      for (int row = 0; row < rows; row++)
	for (int col = 0; col < cols_right; col++)
  	    mj1.set_entry(Deltasize(j+1)-rows+row, (1<<j)-cols_right+col,
 			  Mj1.get_entry(row, col + cols_left));
      
      // central bands, take care of [DS] symmetrization
      InfiniteVector<double, Vector<double>::size_type> v;
      Mj1_t.get_row(cols_left-1, v); // last column of left half
      for (typename InfiniteVector<double, Vector<double>::size_type>::const_iterator it(v.begin());
	   it != v.end(); ++it)
	for (int col = cols_left; col < 1<<(j-1); col++)
	  mj1.set_entry(it.index()+2*(col-cols_left)+2, col, *it);
      
      Mj1_t.get_row(cols_left, v); // first column of right half
      const int offset_right = (Deltasize(j+1)-Deltasize(j0()+1))-(1<<j)+(1<<j0()); // row offset for the right half
      for (typename InfiniteVector<double, Vector<double>::size_type>::const_iterator it(v.begin());
	   it != v.end(); ++it)
	for (int col = 1<<(j-1); col <= (1<<j) - cols_right - 1; col++)
	  mj1.set_entry(it.index()+offset_right+2*(col-(1<<(j-1))), col, *it);
      
      mj1.compress();
    }
  }
  
  template <int d, int dT>
  void
  DSBasis<d, dT>::assemble_Mj1_t(const int j, SparseMatrix<double>& mj1_t) const {
    if (j == j0())
      mj1_t = Mj1_t;
    else {
      SparseMatrix<double> help;
      assemble_Mj1(j, help);
      mj1_t = transpose(help);
    }
  }

  template <int d, int dT>
  void
  DSBasis<d, dT>::assemble_Mj1T(const int j, SparseMatrix<double>& mj1T) const {
    if (j == j0())
      mj1T = Mj1T;
    else {
      mj1T.resize(Deltasize(j+1), 1<<j);
      
      const int rows       = Deltasize(j0()+1);
      const int cols_left  = 1<<(j0()-1);
      const int cols_right = cols_left;
      
      // upper left block
      for (int row = 0; row < rows; row++)
	for (int col = 0; col < cols_left; col++)
	  mj1T.set_entry(row, col, Mj1T.get_entry(row, col));
      
      // lower right block
      for (int row = 0; row < rows; row++)
	for (int col = 0; col < cols_right; col++)
	  mj1T.set_entry(Deltasize(j+1)-rows+row, (1<<j)-cols_right+col,
			 Mj1T.get_entry(row, col + cols_left));
      
      // central bands, take care of [DS] symmetrization
      InfiniteVector<double, Vector<double>::size_type> v;
      Mj1T_t.get_row(cols_left-1, v); // last column of left half
      for (typename InfiniteVector<double, Vector<double>::size_type>::const_iterator it(v.begin());
	   it != v.end(); ++it)
	for (int col = cols_left; col < 1<<(j-1); col++)
	  mj1T.set_entry(it.index()+2*(col-cols_left)+2, col, *it);
      
      Mj1T_t.get_row(cols_left, v); // first column of right half
      const int offset_right = (Deltasize(j+1)-Deltasize(j0()+1))-(1<<j)+(1<<j0()); // row offset for the right half
      for (typename InfiniteVector<double, Vector<double>::size_type>::const_iterator it(v.begin());
	   it != v.end(); ++it)
	for (int col = 1<<(j-1); col <= (1<<j) - cols_right - 1; col++)
	  mj1T.set_entry(it.index()+offset_right+2*(col-(1<<(j-1))), col, *it);
      
      mj1T.compress();
    }
  }
  
  template <int d, int dT>
  void
  DSBasis<d, dT>::assemble_Mj1T_t(const int j, SparseMatrix<double>& mj1T_t) const {
    if (j == j0())
      mj1T_t = Mj1T_t;
    else {
      SparseMatrix<double> help;
      assemble_Mj1T(j, help);
      mj1T_t = transpose(help);
    }
  }

  template <int d, int dT>
  void
  DSBasis<d, dT>::Mj0_get_row(const int j, const Vector<double>::size_type row,
			       InfiniteVector<double, Vector<double>::size_type>& v) const {
    if (j == j0())
      Mj0.get_row(row, v);
    else {
      v.clear();
      const size_t rows_top = (int)ceil(Deltasize(j0()+1)/2.0);
      if (row < rows_top)
	Mj0.get_row(row, v);
      else {
	const size_t bottom = Deltasize(j+1)-Deltasize(j0()+1)/2;
	    typedef Vector<double>::size_type size_type;
	if (row >= bottom)
	  Mj0.get_row(row+rows_top-bottom, v, Deltasize(j)-Deltasize(j0()));
	else
	  Mj0.get_row(rows_top-2+(row-rows_top)%2, v, (row-rows_top)/2+1);
      }
    }
  }
  
  template <int d, int dT>
  void
  DSBasis<d, dT>::Mj0T_get_row(const int j, const Vector<double>::size_type row,
				InfiniteVector<double, Vector<double>::size_type>& v) const {
    if (j == j0())
      Mj0T.get_row(row, v);
    else {
      v.clear();
      const size_t rows_top = (int)ceil(Deltasize(j0()+1)/2.0);
      if (row < rows_top)
	Mj0T.get_row(row, v);
      else {
	const size_t bottom = Deltasize(j+1)-Deltasize(j0()+1)/2;
	if (row >= bottom)
	  Mj0T.get_row(row+rows_top-bottom, v, Deltasize(j)-Deltasize(j0()));
	else
	  Mj0T.get_row(rows_top-2+(row-rows_top)%2, v, (row-rows_top)/2+1);
      }
    }
  }
  
  template <int d, int dT>
  void
  DSBasis<d, dT>::Mj1_get_row(const int j, const Vector<double>::size_type row,
			       InfiniteVector<double, Vector<double>::size_type>& v) const {
    if (j == j0())
      Mj1.get_row(row, v);
    else {
      v.clear();
      
      // Due to the [DS] symmetrization, we have to be a bit careful here.
      int fill_start = 0, fill_end = (1<<j)-1; // first and last column to fill up with interior filter
      const size_type upper_half = Deltasize(j0()+1)/2;
      if (row < upper_half) {
	// read the first half of the corresponding row in Mj1
	const size_type row_j0 = row;
	for (size_type k = 0; k < Mj1.entries_in_row(row_j0) && (int)Mj1.get_nth_index(row_j0,k) < 1<<(j0()-1); k++)
	  v.set_coefficient(Mj1.get_nth_index(row_j0,k),
			    Mj1.get_nth_entry(row_j0,k));
	fill_start = 1<<(j0()-1);
	fill_end   = (1<<j)-1;
      }
      else {
	const size_type bottom_half = Deltasize(j+1)-Deltasize(j0()+1)/2;
	if (row >= bottom_half) {
	  // read the second half of the corresponding row in Mj1
	  const size_type row_j0 = row+Deltasize(j0()+1)-Deltasize(j+1);
	  for (size_type k = 0; k < Mj1.entries_in_row(row_j0); k++)
	    if ((int)Mj1.get_nth_index(row_j0,k) >= 1<<(j0()-1))
	      v.set_coefficient(Mj1.get_nth_index(row_j0,k)+(1<<j)-(1<<j0()),
				Mj1.get_nth_entry(row_j0,k));
	  fill_start = 0;
	  fill_end   = (1<<j)-(1<<(j0()-1))-1;
	}
      }
      
      // Fill in the missing columns from the left half...

      InfiniteVector<double, Vector<double>::size_type> waveTfilter_left;
      Mj1_t.get_row((1<<(j0()-1))-1, waveTfilter_left);
      const int first_row_left = waveTfilter_left.begin().index();  // first nontrivial row in column (1<<(j0-1))-1
      
      // The row ...
      const int last_row_left  = waveTfilter_left.rbegin().index();
      // ... is the last nontrivial row in column (1<<(j0()-1))-1,
      // i.e. the row last_row_left begins at column (1<<(j0()-1))-1, so does the row last_row_left-1.
      // So the row "row" starts at column ...
      const int first_column_left = (1<<(j0()-1))-1+(int)floor(((int)row+1-last_row_left)/2.);
      
      for (int col = first_column_left, filter_row = last_row_left-abs(row-last_row_left)%2;
	   col < (1<<(j-1)) && col <= fill_end && filter_row >= first_row_left; col++, filter_row -= 2)
	if (col >= fill_start)
	  v.set_coefficient(col, waveTfilter_left.get_coefficient(filter_row));
      
      // Analogous strategy for the right half:
      
      InfiniteVector<double, Vector<double>::size_type> waveTfilter_right;
      Mj1_t.get_row((1<<(j0()-1)), waveTfilter_right);
      
      const int offset_right = (Deltasize(j+1)-Deltasize(j0()+1))-(1<<j)+(1<<j0()); // row offset for the right half
      const int first_row_right = waveTfilter_right.begin().index()+offset_right;
      const int last_row_right  = waveTfilter_right.rbegin().index()+offset_right;
      
      // The rows first_row_right and first_row_right+1 end at column 1<<(j-1),
      // so the row "row" ends at column ...
      const int last_column_right = (1<<(j-1))+(int)floor(((int)row-first_row_right)/2.);
      
      for (int col = last_column_right, filter_row = first_row_right-offset_right+abs(row-first_row_right)%2;
	   col >= 1<<(j-1) && col >= fill_start && filter_row <= last_row_right-offset_right; col--, filter_row += 2)
	if (col <= fill_end)
	  v.set_coefficient(col, waveTfilter_right.get_coefficient(filter_row));
    }
  }
  
  template <int d, int dT>
  void
  DSBasis<d, dT>::Mj1T_get_row(const int j, const Vector<double>::size_type row,
				InfiniteVector<double, Vector<double>::size_type>& v) const {
    if (j == j0())
      Mj1T.get_row(row, v);
    else {
      v.clear();
      
      int fill_start = 0, fill_end = (1<<j)-1; // first and last column to fill up with interior filter
      const size_type upper_half = Deltasize(j0()+1)/2;
      if (row < upper_half) {
	// read the first half of the corresponding row in Mj1
	const size_type row_j0 = row;
	for (size_type k = 0; k < Mj1T.entries_in_row(row_j0) && (int)Mj1T.get_nth_index(row_j0,k) < 1<<(j0()-1); k++)
	  v.set_coefficient(Mj1T.get_nth_index(row_j0,k),
			    Mj1T.get_nth_entry(row_j0,k));
	fill_start = 1<<(j0()-1);
	fill_end   = (1<<j)-1;
      }
      else {
	const size_type bottom_half = Deltasize(j+1)-Deltasize(j0()+1)/2;
	if (row >= bottom_half) {
	  // read the second half of the corresponding row in Mj1
	  const size_type row_j0 = row+Deltasize(j0()+1)-Deltasize(j+1);
	  for (size_type k = 0; k < Mj1T.entries_in_row(row_j0); k++)
	    if ((int)Mj1T.get_nth_index(row_j0,k) >= 1<<(j0()-1))
	      v.set_coefficient(Mj1T.get_nth_index(row_j0,k)+(1<<j)-(1<<j0()),
				Mj1T.get_nth_entry(row_j0,k));
	  fill_start = 0;
	  fill_end   = (1<<j)-(1<<(j0()-1))-1;
	}
      }

      // Fill in the missing columns from the left half...

      InfiniteVector<double, Vector<double>::size_type> waveTfilter_left;
      Mj1T_t.get_row((1<<(j0()-1))-1, waveTfilter_left);
      const int first_row_left = waveTfilter_left.begin().index();  // first nontrivial row in column (1<<(j0-1))-1

      // The row ...
      const int last_row_left  = waveTfilter_left.rbegin().index();
      // ... is the last nontrivial row in column (1<<(j0()-1))-1,
      // i.e. the row last_row_left begins at column (1<<(j0()-1))-1, so does the row last_row_left-1.
      // So the row "row" starts at column ...
      const int first_column_left = (1<<(j0()-1))-1+(int)floor(((int)row+1-last_row_left)/2.);
      
      for (int col = first_column_left, filter_row = last_row_left-abs(row-last_row_left)%2;
	   col < (1<<(j-1)) && col <= fill_end && filter_row >= first_row_left; col++, filter_row -= 2)
	if (col >= fill_start)
	  v.set_coefficient(col, waveTfilter_left.get_coefficient(filter_row));
      
      // Analogous strategy for the right half:
      
      InfiniteVector<double, Vector<double>::size_type> waveTfilter_right;
      Mj1T_t.get_row((1<<(j0()-1)), waveTfilter_right);
      
      const int offset_right = (Deltasize(j+1)-Deltasize(j0()+1))-(1<<j)+(1<<j0()); // row offset for the right half
      const int first_row_right = waveTfilter_right.begin().index()+offset_right;
      const int last_row_right  = waveTfilter_right.rbegin().index()+offset_right;
      
      // The rows first_row_right and first_row_right+1 end at column 1<<(j-1),
      // so the row "row" ends at column ...
      const int last_column_right = (1<<(j-1))+(int)floor(((int)row-first_row_right)/2.);
      
      for (int col = last_column_right, filter_row = first_row_right-offset_right+abs(row-first_row_right)%2;
	   col >= 1<<(j-1) && col >= fill_start && filter_row <= last_row_right-offset_right; col--, filter_row += 2)
	if (col <= fill_end)
	  v.set_coefficient(col, waveTfilter_right.get_coefficient(filter_row));
    }
  }
  
  template <int d, int dT>
  void
  DSBasis<d, dT>::Mj0_t_get_row(const int j, const Vector<double>::size_type row,
				 InfiniteVector<double, Vector<double>::size_type>& v) const {
    if (j == j0())
      Mj0_t.get_row(row, v);
    else {
      const size_t rows_top = (int)ceil(Deltasize(j0())/2.0);
      if (row < rows_top)
	Mj0_t.get_row(row, v);
      else {
	const size_t bottom = Deltasize(j)-Deltasize(j0())/2;
	if (row >= bottom)
	  Mj0_t.get_row(row+rows_top-bottom, v, Deltasize(j+1)-Deltasize(j0()+1));
	else
	  Mj0_t.get_row(rows_top-1, v, 2*(row-rows_top)+2);
      }
    }
  }
  
  template <int d, int dT>
  void
  DSBasis<d, dT>::Mj0T_t_get_row(const int j, const Vector<double>::size_type row,
				  InfiniteVector<double, Vector<double>::size_type>& v) const {
    if (j == j0())
      Mj0T_t.get_row(row, v);
    else {
      const size_t rows_top = (int)ceil(Deltasize(j0())/2.0);
      if (row < rows_top)
	Mj0T_t.get_row(row, v);
      else {
	const size_t bottom = Deltasize(j)-Deltasize(j0())/2;
	if (row >= bottom)
	  Mj0T_t.get_row(row+rows_top-bottom, v, Deltasize(j+1)-Deltasize(j0()+1));
	else
	  Mj0T_t.get_row(rows_top-1, v, 2*(row-rows_top)+2);
      }
    }
  }
  
  template <int d, int dT>
  void
  DSBasis<d, dT>::Mj1_t_get_row(const int j, const Vector<double>::size_type row,
				 InfiniteVector<double, Vector<double>::size_type>& v) const {
    if (j == j0())
      Mj1_t.get_row(row, v);
    else {
      const size_t rows_top = 1<<(j0()-1);
      if (row < rows_top)
	Mj1_t.get_row(row, v);
      else {
	const size_t bottom = (1<<j)-(1<<(j0()-1));
	if (row >= bottom)
	  Mj1_t.get_row(row+rows_top-bottom, v, Deltasize(j+1)-Deltasize(j0()+1));
	else {
	  if ((int)row < (1<<(j-1)))
	    Mj1_t.get_row(rows_top-1, v, 2*(row-rows_top)+2);
	  else
	    Mj1_t.get_row(1<<(j0()-1), v, Deltasize(j+1)-Deltasize(j0()+1)+2*((int)row-bottom));
	}
      }
    }
  }
  
  template <int d, int dT>
  void
  DSBasis<d, dT>::Mj1T_t_get_row(const int j, const Vector<double>::size_type row,
				  InfiniteVector<double, Vector<double>::size_type>& v) const {
    if (j == j0())
      Mj1T_t.get_row(row, v);
    else {
      const size_t rows_top = 1<<(j0()-1);
      if (row < rows_top)
	Mj1T_t.get_row(row, v);
      else {
	const size_t bottom = (1<<j)-(1<<(j0()-1));
	if (row >= bottom)
	  Mj1T_t.get_row(row+rows_top-bottom, v, Deltasize(j+1)-Deltasize(j0()+1));
	else {
	  if ((int)row < (1<<(j-1)))
	    Mj1T_t.get_row(rows_top-1, v, 2*(row-rows_top)+2);
	  else
	    Mj1T_t.get_row(1<<(j0()-1), v, Deltasize(j+1)-Deltasize(j0()+1)+2*((int)row-bottom));
	}
      }
    }
  }
}
