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
    assert(std::max(s0,s1) < d);
    assert(std::min(s0,s1) >= d-2);
    
    this->s0 = s0;
    this->s1 = s1;

    setup();
  }

  template <int d, int dT>
  PBasis<d,dT>::PBasis(const bool bc_left, const bool bc_right) {
    this->s0 = bc_left ? 1 : 0;
    this->s1 = bc_right ? 1 : 0;

    setup();
  }

  template <int d, int dT>
  inline
  typename PBasis<d,dT>::Index
  PBasis<d,dT>::first_generator(const int j) const
  {
    assert(j >= j0());
    return Index(j, 0, DeltaLmin(), this);
  }
  
  template <int d, int dT>
  inline
  typename PBasis<d,dT>::Index
  PBasis<d,dT>::last_generator(const int j) const
  {
    assert(j >= j0());
    return Index(j, 0, DeltaRmax(j), this);
  }

  template <int d, int dT>
  inline
  typename PBasis<d,dT>::Index
  PBasis<d,dT>::first_wavelet(const int j) const
  {
    assert(j >= j0());
    return Index(j, 1, Nablamin(), this);
  }
  
  template <int d, int dT>
  inline
  typename PBasis<d,dT>::Index
  PBasis<d,dT>::last_wavelet(const int j) const
  {
    assert(j >= j0());
    return Index(j, 1, Nablamax(j), this);
  }

  template <int d, int dT>
  void
  PBasis<d,dT>::setup() {
    // For simplicity, we do not implement a fully generic setup here
    // but only setup the parameters for several important special cases from [P].
    // Namely, we restrict the setting to the case s0 >= d-2, where no "additional"
    // interior dual generators have to be constructed. In a later version of
    // this class, we will fix this.

#if 1
    // choose j0 s.th. the supports of the dual boundary generators do not overlap
    // (at the end of the setup, j0 may be reduced by one, see below)
    j0_ = (int) ceil(log(ell2T<d,dT>()-ell1T<d,dT>()+std::max(s0,s1)+1.-d)/M_LN2+1);
#else
    // choose j0 s.th. the supports of the primal boundary generators do not overlap
    // the supports of the dual boundary generators from the other endpoint
    // (at the end of the setup, j0 may be reduced by one, see below)
    j0_ = (int) ceil(log(ell2T<d,dT>()-ell1T<d,dT>()+2*std::max(s0,s1)+dT-d+1.)/M_LN2);
#endif
    
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
//     cout << "ML=" << endl << ML_;

    MR_.resize(3*dT+s1-1, dT);
    for (int column = s1; column < d-1; column++)
      for (int row = s1; row < 2*(d-1); row++)
	MR_.set_entry(row-s1, column-s1, ML_0.get_entry(row,column));
    for (int k = d; k <= dT+s1; k++)
      for (int n = 2*k-d; n <= 2*k; n++)
 	MR_.set_entry(n-1-s1,k-1-s1,cdf.a().get_coefficient(MultiIndex<int,1>(-(d/2)+n+d-2*k)));
//     cout << "MR=" << endl << MR_;

    // setup the expansion coefficients of the unbiorthogonalized dual generators
    // w.r.t. the truncated CDF generators, see [DKU, Lemma 3.1]
    Matrix<double> MLTp(3*dT+s0-1, dT);
    for (unsigned int r = 0; r < dT; r++) {
      MLTp.set_entry(r, r, ldexp(1.0, -r));
      for (int m = ellT_l(); m <= 2*ellT_l()+ell1T<d,dT>()-1; m++)
	MLTp.set_entry(dT+m-ellT_l(), r, alpha(m, r)*ldexp(1.0, -r));
      for (int m = 2*ellT_l()+ell1T<d,dT>(); m <= 2*ellT_l()+ell2T<d,dT>()-2; m++)
	MLTp.set_entry(dT+m-ellT_l(), r, betaL(m, r));
    }
//     cout << "MLTp=" << endl << MLTp;
        
    Matrix<double> MRTp(3*dT+s1-1, dT);
    for (unsigned int r = 0; r < dT; r++) {
      MRTp.set_entry(r, r, ldexp(1.0, -r));
      for (int m = ellT_r(); m <= 2*ellT_r()+ell1T<d,dT>()-1; m++)
	MRTp.set_entry(dT+m-ellT_r(), r, alpha(m, r)*ldexp(1.0, -r));
      for (int m = 2*ellT_r()+ell1T<d,dT>(); m <= 2*ellT_r()+ell2T<d,dT>()-2; m++)
	MRTp.set_entry(dT+m-ellT_r(), r, betaR(m, r));
    }
//     cout << "MRTp=" << endl << MRTp;

    // for the biorthogonalization of the generators,
    // compute the gramian matrix of the primal and dual boundary generators
    compute_biorthogonal_boundary_gramian
      <CDFMask_primal<d>,CDFMask_dual<d,dT> >(ML_, MLTp, GammaL);
//     cout << "GammaL=" << endl << GammaL;

    compute_biorthogonal_boundary_gramian
      <CDFMask_primal<d>,CDFMask_dual<d,dT> >(MR_, MRTp, GammaR);
//     cout << "GammaR=" << endl << GammaR;

    setup_CXT();
    setup_CXAT();

    setup_Mj0(ML_, MR_, Mj0); // [DKU, (3.5.1)]
    Mj0.compress(1e-8);
    
    SparseMatrix<double> mj0tp; setup_Mj0Tp(MLTp, MRTp, mj0tp); // [DKU, (3.5.5)]
    mj0tp.compress(1e-8);
    
    // biorthogonalize the generators
    setup_Cj();
    Mj0T = transpose(inv_CjpT) * mj0tp * transpose(CjT); // [DKU, (2.4.3)]
    Mj0T.compress(1e-8);

//     Mj0.scale(M_SQRT2);
//     cout << "Mj0 (without SQRT1_2)=" << endl << Mj0;
//     Mj0.scale(M_SQRT1_2);

//     Mj0T.scale(M_SQRT2);
//     cout << "Mj0T (without SQRT1_2)=" << endl << Mj0T;
//     Mj0T.scale(M_SQRT1_2);

#if 0
    cout << "PBasis(): check biorthogonality of Mj0, Mj0T:" << endl;
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

    // construction of the wavelet basis: initial stable completion, cf. [DKU section 4.1]

    SparseMatrix<double> FF; F(FF);           // [DKU, (4.1.14)]
    SparseMatrix<double> PP; P(ML_, MR_, PP); // [DKU, (4.1.22)]

    SparseMatrix<double> A, H, Hinv;
    GSetup(A, H, Hinv); // [DKU, (4.1.1), (4.1.13)]

#if 0
    SparseMatrix<double> Aold(A); // for the checks below
#endif

    GElim (A, H, Hinv); // elimination [DKU, (4.1.4)ff.]
#if 1
    // choose this branch if you don't want the interior wavelets to be the CDF ones
    // (in that case, the Riesz constants improve, but the condition numbers for
    // the Poisson equation will deteriorate)
    SparseMatrix<double> BB; BT(A, BB); // [DKU, (4.1.13)]
#else
    SparseMatrix<double> BB; double binv = BT(A, BB); // [DKU, (4.1.13)]
//     cout << "b^{-1}=" << binv << endl;
#endif
    
#if 0
    cout << "PBasis(): check properties (4.1.15):" << endl;
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

#if 0
    cout << "PBasis(): check factorization of A:" << endl;
    SparseMatrix<double> testAfact = Aold - Hinv*A;
    cout << "* in infty-norm: " << row_sum_norm(testAfact) << endl;

    cout << "PBasis(): check that H is inverse to Hinv:" << endl;
    SparseMatrix<double> testHinv = H*Hinv;
    for (unsigned int i = 0; i < testHinv.row_dimension(); i++)
      testHinv.set_entry(i, i, testHinv.get_entry(i, i) - 1.0);
    cout << "* in infty-norm: " << row_sum_norm(testHinv) << endl;
#endif

    SparseMatrix<double> mj1ih = PP * Hinv * FF; // [DKU, (4.1.23)]
//     cout << "mj1ih (without SQRT1_2):" << endl << mj1ih;
    mj1ih.scale(M_SQRT1_2); // [P, Prop. 5.6]
//     cout << "PP=" << endl << PP;
//     cout << "Hinv=" << endl << Hinv;

    SparseMatrix<double> PPinv; InvertP(PP, PPinv);

    SparseMatrix<double> help = H * PPinv;
    SparseMatrix<double> gj0ih = transpose(BB) * help; gj0ih.scale(M_SQRT2); // [P, Prop. 5.6]
    SparseMatrix<double> gj1ih = transpose(FF) * help; gj1ih.scale(M_SQRT2); // [DKU, (4.1.24)], [P, Prop. 5.6]

#if 0
    cout << "PBasis(): check initial stable completion:" << endl;
    SparseMatrix<double> mj_initial(Mj0.row_dimension(),
				    Mj0.column_dimension() + mj1ih.column_dimension());
    for (unsigned int i = 0; i < Mj0.row_dimension(); i++)
      for (unsigned int j = 0; j < Mj0.column_dimension(); j++)	{
	const double help = Mj0.get_entry(i, j);
	if (help != 0)
	  mj_initial.set_entry(i, j, help);
      }
    for (unsigned int i = 0; i < mj1ih.row_dimension(); i++)
      for (unsigned int j = 0; j < mj1ih.column_dimension(); j++) {
	const double help = mj1ih.get_entry(i, j);
	if (help != 0)
	  mj_initial.set_entry(i, j+Mj0.column_dimension(), help);
      }
    cout << "Mj=" << endl << mj_initial;

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
    cout << "Gj=" << endl << gj_initial;
    
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
    // (cf. DSBasis, here we have the special case Cj=Cjp=I)
    SparseMatrix<double> I; I.diagonal(Deltasize(j0()+1), 1.0);
    Mj1  = (I - (Mj0*transpose(Mj0T))) * mj1ih;
    Mj1T = transpose(gj1ih);

    Mj0 .compress(1e-8);
    Mj1 .compress(1e-8);
    Mj0T.compress(1e-8);
    Mj1T.compress(1e-8);

#if 0
    // adjust interior wavelets, such that they correspond to the CDF ones, cf. [P, pp. 121-123]
    // (if you want to enable this, binv should have been computed, see some lines above)
    const double scaling_factor = minus1power(d+1)*2*binv;
    Mj1 .scale(scaling_factor);
    Mj1T.scale(1./scaling_factor);
#endif
    
//     Mj1.scale(M_SQRT2);
//     cout << "Mj1 (without SQRT1_2)=" << endl << Mj1;
//     Mj1.scale(M_SQRT1_2);

//     Mj1T.scale(M_SQRT2);
//     cout << "Mj1T (without SQRT1_2)=" << endl << Mj1T;
//     Mj1T.scale(M_SQRT1_2);

#if 0
    cout << "PBasis(): check new stable completion:" << endl;
    
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
    if (d%2 && s0 == s1)
      {
	DS_symmetrization(Mj1, Mj1T);
#if 0
	{
	  cout << "PBasis(): check [DS] symmetrization:" << endl;
	  
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

    // After the overall setup, we try to reduce j0 by one,
    // resulting in much better condition numbers:
    SparseMatrix<double> Mj0_reduced(Deltasize(j0_), Deltasize(j0_-1));
    for (int col = 0; col < Deltasize(j0_-1)/2; col++)
      for (int row = 0; row < Deltasize(j0_); row++)
	Mj0_reduced.set_entry(row, col, Mj0.get_entry(row, col));
    for (int col = Deltasize(j0_-1)/2; col < Deltasize(j0_-1); col++)
      for (int row = 0; row < Deltasize(j0_); row++)
	Mj0_reduced.set_entry(row, col, Mj0.get_entry(row+(1<<j0_), col+(1<<(j0_-1))));

//     cout << "Mj0:" << endl << Mj0;
//     cout << "Mj0 reduced:" << endl << Mj0_reduced;

    SparseMatrix<double> Mj1_reduced(Deltasize(j0_), 1<<(j0_-1));
    for (int col = 0; col < 1<<(j0_-2); col++)
      for (int row = 0; row < Deltasize(j0_); row++)
	Mj1_reduced.set_entry(row, col, Mj1.get_entry(row, col));
    for (int col = 1<<(j0_-2); col < 1<<(j0_-1); col++)
      for (int row = 0; row < Deltasize(j0_); row++)
	Mj1_reduced.set_entry(row, col, Mj1.get_entry(row+(1<<j0_), col+(1<<(j0_-1))));

//     cout << "Mj1:" << endl << Mj1;
//     cout << "Mj1 reduced:" << endl << Mj1_reduced;

    SparseMatrix<double> Mj0T_reduced(Deltasize(j0_), Deltasize(j0_-1));
    for (int col = 0; col < Deltasize(j0_-1)/2; col++)
      for (int row = 0; row < Deltasize(j0_); row++)
	Mj0T_reduced.set_entry(row, col, Mj0T.get_entry(row, col));
    for (int col = Deltasize(j0_-1)/2; col < Deltasize(j0_-1); col++)
      for (int row = 0; row < Deltasize(j0_); row++)
	Mj0T_reduced.set_entry(row, col, Mj0T.get_entry(row+(1<<j0_), col+(1<<(j0_-1))));

//     cout << "Mj0T:" << endl << Mj0T;
//     cout << "Mj0T reduced:" << endl << Mj0T_reduced;

    SparseMatrix<double> Mj1T_reduced(Deltasize(j0_), 1<<(j0_-1));
    for (int col = 0; col < 1<<(j0_-2); col++)
      for (int row = 0; row < Deltasize(j0_); row++)
	Mj1T_reduced.set_entry(row, col, Mj1T.get_entry(row, col));
    for (int col = 1<<(j0_-2); col < 1<<(j0_-1); col++)
      for (int row = 0; row < Deltasize(j0_); row++)
	Mj1T_reduced.set_entry(row, col, Mj1T.get_entry(row+(1<<j0_), col+(1<<(j0_-1))));
    
//     cout << "Mj1T:" << endl << Mj1T;
//     cout << "Mj1T reduced:" << endl << Mj1T_reduced;

    // try to reduce j0:
    SparseMatrix<double> mj_reduced(Mj0_reduced.row_dimension(),
				    Mj0_reduced.column_dimension() + Mj1_reduced.column_dimension());
    for (unsigned int i = 0; i < Mj0_reduced.row_dimension(); i++)
      for (unsigned int j = 0; j < Mj0_reduced.column_dimension(); j++) {
	const double help = Mj0_reduced.get_entry(i, j);
	if (help != 0)
	  mj_reduced.set_entry(i, j, help);
      }
    for (unsigned int i = 0; i < Mj1_reduced.row_dimension(); i++)
      for (unsigned int j = 0; j < Mj1_reduced.column_dimension(); j++) {
	const double help = Mj1_reduced.get_entry(i, j);
	if (help != 0)
	  mj_reduced.set_entry(i, j+Mj0_reduced.column_dimension(), help);
      }
    
    SparseMatrix<double> gj_reduced(Mj0T_reduced.column_dimension() + Mj1T_reduced.column_dimension(),
				    Mj0T_reduced.row_dimension());
    for (unsigned int i = 0; i < Mj0T_reduced.column_dimension(); i++)
      for (unsigned int j = 0; j < Mj0T_reduced.row_dimension(); j++) {
	const double help = Mj0T_reduced.get_entry(j, i);
	if (help != 0)
	  gj_reduced.set_entry(i, j, help);
      }
    for (unsigned int i = 0; i < Mj1T_reduced.column_dimension(); i++)
      for (unsigned int j = 0; j < Mj1T_reduced.row_dimension(); j++) {
	const double help = Mj1T_reduced.get_entry(j, i);
	if (help != 0)
	  gj_reduced.set_entry(i+Mj0T_reduced.column_dimension(), j, help);
      }
    
    SparseMatrix<double> test_reduced = mj_reduced * gj_reduced;
    for (unsigned int i = 0; i < test_reduced.row_dimension(); i++)
      test_reduced.set_entry(i, i, test_reduced.get_entry(i, i) - 1.0);
//     cout << "* for the reduced system, ||Mj*Gj-I||_infty: " << row_sum_norm(test_reduced) << endl;

    if (row_sum_norm(test_reduced) < 1e-10) {
      Mj0  = Mj0_reduced;
      Mj0T = Mj0T_reduced;
      Mj1  = Mj1_reduced;
      Mj1T = Mj1T_reduced;
      j0_--;
    }

    Mj0_t  = transpose(Mj0);
    Mj0T_t = transpose(Mj0T);
    Mj1_t  = transpose(Mj1);
    Mj1T_t = transpose(Mj1T);
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

    Mj0.scale(M_SQRT1_2);
    
//     cout << "Mj0=" << endl << Mj0 << endl;
  }

  template <int d, int dT>
  void
  PBasis<d,dT>::setup_Mj0Tp(const Matrix<double>& MLTp, const Matrix<double>& MRTp, SparseMatrix<double>& Mj0Tp) {
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
    
    int startrow = dT+ellT_l()+ell1T<d,dT>();
    for (int col = dT; col < nj-dT; col++, startrow+=2) {
      int row = startrow;
      for (MultivariateLaurentPolynomial<double, 1>::const_iterator it(cdf.aT().begin());
	   it != cdf.aT().end(); ++it, row++)
	Mj0Tp.set_entry(row, col, *it);
    }
    
    Mj0Tp.scale(M_SQRT1_2);

//     cout << "Mj0Tp=" << endl << Mj0Tp << endl;
  }

  template <int d, int dT>
  void
  PBasis<d,dT>::setup_CXT()
  {
    // IGPMlib reference: I_Mask_Bspline::EvalCL(), ::EvalCR()

    Matrix<double> CLGammaLInv;
    QUDecomposition<double>(GammaL).inverse(CLGammaLInv);
    CLT = transpose(CLGammaLInv);
    
    Matrix<double> CRGammaRInv;
    QUDecomposition<double>(GammaR).inverse(CRGammaRInv);
    CRT = transpose(CRGammaRInv);
    
    QUDecomposition<double>(CLT).inverse(inv_CLT);
    QUDecomposition<double>(CRT).inverse(inv_CRT);
  }

  template <int d, int dT>
  void PBasis<d,dT>::setup_CXAT() {
    // IGPMlib reference: I_Mask_Bspline::EvalCL(), ::EvalCR()

    // setup CLAT <-> Alpha * (CLT)^T
    CLAT.resize(ellT_l()+ell2T<d,dT>()-1, dT);
    for (int i = 1-ell2T<d,dT>(); i <= ellT_l()-1; i++)
      for (unsigned int r = 0; r < dT; r++)
  	CLAT(ell2T<d,dT>()-1+i, r) = alpha(i, r);

//     cout << "CLAT before biorthogonalization:" << endl << CLAT << endl;
    CLAT = CLAT * transpose(CLT);
    CLAT.compress(1e-12);
//     cout << "CLAT after biorthogonalization:" << endl << CLAT << endl;

    // the same for CRAT:
    CRAT.resize(ellT_r()+ell2T<d,dT>()-1, dT);
    for (int i = 1-ell2T<d,dT>(); i <= ellT_r()-1; i++)
      for (unsigned int r = 0; r < dT; r++)
 	CRAT(ell2T<d,dT>()-1+i, r) = alpha(i, r);

//     cout << "CRAT before biorthogonalization:" << endl << CRAT << endl;
    CRAT = CRAT * transpose(CRT);
    CRAT.compress(1e-12);
//     cout << "CRAT after biorthogonalization:" << endl << CRAT << endl;
  }

  template <int d, int dT>
  void PBasis<d,dT>::setup_Cj() {
    // IGPMlib reference: I_Basis_Bspline_s::setup_Cj(), ::put_Mat()

    // [DKU (5.2.5)]

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

#if 0
    cout << "PBasis: testing setup of Cj:" << endl;

    SparseMatrix<double> test1 = CjT * inv_CjT;
    for (unsigned int i = 0; i < test1.row_dimension(); i++)
      test1.set_entry(i, i, test1.get_entry(i, i) - 1.0);
    cout << "* ||CjT*inv_CjT-I||_infty: " << row_sum_norm(test1) << endl;

    SparseMatrix<double> test3 = CjpT * inv_CjpT;
    for (unsigned int i = 0; i < test3.row_dimension(); i++)
      test3.set_entry(i, i, test3.get_entry(i, i) - 1.0);
    cout << "* ||CjpT*inv_CjpT-I||_infty: " << row_sum_norm(test3) << endl;

#endif
  }

  template <int d, int dT>
  void
  PBasis<d, dT>::F(SparseMatrix<double>& FF) {
    // IGPMlib reference: I_Basis_Bspline_s::F()
    
    const int FLow = ell2<d>()-1;      // first column of F_j in Fhat_j [P, p. 113]
    const int FUp  = (1<<j0())+ell1<d>(); // last column of F_j in Fhat_j
    
    // (4.1.14):
    
    FF.resize(Deltasize(j0()+1), 1<<j0());

    for (int r = 0; r < FLow; r++)
      FF.set_entry(r+d-1-s0, r, 1.0);
    
    int i = d-2+ell2<d>()-s0;
    for (int col = FLow; col <= FUp; col++, i+=2)
      FF.set_entry(i, col, 1.0);
    
    i = Deltasize(j0()+1)-d+s1;
    for (int col = (1<<j0())-1; col >= FUp+1; col--, i--)
      FF.set_entry(i, col, 1.0);

//     cout << "F=" << endl << FF << endl;
  }

  template <int d, int dT>
  void
  PBasis<d, dT>::P(const Matrix<double>& ML, const Matrix<double>& MR, SparseMatrix<double>& PP) {
    // IGPMlib reference: I_Basis_Bspline_s::P()
    
    // (4.1.22):

    PP.diagonal(Deltasize(j0()+1), 1.0);
    
    for (unsigned int i = 0; i < ML.row_dimension(); i++)
      for (int k = 0; k < d-1-s0; k++)
 	PP.set_entry(i, k, ML.get_entry(i, k));

    for (unsigned int i = 0; i < MR.row_dimension(); i++)
      for (int k = 0; k < d-1-s1; k++)
 	PP.set_entry(Deltasize(j0()+1)-i-1, Deltasize(j0()+1)-k-1, MR.get_entry(i, k));

//     cout << "P=" << endl << PP << endl;
  }

  template <int d, int dT>
  void
  PBasis<d, dT>::GSetup(SparseMatrix<double>& A, SparseMatrix<double>& H, SparseMatrix<double>& Hinv) {
    // IGPMlib reference: I_Basis_Bspline_s::GSetup()

    // (4.1.13):
    
    const int nj  = Deltasize(j0());
    const int njp = Deltasize(j0()+1);
    A.resize(njp, nj);

    for (int row = 0; row < d-1-s0; row++)
      A.set_entry(row, row, 1.0);
    
    int startrow = d-1-s0;
    for (int col = d-1-s0; col <= nj-d+s1; col++, startrow+=2) {
      int row = startrow;
      for (MultivariateLaurentPolynomial<double, 1>::const_iterator it(cdf.a().begin());
 	   it != cdf.a().end(); ++it, row++) {
  	A.set_entry(row, col, *it);
      }
    }
    
    for (int row = njp-1, col = nj-1; col >= nj-d+s1+1; row--, col--)
      A.set_entry(row, col, 1.0);

    // prepare H, Hinv for elimination process:
    H   .diagonal(njp, 1.0);
    Hinv.diagonal(njp, 1.0);

//     cout << "A=" << endl << A << endl;
  }

  template <int d, int dT>
  void
  PBasis<d, dT>::GElim(SparseMatrix<double>& A, SparseMatrix<double>& H, SparseMatrix<double>& Hinv) {
    // IGPMlib reference: I_Basis_Bspline_s::gelim()
    
    const int firstcol = d-1-s0; // first column of A_j^{(d)} in Ahat_j^{(d)}
    const int lastcol  = Deltasize(j0())-d+s1; // last column
    const int firstrow = d-1-s0; // first row of A_j^{(d)} in Ahat_j^{(d)}
    const int lastrow  = Deltasize(j0()+1)-d+s1; // last row
    
    SparseMatrix<double> help;

//     cout << "A=" << endl << A;
    
    // elimination (4.1.4)ff.:
    for (int i = 1; i <= d; i++) {
      help.diagonal(Deltasize(j0()+1), 1.0);

      // row index of the entry in the first column of A_j^{(i)} (w.r.t. Ahat_j^{(i)})
      const int elimrow = i%2 ? firstrow+(i-1)/2 : lastrow-(int)floor((i-1)/2.);
//       cout << "i%2=" << i%2 << ", elimrow=" << elimrow << endl;

      // factorization [P, p. 112]      
//       const int HhatLow = i%2 ? d-s0-((i+1)/2)%2 : d-s0+1-(d%2)-(i/2)%2;
      const int HhatLow = i%2 ? d-s0-((i+1)/2)%2 : d-s0+1-(d%2)-(i/2)%2;
      const int HhatUp  = i%2
	? Deltasize(j0()+1)-d+s1-abs((d%2)-((i+1)/2)%2)
	: Deltasize(j0()+1)-d+s1-abs((d%2)-(i/2)%2);
      
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
      
//       cout << "Hfactor=" << endl << help;

      A = help * A;
      H = help * H;
      
      A.compress(1e-14);

//       cout << "A=" << endl << A;

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

//     cout << "finally, H=" << endl << H;
//     cout << "finally, Hinv=" << endl << Hinv;
  }


  template <int d, int dT>
  double
  PBasis<d, dT>::BT(const SparseMatrix<double>& A, SparseMatrix<double>& BB) {
    // IGPMlib reference: I_Basis_Bspline_s::Btr()
    
    const int nj  = Deltasize(j0());
    const int njp = Deltasize(j0()+1);
    BB.resize(njp,nj);

    for (int r = 0; r < d-1-s0; r++)
      BB.set_entry(r, r, 1.0);

    const double binv = 1./A.get_entry(d-1-s0+ell2<d>(), d-1-s0);

    int i = d-1-s0+ell2<d>();
    for (int col = d-1-s0; col <= nj-d+s1; col++, i+=2)
      BB.set_entry(i, col, binv);

    for (int row = njp-1, col = nj-1; col >= nj-d+s1+1; row--, col--)
      BB.set_entry(row, col, 1.0);

//     cout << "BThat=" << endl << BB;

    return binv;
  }

  template <int d, int dT>
  void
  PBasis<d, dT>::InvertP(const SparseMatrix<double>& PP, SparseMatrix<double>& PPinv) {
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

#if 0
    cout << "PBasis(): check that PPinv is inverse to PP:" << endl;
    SparseMatrix<double> testPinv = PP*PPinv;
//     cout << "  PP*PPinv=" << endl << testPinv;
    for (unsigned int i = 0; i < testPinv.row_dimension(); i++)
      testPinv.set_entry(i, i, testPinv.get_entry(i, i) - 1.0);
    cout << "* in infty-norm: " << row_sum_norm(testPinv) << endl;
#endif
  }
  
  template <int d, int dT>
  void
  PBasis<d, dT>::DS_symmetrization(SparseMatrix<double>& Mj1, SparseMatrix<double>& Mj1T) {
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
  PBasis<d, dT>::assemble_Mj0(const int j, SparseMatrix<double>& mj0) const {
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
  PBasis<d, dT>::assemble_Mj0_t(const int j, SparseMatrix<double>& mj0_t) const {
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
  PBasis<d, dT>::assemble_Mj0T(const int j, SparseMatrix<double>& mj0T) const {
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
  PBasis<d, dT>::assemble_Mj0T_t(const int j, SparseMatrix<double>& mj0T_t) const {
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
  PBasis<d, dT>::assemble_Mj1(const int j, SparseMatrix<double>& mj1) const {
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
  PBasis<d, dT>::assemble_Mj1_t(const int j, SparseMatrix<double>& mj1_t) const {
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
  PBasis<d, dT>::assemble_Mj1T(const int j, SparseMatrix<double>& mj1T) const {
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
  PBasis<d, dT>::assemble_Mj1T_t(const int j, SparseMatrix<double>& mj1T_t) const {
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
  PBasis<d, dT>::Mj0_get_row(const int j, const Vector<double>::size_type row,
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
  PBasis<d, dT>::Mj0T_get_row(const int j, const Vector<double>::size_type row,
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
  PBasis<d, dT>::Mj1_get_row(const int j, const Vector<double>::size_type row,
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
  PBasis<d, dT>::Mj1T_get_row(const int j, const Vector<double>::size_type row,
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
  PBasis<d, dT>::Mj0_t_get_row(const int j, const Vector<double>::size_type row,
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
  PBasis<d, dT>::Mj0T_t_get_row(const int j, const Vector<double>::size_type row,
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
  PBasis<d, dT>::Mj1_t_get_row(const int j, const Vector<double>::size_type row,
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
  PBasis<d, dT>::Mj1T_t_get_row(const int j, const Vector<double>::size_type row,
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

  template <int d, int dT>
  void
  PBasis<d, dT>::decompose(const InfiniteVector<double, Index>& c,
			    const int jmin,
			    InfiniteVector<double, Index>& v) const {
    v.clear();
    InfiniteVector<double, Index> help;
    for (typename InfiniteVector<double, Index>::const_iterator it(c.begin()), itend(c.end());
	 it != itend; ++it) {
      decompose_1(it.index(), jmin, help); // calls help.clear() first
      v.add(*it, help);
    }
  }
  
  template <int d, int dT>
  void
  PBasis<d, dT>::decompose_t(const InfiniteVector<double, Index>& c,
			      const int jmin,
			      InfiniteVector<double, Index>& v) const {
    v.clear();
    InfiniteVector<double, Index> help;
    for (typename InfiniteVector<double, Index>::const_iterator it(c.begin()), itend(c.end());
	 it != itend; ++it) {
      decompose_t_1(it.index(), jmin, help); // calls help.clear() first
      v.add(*it, help);
    }
  }
  
  template <int d, int dT>
  void
  PBasis<d, dT>::reconstruct(const InfiniteVector<double, Index>& c,
			      const int j,
			      InfiniteVector<double, Index>& v) const {
    v.clear();
    for (typename InfiniteVector<double, Index>::const_iterator it(c.begin()), itend(c.end());
	 it != itend; ++it) {
      InfiniteVector<double, Index> help;
      reconstruct_1(it.index(), j, help);
      v.add(*it, help);
    }
  }

  template <int d, int dT>
  void
  PBasis<d, dT>::reconstruct_t(const InfiniteVector<double, Index>& c,
				const int j,
				InfiniteVector<double, Index>& v) const {
    v.clear();
    for (typename InfiniteVector<double, Index>::const_iterator it(c.begin()), itend(c.end());
	 it != itend; ++it) {
      InfiniteVector<double, Index> help;
      reconstruct_t_1(it.index(), j, help);
      v.add(*it, help);
    }
  }

  template <int d, int dT>
  void
  PBasis<d, dT>::decompose_1(const Index& lambda,
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
	      int fill_start = 0, fill_end = (1<<(lambda.j()-1))-1; // first and last column to fill up with interior filter
	      const size_type upper_half = Deltasize(j0()+1)/2;
	      if (row < upper_half) {
		// read the first half of the corresponding row in Mj1T
		for (size_type k(0); k < Mj1T.entries_in_row(row_j0) && (int)Mj1T.get_nth_index(row_j0,k) < 1<<(j0()-1); k++)
		  c.set_coefficient(Index(lambda.j()-1, 1, Mj1T.get_nth_index(row_j0,k), this),
				    Mj1T.get_nth_entry(row_j0,k));
		fill_start = 1<<(j0()-1);
		fill_end   = (1<<(lambda.j()-1))-1;
	      }
	      else {
		const size_type bottom_half = Deltasize(lambda.j())-Deltasize(j0()+1)/2;
		if (row >= bottom_half) {
		  // read the second half of the corresponding row in Mj1T
		  row_j0 = row+Deltasize(j0()+1)-Deltasize(lambda.j());
		  offset = (1<<(lambda.j()-1))-(1<<j0());
		  for (size_type k(0); k < Mj1T.entries_in_row(row_j0); k++)
		    if ((int)Mj1T.get_nth_index(row_j0,k) >= 1<<(j0()-1))
		      c.set_coefficient(Index(lambda.j()-1, 1, Mj1T.get_nth_index(row_j0,k)+offset, this),
					Mj1T.get_nth_entry(row_j0,k));
		  fill_start = 0;
		  fill_end   = (1<<(lambda.j()-1))-(1<<(j0()-1))-1;
		}
	      }
	      
	      // Fill in the missing columns fron the left half of Mj1T:
	      
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
		   col < (1<<(lambda.j()-2)) && col <= fill_end && filter_row >= first_row_left; col++, filter_row -= 2)
		if (col >= fill_start)
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
		   col >= 1<<(lambda.j()-2) && col >= fill_start && filter_row <= last_row_right-offset_right; col--, filter_row += 2)
		if (col <= fill_end)
		  c.set_coefficient(Index(lambda.j()-1, 1, col, this),
				    Mj1T_t.get_nth_entry(col_right, filter_row+offset_right-first_row_right));
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
	      c.add(Mj0T.get_nth_entry(row_j0,k), dhelp);
	    }
	  }
      }
  }
  
  template <int d, int dT>
  void
  PBasis<d, dT>::decompose_t_1(const Index& lambda,
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
	    }
	    else {
	      // Due to the [DS] symmetrization, we have to be a bit careful here.
	      int fill_start = 0, fill_end = (1<<(lambda.j()-1))-1; // first and last column to fill up with interior filter
	      const size_type upper_half = Deltasize(j0()+1)/2;
	      if (row < upper_half) {
		// read the first half of the corresponding row in Mj1
		for (size_type k(0); k < Mj1.entries_in_row(row_j0) && (int)Mj1.get_nth_index(row_j0,k) < 1<<(j0()-1); k++)
		  c.set_coefficient(Index(lambda.j()-1, 1, Mj1.get_nth_index(row_j0,k), this),
				    Mj1.get_nth_entry(row_j0,k));
		fill_start = 1<<(j0()-1);
		fill_end   = (1<<(lambda.j()-1))-1;
	      }
	      else {
		const size_type bottom_half = Deltasize(lambda.j())-Deltasize(j0()+1)/2;
		if (row >= bottom_half) {
		  // read the second half of the corresponding row in Mj1
		  row_j0 = row+Deltasize(j0()+1)-Deltasize(lambda.j());
		  offset = (1<<(lambda.j()-1))-(1<<j0());
		  for (size_type k(0); k < Mj1.entries_in_row(row_j0); k++)
		    if ((int)Mj1.get_nth_index(row_j0,k) >= 1<<(j0()-1))
		      c.set_coefficient(Index(lambda.j()-1, 1, Mj1.get_nth_index(row_j0,k)+offset, this),
					Mj1.get_nth_entry(row_j0,k));
		  fill_start = 0;
		  fill_end   = (1<<(lambda.j()-1))-(1<<(j0()-1))-1;
		}
	      }
	      
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
		   col < (1<<(lambda.j()-2)) && col <= fill_end && filter_row >= first_row_left; col++, filter_row -= 2)
		if (col >= fill_start)
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
		   col >= 1<<(lambda.j()-2) && col >= fill_start && filter_row <= last_row_right-offset_right; col--, filter_row += 2)
		if (col <= fill_end)
		  c.set_coefficient(Index(lambda.j()-1, 1, col, this),
				    Mj1_t.get_nth_entry(col_right, filter_row+offset_right-first_row_right));
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
	      c.add(Mj0.get_nth_entry(row_j0,k), dhelp);
	    }
	  }
      }
  }
  
  template <int d, int dT>
  void
  PBasis<d, dT>::reconstruct_1(const Index& lambda,
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
	c.add(M.get_nth_entry(row_j0,k), dhelp);
      }
    }
  }

  template <int d, int dT>
  void
  PBasis<d, dT>::reconstruct_t_1(const Index& lambda,
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
	c.add(M.get_nth_entry(row_j0,k), dhelp);
      }
    }
  }

}
