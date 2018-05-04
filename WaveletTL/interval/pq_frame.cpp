// implementation for pq_frame.h
#include <iostream>
#include <cassert>
#include <cmath>
#include <numerics/schoenberg_splines.h>
#include <algebra/triangular_matrix.h>
#include <utils/tiny_tools.h>
#include <interval/boundary_gramian.h>

#include "i_q_index.h"

namespace WaveletTL
{
  template <int d, int dT>
  PQFrame<d,dT>::PQFrame(const int s0, const int s1, const bool boundary_quarks) {
    assert(std::max(s0,s1) < d);
    //assert(std::min(s0,s1) >= d-2);
    
    this->s0 = s0;
    this->s1 = s1;
    this->boundary_quarks = boundary_quarks;
    primal = d;
    dual = dT;
    evaluate_with_pre_computation = false;
    setup();

#ifdef _PRE_COMPUTE_WAVELETS
    pre_compute_wavelets();
#endif

  }

  template <int d, int dT>
  PQFrame<d,dT>::PQFrame(const bool bc_left, const bool bc_right, const bool boundary_quarks) {
    this->s0 = bc_left ? 1 : 0;
    this->s1 = bc_right ? 1 : 0;
    this->boundary_quarks = boundary_quarks;
    primal = d;
    dual = dT;
    evaluate_with_pre_computation = false;

    setup();

#ifdef _PRE_COMPUTE_WAVELETS
    pre_compute_wavelets();
#endif

  }

  template <int d, int dT>
  inline
  typename PQFrame<d,dT>::Index
  PQFrame<d,dT>::first_generator(const int j, const int p) const
  {
    assert(j >= j0());
    return Index(p, j, 0, DeltaLmin(p), this);
     
  }
  
  template <int d, int dT>
  inline
  typename PQFrame<d,dT>::Index
  PQFrame<d,dT>::last_generator(const int j, const int p) const
  {
    assert(j >= j0());
    return Index(p, j, 0, DeltaRmax(j,p), this);
  }

  template <int d, int dT>
  inline
  typename PQFrame<d,dT>::Index
  PQFrame<d,dT>::first_wavelet(const int j, const int p) const
  {
    assert(j >= j0());
    return Index(p, j, 1, Nablamin(p), this);
  }
  
  template <int d, int dT>
  inline
  typename PQFrame<d,dT>::Index
  PQFrame<d,dT>::last_wavelet(const int j, const int p) const
  {
    assert(j >= j0());
    return Index(p, j, 1, Nablamax(j,p), this);
  }
  
  
  template <int d, int dT>
  void PQFrame<d,dT>::setup_additional_duals_L(Matrix<double>& MLTp) {
    // setup gramian matrix from [P], Remark 4.2 using [P], Section 2.5
    Matrix<double> GammaL_tmp;
    compute_biorthogonal_boundary_gramian
      <CDFMask_primal<d>,CDFMask_dual<d,dT> >(ML_, MLTp, GammaL_tmp);

//     cout << "GammaL_tmp = " << endl << GammaL_tmp << endl;

    for (unsigned int k = d-2-s0; k >= 1; k--) {
      // ##############################################
      // [P], Prop. 4.15, Remark 4.16, Section 2.5
      // ##############################################
      Matrix<double> B_k(d+k-1,d+k-1); //correct

      unsigned int offset = d-2-s0;

      // setup B_k
      for (unsigned int n = 1; n <= B_k.row_dimension(); n++) {
	double B_k_ns = 0;
	for (unsigned int s = k+1+s0; s <= d+2*k-1+s0; s++) { //correct
	  for (unsigned int l = 1; l <= 2*(n+s0); l++) {
	    if ((offset+1 <= s) && (s <= MLTp.column_dimension())) {
	      if (l-1 < ML_.column_dimension()) {
		B_k_ns += ML_.get_entry(l-1,n-1)*GammaL_tmp.get_entry(l-1,s-1);
	      }
	    }
	    else if ( ((k+1 <= s) && (s <= offset)) 
		      || (s > MLTp.column_dimension())) {
	      if ((l-1) < ML_.row_dimension()) {
		B_k_ns += ML_.get_entry(l-1,n-1) * (1 ? (l==s) : 0);
	      }
	    }
	  }
	  B_k.set_entry(n-1, s-(k+1+s0), 0.5*B_k_ns);
	  B_k_ns = 0;
	}
      }

//       cout << "B_k = " << endl << B_k << endl;

      Matrix<double> inv_B_k(d+k-1,d+k-1);
      QUDecomposition<double>(B_k).inverse(inv_B_k);

//       cout << "inv_B_k = " << endl << inv_B_k << endl;


      // [P] first equation below (4.19)
      for (unsigned int j = 1; j <= inv_B_k.column_dimension(); j++) {
	// extract k-th column of inv_B_k
	MLTp.set_entry(k-1+j+s0, k-1, inv_B_k.get_entry(j-1,k-1));
      }
    }
  }
  
  template <int d, int dT>
  void PQFrame<d,dT>::setup_additional_duals_R(Matrix<double>& MRTp) {
    // setup gramian matrix from [P], Remark 4.2 using [P], Section 2.5
    Matrix<double> GammaR_tmp;
    compute_biorthogonal_boundary_gramian
      <CDFMask_primal<d>,CDFMask_dual<d,dT> >(MR_, MRTp, GammaR_tmp);


//     cout << "GammaR_tmp = " << endl << GammaR_tmp << endl;

    for (unsigned int k = d-2-s1; k >= 1; k--) {
      // ##############################################
      // [P], Prop. 4.15, Remark 4.16, Section 2.5
      // ##############################################
      Matrix<double> B_k(d+k-1,d+k-1); //correct

      unsigned int offset = d-2-s1;

      // setup B_k
      for (unsigned int n = 1; n <= B_k.row_dimension(); n++) {
	double B_k_ns = 0;
	for (unsigned int s = k+1+s1; s <= d+2*k-1+s1; s++) { //correct
	  for (unsigned int l = 1; l <= 2*(n+s1); l++) {
	    if ((offset+1 <= s) && (s <= MRTp.column_dimension())) {
	      if (l-1 < MR_.column_dimension()) {
		B_k_ns += MR_.get_entry(l-1,n-1)*GammaR_tmp.get_entry(l-1,s-1);
	      }
	    }
	    else if ( ((k+1 <= s) && (s <= offset)) 
		      || (s > MRTp.column_dimension())) {
	      if ((l-1) < MR_.row_dimension()) {
		B_k_ns += MR_.get_entry(l-1,n-1) * (1 ? (l==s) : 0);
	      }
	    }
	  }
	  B_k.set_entry(n-1, s-(k+1+s1), 0.5*B_k_ns);
	  B_k_ns = 0;
	}
      }

//       cout << "B_k = " << endl << B_k << endl;

      Matrix<double> inv_B_k(d+k-1,d+k-1);
      QUDecomposition<double>(B_k).inverse(inv_B_k);

//       cout << "inv_B_k = " << endl << inv_B_k << endl;


      // [P] first equation below (4.19)
      for (unsigned int j = 1; j <= inv_B_k.column_dimension(); j++) {
	// extract k-th column of inv_B_k
	MRTp.set_entry(k-1+j+s1, k-1, inv_B_k.get_entry(j-1,k-1));
      }
    }        
  }


  template <int d, int dT>
  void
  PQFrame<d,dT>::setup() {
    // For simplicity, we do not implement a fully generic setup here
    // but only setup the parameters for several important special cases from [P].
    // Namely, we restrict the setting to the case s0 >= d-2, where no "additional"
    // interior dual generators have to be constructed. In a later version of
    // this class, we will fix this.

#if 1
    // choose j0 s.th. the supports of the dual boundary generators do not overlap
    // (at the end of the setup, j0 may be reduced by one, see below)
    j0_ = (int) ceil(log(ell2T<d,dT>()-ell1T<d,dT>()+std::max(s0,s1)+1.-d)/M_LN2+1);
    //cout << "minLev = " << j0_ << endl;
    //j0_ = 5;
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
 
    const int num_dual_bound_gen_L = (s0 >= (d-2)) ? dT : (dT+d-2-s0);
    const int num_dual_bound_gen_R = (s1 >= (d-2)) ? dT : (dT+d-2-s1);

    const int num_rows_boundref_L = (s0 >= (d-2)) ? (3*dT+s0-1) : (3*dT+2*d-5-s0);
    const int num_rows_boundref_R = (s1 >= (d-2)) ? (3*dT+s1-1) : (3*dT+2*d-5-s1);

    // The total number of (left) boundary generators is always exactly dT
    // (to be able reproduce all polynomials by the dual generators).
    ML_.resize(num_rows_boundref_L, num_dual_bound_gen_L);
    for (int column = s0; column < d-1; column++)
      for (int row = s0; row < 2*(d-1); row++)
	ML_.set_entry(row-s0, column-s0, ML_0.get_entry(row,column));
    for (int k = d; k <= num_dual_bound_gen_L+s0; k++)
      for (int n = 2*k-d; n <= 2*k; n++)
 	ML_.set_entry(n-1-s0,k-1-s0,cdf.a().a(-(d/2)+n+d-2*k));
//     cout << "ML=" << endl << ML_;

    MR_.resize(num_rows_boundref_R, num_dual_bound_gen_R);
    for (int column = s1; column < d-1; column++)
      for (int row = s1; row < 2*(d-1); row++)
	MR_.set_entry(row-s1, column-s1, ML_0.get_entry(row,column));
    for (int k = d; k <= num_dual_bound_gen_R+s1; k++)
      for (int n = 2*k-d; n <= 2*k; n++)
 	MR_.set_entry(n-1-s1,k-1-s1,cdf.a().a(-(d/2)+n+d-2*k));
//     cout << "MR=" << endl << MR_;

    Matrix<double> DTilde;
    Matrix<double> inv_DTilde;
    if (s0 < d-2 || s1 < d-2) {
      // We setup the dual boundary generators for the
      // polynomial reproduction the way it is done in [P], Remark 4.2.
      // We need to do this if 'additional' dual boundary functions have
      // to be constructed.
       
      // setup of matrix \tilde{D}_3
      Matrix<double> DTilde_3(dT, dT);

      // setting up devided differences [0,...,l:t^n]_t, for l = 0,...,dT-1, n fixed
      // this computation can be done in place
      for (int n = 0; n <= dT-1; n++) {
	Vector<double> divided_differences(dT);
	for (int i = 0; i <= dT; i++) {
	  double tmp1 = 0.;
	  double tmp2 = 0.;
	  for (int k = i; k <= dT-1; k++) {
	    if (i == 0) {
	      divided_differences[k] = intpower(k,n);
	    }
	    else {
	      if (k == i) {
		if (k%2 == 0) {
		  tmp1 = divided_differences[k-1];
		  tmp2 = divided_differences[k];
		}
		else {
		  tmp1 = divided_differences[k];
		  tmp2 = divided_differences[k-1];
		}
	      }
	      if (k%2 == 0) {		
		tmp2 = divided_differences[k];
		divided_differences[k] = (divided_differences[k]-tmp1)
		  / ((double) i);
		
	      }
	      else {
		tmp1 = divided_differences[k];
		divided_differences[k] = (divided_differences[k]-tmp2)
		  / ((double) i);
	      }
	    }
	  }
	}
	for (int l = 0; l <= dT-1; l++) {
	  DTilde_3.set_entry(n,l, factorial(l)*divided_differences[l]);
	}
      }
//       cout << "DTilde_3 = " << endl << DTilde_3 << endl;
      
      // setup of matrix \tilde{D}_2
      Matrix<double> DTilde_2(dT, dT);
      for (int i = 0; i <= dT-1; i++) {
	for (int n = 0; n <= dT-1; n++) {
	  double tmp = minus1power(n) * binomial(i,n);
	  if (i-n >= 0)
	    tmp *= intpower(-ell1T<d,dT>()-1,i-n);
	  else 
	    tmp /= intpower(-ell1T<d,dT>()-1,n-i);
	  DTilde_2.set_entry(i,n, tmp);
	}
      }
//       cout << "DTilde_2 = " << endl << DTilde_2 << endl;
       
      // setup of matrix \tilde{D}_1
      Matrix<double> DTilde_1(dT, dT);
      for (int r = 0; r <= dT-1; r++) {
	for (int i = 0; i <= dT-1; i++) {
	  if (r >= i)
	    DTilde_1.set_entry(r, i, binomial(r,i)*alpha(0,r-i));
	}
      }
//       cout << "DTilde_1 = " << endl << DTilde_1 << endl;
       
      // for readability
      DTilde.resize(dT, dT);
      DTilde = DTilde_1 * DTilde_2 * DTilde_3;

//       cout << "DTilde = " << endl << DTilde << endl;

      // now we put together the refinement matrix MLTp
      inv_DTilde.resize(dT, dT);
      QUDecomposition<double>(DTilde).inverse(inv_DTilde);
       
//       cout << "inv_DTilde" << endl << inv_DTilde << endl;
    }

    // now setting up MLTp
    Matrix<double> MLTp(num_rows_boundref_L, num_dual_bound_gen_L);
    if (s0 >= d-2) {
      // setup the expansion coefficients of the unbiorthogonalized dual generators
      // w.r.t. the truncated CDF generators, see [DKU, Lemma 3.1]
      for (int r = 0; r < num_dual_bound_gen_L; r++) {
	MLTp.set_entry(r, r, ldexp(1.0, -r));
	for (int m = ellT_l(); m <= 2*ellT_l()+ell1T<d,dT>()-1; m++)
	  MLTp.set_entry(dT+m-ellT_l(), r, alpha(m, r)*ldexp(1.0, -r));
	for (int m = 2*ellT_l()+ell1T<d,dT>(); m <= 2*ellT_l()+ell2T<d,dT>()-2; m++)
	  MLTp.set_entry(dT+m-ellT_l(), r, betaL(m, r));
      }
    }
    else {
      int helpDim = 3*dT+d-3;
      
      // for simplicity
      Matrix<double> helpMat(helpDim, dT);
      // compute (\tilde{mu}_{n,k}^L)n,k * DTilde^-T * J_dT
      for (int n = 0; n <= helpDim-1; n++) {
	for (int k = 0; k <= dT-1; k++) {
	  double tmp = 0.;
	  for (int i = 0; i <= dT-1; i++) {
	    if ((0 <= n) && (n <= dT-1) && (i == n))
	      tmp += ldexp(1.0, -i) * inv_DTilde.get_entry(dT-1-k,i);// (i,k) -- ^T --> (k,i) 
	    // -- *J_dT from the right --> (dT-1-k,i)
	    else if ((dT <= n) && (n <= 3*dT+d-4)) {
	      // ATTENTION: mistake in [P], Remark 4.2, should be floor(d/2) at \beta_{...} not ceil(d/2) 
	      tmp += betaL(n-ell1<d>()-1, i) * inv_DTilde.get_entry(dT-1-k,i);
	    }
	  }
	  helpMat.set_entry(n, k, tmp);
	}
      }

//       cout << "helpMat" << endl << helpMat << endl;

      unsigned int shift = d-1-s0;
      for (unsigned int n = shift; n <= MLTp.row_dimension(); n++) {
	for (unsigned int k = shift; k <= MLTp.column_dimension(); k++) {
	  if (n <= shift + dT-1) {
	    double tmp = 0.;
	    for (int i = 0; i <= dT-1; i++) {
	      tmp += DTilde.get_entry(i,dT-(n-shift)-1)
		* helpMat.get_entry(i,k-shift);
	    }
	    MLTp.set_entry(n-1, k-1, tmp);
 	  }
	  else {
	    MLTp.set_entry(n-1, k-1, helpMat.get_entry(n-shift,k-shift));
	  }
	}
      }

    }
//     cout << "MLTp=" << endl << MLTp;

    // now setting up MRTp
    Matrix<double> MRTp(num_rows_boundref_R, num_dual_bound_gen_R);     
    if (s1 >= d-2) {
      for (int r = 0; r < num_dual_bound_gen_R; r++) {
	MRTp.set_entry(r, r, ldexp(1.0, -r));
	for (int m = ellT_r(); m <= 2*ellT_r()+ell1T<d,dT>()-1; m++)
	  MRTp.set_entry(dT+m-ellT_r(), r, alpha(m, r)*ldexp(1.0, -r));
	for (int m = 2*ellT_r()+ell1T<d,dT>(); m <= 2*ellT_r()+ell2T<d,dT>()-2; m++)
	  MRTp.set_entry(dT+m-ellT_r(), r, betaR(m, r));
      }
    }
    else {
      int helpDim = 3*dT+d-3;
      
      // for simplicity
      Matrix<double> helpMat(helpDim, dT);
      // compute (\tilde{mu}_{n,k}^L)n,k * DTilde^-T * J_dT
      for (int n = 0; n <= helpDim-1; n++) {
	for (int k = 0; k <= dT-1; k++) {
	  double tmp = 0.;
	  for (int i = 0; i <= dT-1; i++) {
	    if ((0 <= n) && (n <= dT-1) && (i == n))
	      tmp += ldexp(1.0, -i) * inv_DTilde.get_entry(dT-1-k,i);// (i,k) -- ^T --> (k,i) 
	    // -- *J_dT from the right --> (dT-1-k,i)
	    else if ((dT <= n) && (n <= 3*dT+d-4)) {
	      // ATTENTION: mistake in [P], Remark 4.2, should be floor(d/2) at \beta_{...} not ceil(d/2) 
	      tmp += betaR(n-ell1<d>()-1, i) * inv_DTilde.get_entry(dT-1-k,i);
	    }
	  }
	  helpMat.set_entry(n, k, tmp);
	}
      }

//       cout << "helpMat" << endl << helpMat << endl;

      unsigned int shift = d-1-s1;
      for (unsigned int n = shift; n <= MRTp.row_dimension(); n++) {
	for (unsigned int k = shift; k <= MRTp.column_dimension(); k++) {
	  if (n <= shift + dT-1) {
	    double tmp = 0.;
	    for (int i = 0; i <= dT-1; i++) {
	      tmp += DTilde.get_entry(i,dT-(n-shift)-1)
		* helpMat.get_entry(i,k-shift);
	    }
	    MRTp.set_entry(n-1, k-1, tmp);
	  }
	  else {
	    MRTp.set_entry(n-1, k-1, helpMat.get_entry(n-shift,k-shift));
	  }
	}
      }     
    }

//     cout << "MRTp=" << endl << MRTp;
    
    // determine the two-scale relations for the additional
    // dual boundary generators, if necessary
    if (s0 < (d-2)) {
      setup_additional_duals_L(MLTp);
    }
//     cout << "MLTp=" << endl << MLTp;

    if (s1 < (d-2))
      setup_additional_duals_R(MRTp);
//     cout << "MRTp=" << endl << MRTp;

    MLT_BB = MLTp;
    MRT_BB = MRTp;

    // for the biorthogonalization of the generators,
    // compute the gramian matrix of the primal and dual boundary generators
    compute_biorthogonal_boundary_gramian
      <CDFMask_primal<d>,CDFMask_dual<d,dT> >(ML_, MLTp, GammaL);
    //GammaL.compress(1.0e-12);
//     cout << "GammaL=" << endl << GammaL;
 
    compute_biorthogonal_boundary_gramian
      <CDFMask_primal<d>,CDFMask_dual<d,dT> >(MR_, MRTp, GammaR);
//     cout << "GammaR=" << endl << GammaR;
 
    setup_CXT();
    setup_CXAT();
    
    setup_Mj0(ML_, MR_, Mj0); // [DKU, (3.5.1)]

    Mj0.compress(1e-8);
   
    SparseMatrix<double> mj0tp;
    setup_Mj0Tp(MLTp, MRTp, mj0tp); // [DKU, (3.5.5)]
    mj0tp.compress(1e-8);
    
    //cout << "mj0tp = " << endl << mj0tp << endl;
    
    // biorthogonalize the generators
    setup_Cj();
    //cout << "CjT = " << endl << CjT << endl;
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
    //cout << "Mj0^T*Mj0T=" << endl << testbio0 << endl;
    for (unsigned int i = 0; i < testbio0.row_dimension(); i++)
      testbio0.set_entry(i, i, testbio0.get_entry(i, i) - 1.0);
    cout << "* ||Mj0^T*Mj0T-I||_infty: " << row_sum_norm(testbio0) << endl;
    testbio0.compress(1.0e-8);
    //cout << testbio0 << endl;

    testbio0 = transpose(Mj0T) * Mj0;
    //cout << "Mj0T*Mj0^T=" << endl << testbio0 << endl;
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
#if 0
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

#if 1
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
    //if (s0 == s1)
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

//      cout << "Mj1:" << endl << Mj1;
//      cout << "Mj1 reduced:" << endl << Mj1_reduced;

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

    if (row_sum_norm(test_reduced) < 1e-10
	&& (!(d==2 && dT==2 && s0==1 && s1==0)) // re-expand fails here since main band is missing, use SplineBasis
	&& (!(d==2 && dT==2 && s0==0 && s1==1))) { 
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
//    cout << "minLev = " << j0_ << endl;
    
    setup_CVM();
  }

  template <int d, int dT>
  void
  PQFrame<d,dT>::setup_full_collection()
  {
    if (jmax_ == -1 || jmax_ < j0_ || pmax_ == -1) {
      cout << "PQFrame<d,dT>::setup_full_collection(): specify a maximal level of resolution and polynomial degree first!" << endl;
      abort();
    }   
    
    int degrees_of_freedom;
    if (boundary_quarks == true)
        degrees_of_freedom = Deltasize(jmax_+1)*(pmax_+1);
    else
        degrees_of_freedom = Deltasize(jmax_+1)*(pmax_+1)-pmax_*(2*d-2-s0-s1)-pmax_*(jmax_-j0_+1)*(d+dT-2);
    
    cout << "total degrees of freedom between (j0_, 0) and (jmax_, pmax_) is " << degrees_of_freedom << endl;

    cout << "setting up collection of quarklet indices..." << endl;
    full_collection.resize(degrees_of_freedom);
    int k = 0;
    int p = 0;
    for (Index ind = first_generator(j0_); ind <= last_wavelet(jmax_, pmax_); (ind == last_wavelet(jmax_,p))? (++p, ind=first_generator(j0_, p)) : ++ind) {
      full_collection[k] = ind;
//      cout << ind << endl;
      k++;
    }
    cout << "done setting up collection of quarklet indices..." << endl;
    //cout << k << endl;

  }



  template <int d, int dT>
  const double
  PQFrame<d,dT>::alpha(const int m, const unsigned int r) const {
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
	  result += cdf.a().a(k) * help;
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
  PQFrame<d,dT>::betaL(const int m, const unsigned int r) const {
    // [DKU] (3.2.31)
    double result = 0;
    for (int q = (int)ceil((m-ell2T<d,dT>())/2.0); q < ellT_l(); q++)
      result += alpha(q, r) * cdf.aT().a(m-2*q);
    return result;
  }

  template <int d, int dT>
  const double
  PQFrame<d,dT>::betaR(const int m, const unsigned int r) const {
    // [DKU] (3.2.31)
    double result = 0;
    for (int q = (int)ceil((m-ell2T<d,dT>())/2.0); q < ellT_r(); q++)
      result += alpha(q, r) * cdf.aT().a(m-2*q);
    return result;
  }

    template <int d, int dT>
    void
    PQFrame<d,dT>::support(const Index& lambda, int& k1, int& k2) const 
    {
        if (lambda.e() == 0) // generator
        {
            // phi_{j,k}(x) = 2^{j/2} B_{j,k-ell_1}
            k1 = std::max(0, lambda.k() + ell1<d>());
            k2 = std::min(1<<lambda.j(), lambda.k() + ell2<d>());
        }
        else // wavelet
        {
            // cf. [P, p. 125]
            if (lambda.k() < (d+dT)/2-1) 
            {
                // left boundary wavelet
                k1 = 0;

                // hard code special case
                //  	  if ( (d==3 && dT==3) &&  ((basis.get_s0()==1) && (basis.get_s1()==1)) ) {
                // 	    k2 = 8;
                //  	  }
                //  	  else
                //k2 = 2*(d+dT)-2; // overestimate, TODO 
                //k2 = 2*(d+dT)-4; // thats the formula in [P] for free BC
                if (s0 == (d-1) )
                {
                    k2 = 2*(d+dT)-2;
                }
                else
                {
                    k2 = 2*(d+dT)-4;
                }
            } 
            else 
            {
                if ((1<<lambda.j())-lambda.k() <= (d+dT)/2-1) 
                {
                    // right boundary wavelet
                    // hard code special case
                    // 	    if ( (d==3 && dT==3) &&  ((basis.get_s0()==1) && (basis.get_s1()==1)) ) {
                    // 	      k1 = (1<<(lambda.j()+1)) - 8;
                    //  	    }
                    //  	    else
                    //k1 = (1<<(lambda.j()+1))-(2*(d+dT)-2); // overestimate, TODO
                    //k1 = (1<<(lambda.j()+1))-2*(d+dT)+4; // thats the formula in [P] for free BC
                    if (s1 == (d-1) )
                    {
                        k1 = (1<<(lambda.j()+1)) - 2*(d+dT) + 2;
                    }
                    else
                    {
                        k1 = (1<<(lambda.j()+1)) - 2*(d+dT) + 4;
                    }
                    k2 = 1<<(lambda.j()+1);
                } else 
                {
                    // interior wavelet (CDF)
                    // note: despite the fact that the wavelets in the right half of the interval
                    //       are reflected CDF wavelets, their support does not "see" the reflection!
                    // k1 = 2*(lambda.k()-(d+dT)/2+1);
                    k1 = 2*lambda.k()-(d+dT)+2;
                    k2 = k1+2*(d+dT)-2;
                }
            }
        }
    }
  
    template <int d, int dT>
    void
    PQFrame<d,dT>::support(const int j_, const int e_, const int k_, int& k1, int& k2) const 
    {
        if (e_ == 0) // generator
        {
            // phi_{j,k}(x) = 2^{j/2} B_{j,k-ell_1}
            k1 = std::max(0, k_ + ell1<d>());
            k2 = std::min(1<<j_, k_ + ell2<d>());
        }
        else // wavelet
        {
            // cf. [P, p. 125]
            if (k_ < (d+dT)/2-1) {
                // left boundary wavelet
                k1 = 0;
                //k2 = 2*(d+dT)-2; // overestimate, TODO
                //k2 = 2*(d+dT)-4; // thats the formula in [P]
                if (s0 == (d-1) )
                {
                    k2 = 2*(d+dT)-2;
                }
                else
                {
                    k2 = 2*(d+dT)-4;
                }
            } else {
                if ((1<<j_)-k_ <= (d+dT)/2-1) 
                {
                    // right boundary wavelet
                    //k1 = (1<<(j_+1))-(2*(d+dT)-2); // overestimate, TODO
                    //k1 = (1<<(j_+1))-2*(d+dT)+4; // thats the formula in [P]
                    if (s1 == (d-1) )
                    {
                        k1 = (1<<(j_+1)) - 2*(d+dT) + 2;
                    }
                    else
                    {
                        k1 = (1<<(j_+1)) - 2*(d+dT) + 4;
                    }
                    k2 = 1<<(j_+1);
                } else {
                    // interior wavelet (CDF)
                    // note: despite the fact that the wavelets in the right half of the interval
                    //       are reflected CDF wavelets, their support does not "see" the reflection!
                    //k1 = 2*(k_-(d+dT)/2+1);
                    k1 = 2*k_-(d+dT)+2;
                    k2 = k1+2*(d+dT)-2;
                }
            }
        }
    }
  
  template <int d, int dT>
  void
  PQFrame<d,dT>::setup_Mj0(const Matrix<double>& ML, const Matrix<double>& MR, SparseMatrix<double>& Mj0) {
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
    
    int startrow = (s0 >= d-2) ? d+ell_l()+ell1<d>()-2*s0 : d+ell_l()+ell1<d>()-2*s0-(d-2-s0);
    for (int col = d-s0; col < nj-(d-s1); col++, startrow+=2) {
      int row = startrow;
      for (int k = cdf.a().begin(); k <= cdf.a().end(); k++, row++)
	Mj0.set_entry(row, col, cdf.a().a(k));
    }

    Mj0.scale(M_SQRT1_2);
    
    //cout << "Mj0=" << endl << Mj0 << endl;
  }

  template <int d, int dT>
  void
  PQFrame<d,dT>::setup_Mj0Tp(const Matrix<double>& MLTp, const Matrix<double>& MRTp, SparseMatrix<double>& Mj0Tp) {
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

  
    //int startrow = (s0 >= d-2) ? dT+ellT_l()+ell1T<d,dT>() : dT+ellT_l()+ell1T<d,dT>()+2*(d-2-s0)-2;
    int startrow = (s0 >= d-2) ? dT+ellT_l()+ell1T<d,dT>() : dT+ellT_l()+ell1T<d,dT>()+d-2-s0;
    //int startrow = dT+ellT_l()+ell1T<d,dT>()+4;
    //cout << "startrow = " << startrow << endl;
    //for (int col = dT; col < nj-dT; col++, startrow+=2) {

    for (unsigned int col = MLTp.column_dimension(); col < nj-MRTp.column_dimension(); col++, startrow+=2) {
      int row = startrow;
      for (int k = cdf.aT().begin(); k <= cdf.aT().end(); k++, row++)
	Mj0Tp.set_entry(row, col, cdf.aT().a(k));
    }

    Mj0Tp.scale(M_SQRT1_2);

    //cout << "Mj0Tp=" << endl << Mj0Tp << endl;
  }

  template <int d, int dT>
  void
  PQFrame<d,dT>::setup_CXT()
  {
    // IGPMlib reference: I_Mask_Bspline::EvalCL(), ::EvalCR()

    Matrix<double> CLGammaLInv;
    QUDecomposition<double>(GammaL).inverse(CLGammaLInv);
//     cout << "CLGammaLInv = " << endl << CLGammaLInv << endl;
    CLT = transpose(CLGammaLInv);
    
    Matrix<double> CRGammaRInv;
    QUDecomposition<double>(GammaR).inverse(CRGammaRInv);
    CRT = transpose(CRGammaRInv);
    
    QUDecomposition<double>(CLT).inverse(inv_CLT);
    QUDecomposition<double>(CRT).inverse(inv_CRT);
  }

  template <int d, int dT>
  void PQFrame<d,dT>::setup_CXAT() {
    // IGPMlib reference: I_Mask_Bspline::EvalCL(), ::EvalCR()

    if (s0 >= d-2) {
      // setup CLAT <-> Alpha * (CLT)^T
      CLAT.resize(ellT_l()+ell2T<d,dT>()-1, dT);
      for (int i = 1-ell2T<d,dT>(); i <= ellT_l()-1; i++)
	for (unsigned int r = 0; r < dT; r++)
	  CLAT(ell2T<d,dT>()-1+i, r) = alpha(i, r);
      
    }
    else {
      // detrmine row dimension of CLAT first
      int index_last_CDF_used = (MLT_BB.row_dimension()-1) + DeltaLTmin();
      int index_first_CDF_used = 1-ell2T<d,dT>();
      int row_dim_CLAT = index_last_CDF_used-index_first_CDF_used+1;

      CLAT.resize(row_dim_CLAT, dT+d-2-s0);
      Array1D<Vector<double> > CLAT_help(dT+d-2-s0);

      for (unsigned int c = 0; c < CLAT_help.size(); c++) {
	CLAT_help[c].resize(row_dim_CLAT);
      }

      // setup the columns corresponding to the polynomial reproducing
      // dual boundary generators

      // see [P], equation (4.5)
      unsigned int i = dT-1;
      for(unsigned int c = d-2-s0; c < CLAT_help.size(); c++, i--) {
	int n = i;
	int num_rows = (int)(ellT_l()-1-i)-(1-ell2T<d,dT>())+1;
	int s = num_rows-1;
	for (int k = (int)(ellT_l()-1-i); k >= 1-ell2T<d,dT>(); k--, n++, s--) {
	  CLAT_help[c][s] = binomial(n,i);
	}
	//cout << CLAT_help[c] << endl;
      }

      // The matrix CLAT will be set up column-wise. We start with the last column
      // belonging to the right most polynomial reproducing dual left boundary generator.
      // We have to determine the coefficients with respect to the restricted dual CDF 
      // functions on the real line on the level one higher than the level of the boundary
      // functions we are currently dealing with (before biorthogonalization).
      // This is due to the fact that the 'additional' dual boundary generators are only
      // given by their two scale relation.
      
      // loop over all dual left boundary generators
      for(int c = CLAT.column_dimension()-1; c >= 0; c--) {
	// get two scale coefficients
	Vector<double> mask(MLT_BB.row_dimension());
	for (unsigned int r = 0; r < MLT_BB.row_dimension(); r++) {
	  mask[r] = MLT_BB(r,c);
	}
	
	// The following scaling has to be done because of the fact
	// that the columns of CLAT represent the expansions
	// of the dual boundary generators with respect to CDF generators
	// of level j+1 no level j.
	// In [P], eq. (4.5) on the left and on the right we then have 
	// functions of different levels.
	mask.scale(1./M_SQRT2);
	//cout << "mask = " << endl << mask << endl;

	// polynomial reproducing dual boundary generator
	if (d-2-s0 <= c) {
  
	  Vector<double> col_of_CLAT(CLAT.row_dimension());
	  for (unsigned int r = 0; r < mask.size(); r++) {
	    // row corresponding to boundary generator
 	    if ( ((int)r) < dT+d-2-s0 ) {
	      if (mask[r] != 0) {
		col_of_CLAT.add(mask[r], CLAT_help[r]);
	      }
	    }
	    // row corresponding to inner generator
 	    else {
	      if (mask[r] != 0) {
		int tmp = r+DeltaLTmin();
		tmp -= index_first_CDF_used;
		col_of_CLAT[tmp] += mask[r];
	      }
	    }
	  }
	  for (unsigned int r = 0; r < CLAT.row_dimension(); r++) {
	    CLAT(r,c) = col_of_CLAT[r];
	  }
	}
	// 'additional' dual boundary generator
	else {
	  for (unsigned int r = 0; r < mask.size(); r++) {
 	    // row corresponding to boundary generator
 	    if ( ((int)r) < dT+d-2-s0 ) {
	      if (mask[r] != 0) {
		CLAT_help[c].add(mask[r], CLAT_help[r]);
	      }
	    }
 	    // row corresponding to inner generator
 	    else {
	      if (mask[r] != 0) {
		int tmp = r+DeltaLTmin();
		tmp -= index_first_CDF_used;
		CLAT_help[c][tmp] += mask[r];
	      }
 	    }
 	  } 
	  for (unsigned int r = 0; r < CLAT.row_dimension(); r++) {
	    CLAT(r,c) = CLAT_help[c][r];
	  }
	}// end 'additional' dual boundary generator
      }
    }

//     cout << "CLAT before biorthogonalization:" << endl << CLAT << endl;
    CLAT = CLAT * transpose(CLT);
//     CLAT.compress(1e-12);
//     cout << "CLAT after biorthogonalization:" << endl << CLAT << endl;
 

    // now we treat the right end of the interval
    if (s1 >= d-2) {
      // the same for CRAT:
      CRAT.resize(ellT_r()+ell2T<d,dT>()-1, dT);
      for (int i = 1-ell2T<d,dT>(); i <= ellT_r()-1; i++)
	for (unsigned int r = 0; r < dT; r++)
	  CRAT(ell2T<d,dT>()-1+i, r) = alpha(i, r);    
    }
    // see code for left end above for more comments!
    else {
      // determine row dimension of CRAT first
      int offset = ellT_r()-(dT+d-2-s1);
      int index_last_CDF_used = (MRT_BB.row_dimension()-1) + offset;
      int index_first_CDF_used = 1-ell2T<d,dT>();
      int row_dim_CRAT = index_last_CDF_used-index_first_CDF_used+1;

      CRAT.resize(row_dim_CRAT, dT+d-2-s1);
      Array1D<Vector<double> > CRAT_help(dT+d-2-s1);

      for (unsigned int c = 0; c < CRAT_help.size(); c++) {
	CRAT_help[c].resize(row_dim_CRAT);
      }

      // setup the columns corresponding to the polynomial reproducing
      // dual boundary generators

      // see [P], equation (4.5)
      unsigned int i = dT-1;
      for(unsigned int c = d-2-s1; c < CRAT_help.size(); c++, i--) {
	int n = i;
	int num_rows = (int)(ellT_r()-1-i)-(1-ell2T<d,dT>())+1;
	int s = num_rows-1;
	for (int k = (int)(ellT_r()-1-i); k >= 1-ell2T<d,dT>(); k--, n++, s--) {
	  CRAT_help[c][s] = binomial(n,i);
	}
	//cout << CRAT_help[c] << endl;
      }

      // loop over all dual left boundary generators
      for(int c = CRAT.column_dimension()-1; c >= 0; c--) {
	// get two scale coefficients
	Vector<double> mask(MRT_BB.row_dimension());
	for (unsigned int r = 0; r < MRT_BB.row_dimension(); r++) {
	  mask[r] = MRT_BB(r,c);
	}
	
	// The following scaling has to be done because of the fact
	// that the columns of CLAT represent the expansions
	// of the dual boundary generators with respect to CDF generators
	// of level j+1 no level j.
	// In [P], eq. (4.5) on the left and on the right we then have 
	// functions of different levels.
	mask.scale(1./M_SQRT2);
	//cout << "mask = " << endl << mask << endl;

	// polynomial reproducing dual boundary generator
	if (d-2-s1 <= c) {
  
	  Vector<double> col_of_CRAT(CRAT.row_dimension());
	  for (unsigned int r = 0; r < mask.size(); r++) {
	    // row corresponding to boundary generator
 	    if ( ((int)r) < dT+d-2-s1 ) {
	      if (mask[r] != 0) {
		col_of_CRAT.add(mask[r], CRAT_help[r]);
	      }
	    }
	    // row corresponding to inner generator
 	    else {
	      if (mask[r] != 0) {
		int tmp = r+offset;
		tmp -= index_first_CDF_used;
		col_of_CRAT[tmp] += mask[r];
	      }
	    }
	  }
	  for (unsigned int r = 0; r < CRAT.row_dimension(); r++) {
	    CRAT(r,c) = col_of_CRAT[r];
	  }
	}
	// 'additional' dual boundary generator
	else {
	  for (unsigned int r = 0; r < mask.size(); r++) {
 	    // row corresponding to boundary generator
 	    if ( ((int)r) < dT+d-2-s1 ) {
	      if (mask[r] != 0) {
		CRAT_help[c].add(mask[r], CRAT_help[r]);
	      }
	    }
 	    // row corresponding to inner generator
 	    else {
	      if (mask[r] != 0) {
		int tmp = r+offset;
		tmp -= index_first_CDF_used;
		CRAT_help[c][tmp] += mask[r];
	      }
 	    }
 	  } 
	  for (unsigned int r = 0; r < CRAT.row_dimension(); r++) {
	    CRAT(r,c) = CRAT_help[c][r];
	  }
	}// end 'additional' dual boundary generator
      }
    }
  
//     cout << "CRAT before biorthogonalization:" << endl << CRAT << endl;
    CRAT = CRAT * transpose(CRT);
//     CRAT.compress(1e-12);
//     cout << "CRAT after biorthogonalization:" << endl << CRAT << endl;

  }

  template <int d, int dT>
  void PQFrame<d,dT>::setup_Cj() {
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
    //cout << "inv_CjpT= " << endl << inv_CjpT << endl;

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
  void PQFrame<d,dT>::setup_CVM() {
      //for fast implementation, can be changed @PHK
      const int pmax=6;
      CVM.resize(pmax, (dT+1)*(d+dT-2));
      CVM(0,0)= 27;
      
      if(s0==1 && s1==1){
            switch(d){
              case 3: switch(dT){
                      case 3: {
                          
                          Matrix <double> tempmat(pmax, (dT+1)*(d+dT-2),"         0 0.4082482904638632 -0.816496580927726 0.408248290463863            -0.5477225575051661 0.7302967433402217 0.1825741858350556 -0.3651483716701108             -0.5477225575051659 0.7302967433402217 0.1825741858350552 -0.3651483716701112            0.4082482904638641 -0.816496580927726 0.408248290463862 0     "
                                                                         "        -0.5224552221805724 0.8055139979119766 -0.2787671792634939 0.02183116464111837          0.2236067977499873 -0.6708203932499458 0.6708203932499281 -0.2236067977499699          0.2236067977499817 -0.6708203932499401 0.670820393249934 -0.2236067977499756             -0.02183116464265128 0.2787671792660461 -0.8055139979117942 0.5224552221794281     "
                                                                         "       0 0.4082482904639252 -0.8164965809277187 0.4082482904638156           0.547722557504977 -0.730296743340497 -0.1825741858346003 0.3651483716700706          0.5477225575051627 -0.7302967433402219 -0.1825741858350566 0.3651483716701149            0.408248290463864 -0.816496580927726 0.4082482904638622 0     "
                                                                         "       -0.640146461278297 0.6982276160627079 0.2446894817047236 -0.2069245317098367              0.2236067977501682 -0.6708203932500731 0.6708203932497893 -0.2236067977498231            -0.223606796596586 0.6708203920324677 -0.670820394435368 0.2236067989994864             0.2069245144113484 -0.2446894480709285 -0.698227632722497 0.6401464615547835    "
                                                                         "       0 0.4082482904621896 -0.8164965809275739 0.4082482904658408             -0.5477225575012672 0.7302967433437152 0.1825741858292311 -0.3651483716718835                  0.5477225393344308 -0.7302967658399765 -0.1825741454060958 0.36514837414118              0.4082482532810987 -0.8164965800968248 0.4082483293084262 0     "
                                                                         "        -0.5874404691592978 0.4213433865671868 0.6043779008557463 -0.3348295069840767            -0.2236067977542563 0.6708203932589307 -0.6708203932434262 0.2236067977382518                    0.2236044305743063 -0.6708180884643348 0.6708227288772003 -0.2236090723524933             0.3348310956040639 -0.6043812814279892 -0.4213396189204375 0.5874387879634988      ");
                          
//                          Matrix <double> tempmat(pmax, (dT+1)*(d+dT-2)," 0.2025455652183159 -0.6548973275392216 0.6886549217422743 -0.2363031594213687       0.2236067977499787 -0.6708203932499367 0.6708203932499369 -0.223606797749979 0.2236067977499784 -0.6708203932499365 0.6708203932499374 -0.2236067977499794       -0.06398930756904883 0.1706381535174648 -0.7465419216389197 0.6398930756905037"
//                                                                        " -0.2631863804043411 0.683532627964419 -0.647689149490306 0.2097971432366053        0.2236067977499833 -0.6708203932499415 0.6708203932499324 -0.2236067977499744 0.2236067977499784 -0.6708203932499365 0.6708203932499374 -0.2236067977499794        0.05359750604883835 -0.1585898809126361 -0.6145357885456364 0.7709230322243936"
//                                                                        " 0.2167794512514804 -0.6656161436521536 0.6767648011450804 -0.2279281087444072       -0.2236067977499787 0.6708203932499448 -0.6708203932499331 0.2236067977499671 0.2236067977369626 -0.6708203932369206 0.6708203932629532 -0.2236067977629952       0.07505866259329713 -0.4503519751849037 -0.4128226436599776 0.788115956251584"
//                                                                        " -0.2337005659921087 0.6749399235427884 -0.6645436634817194 0.2195947731375184      0.2236067977502918 -0.6708203932501772 0.6708203932496646 -0.2236067977497625 -0.2236067987706151 0.6708203942185108 -0.6708203922553321 0.2236067968074362        0.1040966707159132 -0.5950659374191536 -0.2527858014715041 0.755751117659601"
//                                                                        " 0.2215835947572188 -0.6692635936886583 0.6726040031334517 -0.2249240042020122       0.2236067977492572 -0.6708203932492558 0.6708203932506384 -0.2236067977506397 0.2236068066976706 -0.6708204020865624 0.6708203843577781 -0.2236067889688864       0.07299534619159279 -0.6558770463101348 -0.1527558639038796 0.7356375640224215"
//                                                                        " -0.2262823111496524 0.6720199792759201 -0.6690975041966392 0.222472454458058       -0.2236067977567912 0.6708203932563452 -0.6708203932433455 0.2236067977437161 -0.2236069275296845 0.670820516703335 -0.6708202665819002 0.2236066776140496        0.05914854058002748 -0.6920237805634262 -0.07937587026932126 0.7150552485654124");
                          //cout << "setup: " << endl << tempmat << endl;
                          CVM = tempmat;
                      }
                      }
            }
      }
      
      if(s0==0 && s1==0){
            switch(d){
              case 3: switch(dT){
                      case 3: {
                          //done for symmetric quarklets
                          Matrix <double> tempmat(pmax, (dT+1)*(d+dT-2)," 0.8547043234472845 -0.2849014411490949 0.3988620176087329 -0.170940864689457        0 0.4082482904638632 -0.816496580927726 0.408248290463863     0.4082482904638641 -0.816496580927726 0.408248290463862 0     0.1709408646894492 -0.398862017608727 -0.2849014411490961 0.8547043234472884 "
                                                                         "  -0.1142707226510252 0.5684235947255913 -0.7840729584977717 0.2215094008312083      -0.5224552221805724 0.8055139979119766 -0.2787671792634939 0.02183116464111837     -0.02183116464265128 0.2787671792660461 -0.8055139979117942 0.5224552221794281     0.2215094008295896 -0.7840729584956645 0.5684235947274756 -0.1142707226592489"
                                                                         "  0.9152081163589021 -0.2440554976957071 0.3017169614370021 -0.1086181061173199     0 0.4082482904639252 -0.8164965809277187 0.4082482904638156     0.408248290463864 -0.816496580927726 0.4082482904638622 0       0.1086181061610587 -0.3017169614787351 -0.2440554976909896 0.915208116341211"
                                                                         "  -0.6606713869121563 -0.22074441926776 0.6796712888505161 -0.2298524722294399     -0.640146461278297 0.6982276160627079 0.2446894817047236 -0.2069245317098367     0.2069245144113484 -0.2446894480709285 -0.698227632722497 0.6401464615547835      0.229852470662522 -0.6796713154776119 0.2207444536824991 0.6606713485657464"
                                                                         "  0.929541665392022 -0.2360740737503546 0.2705063410924056 -0.08395024378724425     0 0.4082482904621896 -0.8164965809275739 0.4082482904658408     0.4082482532810987 -0.8164965800968248 0.4082483293084262 0       0.0839502616797303 -0.2705063590741015 -0.2360740721163697 0.9295416589582056"
                                                                         "  -0.8234805555087669 -0.06710692648984101 0.525638845099171 -0.2026826081312094      -0.5874404691592978 0.4213433865671868 0.6043779008557463 -0.3348295069840767     0.3348310956040639 -0.6043812814279892 -0.4213396189204375 0.5874387879634988        0.2026823648370047 -0.525640966362598 0.06710930968033792 0.8234790671411762");
                          //cout << "setup: " << endl << tempmat << endl;
                          
//                          Matrix <double> tempmat(pmax, (dT+1)*(d+dT-2)," -0.8135849383260662 0.4881509629956393 -0.3037383769750648 0.08678239342144715   0.2025455652183159 -0.6548973275392216  0.6886549217422743 -0.2363031594213687 -0.06398930756904883  0.1706381535174648 -0.7465419216389197 0.6398930756905037 0.03848645455033218  0.5388103637046633 -0.7216210228187427  0.4329726136912416 "
//                                                                         " -0.804634462931811  0.4825884114596989 -0.3322361653395849 0.09638950475901535  -0.2631863804043411  0.683532627964419  -0.647689149490306   0.2097971432366053  0.05359750604883835 -0.1585898809126361 -0.6145357885456364 0.7709230322243936 0.02868473974068089 -0.6740913839391539  0.551284841922291   0.4907779690162553"
//                                                                         " -0.8292748489821371 0.4470873098860249 -0.3214083886738285 0.09546103851637334   0.2167794512514804 -0.6656161436521536  0.6767648011450804 -0.2279281087444072  0.07505866259329713 -0.4503519751849037 -0.4128226436599776 0.788115956251584 -0.06123260614795069 -0.4566931886791215  0.2946819175241062  0.8371645398861219"
//                                                                         " -0.8111748791095182 0.4568712663959324 -0.34939052424605   0.1057838504238986   -0.2337005659921087  0.6749399235427884 -0.6645436634817194  0.2195947731375184  0.1040966707159132  -0.5950659374191536 -0.2527858014715041 0.755751117659601 0.01110242632083908  -0.29088348746909    0.03759557501427863 0.9559550750882365"
//                                                                         " -0.810199380138568  0.452767950065765  -0.355958929372318  0.1089559012465069    0.2215835947572188 -0.6692635936886583  0.6726040031334517 -0.2249240042020122  0.07299534619159279 -0.6558770463101348 -0.1527558639038796 0.7356375640224215 0.04944480373238666 -0.2191607550946736 -0.07698606029340599 0.9713891708942911"
//                                                                         " -0.7965210936116699 0.4602213845988072 -0.3745463076927572 0.1160408892733443   -0.2262823111496524  0.6720199792759201 -0.6690975041966392  0.222472454458058   0.05914854058002748 -0.6920237805634262 -0.07937587026932126 0.7150552485654124 0.1530518468523984 -0.1714268707140533 -0.2204061612700494  0.9479499376266967");
//                          //cout << "setup: " << endl << tempmat << endl;
                          CVM = tempmat;
                      }
                      }
            }
      }
      
      if(s0==0 && s1==1){
            switch(d){
              case 3: switch(dT){
                      case 3: {
                          Matrix <double> tempmat(pmax, (dT+1)*(d+dT-2)," 0.8547043234472845 -0.2849014411490949 0.3988620176087329 -0.170940864689457        0 0.4082482904638632 -0.816496580927726 0.408248290463863           -0.5477225575051659 0.7302967433402217 0.1825741858350552 -0.3651483716701112            0.4082482904638641 -0.816496580927726 0.408248290463862 0      "
                                                                         "  -0.1142707226510252 0.5684235947255913 -0.7840729584977717 0.2215094008312083      -0.5224552221805724 0.8055139979119766 -0.2787671792634939 0.02183116464111837            0.2236067977499817 -0.6708203932499401 0.670820393249934 -0.2236067977499756           -0.02183116464265128 0.2787671792660461 -0.8055139979117942 0.5224552221794281     "
                                                                         "  0.9152081163589021 -0.2440554976957071 0.3017169614370021 -0.1086181061173199     0 0.4082482904639252 -0.8164965809277187 0.4082482904638156          0.5477225575051627 -0.7302967433402219 -0.1825741858350566 0.3651483716701149             0.408248290463864 -0.816496580927726 0.4082482904638622 0       "
                                                                         "  -0.6606713869121563 -0.22074441926776 0.6796712888505161 -0.2298524722294399     -0.640146461278297 0.6982276160627079 0.2446894817047236 -0.2069245317098367           -0.223606796596586 0.6708203920324677 -0.670820394435368 0.2236067989994864            0.2069245144113484 -0.2446894480709285 -0.698227632722497 0.6401464615547835      "
                                                                         "  0.929541665392022 -0.2360740737503546 0.2705063410924056 -0.08395024378724425     0 0.4082482904621896 -0.8164965809275739 0.4082482904658408           0.5477225393344308 -0.7302967658399765 -0.1825741454060958 0.36514837414118            0.4082482532810987 -0.8164965800968248 0.4082483293084262 0       "
                                                                         "  -0.8234805555087669 -0.06710692648984101 0.525638845099171 -0.2026826081312094      -0.5874404691592978 0.4213433865671868 0.6043779008557463 -0.3348295069840767         0.2236044305743063 -0.6708180884643348 0.6708227288772003 -0.2236090723524933                 0.3348310956040639 -0.6043812814279892 -0.4213396189204375 0.5874387879634988        ");
                          //cout << "setup: " << endl << tempmat << endl;
                          CVM = tempmat;
                      }
                      }
            }
      }
      
      if(s0==1 && s1==0){
            switch(d){
              case 3: switch(dT){
                      case 3: {
                            Matrix <double> tempmat(pmax, (dT+1)*(d+dT-2),"         0 0.4082482904638632 -0.816496580927726 0.408248290463863      -0.5477225575051661 0.7302967433402217 0.1825741858350556 -0.3651483716701108       0.4082482904638641 -0.816496580927726 0.408248290463862 0     0.1709408646894492 -0.398862017608727 -0.2849014411490961 0.8547043234472884 "
                                                                         "        -0.5224552221805724 0.8055139979119766 -0.2787671792634939 0.02183116464111837            0.2236067977499873 -0.6708203932499458 0.6708203932499281 -0.2236067977499699           -0.02183116464265128 0.2787671792660461 -0.8055139979117942 0.5224552221794281     0.2215094008295896 -0.7840729584956645 0.5684235947274756 -0.1142707226592489"
                                                                         "      0 0.4082482904639252 -0.8164965809277187 0.4082482904638156          0.547722557504977 -0.730296743340497 -0.1825741858346003 0.3651483716700706             0.408248290463864 -0.816496580927726 0.4082482904638622 0       0.1086181061610587 -0.3017169614787351 -0.2440554976909896 0.915208116341211"
                                                                         "      -0.640146461278297 0.6982276160627079 0.2446894817047236 -0.2069245317098367           0.2236067977501682 -0.6708203932500731 0.6708203932497893 -0.2236067977498231           0.2069245144113484 -0.2446894480709285 -0.698227632722497 0.6401464615547835      0.229852470662522 -0.6796713154776119 0.2207444536824991 0.6606713485657464"
                                                                         "       0 0.4082482904621896 -0.8164965809275739 0.4082482904658408          -0.5477225575012672 0.7302967433437152 0.1825741858292311 -0.3651483716718835            0.4082482532810987 -0.8164965800968248 0.4082483293084262 0       0.0839502616797303 -0.2705063590741015 -0.2360740721163697 0.9295416589582056"
                                                                         "        -0.5874404691592978 0.4213433865671868 0.6043779008557463 -0.3348295069840767            -0.2236067977542563 0.6708203932589307 -0.6708203932434262 0.2236067977382518            0.3348310956040639 -0.6043812814279892 -0.4213396189204375 0.5874387879634988        0.2026823648370047 -0.525640966362598 0.06710930968033792 0.8234790671411762");
                          
                          
//                          Matrix <double> tempmat(pmax, (dT+1)*(d+dT-2), "    0.2025455652183159 -0.6548973275392216  0.6886549217422743 -0.2363031594213687   0.2236067977499787 -0.6708203932499367 0.6708203932499369 -0.223606797749979 -0.06398930756904883  0.1706381535174648 -0.7465419216389197 0.6398930756905037 0.03848645455033218  0.5388103637046633 -0.7216210228187427  0.4329726136912416 "
//                                                                         "   -0.2631863804043411  0.683532627964419  -0.647689149490306   0.2097971432366053   0.2236067977499833 -0.6708203932499415 0.6708203932499324 -0.2236067977499744  0.05359750604883835 -0.1585898809126361 -0.6145357885456364 0.7709230322243936 0.02868473974068089 -0.6740913839391539  0.551284841922291   0.4907779690162553"
//                                                                         "    0.2167794512514804 -0.6656161436521536  0.6767648011450804 -0.2279281087444072  -0.2236067977499787 0.6708203932499448 -0.6708203932499331 0.2236067977499671  0.07505866259329713 -0.4503519751849037 -0.4128226436599776 0.788115956251584 -0.06123260614795069 -0.4566931886791215  0.2946819175241062  0.8371645398861219"
//                                                                         "   -0.2337005659921087  0.6749399235427884 -0.6645436634817194  0.2195947731375184   0.2236067977502918 -0.6708203932501772 0.6708203932496646 -0.2236067977497625   0.1040966707159132  -0.5950659374191536 -0.2527858014715041 0.755751117659601 0.01110242632083908  -0.29088348746909    0.03759557501427863 0.9559550750882365"
//                                                                         "    0.2215835947572188 -0.6692635936886583  0.6726040031334517 -0.2249240042020122   0.2236067977492572 -0.6708203932492558 0.6708203932506384 -0.2236067977506397  0.07299534619159279 -0.6558770463101348 -0.1527558639038796 0.7356375640224215 0.04944480373238666 -0.2191607550946736 -0.07698606029340599 0.9713891708942911"
//                                                                         "   -0.2262823111496524  0.6720199792759201 -0.6690975041966392  0.222472454458058   -0.2236067977567912 0.6708203932563452 -0.6708203932433455 0.2236067977437161  0.05914854058002748 -0.6920237805634262 -0.07937587026932126 0.7150552485654124 0.1530518468523984 -0.1714268707140533 -0.2204061612700494  0.9479499376266967");
                          //cout << "setup: " << endl << tempmat << endl;
                          CVM = tempmat;
                      }
                      }
            }
      }
      
      
      
      
  }
  
  template <int d, int dT>
  void
  PQFrame<d, dT>::F(SparseMatrix<double>& FF) {
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
  PQFrame<d, dT>::P(const Matrix<double>& ML, const Matrix<double>& MR, SparseMatrix<double>& PP) {
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
  PQFrame<d, dT>::GSetup(SparseMatrix<double>& A, SparseMatrix<double>& H, SparseMatrix<double>& Hinv) {
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
      for (int k = cdf.a().begin(); k <= cdf.a().end(); k++, row++) {
  	A.set_entry(row, col, cdf.a().a(k));
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
  PQFrame<d, dT>::GElim(SparseMatrix<double>& A, SparseMatrix<double>& H, SparseMatrix<double>& Hinv) {
    // IGPMlib reference: I_Basis_Bspline_s::gelim()
    
//     cout << "A=" << endl << A;

    const int firstcol = d-1-s0; // first column of A_j^{(d)} in Ahat_j^{(d)}
    const int lastcol  = Deltasize(j0())-d+s1; // last column
    const int firstrow = d-1-s0; // first row of A_j^{(d)} in Ahat_j^{(d)}
    const int lastrow  = Deltasize(j0()+1)-d+s1; // last row
    
    SparseMatrix<double> help;

//     cout << "A=" << endl << A;

#if 1
    int incr1 = 0;
    int incr2 = 0;
    // elimination (4.1.4)ff.:
    for (int i = 1; i <= d; i++) {
      help.diagonal(Deltasize(j0()+1), 1.0);

      // row index of the entry in the first column of A_j^{(i)} (w.r.t. Ahat_j^{(i)})
      const int elimrow = i%2 ? firstrow+(i-1)/2 : lastrow-(int)floor((i-1)/2.);
      //cout << "i=" << i << ", i%2=" << i%2 << ", elimrow=" << elimrow << endl;

      // factorization [P, p. 112]      
      const int HhatLow = d-1-s0+incr1;
      const int HhatUp  = Deltasize(j0()+1)-d+s1-incr2;

      if (i%2)
	incr1++;
      else
	incr2++;

      if (i%2) // i odd, elimination from above (4.1.4a)
	{
	  assert(fabs(A.get_entry(elimrow+1, firstcol)) >= 1e-10);
	  const double Uentry = -A.get_entry(elimrow, firstcol) / A.get_entry(elimrow+1, firstcol);
	  
	  // setup elimination matrix 'help' from upper left corner to lower right
 	  for (int k = HhatLow; k+1 <= lastrow; k+= 2)
	    help.set_entry(k, k+1, Uentry);
	}
      else // i even, elimination from below (4.1.4b)
	{
	  assert(fabs(A.get_entry(elimrow-1, lastcol)) >= 1e-10);
	  const double Lentry = -A.get_entry(elimrow, lastcol) / A.get_entry(elimrow-1, lastcol);

	  // setup elimination matrix 'help' from lower right corner to upper left
	  for (int k = HhatUp; k-1 >= firstcol; k-= 2)
	    help.set_entry(k, k-1, Lentry);
	}
      //cout << "Hfactor=" << endl << help;
     

      A = help * A;
      H = help * H;
     
      //cout << "A=" << endl << A;
 
      A.compress(1e-14);


      // invert help
      if (i%2) {
	for (int k = HhatLow; k+1 <= lastrow; k += 2)
	  help.set_entry(k, k+1, -help.get_entry(k, k+1));
      }	else {
	for (int k = HhatUp; k-1 >= firstcol; k -= 2)
	  help.set_entry(k, k-1, -help.get_entry(k, k-1));
      }
      
      Hinv = Hinv * help;
    }
#endif

    // HERE COMES THE DEPRECATED CODE. CAN BE REMOVED AS SOON AS ABOVE CODE IS ACCEPTED.
#if 0
    // elimination (4.1.4)ff.:
    for (int i = 1; i <= d; i++) {
      help.diagonal(Deltasize(j0()+1), 1.0);

      // row index of the entry in the first column of A_j^{(i)} (w.r.t. Ahat_j^{(i)})
      const int elimrow = i%2 ? firstrow+(i-1)/2 : lastrow-(int)floor((i-1)/2.);
//       cout << "i=" << i << ", i%2=" << i%2 << ", elimrow=" << elimrow << endl;

      // factorization [P, p. 112]      
      const int HhatLow = i%2 ? (d-s0)-(((i+1)/2)%2) : d-s0+1-(d%2)-((i/2)%2);

      if (i%2) // i odd, elimination from above (4.1.4a)
	{
	  assert(fabs(A.get_entry(elimrow+1, firstcol)) >= 1e-10);
	  const double Uentry = -A.get_entry(elimrow, firstcol) / A.get_entry(elimrow+1, firstcol);
	  
	  // insert Uentry in Hhat
 	  for (int k = HhatLow; k+1 < Deltasize(j0()+1)-(d-1-s1); k+= 2)
	    help.set_entry(k, k+1, Uentry);
	}
      else // i even, elimination from below (4.1.4b)
	{
	  assert(fabs(A.get_entry(elimrow-1, lastcol)) >= 1e-10);
	  const double Lentry = -A.get_entry(elimrow, lastcol) / A.get_entry(elimrow-1, lastcol);
  
	  // insert Lentry in Hhat
 	  for (int k = HhatLow; k+1 < Deltasize(j0()+1)-(d-1-s1); k+= 2)
	    help.set_entry(k+1, k, Lentry);
	}

//       cout << "A=" << endl << A;
//       cout << "Hfactor=" << endl << help;

      A = help * A;
      H = help * H;
      
      A.compress(1e-14);



      // invert help
      if (i%2) {
	for (int k = HhatLow; k+1 < Deltasize(j0()+1)-(d-1-s1); k += 2)
	  help.set_entry(k, k+1, -help.get_entry(k, k+1));
      }	else {
	for (int k = HhatLow; k+1 < Deltasize(j0()+1)-(d-1-s1); k += 2)
	  help.set_entry(k+1, k, -help.get_entry(k+1, k));
      }
      
      Hinv = Hinv * help;
    }

//     cout << "finally, H=" << endl << H;
//     cout << "finally, Hinv=" << endl << Hinv;
#endif
  }



//   template <int d, int dT>
//   void
//   PBasis<d, dT>::GElim(SparseMatrix<double>& A, SparseMatrix<double>& H, SparseMatrix<double>& Hinv) {
//     // IGPMlib reference: I_Basis_Bspline_s::gelim()
    
//     cout << "s0 = " << s0 << endl;
//     cout << "s1 = " << s1 << endl;

//     const int firstcol = d-1-s0; // first column of A_j^{(d)} in Ahat_j^{(d)}
//     const int lastcol  = Deltasize(j0())-d+s1; // last column
//     const int firstrow = d-1-s0; // first row of A_j^{(d)} in Ahat_j^{(d)}
//     const int lastrow  = Deltasize(j0()+1)-d+s1; // last row
    
//     SparseMatrix<double> help;

// //     cout << "A=" << endl << A;
    
//     // elimination (4.1.4)ff.:
//     for (int i = 1; i <= d; i++) {
//       help.diagonal(Deltasize(j0()+1), 1.0);

//       // row index of the entry in the first column of A_j^{(i)} (w.r.t. Ahat_j^{(i)})
//       const int elimrow = i%2 ? firstrow+(i-1)/2 : lastrow-(int)floor((i-1)/2.);
// //       cout << "i%2=" << i%2 << ", elimrow=" << elimrow << endl;

//       // factorization [P, p. 112]      
// //       const int HhatLow = i%2 ? d-s0-((i+1)/2)%2 : d-s0+1-(d%2)-(i/2)%2;
//       const int HhatLow = i%2 ? d-s0-((i+1)/2)%2 : d-s0+1-(d%2)-(i/2)%2;
//       const int HhatUp  = i%2
// 	? Deltasize(j0()+1)-d+s1-abs((d%2)-((i+1)/2)%2)
// 	: Deltasize(j0()+1)-d+s1-abs((d%2)-(i/2)%2);

//       if (i%2) // i odd, elimination from above (4.1.4a)
// 	{
	  
// 	  assert(fabs(A.get_entry(elimrow+1, firstcol)) >= 1e-10);
// 	  const double Uentry = -A.get_entry(elimrow, firstcol) / A.get_entry(elimrow+1, firstcol);

// 	  // insert Uentry in Hhat
// 	  for (int k = HhatLow; k <= HhatUp; k += 2)
// 	    help.set_entry(k, k+1, Uentry);

// 	}
//       else // i even, elimination from below (4.1.4b)
// 	{
// 	  assert(fabs(A.get_entry(elimrow-1, lastcol)) >= 1e-10);
// 	  const double Lentry = -A.get_entry(elimrow, lastcol) / A.get_entry(elimrow-1, lastcol);
	  
// 	  cout << "333333333333333" << endl;
// 	  // insert Lentry in Hhat	    
// 	  cout << " size = " << help.row_dimension()
// 	       << " " << help.column_dimension() << endl;
// 	  cout << "i = " << i << endl;
// 	  cout << "up = " << HhatUp << endl;
// 	  for (int k = HhatLow; k <= HhatUp; k += 2) {
// 	    cout << " k = " << k << endl;
// 	    help.set_entry(k+1, k, Lentry);
// 	  }
// 	  cout << "44444444444444" << endl;
// 	}
//       cout << "ccccccccccccccccc" << endl;
// //       cout << "Hfactor=" << endl << help;

//       A = help * A;
//       H = help * H;
      
//       A.compress(1e-14);

//       //      cout << "A=" << endl << A;

//       // invert help
//       if (i%2) {
// 	for (int k = HhatLow; k <= HhatUp; k += 2)
// 	  help.set_entry(k, k+1, -help.get_entry(k, k+1));
//       }	else {
// 	for (int k = HhatLow; k <= HhatUp; k += 2)
// 	  help.set_entry(k+1, k, -help.get_entry(k+1, k));
//       }
//       Hinv = Hinv * help;
//     }

//     //     cout << "finally, H=" << endl << H;
// //     cout << "finally, Hinv=" << endl << Hinv;
//   }


  template <int d, int dT>
  double
  PQFrame<d, dT>::BT(const SparseMatrix<double>& A, SparseMatrix<double>& BB) {
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
  PQFrame<d, dT>::InvertP(const SparseMatrix<double>& PP, SparseMatrix<double>& PPinv) {
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
  PQFrame<d, dT>::DS_symmetrization(SparseMatrix<double>& Mj1, SparseMatrix<double>& Mj1T) {
    // IGPMlib reference: I_Basis_Bspline::Modify()
    
    SparseMatrix<double> Hj1(Deltasize(j0()+1), 1<<j0()),
      Hj1T(Deltasize(j0()+1), 1<<j0());
    
    // copy left halves of Mj1, Mj1T, right halves are reflected
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
  PQFrame<d, dT>::assemble_Mj0(const int j, SparseMatrix<double>& mj0) const {
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
  PQFrame<d, dT>::assemble_Mj0_t(const int j, SparseMatrix<double>& mj0_t) const {
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
  PQFrame<d, dT>::assemble_Mj0T(const int j, SparseMatrix<double>& mj0T) const {
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
  PQFrame<d, dT>::assemble_Mj0T_t(const int j, SparseMatrix<double>& mj0T_t) const {
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
  PQFrame<d, dT>::assemble_Mj1(const int j, SparseMatrix<double>& mj1) const {
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
  PQFrame<d, dT>::assemble_Mj1_t(const int j, SparseMatrix<double>& mj1_t) const {
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
  PQFrame<d, dT>::assemble_Mj1T(const int j, SparseMatrix<double>& mj1T) const {
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
  PQFrame<d, dT>::assemble_Mj1T_t(const int j, SparseMatrix<double>& mj1T_t) const {
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
  PQFrame<d, dT>::Mj0_get_row(const int j, const Vector<double>::size_type row,
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
//	typedef Vector<double>::size_type size_type;
	if (row >= bottom)
	  Mj0.get_row(row+rows_top-bottom, v, Deltasize(j)-Deltasize(j0()));
	else
	  Mj0.get_row(rows_top-2+(row-rows_top)%2, v, (row-rows_top)/2+1);
      }
    }
  }
  
  template <int d, int dT>
  void
  PQFrame<d, dT>::Mj0T_get_row(const int j, const Vector<double>::size_type row,
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
  PQFrame<d, dT>::Mj1_get_row(const int j, const Vector<double>::size_type row,
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
  PQFrame<d, dT>::Mj1T_get_row(const int j, const Vector<double>::size_type row,
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
  PQFrame<d, dT>::Mj0_t_get_row(const int j, const Vector<double>::size_type row,
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
  PQFrame<d, dT>::Mj0T_t_get_row(const int j, const Vector<double>::size_type row,
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
  PQFrame<d, dT>::Mj1_t_get_row(const int j, const Vector<double>::size_type row,
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
  PQFrame<d, dT>::Mj1T_t_get_row(const int j, const Vector<double>::size_type row,
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
  PQFrame<d, dT>::decompose(const InfiniteVector<double, Index>& c,
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
  PQFrame<d, dT>::decompose_t(const InfiniteVector<double, Index>& c,
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
  PQFrame<d, dT>::reconstruct(const InfiniteVector<double, Index>& c,
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
  PQFrame<d, dT>::reconstruct_t(const InfiniteVector<double, Index>& c,
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
  PQFrame<d, dT>::decompose_1(const Index& lambda,
			      const int jmin,
			      InfiniteVector<double, Index>& c) const {
    assert(jmin >= j0());
    assert(lambda.j() >= jmin);
    assert(lambda.p() == 0);
    
    c.clear();

    if (lambda.e() == 1) // wavelet
      c.set_coefficient(lambda, 1.0); // true wavelet coefficients don't have to be modified
    else // generator
      {
	if (lambda.j() == jmin || lambda.p() > 0) // generators on the coarsest level don't have to be modified
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
		c.set_coefficient(Index(0, j0(), 1, Mj1T.get_nth_index(row_j0,k), this),
				  Mj1T.get_nth_entry(row_j0,k));
	    } else {
	      // Due to the [DS] symmetrization, we have to be a bit careful here.
	      int fill_start = 0, fill_end = (1<<(lambda.j()-1))-1; // first and last column to fill up with interior filter
	      const size_type upper_half = Deltasize(j0()+1)/2;
	      if (row < upper_half) {
		// read the first half of the corresponding row in Mj1T
		for (size_type k(0); k < Mj1T.entries_in_row(row_j0) && (int)Mj1T.get_nth_index(row_j0,k) < 1<<(j0()-1); k++)
		  c.set_coefficient(Index(0, lambda.j()-1, 1, Mj1T.get_nth_index(row_j0,k), this),
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
		      c.set_coefficient(Index(0, lambda.j()-1, 1, Mj1T.get_nth_index(row_j0,k)+offset, this),
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
		  c.set_coefficient(Index(0, lambda.j()-1, 1, col, this),
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
		  c.set_coefficient(Index(0, lambda.j()-1, 1, col, this),
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
	      decompose_1(Index(lambda.p(), lambda.j()-1, 0, DeltaLmin()+Mj0T.get_nth_index(row_j0,k)+offset, this), jmin, dhelp);
	      c.add(Mj0T.get_nth_entry(row_j0,k), dhelp);
	    }
	  }
      }
  }
  
  template <int d, int dT>
  void
  PQFrame<d, dT>::decompose_t_1(const Index& lambda,
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
	      decompose_t_1(Index(lambda.p(), lambda.j()-1, 0, DeltaLmin()+Mj0.get_nth_index(row_j0,k)+offset, this), jmin, dhelp);
	      c.add(Mj0.get_nth_entry(row_j0,k), dhelp);
	    }
	  }
      }
  }
  
  template <int d, int dT>
  void
  PQFrame<d, dT>::reconstruct_1(const Index& lambda,
			       const int j,
			       InfiniteVector<double, Index>& c) const {
      assert(lambda.p()==0 || lambda.e()==1);
    c.clear();
//    int temp_int1 (DeltaLmin()), temp_int2(Deltasize(lambda.j())), temp_int3(Nablamin());
//    int tempint4(get_s0()), temp_int5(get_s1());
//    int temp_int6(s0), temp_int7(s1);
    
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

      if (lambda.j()+1 >= j) {
	
            if(lambda.e() == 1 && lambda.k() < dT-1  && lambda.p() > 0){//left boundary quarklets
                for (size_type k(0); k < dT+1; k++) { 
                    c.add_coefficient(Index(lambda.p(), lambda.j()+1, 0, lambda.k()+k+DeltaLmin(), this) , CVM(lambda.p()-1,(dT+1)*lambda.k()+k));
                }  
            }         
            else if(lambda.e() == 1 && lambda.k() > Nablamax(lambda.j())-dT+1  && lambda.p() > 0){
                for (size_type k(0); k < dT+1; k++) {            
                    c.add_coefficient(Index(lambda.p(), lambda.j()+1, 0, DeltaRmax(lambda.j()+1,lambda.p())-dT+k-Nablamax(lambda.j())+lambda.k(), this) , CVM(lambda.p()-1, (dT+1)*(lambda.k()-(Nablamax(lambda.j())-dT+1)-1+(d+dT-2)/2)+ k));
                }
                
            }
            else{
                
                for (size_type k(0); k < M.entries_in_row(row_j0); k++) {
                c.add_coefficient(Index(lambda.p(), lambda.j()+1, 0, DeltaLmin()+M.get_nth_index(row_j0,k)+offset, this),
			    M.get_nth_entry(row_j0,k));
                }
            }
      }
      
      
      
      
      
      
      
      else {
	for (size_type k(0); k < M.entries_in_row(row_j0); k++) {
	  InfiniteVector<double, Index> dhelp;
	  reconstruct_1(Index(lambda.p(), lambda.j()+1, 0, DeltaLmin()+M.get_nth_index(row_j0,k)+offset, this), j, dhelp);
	  c.add(M.get_nth_entry(row_j0,k), dhelp);
	}
      }
    }
  }

#if 0
  template <int d, int dT>
  void
  PQFrame<d, dT>::reconstruct_1(const int lamj, const int lame, const int lamk,
			       const int j,
			       InfiniteVector<double, Index>& c) const {
    c.clear();
    
    if (lamj >= j)
      c.set_coefficient(Index(lamj,lame,lamk, this), 1.0);
                        //lambda, 1.0); 
    else {
      // For the reconstruction of psi_lambda, we have to compute
      // the corresponding column of the transformation matrix Mj=(Mj0, Mj1).
      
      // reconstruct by recursion
      typedef Vector<double>::size_type size_type;
      
      const size_type row = (lame == 0 ? lamk - DeltaLmin() : lamk);
      
      const SparseMatrix<double>& M = (lame == 0 ? Mj0_t : Mj1_t);
      
      size_type row_j0 = row;
      size_type offset = 0;
      if (lame == 0) {
	if (lamj > j0()) {
	  const size_type rows_top = (int)ceil(Deltasize(j0())/2.0);
	  if (row >= rows_top) {
	    const size_type bottom = Deltasize(lamj)-Deltasize(j0())/2;
	    if (row >= bottom) {
	      row_j0 = row+rows_top-bottom;
	      offset = Deltasize(lamj+1)-Deltasize(j0()+1);
	    } else {
	      row_j0 = rows_top-1;
	      offset = 2*(row-rows_top)+2;
	    }
	  }
	}
      }
      else {
	if (lamj > j0()) {
	  const size_type rows_top = 1<<(j0()-1);
	  if (row >= rows_top) {
	    const size_type bottom = (1<<lamj)-(1<<(j0()-1));
	    if (row >= bottom) {
	      row_j0 = row+rows_top-bottom;
	      offset = Deltasize(lamj+1)-Deltasize(j0()+1);
	    } else {
	      if ((int)row < (1<<(lamj-1))) {
		row_j0 = rows_top-1;
		offset = 2*(row-rows_top)+2;
	      } else {
		row_j0 = 1<<(j0()-1);
		offset = Deltasize(lamj+1)-Deltasize(j0()+1)+2*((int)row-bottom);
	      }
	    }
	  }
	}
      }
      if (lamj+1 >= j) {
	for (size_type k(0); k < M.entries_in_row(row_j0); k++) {
	  c.add_coefficient(Index(lamj+1, 0, DeltaLmin()+M.get_nth_index(row_j0,k)+offset, this),
			    M.get_nth_entry(row_j0,k));
	}
      } else {
	for (size_type k(0); k < M.entries_in_row(row_j0); k++) {
	  InfiniteVector<double, Index> dhelp;
	  reconstruct_1(Index(lamj+1, 0, DeltaLmin()+M.get_nth_index(row_j0,k)+offset, this), j, dhelp);
	  c.add(M.get_nth_entry(row_j0,k), dhelp);
	}
      }
    }
  }
#endif
  
  template <int d, int dT>
  void
  PQFrame<d, dT>::reconstruct_1(const int lamp, const int lamj, const int lame, const int lamk,
			       const int j,
			       InfiniteVector<double, int>& c) const {
    c.clear();
    assert ((lamj <= j) || (lame == 1)); 
    // ! ( (lamj > j) && (lame == 0) ), i.e., the reconstruction of a generator on a level higher than j will fail (because of a wrong number)
    
    if (lamj >= j)
    {
        if (lame == 0)
        {
            c.set_coefficient(lamk - DeltaLmin(), 1.0);
        }
        else
        {
            c.set_coefficient(Deltasize(lamj) + lamk - Nablamin(), 1.0);
        }
        //c.set_coefficient(lambda, 1.0); 
    }
    else {
      // For the reconstruction of psi_lambda, we have to compute
      // the corresponding column of the transformation matrix Mj=(Mj0, Mj1).
      
      // reconstruct by recursion
      typedef Vector<double>::size_type size_type;
      
      const size_type row = (lame == 0 ? lamk - DeltaLmin() : lamk);

      const SparseMatrix<double>& M = (lame == 0 ? Mj0_t : Mj1_t);
//      bool temp_bool = Mj1_t.empty();
//      int temp_rows(M.row_dimension());
//      int temp_cols(M.column_dimension());
      
      size_type row_j0 = row;
      size_type offset = 0;
      if (lame == 0) {
	if (lamj > j0()) {
	  const size_type rows_top = (int)ceil(Deltasize(j0())/2.0);
	  if (row >= rows_top) {
	    const size_type bottom = Deltasize(lamj)-Deltasize(j0())/2;
	    if (row >= bottom) {
	      row_j0 = row+rows_top-bottom;
	      offset = Deltasize(lamj+1)-Deltasize(j0()+1);
	    } else {
	      row_j0 = rows_top-1;
	      offset = 2*(row-rows_top)+2;
	    }
	  }
	}
      }
      else {
	if (lamj > j0()) {
	  const size_type rows_top = 1<<(j0()-1);
	  if (row >= rows_top) {
	    const size_type bottom = (1<<lamj)-(1<<(j0()-1));
	    if (row >= bottom) {
	      row_j0 = row+rows_top-bottom;
	      offset = Deltasize(lamj+1)-Deltasize(j0()+1);
	    } else {
	      if ((int)row < (1<<(lamj-1))) {
		row_j0 = rows_top-1;
		offset = 2*(row-rows_top)+2;
	      } else {
		row_j0 = 1<<(j0()-1);
		offset = Deltasize(lamj+1)-Deltasize(j0()+1)+2*((int)row-bottom);
	      }
	    }
	  }
	}
      }
      if (lamj+1 >= j) {
          
        if(lame == 1 && lamk < dT-1  && lamp > 0){//left boundary quarklets
            assert (lamp<=6);
            for (size_type k(0); k < dT+1; k++) { 
                c.add_coefficient(lamk+k , CVM(lamp-1,(dT+1)*lamk+k));
            }  
        }  
        else if(lame == 1 && lamk > Nablamax(lamj)-dT+1  && lamp > 0){//right boundary quarklets
            assert (lamp<=6);
            
            for (size_type k(0); k < dT+1; k++) {
                c.add_coefficient(DeltaRmax(lamj+1,lamp)-DeltaLmin(lamp)-dT+k-Nablamax(lamj)+lamk,
                CVM(lamp-1, (dT+1)*(lamk-(Nablamax(lamj)-dT+1)-1+(d+dT-2)/2)+ k));
            }
        }  
        else{
            for (size_type k(0); k < M.entries_in_row(row_j0); k++) {
                c.add_coefficient(M.get_nth_index(row_j0,k)+offset , M.get_nth_entry(row_j0,k));
            } 
        }
	
        
      } else {
	for (size_type k(0); k < M.entries_in_row(row_j0); k++) {
	  InfiniteVector<double, int> dhelp;
          reconstruct_1(lamp, lamj+1, 0, DeltaLmin()+M.get_nth_index(row_j0,k)+offset, j, dhelp);
	  //reconstruct_1(Index(lamj+1, 0, DeltaLmin()+M.get_nth_index(row_j0,k)+offset, this), j, dhelp);
	  c.add(M.get_nth_entry(row_j0,k), dhelp);
	}
      }
    }
  }
  
  
  template <int d, int dT>
  void
  PQFrame<d, dT>::reconstruct_t_1(const Index& lambda,
				  const int j,
				  InfiniteVector<double, Index>& c) const {
      assert(lambda.p()== 0 || lambda.e() == 1);
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

      if (lambda.j()+1 >= j) {
	for (size_type k(0); k < M.entries_in_row(row_j0); k++) {
	  c.add_coefficient(Index(lambda.p(), lambda.j()+1, 0, DeltaLmin()+M.get_nth_index(row_j0,k)+offset, this),
			    M.get_nth_entry(row_j0,k));
	}
      } else {
	for (size_type k(0); k < M.entries_in_row(row_j0); k++) {
	  InfiniteVector<double, Index> dhelp;
	  reconstruct_t_1(Index(lambda.p(), lambda.j()+1, 0, DeltaLmin()+M.get_nth_index(row_j0,k)+offset, this), j, dhelp);
	  c.add(M.get_nth_entry(row_j0,k), dhelp);
	}
      }
    }
  }

  template <int d, int dT>
  inline
  double
  PQFrame<d, dT>::evaluate(const unsigned int derivative, const Index& lambda, const double x) const
  {
    return WaveletTL::evaluate(*this, derivative, lambda, x);
  }
  
  template <int d, int dT>
  inline
  double
  PQFrame<d, dT>::evaluate(const unsigned int derivative, const int p, const int j, const int e, const int k, const double x) const
  {
    return WaveletTL::evaluate(*this, derivative, p, j, e, k, x);
  }

  template <int d, int dT>
  inline
  void
  PQFrame<d,dT>::evaluate
  (const unsigned int derivative,
   const Index& lambda,
   const Array1D<double>& points, Array1D<double>& values) const
  {
    WaveletTL::evaluate(*this, derivative, lambda, points, values);
  }
  
  template <int d, int dT>
  inline
  void
  PQFrame<d,dT>::evaluate
  (const unsigned int derivative,
   const int p_, const int j_, const int e_, const int k_,
   const Array1D<double>& points, Array1D<double>& values) const
  {
    WaveletTL::evaluate(*this, derivative, p_, j_, e_, k_, points, values);
  }
  
  
  //new: evaluate of type SampledMapping<1>
  template <int d, int dT>
  inline
  SampledMapping<1>
  PQFrame<d,dT>::evaluate(const Index& lambda, const bool primal, const int resolution) const
  {
    return WaveletTL::evaluate(*this, lambda, primal, resolution);
  }

   /* Compute the Picewiese expansion of all wavelets and Generatoren for given j, d */
  template <int d, int dT>
  inline
  void
  PQFrame<d,dT>::waveletPP(const int j, Array1D<Piecewise<double> >& wavelets) const
  {
   int i = std::max((1<<j),(1<<d));
   wavelets.resize(i); //-ell1<d>()
   for(int k=0; k<=(i-1); k++){
     Index lambda1(j,1,k,this);
     wavelets[k] = expandAsPP(*this,lambda1);
   }

  }

}
