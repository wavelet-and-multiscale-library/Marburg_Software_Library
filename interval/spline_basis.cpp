// implementation for spline_basis.h

#include <Rd/cdf_utils.h>

namespace WaveletTL
{
  template <int d, int dT>
  SplineBasis<d,dT>::SplineBasis(const char* flavor,
				 const char* options,
				 const int s0, const int s1, const int sT0, const int sT1)
    : SplineBasisData<d,dT>(flavor,options,s0,s1,sT0,sT1)
  {
    if (SplineBasisData<d,dT>::flavor_ == "P") {
      // cf. PBasis<d,dT> ...
      DeltaLmin_       = 1-d-ell1<d>()+s0;
      DeltaRmax_offset = -1-ell1<d>()-s1;
    }
  }

  template <int d, int dT>
  inline
  typename SplineBasis<d,dT>::Index
  SplineBasis<d,dT>::first_generator(const int j) const
  {
    assert(j >= j0());
    return Index(j, 0, DeltaLmin(), this);
  }
  
  template <int d, int dT>
  inline
  typename SplineBasis<d,dT>::Index
  SplineBasis<d,dT>::last_generator(const int j) const
  {
    assert(j >= j0());
    return Index(j, 0, DeltaRmax(j), this);
  }

  template <int d, int dT>
  inline
  typename SplineBasis<d,dT>::Index
  SplineBasis<d,dT>::first_wavelet(const int j) const
  {
    assert(j >= j0());
    return Index(j, 1, Nablamin(), this);
  }
  
  template <int d, int dT>
  inline
  typename SplineBasis<d,dT>::Index
  SplineBasis<d,dT>::last_wavelet(const int j) const
  {
    assert(j >= j0());
    return Index(j, 1, Nablamax(j), this);
  }

  template <int d, int dT>
  void
  SplineBasis<d,dT>::support(const Index& lambda, int& k1, int& k2) const {
    if (lambda.e() == 0) // generator
      {
	// phi_{j,k}(x) = 2^{j/2} B_{j,k-ell_1}
	k1 = std::max(0, lambda.k() + ell1<d>());
	k2 = std::min(1<<lambda.j(), lambda.k() + ell2<d>());
      }
    else // wavelet
      {
	// cf. [P, p. 125]
	if (lambda.k() < (d+dT)/2-1) {
	  // left boundary wavelet
	  k1 = 0;
	  k2 = 2*(d+dT)-2; // overestimate, TODO
	} else {
	  if ((1<<lambda.j())-lambda.k() <= (d+dT)/2-1) {
	    // right boundary wavelet
	    k1 = (1<<(lambda.j()+1))-(2*(d+dT)-2); // overestimate, TODO
	    k2 = 1<<(lambda.j()+1);
	  } else {
	    // interior wavelet (CDF)
	    k1 = 2*(lambda.k()-(d+dT)/2+1);
	    k2 = k1+2*(d+dT)-2;
	  }
	}
      }
  }
  
  template <int d, int dT>
  template <class V>
  void
  SplineBasis<d,dT>::apply_Mj(const int j, const V& x, V& y) const
  {
    SplineBasisData<d,dT>::Mj0_->set_level(j);
    SplineBasisData<d,dT>::Mj1_->set_level(j);

    // decompose x=(x1 x2) appropriately
    SplineBasisData<d,dT>::Mj0_->apply(x, y, 0, 0); // apply Mj0 to first block x1
    SplineBasisData<d,dT>::Mj1_->apply(x, y,        // apply Mj1 to second block x2 and add result
				       SplineBasisData<d,dT>::Mj0_->column_dimension(), 0, true);
  }

  template <int d, int dT>
  template <class V>
  void
  SplineBasis<d,dT>::apply_Mj_transposed(const int j, const V& x, V& y) const
  {
    SplineBasisData<d,dT>::Mj0_->set_level(j);
    SplineBasisData<d,dT>::Mj1_->set_level(j);

    // y=(y1 y2) is a block vector
    SplineBasisData<d,dT>::Mj0_->apply_transposed(x, y, 0, 0); // write into first block y1
    SplineBasisData<d,dT>::Mj1_->apply_transposed(x, y, 0,     // write into second block y2
						  SplineBasisData<d,dT>::Mj0_->column_dimension());
  }

  template <int d, int dT>
  template <class V>
  void
  SplineBasis<d,dT>::apply_Gj(const int j, const V& x, V& y) const
  {
    SplineBasisData<d,dT>::Mj0T_->set_level(j);
    SplineBasisData<d,dT>::Mj1T_->set_level(j);

    // y=(y1 y2) is a block vector
    SplineBasisData<d,dT>::Mj0T_->apply_transposed(x, y, 0, 0); // write into first block y1
    SplineBasisData<d,dT>::Mj1T_->apply_transposed(x, y, 0,     // write into second block y2
						   SplineBasisData<d,dT>::Mj0T_->column_dimension());
  }
  
  template <int d, int dT>
  template <class V>
  void
  SplineBasis<d,dT>::apply_Gj_transposed(const int j, const V& x, V& y) const
  {
    SplineBasisData<d,dT>::Mj0T_->set_level(j);
    SplineBasisData<d,dT>::Mj1T_->set_level(j);
    
    // decompose x=(x1 x2) appropriately
    SplineBasisData<d,dT>::Mj0T_->apply(x, y, 0); // apply Mj0T to first block x1
    SplineBasisData<d,dT>::Mj1T_->apply(x, y,     // apply Mj1T to second block x2 and add result
					SplineBasisData<d,dT>::Mj0T_->column_dimension(), 0, true);
  }
  
  template <int d, int dT>
  template <class V>
  void
  SplineBasis<d,dT>::apply_Tj(const int j, const V& x, V& y) const
  { 
    y = x;
    V z(x);
    apply_Mj(j0(), z, y);
    for (int k = j0()+1; k <= j; k++) {
      apply_Mj(k, y, z);
      y.swap(z);
    }
  }

  template <int d, int dT>
  template <class V>
  void
  SplineBasis<d,dT>::apply_Tj_transposed(const int j, const V& x, V& y) const
  { 
    V z;
    apply_Mj_transposed(j, x, y);
    for (int k = j-1; k >= j0(); k--) {
      z.swap(y);
      apply_Mj_transposed(k, z, y);
      for (int i = Deltasize(k+1); i < Deltasize(j+1); i++)
	y[i] = z[i];
    }
  }
  
  template <int d, int dT>
  void
  SplineBasis<d,dT>::apply_Tj_transposed(const int j,
					 const Vector<double>& x,
					 Vector<double>& y) const
  { 
    Vector<double> z(x.size(), false);
    apply_Mj_transposed(j, x, y);
    for (int k = j-1; k >= j0(); k--) {
      z.swap(y);
      apply_Mj_transposed(k, z, y);
      for (int i = Deltasize(k+1); i < Deltasize(j+1); i++)
	y[i] = z[i];
    }
  }

  template <int d, int dT>
  void
  SplineBasis<d,dT>::apply_Tjinv(const int j, const Vector<double>& x, Vector<double>& y) const
  { 
    // T_j^{-1}=diag(G_{j0},I)*...*diag(G_{j-1},I)*G_j
    Vector<double> z(x.size(), false);
    apply_Gj(j, x, y);
    for (int k = j-1; k >= j0(); k--) {
      z.swap(y);
      apply_Gj(k, z, y);
      for (int i = Deltasize(k+1); i < Deltasize(j+1); i++)
	y[i] = z[i];
    }
  }

}
