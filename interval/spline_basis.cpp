// implementation for spline_basis.h

#include <algebra/sparse_matrix.h>
#include <algebra/vector.h>
#include <numerics/cardinal_splines.h>
#include <numerics/schoenberg_splines.h>
#include <numerics/gauss_data.h>
#include <numerics/quadrature.h>
#include <numerics/iteratsolv.h>

#include <Rd/r_index.h>
#include <Rd/cdf_utils.h>
#include <Rd/cdf_basis.h>
#include <interval/interval_bspline.h>
#include <galerkin/full_gramian.h>

namespace WaveletTL
{
  template <int d, int dT, SplineBasisFlavor flavor, int s0, int s1, int sT0, int sT1>
  const int SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::DeltaLmin() {
    return ell2T<d,dT>()+s0+sT0-dT; // cf. DSBasis<d,dT>
  }

  template <int d, int dT, SplineBasisFlavor flavor, int s0, int s1, int sT0, int sT1>
  const int SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::DeltaRmax_offset() {
    return -(d%2)-(ell2T<d,dT>()+s1+sT1-dT); // cf. DSBasis<d,dT>
  }

  template <int d, int dT, int s0, int s1, int sT0, int sT1>
  const int SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1>::DeltaLmin() {
    return 1-d-ell1<d>()+s0; // cf. PBasis<d,dT>
  }
  
  template <int d, int dT, int s0, int s1, int sT0, int sT1>
  const int SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1>::DeltaRmax_offset() {
    return -1-ell1<d>()-s1; // cf. PBasis<d,dT>
  }
  
  template <int d, int dT, SplineBasisFlavor flavor, int s0, int s1, int sT0, int sT1>
  inline
  typename SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::Index
  SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::first_generator(const int j) const
  {
    assert(j >= j0());
    return Index(j, 0, DeltaLmin());
  }
  
  template <int d, int dT, int s0, int s1, int sT0, int sT1>
  inline
  typename SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1>::Index
  SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1>::first_generator(const int j) const
  {
    assert(j >= j0());
    return Index(j, 0, DeltaLmin());
  }
  
  template <int d, int dT, SplineBasisFlavor flavor, int s0, int s1, int sT0, int sT1>
  inline
  typename SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::Index
  SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::last_generator(const int j) const
  {
    assert(j >= j0());
    return Index(j, 0, DeltaRmax(j));
  }

  template <int d, int dT, int s0, int s1, int sT0, int sT1>
  inline
  typename SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1>::Index
  SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1>::last_generator(const int j) const
  {
    assert(j >= j0());
    return Index(j, 0, DeltaRmax(j));
  }

  template <int d, int dT, SplineBasisFlavor flavor, int s0, int s1, int sT0, int sT1>
  inline
  typename SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::Index
  SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::first_wavelet(const int j) const
  {
    assert(j >= j0());
    return Index(j, 1, Nablamin());
  }
  
  template <int d, int dT, int s0, int s1, int sT0, int sT1>
  inline
  typename SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1>::Index
  SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1>::first_wavelet(const int j) const
  {
    assert(j >= j0());
    return Index(j, 1, Nablamin());
  }
  
  template <int d, int dT, SplineBasisFlavor flavor, int s0, int s1, int sT0, int sT1>
  inline
  typename SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::Index
  SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::last_wavelet(const int j) const
  {
    assert(j >= j0());
    return Index(j, 1, Nablamax(j));
  }

  template <int d, int dT, int s0, int s1, int sT0, int sT1>
  inline
  typename SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1>::Index
  SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1>::last_wavelet(const int j) const
  {
    assert(j >= j0());
    return Index(j, 1, Nablamax(j));
  }

  // generic support routine for all DS-type wavelet bases
  template <int d, int dT, SplineBasisFlavor flavor, int s0, int s1, int sT0, int sT1>
  void
  SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::support(const Index& lambda, int& k1, int& k2) const {
    if (lambda.e() == 0) // generator
      {
 	if (lambda.k() < DeltaLmin()+(int)SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::CLA_.column_dimension())
 	  {
 	    // left boundary generator
 	    k1 = 0;
 	    k2 = SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::CLA_.row_dimension();
 	  }
 	else
 	  {
 	    if (lambda.k() > DeltaRmax(lambda.j())-(int)SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::CRA_.column_dimension())
 	      {
 		// right boundary generator
 		k1 = (1<<lambda.j()) - SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::CRA_.row_dimension();
 		k2 = 1<<lambda.j();
 	      }
 	    else
 	      {
 		k1 = lambda.k() + ell1<d>();
 		k2 = lambda.k() + ell2<d>();
 	      }
 	  }
      }
    else // wavelet
      {
  	// To determine which generators would be necessary to create the
  	// wavelet in question, we mimic a reconstruct_1() call:
 	
 	typedef typename Vector<double>::size_type size_type;
 	std::map<size_type,double> wc, gc;
	wc[lambda.k()-Nablamin()] = 1.0;
	SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::Mj1_.set_level(lambda.j());
	SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::Mj1_.apply(wc, gc, 0, 0);
 	int dummy;
  	support(Index(lambda.j()+1, 0, DeltaLmin()+gc.begin()->first), k1, dummy);
  	support(Index(lambda.j()+1, 0, DeltaLmin()+gc.rbegin()->first), dummy, k2);
      }
  }  

  template <int d, int dT, int s0, int s1, int sT0, int sT1>
  void
  SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1>::support(const Index& lambda, int& k1, int& k2) const {
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
	    // note: despite the fact that the wavelets in the right half of the interval
	    //       are reflected CDF wavelets, their support does not "see" the reflection!
	    k1 = 2*(lambda.k()-(d+dT)/2+1);
	    k2 = k1+2*(d+dT)-2;
	  }
	}
      }
  }
  
  template <int d, int dT, SplineBasisFlavor flavor, int s0, int s1, int sT0, int sT1>
  template <class V>
  void
  SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::apply_Mj(const int j, const V& x, V& y) const
  {
    SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::Mj0_.set_level(j);
    SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::Mj1_.set_level(j);

    // decompose x=(x1 x2) appropriately
    SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::Mj0_.apply
      (x, y, 0, 0); // apply Mj0 to first block x1
    SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::Mj1_.apply
      (x, y,        // apply Mj1 to second block x2 and add result
       SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::Mj0_.column_dimension(), 0, true);
  }

  template <int d, int dT, int s0, int s1, int sT0, int sT1>
  template <class V>
  void
  SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1>::apply_Mj(const int j, const V& x, V& y) const
  {
    SplineBasisData<d,dT,P_construction,s0,s1,sT0,sT1>::Mj0_.set_level(j);
    SplineBasisData<d,dT,P_construction,s0,s1,sT0,sT1>::Mj1_.set_level(j);
    
    // decompose x=(x1 x2) appropriately
    SplineBasisData<d,dT,P_construction,s0,s1,sT0,sT1>::Mj0_.apply
      (x, y, 0, 0); // apply Mj0 to first block x1
    SplineBasisData<d,dT,P_construction,s0,s1,sT0,sT1>::Mj1_.apply
      (x, y,        // apply Mj1 to second block x2 and add result
       SplineBasisData<d,dT,P_construction,s0,s1,sT0,sT1>::Mj0_.column_dimension(), 0, true);
  }
  
  template <int d, int dT, SplineBasisFlavor flavor, int s0, int s1, int sT0, int sT1>
  template <class V>
  void
  SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::apply_Mj_transposed(const int j, const V& x, V& y) const
  {
    SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::Mj0_.set_level(j);
    SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::Mj1_.set_level(j);

    // y=(y1 y2) is a block vector
    SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::Mj0_.apply_transposed
      (x, y, 0, 0); // write into first block y1
    SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::Mj1_.apply_transposed
      (x, y, 0,     // write into second block y2
       SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::Mj0_.column_dimension());
  }

  template <int d, int dT, int s0, int s1, int sT0, int sT1>
  template <class V>
  void
  SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1>::apply_Mj_transposed(const int j, const V& x, V& y) const
  {
    SplineBasisData<d,dT,P_construction,s0,s1,sT0,sT1>::Mj0_.set_level(j);
    SplineBasisData<d,dT,P_construction,s0,s1,sT0,sT1>::Mj1_.set_level(j);
    
    // y=(y1 y2) is a block vector
    SplineBasisData<d,dT,P_construction,s0,s1,sT0,sT1>::Mj0_.apply_transposed
      (x, y, 0, 0); // write into first block y1
    SplineBasisData<d,dT,P_construction,s0,s1,sT0,sT1>::Mj1_.apply_transposed
      (x, y, 0,     // write into second block y2
       SplineBasisData<d,dT,P_construction,s0,s1,sT0,sT1>::Mj0_.column_dimension());
  }

  template <int d, int dT, SplineBasisFlavor flavor, int s0, int s1, int sT0, int sT1>
  template <class V>
  void
  SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::apply_Gj(const int j, const V& x, V& y) const
  {
    SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::Mj0T_.set_level(j);
    SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::Mj1T_.set_level(j);

    // y=(y1 y2) is a block vector
    SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::Mj0T_.apply_transposed
      (x, y, 0, 0); // write into first block y1
    SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::Mj1T_.apply_transposed
      (x, y, 0,     // write into second block y2
       SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::Mj0T_.column_dimension());
  }
  
  template <int d, int dT, int s0, int s1, int sT0, int sT1>
  template <class V>
  void
  SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1>::apply_Gj(const int j, const V& x, V& y) const
  {
    SplineBasisData<d,dT,P_construction,s0,s1,sT0,sT1>::Mj0T_.set_level(j);
    SplineBasisData<d,dT,P_construction,s0,s1,sT0,sT1>::Mj1T_.set_level(j);
    
    // y=(y1 y2) is a block vector
    SplineBasisData<d,dT,P_construction,s0,s1,sT0,sT1>::Mj0T_.apply_transposed
      (x, y, 0, 0); // write into first block y1
    SplineBasisData<d,dT,P_construction,s0,s1,sT0,sT1>::Mj1T_.apply_transposed
      (x, y, 0,     // write into second block y2
       SplineBasisData<d,dT,P_construction,s0,s1,sT0,sT1>::Mj0T_.column_dimension());
  }
  
  template <int d, int dT, SplineBasisFlavor flavor, int s0, int s1, int sT0, int sT1>
  template <class V>
  void
  SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::apply_Gj_transposed(const int j, const V& x, V& y) const
  {
    SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::Mj0T_.set_level(j);
    SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::Mj1T_.set_level(j);
    
    // decompose x=(x1 x2) appropriately
    SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::Mj0T_.apply
      (x, y, 0); // apply Mj0T to first block x1
    SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::Mj1T_.apply
      (x, y,     // apply Mj1T to second block x2 and add result
       SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::Mj0T_.column_dimension(), 0, true);
  }
  
  template <int d, int dT, int s0, int s1, int sT0, int sT1>
  template <class V>
  void
  SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1>::apply_Gj_transposed(const int j, const V& x, V& y) const
  {
    SplineBasisData<d,dT,P_construction,s0,s1,sT0,sT1>::Mj0T_.set_level(j);
    SplineBasisData<d,dT,P_construction,s0,s1,sT0,sT1>::Mj1T_.set_level(j);
    
    // decompose x=(x1 x2) appropriately
    SplineBasisData<d,dT,P_construction,s0,s1,sT0,sT1>::Mj0T_.apply
      (x, y, 0); // apply Mj0T to first block x1
    SplineBasisData<d,dT,P_construction,s0,s1,sT0,sT1>::Mj1T_.apply
      (x, y,     // apply Mj1T to second block x2 and add result
       SplineBasisData<d,dT,P_construction,s0,s1,sT0,sT1>::Mj0T_.column_dimension(), 0, true);
  }
  
  template <int d, int dT, SplineBasisFlavor flavor, int s0, int s1, int sT0, int sT1>
  template <class V>
  void
  SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::apply_Tj(const int j, const V& x, V& y) const
  { 
    y = x;
    V z(x);
    apply_Mj(j0(), z, y);
    for (int k = j0()+1; k <= j; k++) {
      apply_Mj(k, y, z);
      y.swap(z);
    }
  }

  template <int d, int dT, int s0, int s1, int sT0, int sT1>
  template <class V>
  void
  SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1>::apply_Tj(const int j, const V& x, V& y) const
  { 
    y = x;
    V z(x);
    apply_Mj(j0(), z, y);
    for (int k = j0()+1; k <= j; k++) {
      apply_Mj(k, y, z);
      y.swap(z);
    }
  }

  template <int d, int dT, SplineBasisFlavor flavor, int s0, int s1, int sT0, int sT1>
  template <class V>
  void
  SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::apply_Tj_transposed(const int j, const V& x, V& y) const
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
  
  template <int d, int dT, int s0, int s1, int sT0, int sT1>
  template <class V>
  void
  SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1>::apply_Tj_transposed(const int j, const V& x, V& y) const
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
  
  template <int d, int dT, SplineBasisFlavor flavor, int s0, int s1, int sT0, int sT1>
  void
  SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::apply_Tj_transposed(const int j,
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

  template <int d, int dT, int s0, int s1, int sT0, int sT1>
  void
  SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1>::apply_Tj_transposed(const int j,
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

  template <int d, int dT, SplineBasisFlavor flavor, int s0, int s1, int sT0, int sT1>
  void
  SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::apply_Tjinv(const int j, const Vector<double>& x, Vector<double>& y) const
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

  template <int d, int dT, int s0, int s1, int sT0, int sT1>
  void
  SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1>::apply_Tjinv(const int j, const Vector<double>& x, Vector<double>& y) const
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

  //
  //
  // point evaluation subroutines

  template <int d, int dT, SplineBasisFlavor flavor, int s0, int s1, int sT0, int sT1>
  SampledMapping<1>
  SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::evaluate
  (const typename SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::Index& lambda,
   const int resolution) const
  {
    // treat generator case separately
    if (lambda.e() == 0) {
      Grid<1> grid(0, 1, 1<<resolution);
      SampledMapping<1> result(grid); // zero
      Array1D<double> values((1<<resolution)+1);
      for (unsigned int i(0); i < values.size(); i++) {
  	const double x = i*ldexp(1.0, -resolution);
 	values[i] = evaluate(0, lambda, x);
      }     
      return SampledMapping<1>(grid, values);
    }
    
    InfiniteVector<double, typename SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::Index> coeffs;
    coeffs.set_coefficient(lambda, 1.0);
    return evaluate(coeffs, resolution);
  }

  template <int d, int dT, SplineBasisFlavor flavor, int s0, int s1, int sT0, int sT1>
  SampledMapping<1>
  SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::evaluate
  (const InfiniteVector<double, typename SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::Index>& coeffs,
   const int resolution) const
  {
    Grid<1> grid(0, 1, 1<<resolution);
    SampledMapping<1> result(grid); // zero
    if (coeffs.size() > 0) {
      // determine maximal level
      int jmax(0);
      typedef typename SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::Index Index;
      for (typename InfiniteVector<double,Index>::const_iterator it(coeffs.begin()),
  	     itend(coeffs.end()); it != itend; ++it)
  	jmax = std::max(it.index().j()+it.index().e(), jmax);
      
      // insert coefficients into a dense vector
      Vector<double> wcoeffs(Deltasize(jmax));
      for (typename InfiniteVector<double,Index>::const_iterator it(coeffs.begin()),
 	     itend(coeffs.end()); it != itend; ++it) {
 	// determine number of the wavelet
 	typedef typename Vector<double>::size_type size_type;
 	size_type number = 0;
 	if (it.index().e() == 0) {
 	  number = it.index().k()-DeltaLmin();
 	} else {
 	  number = Deltasize(it.index().j())+it.index().k()-Nablamin();
 	}
 	wcoeffs[number] = *it;
      }
      
      // switch to generator representation
      Vector<double> gcoeffs(wcoeffs.size(), false);
      if (jmax == j0())
 	gcoeffs = wcoeffs;
      else
 	apply_Tj(jmax-1, wcoeffs, gcoeffs);
      
      Array1D<double> values((1<<resolution)+1);
      for (unsigned int i(0); i < values.size(); i++) {
 	values[i] = 0;
 	const double x = i*ldexp(1.0, -resolution);
 	for (unsigned int k = 0; k < gcoeffs.size(); k++) {
	  values[i] += gcoeffs[k] * evaluate(0, Index(jmax, 0, DeltaLmin()+k), x);
 	}
      }
      
      return SampledMapping<1>(grid, values);
    }
    
    return result;
  }

  template <int d, int dT, int s0, int s1, int sT0, int sT1>
  SampledMapping<1>
  SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1>::evaluate
  (const InfiniteVector<double, typename SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1>::Index>& coeffs,
   const int resolution) const
  {
    Grid<1> grid(0, 1, 1<<resolution);
    SampledMapping<1> result(grid); // zero
    if (coeffs.size() > 0) {
      // determine maximal level
      int jmax(0);
      typedef typename SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1>::Index Index;
      for (typename InfiniteVector<double,Index>::const_iterator it(coeffs.begin()),
	     itend(coeffs.end()); it != itend; ++it)
	jmax = std::max(it.index().j()+it.index().e(), jmax);
      
      // insert coefficients into a dense vector
      Vector<double> wcoeffs(Deltasize(jmax));
      for (typename InfiniteVector<double,Index>::const_iterator it(coeffs.begin()),
	     itend(coeffs.end()); it != itend; ++it) {
	// determine number of the wavelet
	typedef typename Vector<double>::size_type size_type;
	size_type number = 0;
	if (it.index().e() == 0) {
	  number = it.index().k()-DeltaLmin();
	} else {
	  number = Deltasize(it.index().j())+it.index().k()-Nablamin();
	}
	wcoeffs[number] = *it;
      }
      
      // switch to generator representation
      Vector<double> gcoeffs(wcoeffs.size(), false);
      if (jmax == j0())
	gcoeffs = wcoeffs;
      else
	apply_Tj(jmax-1, wcoeffs, gcoeffs);
      
      Array1D<double> values((1<<resolution)+1);
      for (unsigned int i(0); i < values.size(); i++) {
	values[i] = 0;
	const double x = i*ldexp(1.0, -resolution);
	SchoenbergIntervalBSpline_td<d> sbs(jmax,0);
	for (unsigned int k = 0; k < gcoeffs.size(); k++) {
	  sbs.set_k(DeltaLmin()+k);
	  values[i] += gcoeffs[k] * sbs.value(Point<1>(x));
	}
      }
      
      return SampledMapping<1>(grid, values);
    }
    
    return result;
  }
  
  template <int d, int dT, SplineBasisFlavor flavor, int s0, int s1, int sT0, int sT1>
  double
  SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::evaluate
  (const unsigned int derivative,
   const typename SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::Index& lambda,
   const double x) const
  {
    assert(derivative <= 1); // we only support derivatives up to the first order
    
    double r = 0;
    
    if (lambda.e() == 0) {
      // generator
      if (lambda.k() < DeltaLmin()+(int)SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::CLA_.column_dimension()) {
	// left boundary generator
	for (unsigned int i(0); i < SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::CLA_.row_dimension(); i++) {
	  double help(SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::CLA_.get_entry(i, lambda.k()-DeltaLmin()));
	  if (help != 0)
	    r += help * (derivative == 0
			 ? EvaluateCardinalBSpline_td<d>  (lambda.j(), 1-ell2<d>()+i, x)
			 : EvaluateCardinalBSpline_td_x<d>(lambda.j(), 1-ell2<d>()+i, x));
	}
      }	else {
	if (lambda.k() > DeltaRmax(lambda.j())-(int)SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::CRA_.column_dimension()) {
	  // right boundary generator
	  for (unsigned int i(0); i < SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::CRA_.row_dimension(); i++) {
	    double help(SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::CRA_.get_entry(i, DeltaRmax(lambda.j())-lambda.k()));
	    if (help != 0)
	      r += help * (derivative == 0
			   ? EvaluateCardinalBSpline_td<d>  (lambda.j(), (1<<lambda.j())-(d%2)-(1-ell2<d>()+i), x)
			   : EvaluateCardinalBSpline_td_x<d>(lambda.j(), (1<<lambda.j())-(d%2)-(1-ell2<d>()+i), x));
	  }
	} else {
	  // inner generator
	  r = (derivative == 0
	       ? EvaluateCardinalBSpline_td<d>  (lambda.j(), lambda.k(), x)
	       : EvaluateCardinalBSpline_td_x<d>(lambda.j(), lambda.k(), x));
	}
      }
    } else {
      // wavelet, switch to generator representation
      typedef typename Vector<double>::size_type size_type;
      std::map<size_type,double> wc, gc;
      wc[lambda.k()-Nablamin()] = 1.0;
      SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::Mj1_.set_level(lambda.j());
      SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::Mj1_.apply(wc, gc, 0, 0);
      typedef typename SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::Index Index;
      for (typename std::map<size_type,double>::const_iterator it(gc.begin());
	   it != gc.end(); ++it) {
	r += it->second * evaluate(derivative, Index(lambda.j()+1, 0, DeltaLmin()+it->first), x);
      }
    }
    
    return r;
  }
  
  template <int d, int dT, int s0, int s1, int sT0, int sT1>
  double
  SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1>::evaluate
  (const unsigned int derivative,
   const typename SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1>::Index& lambda,
   const double x) const
  {
    assert(derivative <= 1); // we only support derivatives up to the first order

    double r = 0;

    if (lambda.e() == 0) {
      // generator
      if (lambda.k() > (1<<lambda.j())-ell1<d>()-d) {
	r = (derivative == 0
	     ? MathTL::EvaluateSchoenbergBSpline_td<d>  (lambda.j(),
							 (1<<lambda.j())-d-lambda.k()-2*ell1<d>(),
							 1-x)
	     : -MathTL::EvaluateSchoenbergBSpline_td_x<d>(lambda.j(),
							  (1<<lambda.j())-d-lambda.k()-2*ell1<d>(),
							  1-x));
      } else {
	r = (derivative == 0
	     ? MathTL::EvaluateSchoenbergBSpline_td<d>  (lambda.j(), lambda.k(), x)
	     : MathTL::EvaluateSchoenbergBSpline_td_x<d>(lambda.j(), lambda.k(), x));
      }
    } else {
      // wavelet, switch to generator representation
      typedef typename Vector<double>::size_type size_type;
      std::map<size_type,double> wc, gc;
      wc[lambda.k()-Nablamin()] = 1.0;
      SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::Mj1_.set_level(lambda.j());
      SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::Mj1_.apply(wc, gc, 0, 0);
      typedef typename SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1>::Index Index;
      for (typename std::map<size_type,double>::const_iterator it(gc.begin());
	   it != gc.end(); ++it) {
	r += it->second * evaluate(derivative, Index(lambda.j()+1, 0, DeltaLmin()+it->first), x);
      }
    }
    
    return r;
  }
  
  template <int d, int dT, SplineBasisFlavor flavor, int s0, int s1, int sT0, int sT1>
  void
  SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::evaluate
  (const unsigned int derivative,
   const typename SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::Index& lambda,
   const Array1D<double>& points, Array1D<double>& values) const
  {
    assert(derivative <= 1); // we only support derivatives up to the first order
    
    values.resize(points.size());
    for (unsigned int i(0); i < values.size(); i++)
      values[i] = 0;
    
    if (lambda.e() == 0) {
      // generator
      if (lambda.k() < DeltaLmin()+(int)SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::CLA_.column_dimension()) {
	// left boundary generator
	switch(derivative) {
	case 0:
	  for (unsigned int i(0); i < SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::CLA_.row_dimension(); i++) {
	    const double help(SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::CLA_.get_entry(i, lambda.k()-DeltaLmin()));
	    if (help != 0)
	      for (unsigned int m(0); m < points.size(); m++)
		values[m] += help * EvaluateCardinalBSpline_td<d>  (lambda.j(), 1-ell2<d>()+i, points[m]);
	  }
	  break;
	case 1:
	  for (unsigned int i(0); i < SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::CLA_.row_dimension(); i++) {
	    const double help(SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::CLA_.get_entry(i, lambda.k()-DeltaLmin()));
	    if (help != 0)
	      for (unsigned int m(0); m < points.size(); m++)
		values[m] += help * EvaluateCardinalBSpline_td_x<d>(lambda.j(), 1-ell2<d>()+i, points[m]);
	  }
	  break;
	default:
	  break;
	}
      }	else {
	if (lambda.k() > DeltaRmax(lambda.j())-(int)SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::CRA_.column_dimension()) {
	  // right boundary generator
	  switch(derivative) {
	  case 0:
	    for (unsigned int i(0); i < SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::CRA_.row_dimension(); i++) {
	      const double help(SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::CRA_.get_entry(i, DeltaRmax(lambda.j())-lambda.k()));
	      if (help != 0) {
		const int k = (1<<lambda.j())-ell1<d>()-ell2<d>()-(1-ell2<d>()+i);
		for (unsigned int m(0); m < points.size(); m++)
		  values[m] += help * EvaluateCardinalBSpline_td<d>  (lambda.j(), k, points[m]);
	      }
	    }
	    break;
	  case 1:
	    for (unsigned int i(0); i < SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::CRA_.row_dimension(); i++) {
	      const double help(SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::CRA_.get_entry(i, DeltaRmax(lambda.j())-lambda.k()));
	      if (help != 0) {
		const int k = (1<<lambda.j())-ell1<d>()-ell2<d>()-(1-ell2<d>()+i);
		for (unsigned int m(0); m < points.size(); m++)
		  values[m] += help * EvaluateCardinalBSpline_td_x<d>(lambda.j(), k, points[m]);
	      }
	    }
	    break;
	  default:
	    break;
	  }
	} else {
	  // inner generator
	  switch(derivative) {
	  case 0:
	    for (unsigned int m(0); m < points.size(); m++)
	      values[m] = EvaluateCardinalBSpline_td<d>  (lambda.j(), lambda.k(), points[m]);
	    break;
	  case 1:
	    for (unsigned int m(0); m < points.size(); m++)
	      values[m] = EvaluateCardinalBSpline_td_x<d>(lambda.j(), lambda.k(), points[m]);
	    break;
	  default:
	    break;
	  }
	}
      }
    } else {
      // wavelet, switch to generator representation
      typedef typename Vector<double>::size_type size_type;
      std::map<size_type,double> wc, gc;
      wc[lambda.k()-Nablamin()] = 1.0;
      SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::Mj1_.set_level(lambda.j());
      SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::Mj1_.apply(wc, gc, 0, 0);
      typedef typename SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::Index Index;
      Array1D<double> help(points.size());
      for (typename std::map<size_type,double>::const_iterator it(gc.begin());
	   it != gc.end(); ++it) {
	evaluate(derivative, Index(lambda.j()+1, 0, DeltaLmin()+it->first), points, help);
	for (unsigned int i = 0; i < points.size(); i++)
	  values[i] += it->second * help[i];
      }
    }
  }

//   template <int d, int dT>
//   void
//   SplineBasis<d,dT,DS_construction>::evaluate
//   (const unsigned int derivative,
//    const int j, const int e, const int k,
//    const Array1D<double>& points, Array1D<double>& values) const
//   {
//     assert(derivative <= 1); // we only support derivatives up to the first order
// 
//     values.resize(points.size());
//     for (unsigned int i(0); i < values.size(); i++)
//       values[i] = 0;
//     
//     if (e == 0) {
//       // generator
//       if (k < DeltaLmin()+(int)SplineBasisData<d,dT,DS_construction>::CLA_.column_dimension()) {
// 	// left boundary generator
// 	for (unsigned int i(0); i < SplineBasisData<d,dT,DS_construction>::CLA_.row_dimension(); i++) {
// 	  const double help(SplineBasisData<d,dT,DS_construction>::CLA_.get_entry(i, k-DeltaLmin()));
// 	  if (help != 0)
// 	    for (unsigned int m(0); m < points.size(); m++)
// 	      values[m] += help * (derivative == 0
// 				   ? EvaluateCardinalBSpline_td<d>  (j, 1-ell2<d>()+i, points[m])
// 				   : EvaluateCardinalBSpline_td_x<d>(j, 1-ell2<d>()+i, points[m]));
// 	}
//       }	else {
// 	if (k > DeltaRmax(j)-(int)SplineBasisData<d,dT,DS_construction>::CRA_.column_dimension()) {
// 	  // right boundary generator
// 	  for (unsigned int i(0); i < SplineBasisData<d,dT,DS_construction>::CRA_.row_dimension(); i++) {
// 	    const double help(SplineBasisData<d,dT,DS_construction>::CRA_.get_entry(i, DeltaRmax(j)-k));
//  	    if (help != 0)
// 	      for (unsigned int m(0); m < points.size(); m++)
// 		values[m] += help * (derivative == 0
// 				     ? EvaluateCardinalBSpline_td<d>  (j, (1<<j)-ell1<d>()-ell2<d>()-(1-ell2<d>()+i), points[m])
// 				     : EvaluateCardinalBSpline_td_x<d>(j, (1<<j)-ell1<d>()-ell2<d>()-(1-ell2<d>()+i), points[m]));
// 	  }
// 	} else {
// 	  // inner generator
// 	  for (unsigned int m(0); m < points.size(); m++)
// 	    values[m] = (derivative == 0
// 			 ? EvaluateCardinalBSpline_td<d>  (j, k, points[m])
// 			 : EvaluateCardinalBSpline_td_x<d>(j, k, points[m]));
// 	}
//       }
//     } else {
//       // wavelet, switch to generator representation
//       typedef typename Vector<double>::size_type size_type;
//       size_type number_lambda = Deltasize(j)+k-Nablamin();
//       std::map<size_type,double> wc, gc;
//       wc[number_lambda] = 1.0;
//       apply_Mj(j, wc, gc);
//       typedef typename SplineBasis<d,dT,DS_construction>::Index Index;
//       Array1D<double> help(points.size());
//       for (typename std::map<size_type,double>::const_iterator it(gc.begin());
// 	   it != gc.end(); ++it) {
// 	evaluate(derivative, j+1, 0, DeltaLmin()+it->first, points, help);
// 	for (unsigned int i = 0; i < points.size(); i++)
// 	  values[i] += it->second * help[i];
//       }
//     }
//   }
  
  template <int d, int dT, int s0, int s1, int sT0, int sT1>
  void
  SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1>::evaluate
  (const unsigned int derivative,
   const typename SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1>::Index& lambda,
   const Array1D<double>& points, Array1D<double>& values) const
  {
    assert(derivative <= 1); // we only support derivatives up to the first order

    values.resize(points.size());
    for (unsigned int i(0); i < values.size(); i++)
      values[i] = 0;
    
    if (lambda.e() == 0) {
      // generator
      if (lambda.k() > (1<<lambda.j())-ell1<d>()-d) {
	for (unsigned int m(0); m < points.size(); m++)
	  values[m] = (derivative == 0
		       ? MathTL::EvaluateSchoenbergBSpline_td<d>  (lambda.j(),
								   (1<<lambda.j())-d-lambda.k()-2*ell1<d>(),
								   1-points[m])
		       : -MathTL::EvaluateSchoenbergBSpline_td_x<d>(lambda.j(),
							    (1<<lambda.j())-d-lambda.k()-2*ell1<d>(),
							    1-points[m]));
      } else {
	for (unsigned int m(0); m < points.size(); m++)
	  values[m] = (derivative == 0
		       ? MathTL::EvaluateSchoenbergBSpline_td<d>  (lambda.j(),
								   lambda.k(),
								   points[m])
		       : MathTL::EvaluateSchoenbergBSpline_td_x<d>(lambda.j(),
								   lambda.k(),
								   points[m]));
      }
    } else {
#if 1
      // new version, determine first whether we have an interior wavelet
      typedef typename Vector<double>::size_type size_type;
      const size_type w_number = lambda.k()-Nablamin();
      if (w_number >= SplineBasisData<d,dT,P_construction,s0,s1,sT0,sT1>::Mj1_.ML_column_dimension()
	  && w_number < Nablasize(lambda.j())-SplineBasisData<d,dT,P_construction,s0,s1,sT0,sT1>::Mj1_.MR_column_dimension())
	{
	  // an interior wavelet
	  CDFBasis<d,dT>::evaluate(derivative, RIndex(lambda.j(),lambda.e(),lambda.k()), points, values);
	  if (w_number < (size_type)Nablasize(lambda.j()-1)) {
	    // left half
	    for (unsigned int i = 0; i < values.size(); i++)
	      values[i] *= SplineBasisData<d,dT,P_construction,s0,s1,sT0,sT1>::CDF_factor;
	  } else {
 	    // right half, reflect at 0.5 for d odd (i.e., multiply by (-1)^d)
 	    for (unsigned int i = 0; i < values.size(); i++)
 	      values[i] *= minus1power(d)*SplineBasisData<d,dT,P_construction,s0,s1,sT0,sT1>::CDF_factor;
	  }
	}
      else
	{
	  // boundary wavelet, use old code
	  std::map<size_type,double> wc, gc;
	  wc[w_number] = 1.0;
	  SplineBasisData<d,dT,P_construction,s0,s1,sT0,sT1>::Mj1_.set_level(lambda.j());
	  SplineBasisData<d,dT,P_construction,s0,s1,sT0,sT1>::Mj1_.apply(wc, gc, 0, 0);
	  typedef typename SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1>::Index Index;
	  Array1D<double> help(points.size());
	  for (typename std::map<size_type,double>::const_iterator it(gc.begin());
	       it != gc.end(); ++it) {
	    evaluate(derivative, Index(lambda.j()+1, 0, DeltaLmin()+it->first), points, help);
	    for (unsigned int i = 0; i < points.size(); i++)
	      values[i] += it->second * help[i];
	  }
	}
#else
      // old version, switch to generator representation for all wavelets
      typedef typename Vector<double>::size_type size_type;
      std::map<size_type,double> wc, gc;
      wc[lambda.k()-Nablamin()] = 1.0;
      SplineBasisData<d,dT,P_construction,s0,s1,sT0,sT1>::Mj1_.set_level(lambda.j());
      SplineBasisData<d,dT,P_construction,s0,s1,sT0,sT1>::Mj1_.apply(wc, gc, 0, 0);
      typedef typename SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1>::Index Index;
      Array1D<double> help(points.size());
      for (typename std::map<size_type,double>::const_iterator it(gc.begin());
	   it != gc.end(); ++it) {
	evaluate(derivative, Index(lambda.j()+1, 0, DeltaLmin()+it->first), points, help);
	for (unsigned int i = 0; i < points.size(); i++)
	  values[i] += it->second * help[i];
      }
#endif
    }
  }
  
  template <int d, int dT, SplineBasisFlavor flavor, int s0, int s1, int sT0, int sT1>
  void
  SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::evaluate
  (const typename SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::Index& lambda,
   const Array1D<double>& points, Array1D<double>& funcvalues, Array1D<double>& dervalues) const
  {
    const unsigned int npoints(points.size());
    funcvalues.resize(npoints);
    dervalues.resize(npoints);
    for (unsigned int i(0); i < npoints; i++) {
      funcvalues[i] = 0;
      dervalues[i] = 0;
    }

    if (lambda.e() == 0) {
      // generator
      if (lambda.k() < DeltaLmin()+(int)SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::CLA_.column_dimension()) {
	// left boundary generator
	for (unsigned int i(0); i < SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::CLA_.row_dimension(); i++) {
	  const double help(SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::CLA_.get_entry(i, lambda.k()-DeltaLmin()));
	  if (help != 0)
	    for (unsigned int m(0); m < npoints; m++) {
	      funcvalues[m] += help * EvaluateCardinalBSpline_td<d>  (lambda.j(), 1-ell2<d>()+i, points[m]);
	      dervalues[m]  += help * EvaluateCardinalBSpline_td_x<d>(lambda.j(), 1-ell2<d>()+i, points[m]);
	    }
	}
      }	else {
	if (lambda.k() > DeltaRmax(lambda.j())-(int)SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::CRA_.column_dimension()) {
	  // right boundary generator
	  for (unsigned int i(0); i < SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::CRA_.row_dimension(); i++) {
	    const double help(SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::CRA_.get_entry(i, DeltaRmax(lambda.j())-lambda.k()));
	    if (help != 0) {
	      const int k = (1<<lambda.j())-(d%2)-(1-ell2<d>()+i);
	      for (unsigned int m(0); m < npoints; m++) {
		funcvalues[m] += help * EvaluateCardinalBSpline_td<d>  (lambda.j(), k, points[m]);
		dervalues[m]  += help * EvaluateCardinalBSpline_td_x<d>(lambda.j(), k, points[m]);
	      }
	    }
	  }
	} else {
	  // inner generator
	  for (unsigned int m(0); m < npoints; m++) {
	    funcvalues[m] = EvaluateCardinalBSpline_td<d>  (lambda.j(), lambda.k(), points[m]);
	    dervalues[m]  = EvaluateCardinalBSpline_td_x<d>(lambda.j(), lambda.k(), points[m]);
	  }
	}
      }
    } else {
      // wavelet, switch to generator representation
      typedef typename Vector<double>::size_type size_type;
      std::map<size_type,double> wc, gc;
      wc[lambda.k()-Nablamin()] = 1.0;
      SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::Mj1_.set_level(lambda.j());
      SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::Mj1_.apply(wc, gc, 0, 0);
      typedef typename SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::Index Index;
      Array1D<double> help1, help2;
      for (typename std::map<size_type,double>::const_iterator it(gc.begin());
	   it != gc.end(); ++it) {
	evaluate(Index(lambda.j()+1, 0, DeltaLmin()+it->first), points, help1, help2);
	for (unsigned int i = 0; i < npoints; i++) {
	  funcvalues[i] += it->second * help1[i];
	  dervalues[i]  += it->second * help2[i];
	}
      }
    }
  }
  
  template <int d, int dT, int s0, int s1, int sT0, int sT1>
  void
  SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1>::evaluate
  (const typename SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1>::Index& lambda,
   const Array1D<double>& points, Array1D<double>& funcvalues, Array1D<double>& dervalues) const
  {
    const unsigned int npoints(points.size());
    funcvalues.resize(npoints);
    dervalues.resize(npoints);
    for (unsigned int i(0); i < npoints; i++) {
      funcvalues[i] = 0;
      dervalues[i] = 0;
    }
    
    if (lambda.e() == 0) {
      // generator
      if (lambda.k() > (1<<lambda.j())-ell1<d>()-d) {
	for (unsigned int m(0); m < npoints; m++) {
	  funcvalues[m] = MathTL::EvaluateSchoenbergBSpline_td<d>  (lambda.j(),
								    (1<<lambda.j())-d-lambda.k()-2*ell1<d>(),
								    1-points[m]);
	  dervalues[m]  = -MathTL::EvaluateSchoenbergBSpline_td_x<d>(lambda.j(),
								     (1<<lambda.j())-d-lambda.k()-2*ell1<d>(),
								     1-points[m]);
	}
      } else {
	for (unsigned int m(0); m < npoints; m++) {
	  funcvalues[m] = MathTL::EvaluateSchoenbergBSpline_td<d>  (lambda.j(),
								    lambda.k(),
								    points[m]);
	  dervalues[m]  = MathTL::EvaluateSchoenbergBSpline_td_x<d>(lambda.j(),
								    lambda.k(), 
								    points[m]);
	}
      }
    } else {
      // old version, switch to generator representation for all wavelets
      typedef typename Vector<double>::size_type size_type;
      std::map<size_type,double> wc, gc;
      wc[lambda.k()-Nablamin()] = 1.0;
      SplineBasisData<d,dT,P_construction,s0,s1,sT0,sT1>::Mj1_.set_level(lambda.j());
      SplineBasisData<d,dT,P_construction,s0,s1,sT0,sT1>::Mj1_.apply(wc, gc, 0, 0);
      typedef typename SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1>::Index Index;
      Array1D<double> help1, help2;
      for (typename std::map<size_type,double>::const_iterator it(gc.begin());
	   it != gc.end(); ++it) {
	evaluate(Index(lambda.j()+1, 0, DeltaLmin()+it->first), points, help1, help2);
	for (unsigned int i = 0; i < npoints; i++) {
	  funcvalues[i] += it->second * help1[i];
	  dervalues[i]  += it->second * help2[i];
	}
      }
    }
  }
  
  //
  //
  // expansion subroutines

  template <int d, int dT, int s0, int s1, int sT0, int sT1>
  void
  SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1>::expand
  (const Function<1>* f,
   const bool primal,
   const int jmax,
   InfiniteVector<double, typename SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1>::Index>& coeffs) const
  {
    Vector<double> coeffs_vector;
    expand(f, primal, jmax, coeffs_vector);
    typedef typename Vector<double>::size_type size_type;
    typedef typename SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1>::Index Index;
    size_type i(0);
    for (Index lambda(first_generator(j0()));
	 i < coeffs_vector.size(); ++lambda, i++)
      {
	const double coeff = coeffs_vector[i];
	if (fabs(coeff)>1e-15)
	  coeffs.set_coefficient(lambda, coeff);
      } 
  }
  
  template <int d, int dT, SplineBasisFlavor flavor, int s0, int s1, int sT0, int sT1>
  double
  SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::integrate
  (const Function<1>* f,
   const typename SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::Index& lambda) const
  {
    double r = 0;
    
    const int j = lambda.j()+lambda.e();
    int k1, k2;
    support(lambda, k1, k2);
    
    // setup Gauss points and weights for a composite quadrature formula:
    const unsigned int N_Gauss = 5;
    const double h = ldexp(1.0, -j);
    Array1D<double> gauss_points (N_Gauss*(k2-k1));
    for (int patch = k1; patch < k2; patch++) // refers to 2^{-j}[patch,patch+1]
      for (unsigned int n = 0; n < N_Gauss; n++)
	gauss_points[(patch-k1)*N_Gauss+n] = h*(2*patch+1+GaussPoints[N_Gauss-1][n])/2;
    
    // add all integral shares
    for (unsigned int n = 0; n < N_Gauss; n++)
      {
	const double gauss_weight = GaussWeights[N_Gauss-1][n] * h;
	for (int patch = k1; patch < k2; patch++)
	  {
	    const double t = gauss_points[(patch-k1)*N_Gauss+n];
	    
	    const double ft = f->value(Point<1>(t));
	    if (ft != 0)
	      r += ft
		* evaluate(0, lambda, t)
		* gauss_weight;
	  }
      }
    
    return r;
  }
  
  template <int d, int dT, SplineBasisFlavor flavor, int s0, int s1, int sT0, int sT1>
  double
  SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::integrate
  (const typename SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::Index& lambda,
   const typename SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::Index& mu) const
  {
    double r = 0;
    
    // First we compute the support intersection of \psi_\lambda and \psi_\mu:
    typedef typename SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::Support Support;
    Support supp;
    
    if (intersect_supports(*this, lambda, mu, supp))
      {
 	// Set up Gauss points and weights for a composite quadrature formula:
 	const unsigned int N_Gauss = d;
 	const double h = ldexp(1.0, -supp.j);
 	Array1D<double> gauss_points (N_Gauss*(supp.k2-supp.k1)), func1values, func2values;
 	for (int patch = supp.k1, id = 0; patch < supp.k2; patch++) // refers to 2^{-j}[patch,patch+1]
 	  for (unsigned int n = 0; n < N_Gauss; n++, id++)
 	    gauss_points[id] = h*(2*patch+1+GaussPoints[N_Gauss-1][n])/2.;
	
 	// - compute point values of the integrands
  	evaluate(0, lambda, gauss_points, func1values);
 	evaluate(0, mu, gauss_points, func2values);
	
 	// - add all integral shares
 	for (int patch = supp.k1, id = 0; patch < supp.k2; patch++)
 	  for (unsigned int n = 0; n < N_Gauss; n++, id++) {
 	    const double gauss_weight = GaussWeights[N_Gauss-1][n] * h;
	    r += func1values[id] * func2values[id] * gauss_weight;
 	  }
      }
    
    return r;
  }

  // generic expansion routine for DS-type wavelet bases
  template <int d, int dT, SplineBasisFlavor flavor, int s0, int s1, int sT0, int sT1>
  void
  SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::expand
  (const Function<1>* f,
   const bool primal,
   const int jmax,
   InfiniteVector<double, typename SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::Index>& coeffs) const
  {
    typedef typename SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::Index Index;
    
    for (Index lambda = first_generator(j0());;++lambda)
      {
 	coeffs.set_coefficient(lambda, integrate(f, lambda));
 	if (lambda == last_wavelet(jmax))
 	  break;
      }
    
    if (!primal) {
      // setup active index set
      std::set<Index> Lambda;
      for (Index lambda = first_generator(j0());; ++lambda) {
 	Lambda.insert(lambda);
	if (lambda == last_wavelet(jmax)) break;
      }
      
      // setup Gramian A_Lambda
      SparseMatrix<double> A_Lambda(Lambda.size());
      typedef typename SparseMatrix<double>::size_type size_type;     
      size_type row = 0;
      for (typename std::set<Index>::const_iterator it1(Lambda.begin()), itend(Lambda.end());
	   it1 != itend; ++it1, ++row)
	{
	  std::list<size_type> indices;
	  std::list<double> entries;
	  
	  size_type column = 0;
	  for (typename std::set<Index>::const_iterator it2(Lambda.begin());
	       it2 != itend; ++it2, ++column)
	    {
	      double entry = integrate(*it2, *it1);
	      
	      if (entry != 0) {
		indices.push_back(column);
		entries.push_back(entry);
	      }
	    }
	  A_Lambda.set_row(row, indices, entries);
	} 
 
      // solve A_Lambda*x = b
      Vector<double> b(Lambda.size());
      row = 0;
      for (typename std::set<Index>::const_iterator it(Lambda.begin()), itend(Lambda.end());
	   it != itend; ++it, ++row)
	b[row] = coeffs.get_coefficient(*it);
      
      Vector<double> x(b);
      unsigned int iterations;
      CG(A_Lambda, b, x, 1e-15, 500, iterations);
      
      coeffs.clear();
      row = 0;
      for (typename std::set<Index>::const_iterator it(Lambda.begin()), itend(Lambda.end());
	   it != itend; ++it, ++row)
	coeffs.set_coefficient(*it, x[row]);
    }
  }
  
  template <int d, int dT, int s0, int s1, int sT0, int sT1>
  void
  SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1>::expand
  (const Function<1>* f,
   const bool primal,
   const int jmax,
   Vector<double>& coeffs) const
  {
    assert(jmax >= j0());
	   
    coeffs.resize(Deltasize(jmax+1));

    // 1. compute integrals w.r.t. the primal generators on level jmax
    Vector<double> coeffs_phijk(coeffs.size(), false);
    SimpsonRule simpson;
    CompositeRule<1> composite(simpson, 24); // should be sufficient for many cases
    SchoenbergIntervalBSpline_td<d> sbs(jmax+1,0);
    for (int k = DeltaLmin(); k <= DeltaRmax(jmax+1); k++) {
      sbs.set_k(k);
      ProductFunction<1> integrand(f, &sbs);
      coeffs_phijk[k-DeltaLmin()]
	= composite.integrate(integrand,
			      Point<1>(std::max(0.0, (k+ell1<d>())*ldexp(1.0, -jmax-1))),
			      Point<1>(std::min(1.0, (k+ell2<d>())*ldexp(1.0, -jmax-1))));
    }
    // 2. transform rhs into that of psi_{j,k} basis: apply T_{j-1}^T
    Vector<double> rhs(coeffs.size(), false);
    apply_Tj_transposed(jmax, coeffs_phijk, rhs);
    
    if (!primal) {
      FullGramian<d,dT,s0,s1,sT0,sT1> G(*this);
      G.set_level(jmax+1);
      unsigned int iterations;
      CG(G, rhs, coeffs, 1e-15, 500, iterations);
    } else {
      coeffs.swap(rhs);
    }
  }

  template <int d, int dT, SplineBasisFlavor flavor, int s0, int s1, int sT0, int sT1>
  void
  SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::setup_full_collection()
  {
    if (jmax_ == -1 || jmax_ < j0()) {
      cout << "SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::setup_full_collection(): specify a mximal level of resolution first!" << endl;
      abort();
    }   

    int degrees_of_freedom = Deltasize(jmax_+1);
    cout << "total degrees of freedom between j0_ and jmax_ is " << degrees_of_freedom << endl;

    cout << "setting up collection of wavelet indices..." << endl;
    full_collection.resize(degrees_of_freedom);
    int k = 0;
    for (Index ind = first_generator(j0()); ind <= last_wavelet(jmax_); ++ind) {
      full_collection[k] = ind;
      k++;
    }
    cout << "done setting up collection of wavelet indices..." << endl;
  }

  template <int d, int dT, int s0, int s1, int sT0, int sT1>
  void
  SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1>::setup_full_collection()
  {
    if (jmax_ == -1 || jmax_ < j0()) {
      cout << "SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::setup_full_collection(): specify a mximal level of resolution first!" << endl;
      abort();
    }   

    int degrees_of_freedom = Deltasize(jmax_+1);
    cout << "total degrees of freedom between j0_ and jmax_ is " << degrees_of_freedom << endl;

    cout << "setting up collection of wavelet indices..." << endl;
    full_collection.resize(degrees_of_freedom);
    int k = 0;
    for (Index ind = first_generator(j0()); ind <= last_wavelet(jmax_); ++ind) {
      full_collection[k] = ind;
      k++;
    }
    cout << "done setting up collection of wavelet indices..." << endl;
  }


}
