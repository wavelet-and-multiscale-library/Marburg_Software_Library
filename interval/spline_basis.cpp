// implementation for spline_basis.h

#include <numerics/cardinal_splines.h>
#include <numerics/schoenberg_splines.h>

#include <Rd/r_index.h>
#include <Rd/cdf_utils.h>
#include <interval/interval_bspline.h>

namespace WaveletTL
{
  template <int d, int dT, SplineBasisFlavor flavor>
  SplineBasis<d,dT,flavor>::SplineBasis
  (const char* options,
   const int s0, const int s1, const int sT0, const int sT1)
    : SplineBasisData<d,dT,flavor>(options,s0,s1,sT0,sT1)
  {
    assert(flavor == P_construction);
    
    // cf. PBasis<d,dT> ...
    DeltaLmin_       = 1-d-ell1<d>()+s0;
    DeltaRmax_offset = -1-ell1<d>()-s1;
  }

  template <int d, int dT>
  SplineBasis<d,dT,DS_construction>::SplineBasis
  (const char* options,
   const int s0, const int s1, const int sT0, const int sT1)
    : SplineBasisData<d,dT,DS_construction>(options,s0,s1,sT0,sT1)
  {
    // cf. DSBasis<d,dT> ...
    DeltaLmin_       = ell2T<d,dT>()+s0+sT0-dT;
    DeltaRmax_offset = -(d%2)-(ell2T<d,dT>()+s1+sT1-dT);
  }
  
  template <int d, int dT, SplineBasisFlavor flavor>
  inline
  typename SplineBasis<d,dT,flavor>::Index
  SplineBasis<d,dT,flavor>::first_generator(const int j) const
  {
    assert(j >= j0());
    return Index(j, 0, DeltaLmin(), this);
  }
  
  template <int d, int dT, SplineBasisFlavor flavor>
  inline
  typename SplineBasis<d,dT,flavor>::Index
  SplineBasis<d,dT,flavor>::last_generator(const int j) const
  {
    assert(j >= j0());
    return Index(j, 0, DeltaRmax(j), this);
  }

  template <int d, int dT, SplineBasisFlavor flavor>
  inline
  typename SplineBasis<d,dT,flavor>::Index
  SplineBasis<d,dT,flavor>::first_wavelet(const int j) const
  {
    assert(j >= j0());
    return Index(j, 1, Nablamin(), this);
  }
  
  template <int d, int dT, SplineBasisFlavor flavor>
  inline
  typename SplineBasis<d,dT,flavor>::Index
  SplineBasis<d,dT,flavor>::last_wavelet(const int j) const
  {
    assert(j >= j0());
    return Index(j, 1, Nablamax(j), this);
  }

  template <int d, int dT, SplineBasisFlavor flavor>
  void
  SplineBasis<d,dT,flavor>::support(const Index& lambda, int& k1, int& k2) const {
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
  
  template <int d, int dT, SplineBasisFlavor flavor>
  template <class V>
  void
  SplineBasis<d,dT,flavor>::apply_Mj(const int j, const V& x, V& y) const
  {
    SplineBasisData<d,dT,flavor>::Mj0_.set_level(j);
    SplineBasisData<d,dT,flavor>::Mj1_.set_level(j);

    // decompose x=(x1 x2) appropriately
    SplineBasisData<d,dT,flavor>::Mj0_.apply(x, y, 0, 0); // apply Mj0 to first block x1
    SplineBasisData<d,dT,flavor>::Mj1_.apply(x, y,        // apply Mj1 to second block x2 and add result
					     SplineBasisData<d,dT,flavor>::Mj0_.column_dimension(), 0, true);
  }

  template <int d, int dT>
  template <class V>
  void
  SplineBasis<d,dT,DS_construction>::apply_Mj(const int j, const V& x, V& y) const
  {
    SplineBasisData<d,dT,DS_construction>::Mj0_.set_level(j);
    SplineBasisData<d,dT,DS_construction>::Mj1_.set_level(j);

    // decompose x=(x1 x2) appropriately
    SplineBasisData<d,dT,DS_construction>::Mj0_.apply(x, y, 0, 0); // apply Mj0 to first block x1
    SplineBasisData<d,dT,DS_construction>::Mj1_.apply(x, y,        // apply Mj1 to second block x2 and add result
						      SplineBasisData<d,dT,DS_construction>::Mj0_.column_dimension(), 0, true);
  }
  
  template <int d, int dT, SplineBasisFlavor flavor>
  template <class V>
  void
  SplineBasis<d,dT,flavor>::apply_Mj_transposed(const int j, const V& x, V& y) const
  {
    SplineBasisData<d,dT,flavor>::Mj0_.set_level(j);
    SplineBasisData<d,dT,flavor>::Mj1_.set_level(j);

    // y=(y1 y2) is a block vector
    SplineBasisData<d,dT,flavor>::Mj0_.apply_transposed(x, y, 0, 0); // write into first block y1
    SplineBasisData<d,dT,flavor>::Mj1_.apply_transposed(x, y, 0,     // write into second block y2
							SplineBasisData<d,dT,flavor>::Mj0_.column_dimension());
  }

  template <int d, int dT>
  template <class V>
  void
  SplineBasis<d,dT,DS_construction>::apply_Mj_transposed(const int j, const V& x, V& y) const
  {
    SplineBasisData<d,dT,DS_construction>::Mj0_.set_level(j);
    SplineBasisData<d,dT,DS_construction>::Mj1_.set_level(j);

    // y=(y1 y2) is a block vector
    SplineBasisData<d,dT,DS_construction>::Mj0_.apply_transposed(x, y, 0, 0); // write into first block y1
    SplineBasisData<d,dT,DS_construction>::Mj1_.apply_transposed(x, y, 0,     // write into second block y2
								 SplineBasisData<d,dT,DS_construction>::Mj0_.column_dimension());
  }

  template <int d, int dT, SplineBasisFlavor flavor>
  template <class V>
  void
  SplineBasis<d,dT,flavor>::apply_Gj(const int j, const V& x, V& y) const
  {
    SplineBasisData<d,dT,flavor>::Mj0T_.set_level(j);
    SplineBasisData<d,dT,flavor>::Mj1T_.set_level(j);

    // y=(y1 y2) is a block vector
    SplineBasisData<d,dT,flavor>::Mj0T_.apply_transposed(x, y, 0, 0); // write into first block y1
    SplineBasisData<d,dT,flavor>::Mj1T_.apply_transposed(x, y, 0,     // write into second block y2
							 SplineBasisData<d,dT,flavor>::Mj0T_.column_dimension());
  }
  
  template <int d, int dT>
  template <class V>
  void
  SplineBasis<d,dT,DS_construction>::apply_Gj(const int j, const V& x, V& y) const
  {
    SplineBasisData<d,dT,DS_construction>::Mj0T_.set_level(j);
    SplineBasisData<d,dT,DS_construction>::Mj1T_.set_level(j);

    // y=(y1 y2) is a block vector
    SplineBasisData<d,dT,DS_construction>::Mj0T_.apply_transposed(x, y, 0, 0); // write into first block y1
    SplineBasisData<d,dT,DS_construction>::Mj1T_.apply_transposed(x, y, 0,     // write into second block y2
								  SplineBasisData<d,dT,DS_construction>::Mj0T_.column_dimension());
  }
  
  template <int d, int dT, SplineBasisFlavor flavor>
  template <class V>
  void
  SplineBasis<d,dT,flavor>::apply_Gj_transposed(const int j, const V& x, V& y) const
  {
    SplineBasisData<d,dT,flavor>::Mj0T_.set_level(j);
    SplineBasisData<d,dT,flavor>::Mj1T_.set_level(j);
    
    // decompose x=(x1 x2) appropriately
    SplineBasisData<d,dT,flavor>::Mj0T_.apply(x, y, 0); // apply Mj0T to first block x1
    SplineBasisData<d,dT,flavor>::Mj1T_.apply(x, y,     // apply Mj1T to second block x2 and add result
					      SplineBasisData<d,dT,flavor>::Mj0T_.column_dimension(), 0, true);
  }
  
  template <int d, int dT>
  template <class V>
  void
  SplineBasis<d,dT,DS_construction>::apply_Gj_transposed(const int j, const V& x, V& y) const
  {
    SplineBasisData<d,dT,DS_construction>::Mj0T_.set_level(j);
    SplineBasisData<d,dT,DS_construction>::Mj1T_.set_level(j);
    
    // decompose x=(x1 x2) appropriately
    SplineBasisData<d,dT,DS_construction>::Mj0T_.apply(x, y, 0); // apply Mj0T to first block x1
    SplineBasisData<d,dT,DS_construction>::Mj1T_.apply(x, y,     // apply Mj1T to second block x2 and add result
						       SplineBasisData<d,dT,DS_construction>::Mj0T_.column_dimension(), 0, true);
  }
  
  template <int d, int dT, SplineBasisFlavor flavor>
  template <class V>
  void
  SplineBasis<d,dT,flavor>::apply_Tj(const int j, const V& x, V& y) const
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
  SplineBasis<d,dT,DS_construction>::apply_Tj(const int j, const V& x, V& y) const
  { 
    y = x;
    V z(x);
    apply_Mj(j0(), z, y);
    for (int k = j0()+1; k <= j; k++) {
      apply_Mj(k, y, z);
      y.swap(z);
    }
  }

  template <int d, int dT, SplineBasisFlavor flavor>
  template <class V>
  void
  SplineBasis<d,dT,flavor>::apply_Tj_transposed(const int j, const V& x, V& y) const
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
  template <class V>
  void
  SplineBasis<d,dT,DS_construction>::apply_Tj_transposed(const int j, const V& x, V& y) const
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
  
  template <int d, int dT, SplineBasisFlavor flavor>
  void
  SplineBasis<d,dT,flavor>::apply_Tj_transposed(const int j,
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
  SplineBasis<d,dT,DS_construction>::apply_Tj_transposed(const int j,
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

  template <int d, int dT, SplineBasisFlavor flavor>
  void
  SplineBasis<d,dT,flavor>::apply_Tjinv(const int j, const Vector<double>& x, Vector<double>& y) const
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

  template <int d, int dT>
  void
  SplineBasis<d,dT,DS_construction>::apply_Tjinv(const int j, const Vector<double>& x, Vector<double>& y) const
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

  template <int d, int dT, SplineBasisFlavor flavor>
  SampledMapping<1>
  SplineBasis<d,dT,flavor>::evaluate
  (const InfiniteVector<double, typename SplineBasis<d,dT,flavor>::Index>& coeffs,
   const int resolution) const
  {
    assert(flavor == P_construction); // the other cases are done in a template specialization
    
    Grid<1> grid(0, 1, 1<<resolution);
    SampledMapping<1> result(grid); // zero
    if (coeffs.size() > 0) {
      // determine maximal level
      int jmax(0);
      typedef typename SplineBasis<d,dT,flavor>::Index Index;
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
  
  template <int d, int dT>
  SampledMapping<1>
  SplineBasis<d,dT,DS_construction>::evaluate
  (const InfiniteVector<double, typename SplineBasis<d,dT,DS_construction>::Index>& coeffs,
   const int resolution) const
  {
    Grid<1> grid(0, 1, 1<<resolution);
    SampledMapping<1> result(grid); // zero
    if (coeffs.size() > 0) {
      // determine maximal level
      int jmax(0);
      typedef typename SplineBasis<d,dT,DS_construction>::Index Index;
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
	  values[i] += gcoeffs[k] * evaluate(0, Index(jmax, 0, DeltaLmin()+k, this), x);
 	}
      }
      
      return SampledMapping<1>(grid, values);
    }
    
    return result;
  }
  
  template <int d, int dT, SplineBasisFlavor flavor>
  double
  SplineBasis<d,dT,flavor>::evaluate
  (const unsigned int derivative,
   const typename SplineBasis<d,dT,flavor>::Index& lambda,
   const double x) const
  {
    assert(flavor == P_construction); // the other cases are done in a template specialization
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
      size_type number_lambda = Deltasize(lambda.j())+lambda.k()-Nablamin();
      std::map<size_type,double> wc, gc;
      wc[number_lambda] = 1.0;
      apply_Mj(lambda.j(), wc, gc);
      typedef typename SplineBasis<d,dT,flavor>::Index Index;
      for (typename std::map<size_type,double>::const_iterator it(gc.begin());
	   it != gc.end(); ++it) {
	r += it->second * evaluate(derivative, Index(lambda.j()+1, 0, DeltaLmin()+it->first, this), x);
      }
    }
    
    return r;
  }
  
  template <int d, int dT>
  double
  SplineBasis<d,dT,DS_construction>::evaluate
  (const unsigned int derivative,
   const typename SplineBasis<d,dT,DS_construction>::Index& lambda,
   const double x) const
  {
    assert(derivative <= 1); // we only support derivatives up to the first order
    
    double r = 0;
    
    if (lambda.e() == 0) {
      // generator
      if (lambda.k() < DeltaLmin()+(int)SplineBasisData<d,dT,DS_construction>::CLA_.column_dimension()) {
	// left boundary generator
	for (unsigned int i(0); i < SplineBasisData<d,dT,DS_construction>::CLA_.row_dimension(); i++) {
	  double help(SplineBasisData<d,dT,DS_construction>::CLA_.get_entry(i, lambda.k()-DeltaLmin()));
	  if (help != 0)
	    r += help * (derivative == 0
			 ? EvaluateCardinalBSpline_td<d>  (lambda.j(), 1-ell2<d>()+i, x)
			 : EvaluateCardinalBSpline_td_x<d>(lambda.j(), 1-ell2<d>()+i, x));
	}
      }	else {
	if (lambda.k() > DeltaRmax(lambda.j())-(int)SplineBasisData<d,dT,DS_construction>::CRA_.column_dimension()) {
	  // right boundary generator
	  for (unsigned int i(0); i < SplineBasisData<d,dT,DS_construction>::CRA_.row_dimension(); i++) {
	    double help(SplineBasisData<d,dT,DS_construction>::CRA_.get_entry(i, DeltaRmax(lambda.j())-lambda.k()));
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
      size_type number_lambda = Deltasize(lambda.j())+lambda.k()-Nablamin();
      std::map<size_type,double> wc, gc;
      wc[number_lambda] = 1.0;
      apply_Mj(lambda.j(), wc, gc);
      typedef typename SplineBasis<d,dT,DS_construction>::Index Index;
      for (typename std::map<size_type,double>::const_iterator it(gc.begin());
	   it != gc.end(); ++it) {
	r += it->second * evaluate(derivative, Index(lambda.j()+1, 0, DeltaLmin()+it->first, this), x);
      }
    }
    
    return r;
  }

  template <int d, int dT, SplineBasisFlavor flavor>
  void
  SplineBasis<d,dT,flavor>::evaluate
  (const unsigned int derivative,
   const typename SplineBasis<d,dT,flavor>::Index& lambda,
   const Array1D<double>& points, Array1D<double>& values) const
  {
    assert(flavor == P_construction); // the other cases are done in a template specialization
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
      // wavelet, switch to generator representation
      typedef typename Vector<double>::size_type size_type;
      size_type number_lambda = Deltasize(lambda.j())+lambda.k()-Nablamin();
      std::map<size_type,double> wc, gc;
      wc[number_lambda] = 1.0;
      apply_Mj(lambda.j(), wc, gc);
      typedef typename SplineBasis<d,dT,flavor>::Index Index;
      Array1D<double> help(points.size());
      for (typename std::map<size_type,double>::const_iterator it(gc.begin());
	   it != gc.end(); ++it) {
	evaluate(derivative, Index(lambda.j()+1, 0, DeltaLmin()+it->first, this), points, help);
	for (unsigned int i = 0; i < points.size(); i++)
	  values[i] += it->second * help[i];
      }
    }
  }
  
  template <int d, int dT>
  void
  SplineBasis<d,dT,DS_construction>::evaluate
  (const unsigned int derivative,
   const typename SplineBasis<d,dT,DS_construction>::Index& lambda,
   const Array1D<double>& points, Array1D<double>& values) const
  {
    assert(derivative <= 1); // we only support derivatives up to the first order

    values.resize(points.size());
    for (unsigned int i(0); i < values.size(); i++)
      values[i] = 0;
    
    if (lambda.e() == 0) {
      // generator
      if (lambda.k() < DeltaLmin()+(int)SplineBasisData<d,dT,DS_construction>::CLA_.column_dimension()) {
	// left boundary generator
	for (unsigned int i(0); i < SplineBasisData<d,dT,DS_construction>::CLA_.row_dimension(); i++) {
	  const double help(SplineBasisData<d,dT,DS_construction>::CLA_.get_entry(i, lambda.k()-DeltaLmin()));
	  // 	  if (help != 0)
	  for (unsigned int m(0); m < points.size(); m++)
	    values[m] += help * (derivative == 0
				 ? EvaluateCardinalBSpline_td<d>  (lambda.j(), 1-ell2<d>()+i, points[m])
				 : EvaluateCardinalBSpline_td_x<d>(lambda.j(), 1-ell2<d>()+i, points[m]));
	}
      }	else {
	if (lambda.k() > DeltaRmax(lambda.j())-(int)SplineBasisData<d,dT,DS_construction>::CRA_.column_dimension()) {
	  // right boundary generator
	  for (unsigned int i(0); i < SplineBasisData<d,dT,DS_construction>::CRA_.row_dimension(); i++) {
	    const double help(SplineBasisData<d,dT,DS_construction>::CRA_.get_entry(i, DeltaRmax(lambda.j())-lambda.k()));
	    // 	    if (help != 0)
	    for (unsigned int m(0); m < points.size(); m++)
	      values[m] += help * (derivative == 0
				   ? EvaluateCardinalBSpline_td<d>  (lambda.j(), (1<<lambda.j())-ell1<d>()-ell2<d>()-(1-ell2<d>()+i), points[m])
				   : EvaluateCardinalBSpline_td_x<d>(lambda.j(), (1<<lambda.j())-ell1<d>()-ell2<d>()-(1-ell2<d>()+i), points[m]));
	  }
	} else {
	  // inner generator
	  for (unsigned int m(0); m < points.size(); m++)
	    values[m] = (derivative == 0
			 ? EvaluateCardinalBSpline_td<d>  (lambda.j(), lambda.k(), points[m])
			 : EvaluateCardinalBSpline_td_x<d>(lambda.j(), lambda.k(), points[m]));
	}
      }
    } else {
      // wavelet, switch to generator representation
      typedef typename Vector<double>::size_type size_type;
      size_type number_lambda = Deltasize(lambda.j())+lambda.k()-Nablamin();
      std::map<size_type,double> wc, gc;
      wc[number_lambda] = 1.0;
      apply_Mj(lambda.j(), wc, gc);
      typedef typename SplineBasis<d,dT,DS_construction>::Index Index;
      Array1D<double> help(points.size());
      for (typename std::map<size_type,double>::const_iterator it(gc.begin());
	   it != gc.end(); ++it) {
	evaluate(derivative, Index(lambda.j()+1, 0, DeltaLmin()+it->first, this), points, help);
	for (unsigned int i = 0; i < points.size(); i++)
	  values[i] += it->second * help[i];
      }
    }
  }
  
  template <int d, int dT, SplineBasisFlavor flavor>
  void
  SplineBasis<d,dT,flavor>::evaluate
  (const typename SplineBasis<d,dT,flavor>::Index& lambda,
   const Array1D<double>& points, Array1D<double>& funcvalues, Array1D<double>& dervalues) const
  {
    assert(flavor == P_construction); // the other cases are done in a template specialization

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
      // wavelet, switch to generator representation
      typedef typename Vector<double>::size_type size_type;
      size_type number_lambda = Deltasize(lambda.j())+lambda.k()-Nablamin();
      std::map<size_type,double> wc, gc;
      wc[number_lambda] = 1.0;
      apply_Mj(lambda.j(), wc, gc);
      typedef typename SplineBasis<d,dT,flavor>::Index Index;
      Array1D<double> help1, help2;
      for (typename std::map<size_type,double>::const_iterator it(gc.begin());
	   it != gc.end(); ++it) {
	evaluate(Index(lambda.j()+1, 0, DeltaLmin()+it->first, this), points, help1, help2);
	for (unsigned int i = 0; i < npoints; i++) {
	  funcvalues[i] += it->second * help1[i];
	  dervalues[i]  += it->second * help2[i];
	}
      }
    }
  }
  
  template <int d, int dT>
  void
  SplineBasis<d,dT,DS_construction>::evaluate
  (const typename SplineBasis<d,dT,DS_construction>::Index& lambda,
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
      if (lambda.k() < DeltaLmin()+(int)SplineBasisData<d,dT,DS_construction>::CLA_.column_dimension()) {
	// left boundary generator
	for (unsigned int i(0); i < SplineBasisData<d,dT,DS_construction>::CLA_.row_dimension(); i++) {
	  const double help(SplineBasisData<d,dT,DS_construction>::CLA_.get_entry(i, lambda.k()-DeltaLmin()));
	  // 	  if (help != 0)
	  for (unsigned int m(0); m < npoints; m++) {
	    funcvalues[m] += help * EvaluateCardinalBSpline_td<d>  (lambda.j(), 1-ell2<d>()+i, points[m]);
	    dervalues[m]  += help * EvaluateCardinalBSpline_td_x<d>(lambda.j(), 1-ell2<d>()+i, points[m]);
	  }
	}
      }	else {
	if (lambda.k() > DeltaRmax(lambda.j())-(int)SplineBasisData<d,dT,DS_construction>::CRA_.column_dimension()) {
	  // right boundary generator
	  for (unsigned int i(0); i < SplineBasisData<d,dT,DS_construction>::CRA_.row_dimension(); i++) {
	    const double help(SplineBasisData<d,dT,DS_construction>::CRA_.get_entry(i, DeltaRmax(lambda.j())-lambda.k()));
	    // 	    if (help != 0)
	    for (unsigned int m(0); m < npoints; m++) {
	      funcvalues[m] += help * EvaluateCardinalBSpline_td<d>  (lambda.j(), (1<<lambda.j())-(d%2)-(1-ell2<d>()+i), points[m]);
	      dervalues[m]  += help * EvaluateCardinalBSpline_td_x<d>(lambda.j(), (1<<lambda.j())-(d%2)-(1-ell2<d>()+i), points[m]);
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
      size_type number_lambda = Deltasize(lambda.j())+lambda.k()-Nablamin();
      std::map<size_type,double> wc, gc;
      wc[number_lambda] = 1.0;
      apply_Mj(lambda.j(), wc, gc);
      typedef typename SplineBasis<d,dT,DS_construction>::Index Index;
      Array1D<double> help1, help2;
      for (typename std::map<size_type,double>::const_iterator it(gc.begin());
	   it != gc.end(); ++it) {
	evaluate(Index(lambda.j()+1, 0, DeltaLmin()+it->first, this), points, help1, help2);
	for (unsigned int i = 0; i < npoints; i++) {
	  funcvalues[i] += it->second * help1[i];
	  dervalues[i]  += it->second * help2[i];
	}
      }
    }
  }
  

}
