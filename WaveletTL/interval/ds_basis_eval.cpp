// implementation for ds_basis.h, point evaluation on a dyadic grid

#include <Rd/r_index.h>
#include <Rd/cdf_utils.h>
#include <numerics/cardinal_splines.h>

namespace WaveletTL
{
  template <int d, int dT, DSBiorthogonalizationMethod BIO>
  SampledMapping<1>
  DSBasis<d,dT,BIO>::evaluate
  (const typename DSBasis<d,dT,BIO>::Index& lambda,
   const bool primal,
   const int resolution) const
  {
    if (lambda.e() == 0) { // generator
      if (primal) {
	const Matrix<double>& CLA = get_CLA();
	if (lambda.k() < DeltaLmin()+(int)CLA.column_dimension()) {
	  // left boundary generator
	  InfiniteVector<double, RIndex> coeffs;
	  for(unsigned int i(0); i < CLA.row_dimension(); i++) {
	    double v(CLA(i, lambda.k()-DeltaLmin()));
	    if (v != 0)
	      coeffs.set_coefficient(RIndex(lambda.j(), 0, 1-ell2<d>()+i), v);
	  }
	  return get_CDF_basis().evaluate(0, coeffs, primal, 0, 1, resolution);
	} else {
	  const Matrix<double>& CRA = get_CRA();
	  if (lambda.k() > DeltaRmax(lambda.j())-(int)CRA.column_dimension()) {
	    // right boundary generator
	    InfiniteVector<double, RIndex> coeffs;
	    for (unsigned int i(0); i < CRA.row_dimension(); i++) {
	      double v(CRA(i, DeltaRmax(lambda.j())-lambda.k()));
	      if (v != 0)
		coeffs.set_coefficient(RIndex(lambda.j(), 0, (1<<lambda.j())-ell1<d>()-ell2<d>()-(1-ell2<d>()+i)), v);
	    }
	    return get_CDF_basis().evaluate(0, coeffs, primal, 0, 1, resolution);
	  } else {
	    // inner generator
	    return get_CDF_basis().evaluate(0, RIndex(lambda.j(), 0, lambda.k()),
						  primal, 0, 1, resolution);
	  }
	}
      } else {
	// dual
	const Matrix<double>& CLAT = get_CLAT();
	if (lambda.k() < DeltaLTmin()+(int)CLAT.column_dimension()) {
	  // left boundary generator
	  InfiniteVector<double, RIndex> coeffs;
	  for (unsigned int i(0); i < CLAT.row_dimension(); i++) {
	    double v(CLAT(i, lambda.k()-DeltaLTmin()));
	    if (v != 0)
	      coeffs.set_coefficient(RIndex(lambda.j(), 0, 1-ell2T<d,dT>()+i), v);
	  }
	  return get_CDF_basis().evaluate(0, coeffs, primal, 0, 1, resolution);
	} else {
	  const Matrix<double>& CRAT = get_CRAT();
	  if (lambda.k() > DeltaRTmax(lambda.j())-(int)CRAT.column_dimension()) {
	    // right boundary generator
	    InfiniteVector<double, RIndex> coeffs;
	    for (unsigned int i(0); i < CRAT.row_dimension(); i++) {
	      double v(CRAT(i, DeltaRTmax(lambda.j())-lambda.k()));
	      if (v != 0)
		coeffs.set_coefficient(RIndex(lambda.j(), 0, (1<<lambda.j())-ell1<d>()-ell2<d>()-(1-ell2T<d,dT>()+i)), v);
	    }
	    return get_CDF_basis().evaluate(0, coeffs, primal, 0, 1, resolution);
	  } else {
	    // inner generator
	    return get_CDF_basis().evaluate(0, RIndex(lambda.j(), 0, lambda.k()),
						  primal, 0, 1, resolution);
	  }
	}
      }
    } else { // wavelet
      InfiniteVector<double, typename DSBasis<d,dT,BIO>::Index> gcoeffs;
      if (primal)
	reconstruct_1(lambda, lambda.j()+1, gcoeffs);
      else
	reconstruct_t_1(lambda, lambda.j()+1, gcoeffs);

      return evaluate(gcoeffs, primal, resolution);
    }
    
    return SampledMapping<1>(); // dummy return for the compiler
  }

  template <int d, int dT, DSBiorthogonalizationMethod BIO>
  SampledMapping<1>
  DSBasis<d,dT,BIO>::evaluate
  (const InfiniteVector<double, typename DSBasis<d,dT,BIO>::Index>& coeffs,
   const bool primal,
   const int resolution) const
  {
    SampledMapping<1> result(Grid<1>(0, 1, 1<<resolution)); // zero
   
    if (coeffs.size() > 0) {
      // determine maximal level
      int jmax(0);
      typedef typename DSBasis<d,dT,BIO>::Index Index;
      for (typename InfiniteVector<double,Index>::const_iterator it(coeffs.begin()),
	     itend(coeffs.end()); it != itend; ++it)
	jmax = std::max(it.index().j()+it.index().e(), jmax);

      // switch to generator representation
      InfiniteVector<double,Index> gcoeffs;
      if (primal)
	reconstruct(coeffs,jmax,gcoeffs);
      else
	reconstruct_t(coeffs,jmax,gcoeffs);

      for (typename InfiniteVector<double,Index>::const_iterator it(gcoeffs.begin()),
	     itend(gcoeffs.end()); it != itend; ++it)
	result.add(*it, evaluate(it.index(), primal, resolution));
    }
    
    return result;
  }

  template <int d, int dT, DSBiorthogonalizationMethod BIO>
  inline
  double
  DSBasis<d,dT,BIO>::evaluate(const unsigned int derivative,
			      const Index& lambda, const double x) const
  {
    return WaveletTL::evaluate(*this, derivative, lambda, x);
  }

  template <int d, int dT, DSBiorthogonalizationMethod BIO>
  inline
  void
  DSBasis<d,dT,BIO>::evaluate
  (const unsigned int derivative,
   const Index& lambda,
   const Array1D<double>& points, Array1D<double>& values) const
  {
    WaveletTL::evaluate(*this, derivative, lambda, points, values);
  }


}
