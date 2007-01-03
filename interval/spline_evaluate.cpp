// implementation for spline_evaluate.h

namespace WaveletTL
{
  template <int d, int dT, SplineBasisFlavor flavor>
  inline
  SampledMapping<1>
  evaluate(const SplineBasis<d,dT,flavor>& basis,
	   const InfiniteVector<double, typename SplineBasis<d,dT,flavor>::Index>& coeffs,
	   const int resolution)
  {
    return basis.evaluate(coeffs, resolution);
  }


  template <int d, int dT, SplineBasisFlavor flavor>
  inline
  double
  evaluate(const SplineBasis<d,dT,flavor>& basis, const unsigned int derivative,
	   const typename SplineBasis<d,dT,flavor>::Index& lambda,
	   const double x)
  {
    return basis.evaluate(derivative, lambda, x);
  }
  
  template <int d, int dT, SplineBasisFlavor flavor>
  inline
  void
  evaluate(const SplineBasis<d,dT,flavor>& basis, const unsigned int derivative,
	   const typename SplineBasis<d,dT,flavor>::Index& lambda,
	   const Array1D<double>& points, Array1D<double>& values)
  {
    basis.evaluate(derivative, lambda, points, values);
  }
  
  template <int d, int dT, SplineBasisFlavor flavor>
  inline
  void
  evaluate(const SplineBasis<d,dT,flavor>& basis,
	   const typename SplineBasis<d,dT,flavor>::Index& lambda,
	   const Array1D<double>& points, Array1D<double>& funcvalues, Array1D<double>& dervalues)
  {
    basis.evaluate(lambda, points, funcvalues, dervalues);
  }
  
}
