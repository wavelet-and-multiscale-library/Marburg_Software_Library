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
    // universal routine for all values of "flavor"
    return basis.evaluate(coeffs, resolution);
  }

  template <int d, int dT, SplineBasisFlavor flavor>
  inline
  SampledMapping<1>
  evaluate(const SplineBasis<d,dT,flavor>& basis,
	   const typename SplineBasis<d,dT,flavor>::Index& lambda,
	   const int resolution)
  {
    // universal routine for all values of "flavor"
    return basis.evaluate(lambda, resolution);
  }

  template <int d, int dT, SplineBasisFlavor flavor>
  inline
  double
  evaluate(const SplineBasis<d,dT,flavor>& basis, const unsigned int derivative,
	   const typename SplineBasis<d,dT,flavor>::Index& lambda,
	   const double x)
  {
    // universal routine for all values of "flavor"
    return basis.evaluate(derivative, lambda, x);
  }
  
  template <int d, int dT, SplineBasisFlavor flavor>
  inline
  void
  evaluate(const SplineBasis<d,dT,flavor>& basis, const unsigned int derivative,
	   const typename SplineBasis<d,dT,flavor>::Index& lambda,
	   const Array1D<double>& points, Array1D<double>& values)
  {
    // universal routine for all values of "flavor"
    basis.evaluate(derivative, lambda, points, values);
  }
  
  template <int d, int dT, SplineBasisFlavor flavor>
  inline
  void
  evaluate(const SplineBasis<d,dT,flavor>& basis,
	   const typename SplineBasis<d,dT,flavor>::Index& lambda,
	   const Array1D<double>& points, Array1D<double>& funcvalues, Array1D<double>& dervalues)
  {
    // universal routine for all values of "flavor"
    basis.evaluate(lambda, points, funcvalues, dervalues);
  }
  
}
