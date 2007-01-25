// implementation for a_evaluate.h

#include <geometry/grid.h>
#include <utils/array1d.h>

using MathTL::Grid;

namespace WaveletTL
{
  template <int n>
  SampledMapping<1>
  evaluate(const ABasis<n>& basis,
           const typename ABasis<n>::Index& lambda,
           const int resolution)
  {
    if (lambda.e() == 0) {
      // generator
      Grid<1> grid(0, 1, 1<<resolution);
      MathTL::Array1D<double> values((1<<resolution)+1);
      for (unsigned int i(0); i < values.size(); i++)
        values[i] = evaluate(basis, 0, lambda, i*ldexp(1.0,-resolution));
      return SampledMapping<1>(grid, values);
    }
    else {
      // wavelet
      Grid<1> grid(0, 1, 1<<resolution);
      MathTL::Array1D<double> values((1<<resolution)+1);
      for (unsigned int i(0); i < values.size(); i++)
        values[i] = evaluate(basis, 0, lambda, i*ldexp(1.0,-resolution));
      return SampledMapping<1>(grid, values);
    }
    
    return SampledMapping<1>(); // dummy return for the compiler
  }
  
  template <int n>
  SampledMapping<1>
  evaluate(const ABasis<n>& basis,
           const InfiniteVector<double, typename ABasis<n>::Index>& coeffs,
           const int resolution)
  {
    SampledMapping<1> result(Grid<1>(0, 1, 1<<resolution)); // zero
   
    if (coeffs.size() > 0) {
      // determine maximal level
      int jmax(0);
      typedef typename ABasis<n>::Index Index;
      for (typename InfiniteVector<double,Index>::const_iterator it(coeffs.begin()),
             itend(coeffs.end()); it != itend; ++it)
        jmax = std::max(it.index().j()+it.index().e(), jmax);

      // switch to generator representation
      InfiniteVector<double,Index> gcoeffs;
      basis.reconstruct(coeffs,jmax,gcoeffs);
      for (typename InfiniteVector<double,Index>::const_iterator it(gcoeffs.begin()),
             itend(gcoeffs.end()); it != itend; ++it)
        result.add(*it, evaluate(basis, it.index(), resolution));
    }
    
    return result;
  }

  template <int n>
  inline
  double evaluate(const ABasis<n>& basis, const unsigned int derivative,
                  const typename ABasis<n>::Index& lambda,
                  const double x)
  {
//    assert(derivative <= 1); // currently we only support derivatives up to the first order

    Piecewise<double> fkt;
    unsigned int i;

    fkt = basis.get_function(lambda);
    for (i = 0; i < derivative; i++)
      fkt.differentiate();

    return fkt.value(x);
  }
 
  template <int n>
  void
  evaluate(const ABasis<n>& basis, const unsigned int derivative,
           const typename ABasis<n>::Index& lambda,
           const Array1D<double>& points, Array1D<double>& values)
  {
//    assert(derivative <= 1); // we only support derivatives up to the first order

    Piecewise<double> fkt;
    const unsigned int npoints(points.size()); // number of points
    unsigned int i;

    fkt = basis.get_function(lambda);
    for (i = 0; i < derivative; i++)
      fkt.differentiate();

    for (i = 0; i < npoints; i++)
      values[i] = fkt.value(points[i]);
  }
  
  template <int n>
  inline
  void evaluate(const ABasis<n>& basis,
                const typename ABasis<n>::Index& lambda,
                const Array1D<double>& points, Array1D<double>& funcvalues, Array1D<double>& dervalues)
  {
    evaluate(basis, 0, lambda, points, funcvalues);
    evaluate(basis, 1, lambda, points, dervalues);
  }

}
