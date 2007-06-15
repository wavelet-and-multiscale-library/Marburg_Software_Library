// implementation for interval_evaluate.h

#include <geometry/grid.h>
#include <utils/array1d.h>

using MathTL::Grid;

namespace WaveletTL
{
  template <class IBasis>
  SampledMapping<1>
  evaluate(const IBasis& basis,
           const typename IBasis::Index& lambda,
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

  template <class IBasis>
  SampledMapping<1>
  evaluate(const IBasis& basis,
           const typename IBasis::Index& lambda,
           const bool primal,
           const int resolution)
  {
    assert(primal);
    return evaluate(basis, lambda, resolution);
  }
  
  template <class IBasis>
  SampledMapping<1>
  evaluate(const IBasis& basis,
           const InfiniteVector<double, typename IBasis::Index>& coeffs,
           const int resolution)
  {
    SampledMapping<1> result(Grid<1>(0, 1, 1<<resolution)); // zero
   
    if (coeffs.size() > 0) {
      // determine maximal level
      int jmax(0);
      typedef typename IBasis::Index Index;
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

  template <class IBasis>
  inline
  double evaluate(const IBasis& basis, const unsigned int derivative,
                  const typename IBasis::Index& lambda,
                  const double x)
  {
    return basis.primal_evaluate(derivative, lambda, x);
  }
 
  template <class IBasis>
  void
  evaluate(const IBasis& basis, const unsigned int derivative,
           const typename IBasis::Index& lambda,
           const Array1D<double>& points, Array1D<double>& values)
  {
#if 1 // new code, calls method of basis class
    basis.primal_evaluate(derivative, lambda, points, values);
#else // old code, brute force solution, loops through points and calls evaluate method for each point
    const unsigned int npoints(points.size()); // number of points
    unsigned int i;

    for (i = 0; i < npoints; i++)
      values[i] = basis.primal_evaluate(derivative, lambda, points[i]);
#endif
  }
  
  template <class IBasis>
  inline
  void evaluate(const IBasis& basis,
                const typename IBasis::Index& lambda,
                const Array1D<double>& points, Array1D<double>& funcvalues, Array1D<double>& dervalues)
  {
    evaluate(basis, 0, lambda, points, funcvalues);
    evaluate(basis, 1, lambda, points, dervalues);
  }
}
