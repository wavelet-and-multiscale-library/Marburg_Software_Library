// implementation for spline_evaluate.h

#include <interval/interval_bspline.h>

namespace WaveletTL
{
  template <int d, int dT>
  SampledMapping<1>
  evaluate(const SplineBasis<d,dT>& basis,
	   const InfiniteVector<double, typename SplineBasis<d,dT>::Index>& coeffs,
	   const int resolution)
  {
    Grid<1> grid(0, 1, 1<<resolution);
    SampledMapping<1> result(grid); // zero
    if (coeffs.size() > 0) {
      // determine maximal level
      int jmax(0);
      typedef typename SplineBasis<d,dT>::Index Index;
      for (typename InfiniteVector<double,Index>::const_iterator it(coeffs.begin()),
	     itend(coeffs.end()); it != itend; ++it)
	jmax = std::max(it.index().j()+it.index().e(), jmax);

      // insert coefficients into a dense vector
      Vector<double> wcoeffs(basis.Deltasize(jmax));
      for (typename InfiniteVector<double,Index>::const_iterator it(coeffs.begin()),
	     itend(coeffs.end()); it != itend; ++it) {
	// determine number of the wavelet
	typedef typename Vector<double>::size_type size_type;
	size_type number = 0;
	if (it.index().e() == 0) {
	  number = it.index().k()-basis.DeltaLmin();
	} else {
	  number = basis.Deltasize(it.index().j())+it.index().k()-basis.Nablamin();
	}
	wcoeffs[number] = *it;
      }
      
      // switch to generator representation
      Vector<double> gcoeffs(wcoeffs.size(), false);
      if (jmax == basis.j0())
	gcoeffs = wcoeffs;
      else
	basis.apply_Tj(jmax-1, wcoeffs, gcoeffs);
      
      Array1D<double> values((1<<resolution)+1);
      for (unsigned int i(0); i < values.size(); i++) {
	values[i] = 0;
	const double x = i*ldexp(1.0, -resolution);
	SchoenbergIntervalBSpline_td<d> sbs(jmax,0);
	for (unsigned int k = 0; k < gcoeffs.size(); k++) {
	  sbs.set_k(basis.DeltaLmin()+k);
	  values[i] += gcoeffs[k] * sbs.value(Point<1>(x));
	}
      }
      
      return SampledMapping<1>(grid, values);
    }
    
    return result;
  }


}
