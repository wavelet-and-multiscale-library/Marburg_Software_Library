// implementation for ldomain_evaluate.h

namespace WaveletTL
{
  template <class IBASIS>
  Array1D<SampledMapping<2> >
  evaluate(const LDomainBasis<IBASIS>& basis,
	   const typename LDomainBasis<IBASIS>::Index& lambda,
	   const bool primal,
	   const int resolution)
  {
    Array1D<SampledMapping<2> > r(3);

    Grid<2> grid0(Point<2>(0), Point<2>(1), 1<<resolution);
    FixedArray1D<Array1D<double>,2> values;
    values[0].resize((1<<resolution)+1);
    values[1].resize((1<<resolution)+1);
    
    r[0] = SampledMapping<2>(Point<2>(-1, 0), Point<2>(0,1), values);
    r[1] = SampledMapping<2>(Point<2>(-1,-1), Point<2>(0,0), values);
    r[2] = SampledMapping<2>(Point<2>( 0,-1), Point<2>(1,0), values);

    return r;
  }
    
//   template <class IBASIS>
//   Array1D<SampledMapping<2> >
//   evaluate(const LDomainBasis<IBASIS>& basis,
// 	   const InfiniteVector<double, typename LDomainBasis<IBASIS>::Index>& coeffs,
// 	   const bool primal,
// 	   const int resolution)
//   {
//     Grid<2> grid(Point<2>(0), Point<2>(1), 1<<resolution);
//     SampledMapping<2> result(grid); // zero
      
//     typedef typename LDomainBasis<IBASIS,DIM>::Index Index;
//     for (typename InfiniteVector<double,Index>::const_iterator it(coeffs.begin()),
// 	   itend(coeffs.end()); it != itend; ++it)
//       result.add(*it, evaluate(basis, it.index(), primal, resolution));
      
//     return result;
//   }

}
