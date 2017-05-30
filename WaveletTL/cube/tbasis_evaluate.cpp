// implementation for tbasis_evaluate.h

// using namespace MathTL;

namespace WaveletTL
{
  template <class IBASIS, unsigned int DIM>
  SampledMapping<DIM>
  evaluate(const TensorBasis<IBASIS,DIM>& basis,
	   const typename TensorBasis<IBASIS,DIM>::Index& lambda,
	   const bool primal,
	   const int resolution)
  {
    FixedArray1D<Array1D<double>,DIM> values; // point values of the factors within psi_lambda
    for (unsigned int i = 0; i < DIM; i++)
    {
        values[i] = evaluate(*(basis.bases()[i]),
			   typename IBASIS::Index(lambda.j()[i],
						  lambda.e()[i],
						  lambda.k()[i],
						  basis.bases()[i]),
			   primal,
			   resolution).values();
    }
    SampledMapping<DIM> result(Point<DIM>(0), Point<DIM>(1), values);
    return result; // gcc 2.95 does not like these two lines melted into one
  }
  
  template <class IBASIS, unsigned int DIM>
  SampledMapping<DIM>
  evaluate(const TensorBasis<IBASIS,DIM>& basis,
	   const int& lambdanum,
	   const bool primal,
	   const int resolution)
  {
    typedef typename TensorBasis<IBASIS,DIM>::Index Index;
    Index lambda(basis.get_wavelet(lambdanum));
    FixedArray1D<Array1D<double>,DIM> values; // point values of the factors within psi_lambda
    for (unsigned int i = 0; i < DIM; i++)
    {
        values[i] = evaluate(*(basis.bases()[i]),
			   typename IBASIS::Index(lambda.j()[i],
						  lambda.e()[i],
						  lambda.k()[i],
						  basis.bases()[i]),
			   primal,
			   resolution).values();
    }
    SampledMapping<DIM> result(Point<DIM>(0), Point<DIM>(1), values);
    return result; // gcc 2.95 does not like these two lines melted into one
  }


  template <class IBASIS>
  SampledMapping<1>
  evaluate(const TensorBasis<IBASIS,1>& basis,
	   const InfiniteVector<double, typename TensorBasis<IBASIS,1>::Index>& coeffs,
	   const bool primal,
	   const int resolution)
  {
    Grid<1> grid(0.0, 1.0, 1<<resolution);
    SampledMapping<1> result(grid); // zero

    typedef typename TensorBasis<IBASIS,1>::Index Index;
    for (typename InfiniteVector<double,Index>::const_iterator it(coeffs.begin()),
	   itend(coeffs.end()); it != itend; ++it)
           {
               result.add(*it, evaluate(basis, it.index(), primal, resolution));
           }


    return result;
  }
  
  template <class IBASIS>
  SampledMapping<1>
  evaluate(const TensorBasis<IBASIS,1>& basis,
	   const InfiniteVector<double, int>& coeffs,
	   const bool primal,
	   const int resolution)
  {
    Grid<1> grid(0.0, 1.0, 1<<resolution);
    SampledMapping<1> result(grid); // zero

    typedef typename TensorBasis<IBASIS,1>::Index Index;
    for (typename InfiniteVector<double,Index>::const_iterator it(coeffs.begin()),
	   itend(coeffs.end()); it != itend; ++it)
           {
               result.add(*it, evaluate(basis, it.index(), primal, resolution));
           }


    return result;
  }
  

  template <class IBASIS, unsigned int DIM>
  SampledMapping<DIM>
  evaluate(const TensorBasis<IBASIS,DIM>& basis,
	   const InfiniteVector<double, typename TensorBasis<IBASIS,DIM>::Index>& coeffs,
	   const bool primal,
	   const int resolution)
  {
    Grid<DIM> grid(Point<DIM>(0), Point<DIM>(1), 1<<resolution);
    SampledMapping<DIM> result(grid); // zero

    typedef typename TensorBasis<IBASIS,DIM>::Index Index;
    for (typename InfiniteVector<double,Index>::const_iterator it(coeffs.begin()),
	   itend(coeffs.end()); it != itend; ++it)
           {
               result.add(*it, evaluate(basis, it.index(), primal, resolution));
           }
      

    return result;
  }
  
  template <class IBASIS, unsigned int DIM>
  SampledMapping<DIM>
  evaluate(const TensorBasis<IBASIS,DIM>& basis,
	   const InfiniteVector<double, int>& coeffs,
	   const bool primal,
	   const int resolution)
  {
    Grid<DIM> grid(Point<DIM>(0), Point<DIM>(1), 1<<resolution);
    SampledMapping<DIM> result(grid); // zero

//    typedef typename TensorBasis<IBASIS,DIM>::Index Index;
    for (typename InfiniteVector<double,int>::const_iterator it(coeffs.begin()),
	   itend(coeffs.end()); it != itend; ++it)
           {
               result.add(*it, evaluate(basis, it.index(), primal, resolution));
           }
      

    return result;
  }
}

