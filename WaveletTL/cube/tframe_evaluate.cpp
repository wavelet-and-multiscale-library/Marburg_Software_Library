// implementation for tframe_evaluate.h

// using namespace MathTL;

namespace WaveletTL
{
  template <class IFRAME, unsigned int DIM>
  SampledMapping<DIM>
  evaluate(const TensorFrame<IFRAME,DIM>& frame,
	   const typename TensorFrame<IFRAME,DIM>::Index& lambda,
	   const bool primal,
	   const int resolution)
  {
    FixedArray1D<Array1D<double>,DIM> values; // point values of the factors within psi_lambda
    for (unsigned int i = 0; i < DIM; i++)
    {
        values[i] = evaluate(*(frame.frames()[i]),
			   typename IFRAME::Index(lambda.p()[i],
                                                  lambda.j()[i],
						  lambda.e()[i],
						  lambda.k()[i],
						  frame.frames()[i]),
			   primal,
			   resolution).values();
    }
    SampledMapping<DIM> result(Point<DIM>(0), Point<DIM>(1), values);
    return result; // gcc 2.95 does not like these two lines melted into one
  }
  
  template <class IFRAME, unsigned int DIM>
  SampledMapping<DIM>
  evaluate(const TensorFrame<IFRAME,DIM>& frame,
	   const int& lambdanum,
	   const bool primal,
	   const int resolution)
  {
    typedef typename TensorFrame<IFRAME,DIM>::Index Index;
    Index lambda(frame.get_quarklet(lambdanum));
    FixedArray1D<Array1D<double>,DIM> values; // point values of the factors within psi_lambda
    for (unsigned int i = 0; i < DIM; i++)
    {
        values[i] = evaluate(*(frame.frames()[i]),
			   typename IFRAME::Index(lambda.p()[i],
                                                  lambda.j()[i],
						  lambda.e()[i],
						  lambda.k()[i],
						  frame.frames()[i]),
			   primal,
			   resolution).values();
    }
    SampledMapping<DIM> result(Point<DIM>(0), Point<DIM>(1), values);
    return result; // gcc 2.95 does not like these two lines melted into one
  }


  template <class IFRAME>
  SampledMapping<1>
  evaluate(const TensorFrame<IFRAME,1>& frame,
	   const InfiniteVector<double, typename TensorFrame<IFRAME,1>::Index>& coeffs,
	   const bool primal,
	   const int resolution)
  {
    Grid<1> grid(0.0, 1.0, 1<<resolution);
    SampledMapping<1> result(grid); // zero

    typedef typename TensorFrame<IFRAME,1>::Index Index;
    for (typename InfiniteVector<double,Index>::const_iterator it(coeffs.begin()),
	   itend(coeffs.end()); it != itend; ++it)
           {
               result.add(*it, evaluate(frame, it.index(), primal, resolution));
           }


    return result;
  }
  
  template <class IFRAME>
  SampledMapping<1>
  evaluate(const TensorFrame<IFRAME,1>& frame,
	   const InfiniteVector<double, int>& coeffs,
	   const bool primal,
	   const int resolution)
  {
    Grid<1> grid(0.0, 1.0, 1<<resolution);
    SampledMapping<1> result(grid); // zero

    typedef typename TensorFrame<IFRAME,1>::Index Index;
    for (typename InfiniteVector<double,Index>::const_iterator it(coeffs.begin()),
	   itend(coeffs.end()); it != itend; ++it)
           {
               result.add(*it, evaluate(frame, it.index(), primal, resolution));
           }


    return result;
  }
  

  template <class IFRAME, unsigned int DIM>
  SampledMapping<DIM>
  evaluate(const TensorFrame<IFRAME,DIM>& frame,
	   const InfiniteVector<double, typename TensorFrame<IFRAME,DIM>::Index>& coeffs,
	   const bool primal,
	   const int resolution)
  {
    Grid<DIM> grid(Point<DIM>(0), Point<DIM>(1), 1<<resolution);
    SampledMapping<DIM> result(grid); // zero

    typedef typename TensorFrame<IFRAME,DIM>::Index Index;
    for (typename InfiniteVector<double,Index>::const_iterator it(coeffs.begin()),
	   itend(coeffs.end()); it != itend; ++it)
           {
               result.add(*it, evaluate(frame, it.index(), primal, resolution));
           }
      

    return result;
  }
  
  template <class IFRAME, unsigned int DIM>
  SampledMapping<DIM>
  evaluate(const TensorFrame<IFRAME,DIM>& frame,
	   const InfiniteVector<double, int>& coeffs,
	   const bool primal,
	   const int resolution)
  {
    Grid<DIM> grid(Point<DIM>(0), Point<DIM>(1), 1<<resolution);
    SampledMapping<DIM> result(grid); // zero

//    typedef typename TensorFrame<IFRAME,DIM>::Index Index;
    for (typename InfiniteVector<double,int>::const_iterator it(coeffs.begin()),
	   itend(coeffs.end()); it != itend; ++it)
           {
               result.add(*it, evaluate(frame, it.index(), primal, resolution));
           }
      

    return result;
  }
}

