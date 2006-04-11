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

    typedef typename LDomainBasis<IBASIS>::Index Index;

    typename Index::type_type zero;
    if (lambda.e() == zero) {
      // psi_lambda is a generator. Here we know that per patch,
      // we can resort to a single tensor product.
      
      FixedArray1D<Array1D<double>,2> values;
      values[0].resize((1<<resolution)+1); // values in x-direction
      values[1].resize((1<<resolution)+1); // values in y-direction

      // patch 0
      switch (lambda.p()) {
      case 0:
	// psi_lambda completely lives on patch 0
	values[0] = evaluate(basis.basis00(),
			     typename IBASIS::Index(lambda.j(),
						    0,
						    lambda.k()[0],
						    &basis.basis00()),
			     primal,
			     resolution).values();
	values[1] = evaluate(basis.basis10(),
			     typename IBASIS::Index(lambda.j(),
						    0,
						    lambda.k()[1],
						    &basis.basis10()),
			     primal,
			     resolution).values();
	r[0] = SampledMapping<2>(Point<2>(-1, 0), Point<2>(0,1), values);
	
	for (int i = 0; i <= 1<<resolution; i++) {
	  values[0][i] = values[1][i] = 0;
	}
	r[1] = SampledMapping<2>(Point<2>(-1,-1), Point<2>(0,0), values);
	r[2] = SampledMapping<2>(Point<2>( 0,-1), Point<2>(1,0), values);
	break;
      case 1:
	// psi_lambda completely lives on patch 1
	values[0] = evaluate(basis.basis01(),
			     typename IBASIS::Index(lambda.j(),
						    0,
						    lambda.k()[0],
						    &basis.basis01()),
			     primal,
			     resolution).values();
	values[1] = evaluate(basis.basis01(),
			     typename IBASIS::Index(lambda.j(),
						    0,
						    lambda.k()[1],
						    &basis.basis01()),
			     primal,
			     resolution).values();
	r[1] = SampledMapping<2>(Point<2>(-1,-1), Point<2>(0,0), values);

	for (int i = 0; i <= 1<<resolution; i++) {
	  values[0][i] = values[1][i] = 0;
	}
	r[0] = SampledMapping<2>(Point<2>(-1, 0), Point<2>(0,1), values);
	r[2] = SampledMapping<2>(Point<2>( 0,-1), Point<2>(1,0), values);
	break;
      case 2:
	// psi_lambda completely lives on patch 2
	values[0] = evaluate(basis.basis10(),
			     typename IBASIS::Index(lambda.j(),
						    0,
						    lambda.k()[0],
						    &basis.basis10()),
			     primal,
			     resolution).values();
	values[1] = evaluate(basis.basis00(),
			     typename IBASIS::Index(lambda.j(),
						    0,
						    lambda.k()[1],
						    &basis.basis00()),
			     primal,
			     resolution).values();
	r[2] = SampledMapping<2>(Point<2>( 0,-1), Point<2>(1,0), values);

	for (int i = 0; i <= 1<<resolution; i++) {
	  values[0][i] = values[1][i] = 0;
	}
	r[0] = SampledMapping<2>(Point<2>(-1, 0), Point<2>(0,1), values);
	r[1] = SampledMapping<2>(Point<2>(-1,-1), Point<2>(0,0), values);
	break;
      case 3:
	// psi_lambda lives on patches 0 and 1
	values[0] = evaluate(basis.basis00(),
			     typename IBASIS::Index(lambda.j(),
						    0,
						    lambda.k()[0],
						    &basis.basis00()),
			     primal,
			     resolution).values();
	values[1] = evaluate(basis.basis10(),
			     typename IBASIS::Index(lambda.j(),
						    0,
						    basis.basis10().DeltaLmin(),
						    &basis.basis10()),
			     primal,
			     resolution).values();
	for (int i = 0; i <= 1<<resolution; i++) values[1][i] *= M_SQRT1_2;
	r[0] = SampledMapping<2>(Point<2>(-1, 0), Point<2>(0,1), values);

	values[0] = evaluate(basis.basis01(),
			     typename IBASIS::Index(lambda.j(),
						    0,
						    lambda.k()[0],
						    &basis.basis01()),
			     primal,
			     resolution).values();
	values[1] = evaluate(basis.basis01(),
			     typename IBASIS::Index(lambda.j(),
						    0,
						    basis.basis01().DeltaRmax(lambda.j()),
						    &basis.basis01()),
			     primal,
			     resolution).values();
	for (int i = 0; i <= 1<<resolution; i++) values[1][i] *= M_SQRT1_2;
	r[1] = SampledMapping<2>(Point<2>(-1,-1), Point<2>(0,0), values);

	for (int i = 0; i <= 1<<resolution; i++) {
	  values[0][i] = values[1][i] = 0;
	}
	r[2] = SampledMapping<2>(Point<2>( 0,-1), Point<2>(1,0), values);
	break;
      case 4:
	// psi_lambda lives on patches 1 and 2
	values[0] = evaluate(basis.basis01(),
			     typename IBASIS::Index(lambda.j(),
						    0,
						    basis.basis01().DeltaRmax(lambda.j()),
						    &basis.basis01()),
			     primal,
			     resolution).values();
	values[1] = evaluate(basis.basis01(),
			     typename IBASIS::Index(lambda.j(),
						    0,
						    lambda.k()[1],
						    &basis.basis01()),
			     primal,
			     resolution).values();
	for (int i = 0; i <= 1<<resolution; i++) values[1][i] *= M_SQRT1_2;
	r[1] = SampledMapping<2>(Point<2>(-1,-1), Point<2>(0,0), values);

	values[0] = evaluate(basis.basis10(),
			     typename IBASIS::Index(lambda.j(),
						    0,
						    basis.basis10().DeltaLmin(),
						    &basis.basis10()),
			     primal,
			     resolution).values();
	values[1] = evaluate(basis.basis00(),
			     typename IBASIS::Index(lambda.j(),
						    0,
						    lambda.k()[1],
						    &basis.basis00()),
			     primal,
			     resolution).values();
	for (int i = 0; i <= 1<<resolution; i++) values[1][i] *= M_SQRT1_2;
	r[2] = SampledMapping<2>(Point<2>( 0,-1), Point<2>(1,0), values);

	for (int i = 0; i <= 1<<resolution; i++) {
	  values[0][i] = values[1][i] = 0;
	}
	r[0] = SampledMapping<2>(Point<2>(-1, 0), Point<2>(0,1), values);
	break;
      }
    } else {
      // psi_lambda is a wavelet. We leave this delicate case to the wavelet transform routines...
      InfiniteVector<double, Index> gcoeffs;
      if (primal)
 	basis.reconstruct_1(lambda, lambda.j()+1, gcoeffs);
      else
 	basis.reconstruct_t_1(lambda, lambda.j()+1, gcoeffs);
      
      return evaluate(basis, gcoeffs, primal, resolution);
    }
    
    return r;
  }
    
  template <class IBASIS>
  Array1D<SampledMapping<2> >
  evaluate(const LDomainBasis<IBASIS>& basis,
 	   const InfiniteVector<double, typename LDomainBasis<IBASIS>::Index>& coeffs,
 	   const bool primal,
 	   const int resolution)
  {
    // first prepare an empty plot
    FixedArray1D<Array1D<double>,2> values;
    values[0].resize((1<<resolution)+1);
    values[1].resize((1<<resolution)+1);
    for (int i = 0; i <= 1<<resolution; i++) {
      values[0][i] = values[1][i] = 0;
    }
    Array1D<SampledMapping<2> > result(3);
    result[0] = SampledMapping<2>(Point<2>(-1, 0), Point<2>(0,1), values);
    result[1] = SampledMapping<2>(Point<2>(-1,-1), Point<2>(0,0), values);
    result[2] = SampledMapping<2>(Point<2>( 0,-1), Point<2>(1,0), values);
          
    // add all plots of the single functions
    typedef typename LDomainBasis<IBASIS>::Index Index;
    for (typename InfiniteVector<double,Index>::const_iterator it(coeffs.begin()),
 	   itend(coeffs.end()); it != itend; ++it) {
      Array1D<SampledMapping<2> > temp(evaluate(basis, it.index(), primal, resolution));
      result[0].add(*it, temp[0]);
      result[1].add(*it, temp[1]);
      result[2].add(*it, temp[2]);
    }
    
    return result;
  }

}
