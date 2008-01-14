// implementation for ring_basis.h

namespace WaveletTL
{
  template <int d, int dt, int s0, int s1>
  RingBasis<d,dt,s0,s1>::RingBasis(const double r0, const double r1)
    : r0_(r0), r1_(r1)
  {
  }

  template <int d, int dt, int s0, int s1>
  typename RingBasis<d,dt,s0,s1>::Index
  RingBasis<d,dt,s0,s1>::first_generator(const int j)
  {
    assert(j >= j0());

    typename Index::type_type e;
    typename Index::translation_type k
      (PeriodicBasis<CDFBasis<d,dt> >::DeltaLmin(),
       SplineBasis<d,dt,P_construction,s0,s1,0,0>::DeltaLmin());
    
    return Index(j, e, k);
  }
  
  template <int d, int dt, int s0, int s1>
  typename RingBasis<d,dt,s0,s1>::Index
  RingBasis<d,dt,s0,s1>::last_generator(const int j)
  {
    assert(j >= j0());

    typename Index::type_type e;
    typename Index::translation_type k
      (PeriodicBasis<CDFBasis<d,dt> >::DeltaRmax(j),
       SplineBasis<d,dt,P_construction,s0,s1,0,0>::DeltaRmax(j));
    
    return Index(j, e, k);
  }

  template <int d, int dt, int s0, int s1>
  typename RingBasis<d,dt,s0,s1>::Index
  RingBasis<d,dt,s0,s1>::first_wavelet(const int j)
  {
    assert(j >= j0());
    
    typename Index::type_type e(0, 1);
    typename Index::translation_type k
      (PeriodicBasis<CDFBasis<d,dt> >::DeltaLmin(),
       SplineBasis<d,dt,P_construction,s0,s1,0,0>::Nablamin());

    return Index(j, e, k);
  }

  template <int d, int dt, int s0, int s1>
  typename RingBasis<d,dt,s0,s1>::Index
  RingBasis<d,dt,s0,s1>::last_wavelet(const int j)
  {
    assert(j >= j0());
    
    typename Index::type_type e(1, 1);
    typename Index::translation_type k
      (PeriodicBasis<CDFBasis<d,dt> >::Nablamax(j),
       SplineBasis<d,dt,P_construction,s0,s1,0,0>::Nablamax(j));
    
    return Index(j, e, k);
  }
  
  
  //
  //
  // point evaluation subroutines

  template <int d, int dt, int s0, int s1>
  SampledMapping<2>
  RingBasis<d,dt,s0,s1>::evaluate(const typename RingBasis<d,dt,s0,s1>::Index& lambda,
				  const int resolution) const
  {
    FixedArray1D<Array1D<double>,2> values; // point values of the factors within psi_lambda

    typedef PeriodicBasis<CDFBasis<d,dt> >             Basis0;
    typedef SplineBasis<d,dt,P_construction,s0,s1,0,0> Basis1;

    values[0] = basis0.evaluate(typename Basis0::Index(lambda.j(), lambda.e()[0], lambda.k()[0]),
				resolution).values();
    values[1] = basis1.evaluate(typename Basis1::Index(lambda.j(), lambda.e()[1], lambda.k()[1]),
				resolution).values();

    SampledMapping<2> result(Point<2>(0), Point<2>(1), values);
    return result; // gcc 2.95 does not like these two lines melted into one
  }
    
  template <int d, int dt, int s0, int s1>
  SampledMapping<2>
  RingBasis<d,dt,s0,s1>::evaluate(const InfiniteVector<double, typename RingBasis<d,dt,s0,s1>::Index>& coeffs,
 				  const int resolution) const
  {
    Grid<2> grid(Point<2>(0), Point<2>(1), 1<<resolution);
    SampledMapping<2> result(grid); // zero
    
    typedef typename RingBasis<d,dt,s0,s1>::Index Index;
    for (typename InfiniteVector<double,Index>::const_iterator it(coeffs.begin()),
  	   itend(coeffs.end()); it != itend; ++it)
      result.add(*it, evaluate(it.index(), resolution));
    
    return result;
  }

  //
  //
  // expansion subroutines

  template <int d, int dt, int s0, int s1>
  void
  RingBasis<d,dt,s0,s1>::expand(const Function<2>* f,
				const bool primal,
				const int jmax,
				InfiniteVector<double, Index>& coeffs) const
  {
    typedef typename RingBasis<d,dt,s0,s1>::Index Index;

    assert(primal); // for the moment

    coeffs.clear();
    for (Index lambda = first_generator(j0());;++lambda) {
      coeffs.set_coefficient(lambda, integrate(f, lambda));
      if (lambda == last_wavelet(jmax))
	break;
    }
    
  }

  template <int d, int dt, int s0, int s1>
  double
  RingBasis<d,dt,s0,s1>::integrate
  (const Function<2>* f,
   const typename RingBasis<d,dt,s0,s1>::Index& lambda) const
  {
    // We have to compute
    //   int_R f(x) psi_lambda(x) dx = int_{r=r0}^r1 int_0^{2*pi} f(x) psi_lambda(x) dphi r dr
    // with x=x(r,phi)=(r*cos(phi), r*sin(phi)).

    double r = 0;

    typedef PeriodicBasis<CDFBasis<d,dt> >             Basis0;
    typedef SplineBasis<d,dt,P_construction,s0,s1,0,0> Basis1;

    // first compute supp(psi_lambda)
    const unsigned int jplus = multi_degree(lambda.e()) > 0 ? 1 : 0;
    int j = lambda.j() + jplus;
    int k1[2], k2[2], length0;
    basis0.support(typename Basis0::Index(lambda.j(), lambda.e()[0], lambda.k()[0]), k1[0], k2[0]);
    if (lambda.e()[0] == 0 && jplus > 0) {
      k1[0] *= 2;
      k2[0] *= 2;
    }
    length0 = (k2[0] > k1[0] ? k2[0]-k1[0] : k2[0]+(1<<j)-k1[0]); // number of subintervals
    basis1.support(typename Basis1::Index(lambda.j(), lambda.e()[1], lambda.k()[1]), k1[1], k2[1]);
    if (lambda.e()[1] == 0 && jplus > 0) {
      k1[1] *= 2;
      k2[1] *= 2;
    }

    cout << "RingBasis::integrate() called with lambda=" << lambda
	 << ", j=" << j
	 << ", k1[0]=" << k1[0]
	 << ", k2[0]=" << k2[0]
	 << ", k1[1]=" << k1[1]
	 << ", k2[1]=" << k2[1]
	 << endl;

    r = 42;

    return r;
  }



//   template <class IBASIS, unsigned int DIM>
//   double
//   CubeBasis<IBASIS,DIM>::integrate
//   (const Function<DIM>* f,
//    const Index& lambda) const
//   {
//     // f(v) = \int_0^1 g(t)v(t) dt
    
//     double r = 0;
    
//     // first compute supp(psi_lambda)
//     Support supp;
//     support(lambda, supp);
    
//     // setup Gauss points and weights for a composite quadrature formula:
//     const int N_Gauss = 5;
//     const double h = ldexp(1.0, -supp.j); // granularity for the quadrature
//     FixedArray1D<Array1D<double>,DIM> gauss_points, gauss_weights, v_values;
//     for (unsigned int i = 0; i < DIM; i++) {
//       gauss_points[i].resize(N_Gauss*(supp.b[i]-supp.a[i]));
//       gauss_weights[i].resize(N_Gauss*(supp.b[i]-supp.a[i]));
//       for (int patch = supp.a[i]; patch < supp.b[i]; patch++)
// 	for (int n = 0; n < N_Gauss; n++) {
// 	  gauss_points[i][(patch-supp.a[i])*N_Gauss+n]
// 	    = h*(2*patch+1+GaussPoints[N_Gauss-1][n])/2.;
// 	  gauss_weights[i][(patch-supp.a[i])*N_Gauss+n]
// 	    = h*GaussWeights[N_Gauss-1][n];
// 	}
//     }

//     // compute the point values of the integrand (where we use that it is a tensor product)
//     for (unsigned int i = 0; i < DIM; i++)
//       bases()[i]->evaluate(0,
// 			   typename IBASIS::Index(lambda.j(),
// 						  lambda.e()[i],
// 						  lambda.k()[i],
// 						  bases()[i]),
// 			   gauss_points[i], v_values[i]);
    
//     // iterate over all points and sum up the integral shares
//     int index[DIM]; // current multiindex for the point values
//     for (unsigned int i = 0; i < DIM; i++)
//       index[i] = 0;
    
//     Point<DIM> x;
//     while (true) {
//       for (unsigned int i = 0; i < DIM; i++)
// 	x[i] = gauss_points[i][index[i]];
//       double share = f->value(x);
//       for (unsigned int i = 0; i < DIM; i++)
// 	share *= gauss_weights[i][index[i]] * v_values[i][index[i]];
//       r += share;

//       // "++index"
//       bool exit = false;
//       for (unsigned int i = 0; i < DIM; i++) {
// 	if (index[i] == N_Gauss*(supp.b[i]-supp.a[i])-1) {
// 	  index[i] = 0;
// 	  exit = (i == DIM-1);
// 	} else {
// 	  index[i]++;
// 	  break;
// 	}
//       }
//       if (exit) break;
//     }
    
//     return r;
//   }



}
