// implementation for ring_basis.h

#include <set>
#include <numerics/quadrature.h>
#include <numerics/iteratsolv.h>
#include <algebra/sparse_matrix.h>

namespace WaveletTL
{
  template <int d, int dt, int s0, int s1>
  RingBasis<d,dt,s0,s1>::RingBasis(const double r0, const double r1)
    : r0_(r0), r1_(r1), chart_(r0,r1)
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
    // multiply by the normalization factor r^{-1/2}
    for (int i = 0; i <= 1<<resolution; i++) {
      const double r = r0_+i*ldexp(1.0, -resolution)*(r1_-r0_);
      values[1][i] /= sqrt(r);
    }
    
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

    coeffs.clear();
    for (Index lambda = first_generator(j0());;++lambda) {
      double integral = integrate(f, lambda);
      if (fabs(integral) < 1e-15) integral = 0;
      coeffs.set_coefficient(lambda, integral);
      if (lambda == last_wavelet(jmax))
	break;
    }

    if (!primal) {
      // In the following, we setup the Gramian matrix A_Lambda, using all wavelets
      // up to the maximal level jmax. By the tensor product structure of the basis,
      // the entries of A_Lambda can be written as product of 1D integrals:
      //   int_R psi_lambda(x) psi_mu(x) dx
      //   = 2*pi * int_{r0}^{r1} int_0^1 psi_lambda(x) psi_mu(x) dphi r dr
      //   = 2*pi*(r1-r0) * int_0^1 int_0^1 psi_lambda(x) psi_mu(x) dphi (r0+s*(r1-r0)) ds
      // with x=x(s,phi)=r(s)*(cos(2*pi*phi),sin(2*pi*phi)), r(s)=r0+s*(r1-r0).
      // Although A_Lambda is not a Kronecker product of the 1D Gramians (since the
      // order of 2D indices is defined somehow diffently), it does count to
      // setup the 1D Gramians first.

      typedef PeriodicBasis<CDFBasis<d,dt> >             Basis0;
      typedef SplineBasis<d,dt,P_construction,s0,s1,0,0> Basis1;
      typedef typename Basis0::Index Index0;
      typedef typename Basis1::Index Index1;
      typedef typename SparseMatrix<double>::size_type size_type;     

      // setup active wavelet index set
      std::set<Index> Lambda;
      for (Index lambda = first_generator(j0());; ++lambda) {
	Lambda.insert(lambda);
	if (lambda == last_wavelet(jmax)) break;
      }
//       cout << "active index set:" << endl;
//       for (typename std::set<Index>::const_iterator it(Lambda.begin()), itend(Lambda.end());
//   	   it != itend; ++it)
//   	cout << *it << endl;
      
      // setup active wavelet index set in angular direction
      std::set<Index0> Lambda0;
      for (Index0 lambda = basis0.first_generator(j0());; ++lambda) {
	Lambda0.insert(lambda);
	if (lambda == basis0.last_wavelet(jmax)) break;
      }
//       cout << "active index set in angular direction:" << endl;
//       for (typename std::set<Index0>::const_iterator it(Lambda0.begin()), itend(Lambda0.end());
//  	   it != itend; ++it)
//  	cout << *it << endl;

      // setup angular Gramian
      SparseMatrix<double> A_Lambda0(Lambda0.size());
      size_type row = 0;
      for (typename std::set<Index0>::const_iterator it1(Lambda0.begin()), itend(Lambda0.end());
 	   it1 != itend; ++it1, ++row)
 	{
 	  std::list<size_type> indices;
 	  std::list<double> entries;
	  
 	  size_type column = 0;
 	  for (typename std::set<Index0>::const_iterator it2(Lambda0.begin());
 	       it2 != itend; ++it2, ++column)
 	    {
	      double entry = basis0.integrate(*it2, *it1);
	      
 	      if (fabs(entry) > 1e-15) {
 		indices.push_back(column);
 		entries.push_back(2*M_PI*entry);
 	      }
 	    }
 	  A_Lambda0.set_row(row, indices, entries);
 	} 
//       cout << "A_Lambda0=" << endl << A_Lambda0;

      // setup active wavelet index set in radial direction
      std::set<Index1> Lambda1;
      for (Index1 lambda = basis1.first_generator(j0());; ++lambda) {
	Lambda1.insert(lambda);
	if (lambda == basis1.last_wavelet(jmax)) break;
      }
//       cout << "active index set in radial direction:" << endl;
//       for (typename std::set<Index1>::const_iterator it(Lambda1.begin()), itend(Lambda1.end());
//   	   it != itend; ++it)
//   	cout << *it << endl;

      // setup radial Gramian
      SparseMatrix<double> A_Lambda1(Lambda1.size());
      row = 0;
      for (typename std::set<Index1>::const_iterator it1(Lambda1.begin()), itend(Lambda1.end());
 	   it1 != itend; ++it1, ++row)
 	{
 	  std::list<size_type> indices;
 	  std::list<double> entries;
	  
 	  size_type column = 0;
 	  for (typename std::set<Index1>::const_iterator it2(Lambda1.begin());
 	       it2 != itend; ++it2, ++column)
 	    {
	      double entry = integrate_radial(*it2, *it1);
	      
 	      if (fabs(entry) > 1e-15) {
 		indices.push_back(column);
 		entries.push_back(entry);
 	      }
 	    }
 	  A_Lambda1.set_row(row, indices, entries);
 	} 
//       cout << "A_Lambda1=" << endl << A_Lambda1;
      
      // setup global Gramian A_Lambda
      SparseMatrix<double> A_Lambda(Lambda.size());
      row = 0;
      for (typename std::set<Index>::const_iterator it1(Lambda.begin()), itend(Lambda.end());
	   it1 != itend; ++it1, ++row)
	{
	  std::list<size_type> indices;
	  std::list<double> entries;
	  
	  size_type column = 0;
	  for (typename std::set<Index>::const_iterator it2(Lambda.begin());
	       it2 != itend; ++it2, ++column)
	    {
	      double entry
		= A_Lambda0.get_entry(Index0(it2->j(), it2->e()[0], it2->k()[0]).number(),
				      Index0(it1->j(), it1->e()[0], it1->k()[0]).number())
		* A_Lambda1.get_entry(Index1(it2->j(), it2->e()[1], it2->k()[1]).number(),
				      Index1(it1->j(), it1->e()[1], it1->k()[1]).number());
	      
	      if (fabs(entry) > 1e-15) {
		indices.push_back(column);
		entries.push_back(entry);
	      }
	    }
	  A_Lambda.set_row(row, indices, entries);
	} 

      // solve A_Lambda*x = b
      Vector<double> b(A_Lambda.row_dimension());
      row = 0;
      for (Index lambda = first_generator(j0());; ++lambda, ++row) {
 	b[row] = coeffs.get_coefficient(lambda);
	if (lambda == last_wavelet(jmax)) break;
      }
      
//       cout << "b=" << b << endl;
      
      Vector<double> x(b);
      unsigned int iterations;
      CG(A_Lambda, b, x, 1e-15, 500, iterations);
      
      coeffs.clear();
      row = 0;
      for (typename std::set<Index>::const_iterator it(Lambda.begin()), itend(Lambda.end());
	   it != itend; ++it, ++row) {
	if (fabs(x[row]) > 1e-15)
	  coeffs.set_coefficient(*it, x[row]);
      }
    }
  }

  template <int d, int dt, int s0, int s1>
  double
  RingBasis<d,dt,s0,s1>::integrate
  (const Function<2>* f,
   const typename RingBasis<d,dt,s0,s1>::Index& lambda) const
  {
    // We have to compute
    //   int_R f(x) psi_lambda(x) dx
    //   = 2*pi * int_{r0}^{r1} int_0^1 f(x) psi_lambda(x) dphi r dr
    //   = 2*pi*(r1-r0) * int_0^1 int_0^1 f(x) psi_lambda(x) dphi (r0+s*(r1-r0)) ds
    // with x=x(s,phi)=r(s)*(cos(2*pi*phi),sin(2*pi*phi)), r(s)=r0+s*(r1-r0).
    
    double r = 0;

    typedef PeriodicBasis<CDFBasis<d,dt> >             Basis0;
    typedef SplineBasis<d,dt,P_construction,s0,s1,0,0> Basis1;

    // first compute supp(psi_lambda)
    const unsigned int jplus = multi_degree(lambda.e()) > 0 ? 1 : 0;
    int j = lambda.j() + jplus;
    int k1[2], k2[2], length[2];
    basis0.support(typename Basis0::Index(lambda.j(), lambda.e()[0], lambda.k()[0]), k1[0], k2[0]);
    basis1.support(typename Basis1::Index(lambda.j(), lambda.e()[1], lambda.k()[1]), k1[1], k2[1]);
    if (jplus > 0) {
      if (lambda.e()[0] == 0) {
	k1[0] *= 2;
	k2[0] *= 2;	
      } else {
	k1[1] *= 2;
	k2[1] *= 2;
      }
    }
    length[0] = (k2[0] > k1[0] ? k2[0]-k1[0] : k2[0]+(1<<j)-k1[0]); // number of "angular" subintervals
    length[1] = k2[1]-k1[1];

//     cout << "RingBasis::integrate() called with lambda=" << lambda
//    	 << ", j=" << j
//  	 << ", k1[0]=" << k1[0]
//  	 << ", k2[0]=" << k2[0]
//  	 << ", k1[1]=" << k1[1]
//  	 << ", k2[1]=" << k2[1]
//    	 << endl;

    // setup Gauss points and weights for a composite quadrature formula:
    const int N_Gauss = 5;
    const double h = 1.0/(1<<j); // granularity for the quadrature
    FixedArray1D<Array1D<double>,2> gauss_points, gauss_weights, v_values;
    for (int i = 0; i < 2; i++) {
      const unsigned int lengthi = N_Gauss*length[i];
      gauss_points[i].resize(lengthi);
      gauss_weights[i].resize(lengthi);
      v_values[i].resize(lengthi);
    }
    
    // angular direction
    int k = k1[0];
    for (int patch = 0; patch < length[0]; patch++, k = dyadic_modulo(++k,j)) // work on 2^{-j}[k,k+1]
      for (int n = 0; n < N_Gauss; n++) {
 	gauss_points[0][patch*N_Gauss+n] = h*(2*k+1+GaussPoints[N_Gauss-1][n])/2;
	gauss_weights[0][patch*N_Gauss+n] = h*GaussWeights[N_Gauss-1][n];
      }

//     cout << "angular Gauss points: " << gauss_points[0] << endl;
//     cout << "angular Gauss weights: " << gauss_weights[0] << endl;

    // radial direction
    for (int patch = k1[1]; patch < k2[1]; patch++)
      for (int n = 0; n < N_Gauss; n++) {
	gauss_points[1][(patch-k1[1])*N_Gauss+n]
	  = h*(2*patch+1+GaussPoints[N_Gauss-1][n])/2.;
	gauss_weights[1][(patch-k1[1])*N_Gauss+n]
	  = h*GaussWeights[N_Gauss-1][n];
      }
    
//     cout << "radial Gauss points: " << gauss_points[1] << endl;
//     cout << "radial Gauss weights: " << gauss_weights[1] << endl;
    
    // compute the point values of the wavelet (where we use that it is a tensor product)
    v_values[0].resize(gauss_points[0].size());
    v_values[1].resize(gauss_points[1].size());
    basis0.evaluate(0, typename Basis0::Index(lambda.j(), lambda.e()[0], lambda.k()[0]),
		    gauss_points[0], v_values[0]);
    basis1.evaluate(0, typename Basis1::Index(lambda.j(), lambda.e()[1], lambda.k()[1]),
		    gauss_points[1], v_values[1]);

//     cout << "angular point values of psi_lambda: " << v_values[0] << endl;
//     cout << "radial point values of psi_lambda: " << v_values[1] << endl;

    // iterate over all points, evaluate f, and sum up the integral shares
    int index[2]; // current multiindex for the point values
    for (unsigned int i = 0; i < 2; i++)
      index[i] = 0;
    
    Point<2> x;
    Point<2> x_ring;
    
    while (true) {
      // read (r,s)
      for (unsigned int i = 0; i < 2; i++)
 	x[i] = gauss_points[i][index[i]];
      
      // compute x=x(r,s)=r(s)*(cos(2*pi*phi),sin(2*pi*phi))
      const double arg = 2*M_PI*x[0];
      const double rofs = r0_+x[1]*(r1_-r0_);
      x_ring[0] = rofs*cos(arg);
      x_ring[1] = rofs*sin(arg);
      
      double share = 2*M_PI * f->value(x_ring) * rofs; // 2*pi*f(x)*r(s)
      for (unsigned int i = 0; i < 2; i++)
 	share *= gauss_weights[i][index[i]] * v_values[i][index[i]];
      r += share;
      
      // "++index"
      bool exit = false;
      for (unsigned int i = 0; i < 2; i++) {
 	if (index[i] == N_Gauss*length[i]-1) {
 	  index[i] = 0;
 	  exit = (i == 1);
 	} else {
 	  index[i]++;
 	  break;
 	}
      }
      if (exit) break;
    }

    return r;
  }

  template <int d, int dt, int s0, int s1>
  double
  RingBasis<d,dt,s0,s1>::integrate_radial
  (const typename SplineBasis<d,dt,P_construction,s0,s1,0,0>::Index& lambda,
   const typename SplineBasis<d,dt,P_construction,s0,s1,0,0>::Index& mu) const
  {
    double r = 0;
    
    // First we compute the support intersection of \psi_\lambda and \psi_\mu:
    typedef typename SplineBasis<d,dt,P_construction,s0,s1,0,0>::Support Support;
    Support supp;
    if (basis1.intersect_supports(lambda, mu, supp))
      {
  	// Set up Gauss points and weights for a composite quadrature formula:
  	const unsigned int N_Gauss = d+2;
  	const double h = 1.0/(1<<supp.j);
  	Array1D<double> gauss_points (N_Gauss*(supp.k2-supp.k1)), func1values, func2values;
  	for (int patch = supp.k1, id = 0; patch < supp.k2; patch++) // refers to 2^{-j}[patch,patch+1]
  	  for (unsigned int n = 0; n < N_Gauss; n++, id++)
  	    gauss_points[id] = h*(2*patch+1+GaussPoints[N_Gauss-1][n])/2.;
	
  	// - compute point values of the wavelets
   	basis1.evaluate(0, lambda, gauss_points, func1values);
  	basis1.evaluate(0, mu, gauss_points, func2values);
	
  	// - add all integral shares
  	for (int patch = supp.k1, id = 0; patch < supp.k2; patch++)
  	  for (unsigned int n = 0; n < N_Gauss; n++, id++) {
  	    const double gauss_weight = GaussWeights[N_Gauss-1][n] * h;
 	    r += (r1_-r0_) * func1values[id] * func2values[id]
	      * (r0_+gauss_points[id]*(r1_-r0_)) * gauss_weight;
  	  }
      }
    
    return r;
  }
  
}
