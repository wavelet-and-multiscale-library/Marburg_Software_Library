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

  template <int d, int dt, int s0, int s1>
  inline int
  RingBasis<d,dt,s0,s1>::Deltasize(const int j)
  {
    assert(j >= j0());

    typedef PeriodicBasis<CDFBasis<d,dt> >             Basis0;
    typedef SplineBasis<d,dt,P_construction,s0,s1,0,0> Basis1;

    return Basis0::Deltasize(j) * Basis1::Deltasize(j);
  }
  
  template <int d, int dt, int s0, int s1>
  inline int
  RingBasis<d,dt,s0,s1>::Nabla01size(const int j)
  {
    assert(j >= j0());

    return 1<<(2*j);
  }

  template <int d, int dt, int s0, int s1>
  inline int
  RingBasis<d,dt,s0,s1>::Nabla10size(const int j)
  {
    assert(j >= j0());

    return 1<<(2*j);
  }

  template <int d, int dt, int s0, int s1>
  inline int
  RingBasis<d,dt,s0,s1>::Nabla11size(const int j)
  {
    assert(j >= j0());

    return 1<<(2*j);
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
    // adjust "radial" values by the normalization factor
    const double help = sqrt(2*M_PI*(r1_-r0_));
    for (unsigned int i = 0; i < values[1].size(); i++) {
      values[1][i] /= help*sqrt(r0_+i*ldexp(1.0, -resolution)*(r1_-r0_));
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
      // Note that A_Lambda is not exactly a Kronecker product of the 1D Gramians,
      // since the order of 2D indices is defined somehow diffently.
      // Note also that due to the normalization of the mapped
      // wavelets, it suffices to compute the (1D) integrals
      //   int_0^1 int_0^1 psi^0_lambda(phi,s) psi^0_mu(phi,s) dphi ds,
      // since both the 2*pi*(r1-r0) and the r(s) completely cancel out.
      
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
      
      // setup global Gramian A_Lambda
      SparseMatrix<double> A_Lambda(Lambda.size());
      size_type row = 0;
      for (typename std::set<Index>::const_iterator it1(Lambda.begin()), itend(Lambda.end());
	   it1 != itend; ++it1, ++row)
	{
	  std::list<size_type> indices;
	  std::list<double> entries;
	  
	  size_type column = 0;
	  for (typename std::set<Index>::const_iterator it2(Lambda.begin());
	       it2 != itend; ++it2, ++column)
	    {
	      double entry = integrate(*it2, *it1);
	      
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
    //   = 2*pi*(r1-r0) * int_0^1 int_0^1 f(x(s,phi)) psi_lambda(x(s,phi)) dphi (r0+s*(r1-r0)) ds
    // with x=x(s,phi)=r(s)*(cos(2*pi*phi),sin(2*pi*phi)), r(s)=r0+s*(r1-r0).
    // Note that due to the normalization of the mapped
    // wavelets, it suffices to compute the integrals
    //   sqrt(2*pi*(r1-r0)) int_0^1 int_0^1 f(x(s,phi)) psi^0_lambda(s,phi) dphi sqrt(r0+s*(r1-r0)) ds.

    double r = 0;

    typedef PeriodicBasis<CDFBasis<d,dt> >             Basis0;
    typedef SplineBasis<d,dt,P_construction,s0,s1,0,0> Basis1;

    // first compute supp(psi_lambda)
    const unsigned int jplus = multi_degree(lambda.e()) > 0 ? 1 : 0;
    int j = lambda.j() + jplus;
    int k1[2], k2[2], length[2];
    basis0.support(typename Basis0::Index(lambda.j(), lambda.e()[0], lambda.k()[0]), k1[0], k2[0]);
    basis1.support(typename Basis1::Index(lambda.j(), lambda.e()[1], lambda.k()[1]), k1[1], k2[1]);
    if (jplus > 0) { // in case of a wavelet, adjust the granularity of eventual generator factors
      for (int i = 0; i < 2; i++) {
	if (lambda.e()[i] == 0) {
	  k1[i] *= 2;
	  k2[i] *= 2;	
	}
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
      // read (phi,s)
      for (unsigned int i = 0; i < 2; i++)
 	x[i] = gauss_points[i][index[i]];
      
      // map x into the ring
      chart_.map_point(x, x_ring);
      
      // compute integral share
      double share = f->value(x_ring) * chart_.Gram_factor(x); // sqrt(2*pi*(r1-r0))*f(x)*sqrt(r(s))
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
  inline
  double
  RingBasis<d,dt,s0,s1>::integrate(const Index& lambda,
				   const Index& mu) const
  {
    typedef PeriodicBasis<CDFBasis<d,dt> >             Basis0;
    typedef SplineBasis<d,dt,P_construction,s0,s1,0,0> Basis1;
    typedef typename Basis0::Index Index0;
    typedef typename Basis1::Index Index1;
    
    return
      basis0.integrate(Index0(lambda.j(), lambda.e()[0], lambda.k()[0]),
		       Index0(mu.j(), mu.e()[0], mu.k()[0]))
      * basis1.integrate(Index1(lambda.j(), lambda.e()[1], lambda.k()[1]),
			 Index1(mu.j(), mu.e()[1], mu.k()[1]));
  }
}
