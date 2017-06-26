// implementation of PeriodicFrame2 methods

#include <cmath>
#include <algebra/sparse_matrix.h>
#include <geometry/point.h>
#include <geometry/grid.h>
#include <utils/tiny_tools.h>
#include <numerics/gauss_data.h>
#include <numerics/iteratsolv.h>

using MathTL::Grid;
using MathTL::SampledMapping;
using MathTL::Point;
using MathTL::SparseMatrix;

namespace WaveletTL
{
    
  template <class RFRAME>
  PeriodicFrame<RFRAME>::PeriodicFrame()
    : r_frame(),
      j0_(j0()),
      Mj0_(j0(), r_frame.offset_a, r_frame.band_a, M_SQRT1_2),
      Mj1_(j0(), r_frame.offset_b, r_frame.band_b, M_SQRT1_2),
      Mj0T_(j0(), r_frame.offset_aT, r_frame.band_aT, M_SQRT1_2),
      Mj1T_(j0(), r_frame.offset_bT, r_frame.band_bT, M_SQRT1_2)
      
  {
//     cout << "Mj0=" << endl << Mj0_ << endl;
//     cout << "Mj0T=" << endl << Mj0T_ << endl;
//     cout << "Mj1=" << endl << Mj1_ << endl;
//     cout << "Mj1T=" << endl << Mj1T_ << endl;
      
  }
  
  
  
  template <class RFRAME>
  template <class V>
  inline
  void
  PeriodicFrame<RFRAME>::apply_Mj0(const int j, const V& x, V& y,
				   const size_type x_offset, const size_type y_offset,
				   const bool add_to) const
  {
    Mj0_.set_level(j);
    Mj0_.apply(x, y, x_offset, y_offset, add_to);
  }

  template <class RFRAME>
  template <class V>
  inline
  void
  PeriodicFrame<RFRAME>::apply_Mj1(const int j, const V& x, V& y,
				   const size_type x_offset, const size_type y_offset,
				   const bool add_to) const
  {
    Mj1_.set_level(j);
    Mj1_.apply(x, y, x_offset, y_offset, add_to);
  }
  
  
  
  template <class RFRAME>
  template <class V>
  void
  PeriodicFrame<RFRAME>::apply_Mj(const int j, const V& x, V& y) const
  {
    Mj0_.set_level(j);
    Mj1_.set_level(j);
    
    // decompose x=(x1 x2) appropriately
    Mj0_.apply(x, y, 0, 0, false);                      // apply Mj0 to first block x1
    Mj1_.apply(x, y, Mj0_.column_dimension(), 0, true); // apply Mj1 to second block x2 and add result
  }
  
  
  template <class RFRAME>
  template <class V>
  void
  PeriodicFrame<RFRAME>::apply_Tj(const int j, const V& x, V& y) const
  { 
    y = x;
    V z(x);
    apply_Mj(j0(), z, y);
    for (int k = j0()+1; k <= j; k++) {
      apply_Mj(k, y, z);
      y.swap(z);
    }
  }

  template <class RFRAME>
  void
  PeriodicFrame<RFRAME>::support(const Index& lambda, int& k1, int& k2)
  {
    if (lambda.e() == 0) // generator
      {
	// For the generators on the real line, we have
        //   \supp\phi_{j,k} = 2^{-j}[ell1+k,ell2+k]
	k1 = dyadic_modulo(RFRAME::primal_mask::begin()+lambda.k(), lambda.j());
	k2 = dyadic_modulo(RFRAME::primal_mask::end()  +lambda.k(), lambda.j());
      }
    else // wavelet
      {
	// For the wavelets on the real line, we have
        //   \supp\phi_{j,k} = 2^{-(j+1)}[ell1+1-ell2T+2*k,ell2+1-ell1T+2*k]
	k1 = dyadic_modulo(RFRAME::primal_mask::begin()+1-RFRAME::dual_mask::end()+2*lambda.k(), lambda.j()+1);
	k2 = dyadic_modulo(RFRAME::primal_mask::end()+1-RFRAME::dual_mask::begin()+2*lambda.k(), lambda.j()+1);
      }
      }

  template <class RFRAME>
  bool
  PeriodicFrame<RFRAME>::intersect_supports(const Index& lambda, const Index& mu)
  {
    // A note on the strategy:
    // Both supp(psi_lambda) and supp(psi_mu) are subintervals of the circle.
    // In order to decide whether they intersect or not, we compute the distance
    // between the respective midpoints.

    // compute dyadic coordinate of the midpoint of supp(psi_lambda)
    const int j_lambda = lambda.j()+1; // independent from lambda.e()!

    int k1_lambda, k2_lambda;
    if (lambda.e() == 0) {
      k1_lambda = 2*dyadic_modulo(RFRAME::primal_mask::begin()+lambda.k(), lambda.j());
      k2_lambda = 2*dyadic_modulo(RFRAME::primal_mask::end()  +lambda.k(), lambda.j());
    } else {
      k1_lambda = dyadic_modulo(RFRAME::primal_mask::begin()+1-RFRAME::dual_mask::end()+2*lambda.k(), j_lambda);
      k2_lambda = dyadic_modulo(RFRAME::primal_mask::end()+1-RFRAME::dual_mask::begin()+2*lambda.k(), j_lambda);
    }
    int length_lambda = k2_lambda-k1_lambda+(k2_lambda>k1_lambda ? 0 : 1<<j_lambda);
    int mid_lambda = dyadic_modulo((k2_lambda+k1_lambda+(k2_lambda>k1_lambda ? 0 : 1<<j_lambda))/2, j_lambda);

    // do the same for supp(psi_mu)
    const int j_mu = mu.j()+1;
    int k1_mu, k2_mu;
    if (mu.e() == 0) {
      k1_mu = 2*dyadic_modulo(RFRAME::primal_mask::begin()+mu.k(), mu.j());
      k2_mu = 2*dyadic_modulo(RFRAME::primal_mask::end()  +mu.k(), mu.j());
    } else {
      k1_mu = dyadic_modulo(RFRAME::primal_mask::begin()+1-RFRAME::dual_mask::end()+2*mu.k(), mu.j()+1);
      k2_mu = dyadic_modulo(RFRAME::primal_mask::end()+1-RFRAME::dual_mask::begin()+2*mu.k(), mu.j()+1);
    }
    int length_mu = k2_mu-k1_mu+(k2_mu>k1_mu ? 0 : 1<<j_mu);
    int mid_mu = dyadic_modulo((k2_mu+k1_mu+(k2_mu>k1_mu ? 0 : 1<<j_mu))/2, j_mu);
    
    const int j = std::max(j_lambda, j_mu);
    mid_lambda    <<= (j-j_lambda);
    length_lambda <<= (j-j_lambda);
    mid_mu        <<= (j-j_mu);
    length_mu     <<= (j-j_mu);

    const int dist_midpoints =
      mid_lambda > mid_mu
      ? std::min(mid_lambda-mid_mu, mid_mu+(1<<j)-mid_lambda)
      : std::min(mid_mu-mid_lambda, mid_lambda+(1<<j)-mid_mu);
    
    return (2*dist_midpoints < length_lambda+length_mu);
  }

  template <class RFRAME>
  inline
  typename PeriodicFrame<RFRAME>::Index
  PeriodicFrame<RFRAME>::first_generator(const int j, const int p)
  {
    assert(j >= j0());
    return Index(p, j, 0, PeriodicFrame<RFRAME>::DeltaLmin());
  }
  
  template <class RFRAME>
  inline
  typename PeriodicFrame<RFRAME>::Index
  PeriodicFrame<RFRAME>::last_generator(const int j, const int p)
  {
    assert(j >= j0());
    return Index(p, j, 0, PeriodicFrame<RFRAME>::DeltaRmax(j));
  }

  template <class RFRAME>
  inline
  typename PeriodicFrame<RFRAME>::Index
  PeriodicFrame<RFRAME>::first_wavelet(const int j, const int p)
  {
    assert(j >= j0());
    return Index(p, j, 1, PeriodicFrame<RFRAME>::Nablamin());
  }
  
  template <class RFRAME>
  inline
  typename PeriodicFrame<RFRAME>::Index
  PeriodicFrame<RFRAME>::last_wavelet(const int j, const int p)
  {
    assert(j >= j0());
    return Index(p, j, 1, PeriodicFrame<RFRAME>::Nablamax(j));
  }

  template <class RFRAME>
  inline
  typename PeriodicFrame<RFRAME>::Index
  PeriodicFrame<RFRAME>::first_index(const int j, const int e)
  {
    return (e == 0 ? first_generator(j) : first_wavelet(j));
  }
  
  template <class RFRAME>
  inline
  typename PeriodicFrame<RFRAME>::Index
  PeriodicFrame<RFRAME>::last_index(const int j, const int e)
  {
    return (e == 0 ? last_generator(j) : last_wavelet(j));
    }

  template <class RFRAME>
  SampledMapping<1>
  PeriodicFrame<RFRAME>::evaluate
  (const typename PeriodicFrame<RFRAME>::Index& lambda,
   const int resolution, const bool normalization) const
  {
    Grid<1> grid(0, 1, 1<<resolution);
    Array1D<double> values((1<<resolution)+1);
    for (unsigned int i(0); i < values.size(); i++) {
      const double x = i*ldexp(1.0, -resolution);
      values[i] = evaluate(0, lambda, x, normalization);
    }     
    return SampledMapping<1>(grid, values);
  }
  
  template <class RFRAME>
    SampledMapping<1>
    PeriodicFrame<RFRAME>::evaluate
    (const InfiniteVector<double, Index>& coeffs,
            const int resolution,
            const int derivative, const bool normalization) const {
        Grid<1> grid(0, 1, 1 << resolution);
        SampledMapping<1> result(grid); // zero
        if (coeffs.size() > 0) {
            // determine maximal level
            
            int number(0);
            
            //typedef typename SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::Index Index;
            for (typename InfiniteVector<double, Index>::const_iterator it(coeffs.begin()),
                    itend(coeffs.end()); it != itend; ++it){
                number = std:: max((1<<(jmax_+1)) * it.index().p() + it.index().number(), number);
            }
            //cout << "nummer: " << number << endl;@PHK

            // insert coefficients into a dense vector
            //Vector<double> wcoeffs((1<<(jmax_+1)) * pmax + Deltasize(jmax));
            Vector<double> wcoeffs(number+1);
            for (typename InfiniteVector<double, Index>::const_iterator it(coeffs.begin()),
                    itend(coeffs.end()); it != itend; ++it) {
                // determine number of the wavelet
                typedef typename Vector<double>::size_type size_type;
                size_type number = 0;
                if (it.index().e() == 0) {
                    number = (1<<(jmax_+1)) * it.index().p() + it.index().k() - DeltaLmin();
                } else {
                    number = (1<<(jmax_+1)) * it.index().p() + Deltasize(it.index().j()) + it.index().k() - Nablamin();
                }
                wcoeffs[number] = *it;
            }
            
            //cout << wcoeffs << endl;@PHK
            
            

//            // switch to generator representation
//            Vector<double> gcoeffs(wcoeffs.size(), false);
//            if (jmax == j0()){
//                gcoeffs = wcoeffs;
//            cout << "ja" << endl;
//            }
//            else
//                apply_Tj(jmax - 1, wcoeffs, gcoeffs);

            Array1D<double> values((1 << resolution) + 1);
            for (unsigned int i(0); i < values.size(); i++) {
                values[i] = 0;
                const double x = i * ldexp(1.0, -resolution);
                for (unsigned int k = 0; k < wcoeffs.size(); k++) {
                    if(wcoeffs[k] != 0)
                    values[i] += wcoeffs[k] * evaluate(derivative, *get_quarklet(k), x, normalization);
                }
            }

            return SampledMapping<1>(grid, values);
        }

        return result;
    }

  
  template <class RFRAME>
  inline
  double
  PeriodicFrame<RFRAME>::evaluate
  (const unsigned int derivative, const Index& lambda, const double x, const bool normalization) const
  {
    
      if(lambda.e() == 0){
        if(derivative == 0 && normalization == 1)
              return r_frame.evaluate(derivative, lambda,
            		    x-floor(x-ldexp(1.0,-lambda.j())
            			    *(RFRAME::primal_mask::begin()
            			       +lambda.k())))
                      - r_frame.integrate(lambda);
        else
               return r_frame.evaluate(derivative, lambda,
                                x-floor(x-ldexp(1.0,-lambda.j())
				    *(RFRAME::primal_mask::begin()
				      +lambda.k())));
      }
      
          
      else
      return r_frame.evaluate(derivative, lambda,
			    x-floor(x-ldexp(1.0,-lambda.j())
				    *((lambda.e() == 0
				       ? RFRAME::primal_mask::begin()
				       : (RFRAME::primal_mask::begin()+1-RFRAME::dual_mask::end())/2)
				      +lambda.k())));
  }
  
  template <class RFRAME>
  inline
  void
  PeriodicFrame<RFRAME>::evaluate(const unsigned int derivative,
				  const Index& lambda,
				  const Array1D<double>& points,
				  Array1D<double>& values, const bool normalization) const
  {
    values.resize(points.size());
    for (unsigned int i = 0; i < points.size(); i++)
      values[i] = evaluate(derivative, lambda, points[i], normalization);
  }
    
  template <class RFRAME>
  inline
  void
  PeriodicFrame<RFRAME>::evaluate(const Index& lambda,
				  const Array1D<double>& points,
				  Array1D<double>& funcvalues,
				  Array1D<double>& dervalues, const bool normalization) const
  {
    funcvalues.resize(points.size());
    dervalues.resize(points.size());
    evaluate(0, lambda, points, funcvalues, normalization);
    evaluate(1, lambda, points, dervalues );
  }

  template <class RFRAME>
  void
  PeriodicFrame<RFRAME>::expand(const Function<1>* f,
				const bool primal,
				const int jmax,
                                const int pmax, 
				InfiniteVector<double, Index>& coeffs, const bool normalization) const
  {
      int p = 0;
    for (Index lambda = first_generator(j0());;)
      {
	coeffs.set_coefficient(lambda, integrate(f, lambda, normalization));
        
            
	if (lambda == last_wavelet(jmax,pmax))
	  break;
        if (lambda == last_wavelet(jmax,p)){
            ++p;
            lambda = first_generator(j0(), p);
        }
        else
            ++lambda;
        
      }

    if (!primal) {
      // setup active index set
      std::set<Index> Lambda;
      p = 0;
      for (Index lambda = first_generator(j0());;) {
 	Lambda.insert(lambda);
	if (lambda == last_wavelet(jmax,pmax))
	  break;
        if (lambda == last_wavelet(jmax,p)){
            ++p;
            lambda = first_generator(j0(), p);
        }
        else
            ++lambda;
        
      }
      
      // setup Gramian A_Lambda
      SparseMatrix<double> A_Lambda(Lambda.size());
      typedef typename SparseMatrix<double>::size_type size_type;     
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
	      double entry = integrate(0, *it2, *it1, normalization);
	      
	      if (entry != 0) {
		indices.push_back(column);
		entries.push_back(entry);
	      }
	    }
	  A_Lambda.set_row(row, indices, entries);
	} 
      
      // solve A_Lambda*x = b
      Vector<double> b(Lambda.size());
      row = 0;
      for (typename std::set<Index>::const_iterator it(Lambda.begin()), itend(Lambda.end());
	   it != itend; ++it, ++row)
	b[row] = coeffs.get_coefficient(*it);
      
      Vector<double> x(b);
      unsigned int iterations;
      CG(A_Lambda, b, x, 1e-15, 500, iterations);
      cout << iterations << endl;
      
      coeffs.clear();
      row = 0;
      for (typename std::set<Index>::const_iterator it(Lambda.begin()), itend(Lambda.end());
	   it != itend; ++it, ++row)
	coeffs.set_coefficient(*it, x[row]);
    }
//    InfiniteVector<double, Index> helpcoeffs;
//    coeffs.COARSE(1e-6,helpcoeffs);
//    coeffs = helpcoeffs;
  }
//
//  template <class RBASIS>
//  void
//  PeriodicBasis<RBASIS>::expand(const Function<1>* f,
//				const bool primal,
//				const int jmax,
//				Vector<double>& coeffs) const
//  {
//    assert(primal);
//
//    coeffs.resize(Deltasize(jmax+1));
//
//    int id = 0;
//    for (Index lambda = first_generator(j0());;++lambda, ++id)
//      {
//	coeffs[id] = integrate(f, lambda);
//	if (lambda == last_wavelet(jmax))
//	  break;
//      }   
//  }
  
  template <class RFRAME>
  double
  PeriodicFrame<RFRAME>::integrate(const Function<1>* f,
				   const Index& lambda, const bool& normalization) const
  {
    double r = 0;

    // first we compute the support of psi_lambda
    const int j = lambda.j()+lambda.e();
    int k1, k2;
    support(lambda, k1, k2); // note: k2 may be less or equal to k1 in the case of overlap
    const int length = (k2 > k1 ? k2-k1 : k2+(1<<j)-k1); // number of subintervals
    
    // setup Gauss points and weights for a composite quadrature formula:
    unsigned int N_Gauss = 5;
    const double h = 1.0/(1<<j);

    Array1D<double> gauss_points (N_Gauss*length);
    int k = k1;
    for (int patch = 0; patch < length; patch++, k = dyadic_modulo(++k,j)) // work on 2^{-j}[k,k+1]
      for (unsigned int n = 0; n < N_Gauss; n++)
 	gauss_points[patch*N_Gauss+n] = h*(2*k+1+GaussPoints[N_Gauss-1][n])/2;
    
    // add all integral shares
    for (unsigned int n = 0; n < N_Gauss; n++)
      {
 	const double gauss_weight = GaussWeights[N_Gauss-1][n] * h;
 	for (int patch = k1; patch < k1+length; patch++)
 	  {
 	    const double t = gauss_points[(patch-k1)*N_Gauss+n];
	    
 	    const double ft = f->value(MathTL::Point<1>(t));
 	    if (ft != 0)
 	      r += ft
 		* evaluate(0, lambda, t,0)
 		* gauss_weight;
 	  }
      }
    
      if(lambda.e() == 0 && normalization == 1){
        N_Gauss = 9;
        Array1D<double> gauss_points2 (N_Gauss);
            
            for (unsigned int n = 0; n < N_Gauss; n++){
                    gauss_points2[n] = 0.5+GaussPoints[N_Gauss-1][n]/2;
            }
                    
            
            
//                    cout << h*k1 << ", " << h*k2 << endl;
//                    cout << gauss_points2 << endl;
            
                    const double pvalue = evaluate(0, lambda, h*k1);
//                    cout << pvalue << endl;
            // - add all integral shares
                for (unsigned int n = 0; n < N_Gauss; n++) {
                        const double t = gauss_points2[n];
                        const double gauss_weight = GaussWeights[N_Gauss - 1][n];

                        const double ft = f->value(MathTL::Point<1>(t));
                        //cout << "(" << t << ", " << gt << ")" <<endl;
                        if (ft != 0)
                            r += ft
                                * pvalue
                                * gauss_weight;
                        //cout << r << endl;
                }
            
       
       
        }
    
    return r;
    }
  
  

  template <class RFRAME>
  double
  PeriodicFrame<RFRAME>::integrate(const unsigned int& derivative,
				   const Index& lambda,
				   const Index& mu, const bool& normalization) const
  {
    double r = 0;
    
    if (intersect_supports(lambda, mu))
      {
	// first we determine the support over which to integrate
	int j, k1, k2, length;
	if (lambda.j()+lambda.e() >= mu.j()+mu.e()) {
	  j = lambda.j()+lambda.e();
	  support(lambda, k1, k2); // note: k2 may be less or equal to k1 in the case of overlap
	} else {
	  j = mu.j()+mu.e();
	  support(mu, k1, k2); // note: k2 may be less or equal to k1 in the case of overlap
	}
	length = (k2 > k1 ? k2-k1 : k2+(1<<j)-k1); // number of subintervals
	
	// setup Gauss points and weights for a composite quadrature formula:
 	const unsigned int N_Gauss = primal_polynomial_degree()+(lambda.p()+mu.p())/2;//@PHK changed to p-setting
	const double h = 1.0/(1<<j);
	
	Array1D<double> gauss_points (N_Gauss*(length)), func1values, func2values;
	int k = k1;
	for (int patch = 0; patch < length; patch++, k = dyadic_modulo(++k,j)) // work on 2^{-j}[k,k+1]
	  for (unsigned int n = 0; n < N_Gauss; n++)
	    gauss_points[patch*N_Gauss+n] = h*(2*k+1+GaussPoints[N_Gauss-1][n])/2;

  	// - compute point values of the integrands PROBLEM WITH GRAMIAN @PHK stimmt nicht, da nicht das gesamte Intervall berücksichtigt wird
//   	evaluate(derivative, lambda, gauss_points, func1values, normalization);
//  	evaluate(derivative, mu, gauss_points, func2values, normalization);
        
        evaluate(derivative, lambda, gauss_points, func1values, 0);
  	evaluate(derivative, mu, gauss_points, func2values, 0);
	
  	// - add all integral shares
	for (int patch = k1, id = 0; patch < k1+length; patch++)
	  for (unsigned int n = 0; n < N_Gauss; n++, id++) {
	    r += func1values[id] * func2values[id] * GaussWeights[N_Gauss-1][n] * h;
	  }

	return r;
      }

    return 0;
  } 

  
  template <class RFRAME>
  void
  PeriodicFrame<RFRAME>::setup_full_collection()
  {
    if (jmax_ == -1 || jmax_ < j0_ || pmax_ == -1) {
      cout << "PeriodicFrame<RFRAME>::setup_full_collection(): specify a maximal level of resolution first!" << endl;
      abort();
    }   

    int degrees_of_freedom = Deltasize(jmax_+1)*(pmax_+1);
    cout << "total degrees of freedom between (j0_, 0) and (jmax_, pmax_) is " << degrees_of_freedom << endl;
    //cout << "j0_:" << j0_ << endl;
    cout << "setting up collection of wavelet indices..." << endl;
    
    //cout << "jmax= " << jmax_ << endl;
    full_collection.resize(degrees_of_freedom);
    int k = 0;
    int p = 0;
    
    for (Index ind = first_generator(j0_); ind <= last_wavelet(jmax_, pmax_); (ind == last_wavelet(jmax_,p))? (++p, ind=first_generator(j0_, p)) : ++ind) {
      full_collection[k] = ind;
      //cout << ind << endl;
      k++;
    }
    cout << "done setting up collection of wavelet indices..." << endl;
    

  }
}
