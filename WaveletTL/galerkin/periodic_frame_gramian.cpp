// implementation for periodic_gramian.h

#include <cmath>
#include <algorithm>
#include <list>
#include <utils/array1d.h>
#include <algebra/vector.h>
#include <algebra/sparse_matrix.h>
#include <numerics/eigenvalues.h>
#include <numerics/gauss_data.h>


namespace WaveletTL
{
  template <class RFRAME>
  PeriodicFrameIntervalGramian<RFRAME>::PeriodicFrameIntervalGramian
  (const PeriodicLaplacianProblem& plp, const PeriodicFrame<RFRAME>& frame)
    : basis_(frame), plp_(plp),
      normA(1.0),
      normAinv(ldexp(1.0, 2*(RFRAME::primal_polynomial_degree()-1))) // lower bound from [Bittner]
  {
      precompute_rhs();
      //const int jmax = 12;
      //basis_.set_jmax(jmax);
      //cout << "TESTTEST" << endl;
  }
  
  template <class RFRAME>
  inline
  double
  PeriodicFrameIntervalGramian<RFRAME>::a(const typename WaveletBasis::Index& lambda,
				     const typename WaveletBasis::Index& mu) const
  {
    return basis_.integrate(0, lambda, mu);
  }
  
  
  template <class RFRAME>
  inline
  void
  PeriodicFrameIntervalGramian<RFRAME>::RHS(const double eta,
			     InfiniteVector<double, Index>& coeffs) const
  {
    if (!rhs_precomputed) precompute_rhs();

    coeffs.clear();
    double coarsenorm(0);
    double bound(fnorm_sqr - eta*eta);
    //typedef typename RBASIS::Index Index;
    typename Array1D<std::pair<Index, double> >::const_iterator it(fcoeffs.begin());
    do {
      coarsenorm += it->second * it->second;
      coeffs.set_coefficient(it->first, it->second);
      ++it;
    } while (it != fcoeffs.end() && coarsenorm < bound);
  }
  
  
  template <class RFRAME>
  void
  PeriodicFrameIntervalGramian<RFRAME>::precompute_rhs() const
  {
    //typedef typename WaveletBasis::Index Index;
    
    // precompute the right-hand side on a fine level
    InfiniteVector<double,Index> fhelp;
    const int j0   = basis().j0();
    const int jmax = basis_.get_jmax_();
    const int pmax = basis_.get_pmax_();
    int p= 0;
    for (Index lambda(basis_.first_generator(j0));; )
      {
	const double coeff = f(lambda)/D(lambda);
        //cout << lambda << ", " << f(lambda) << ", " << D(lambda) << endl;;
	if (fabs(coeff)>1e-15)
	  fhelp.set_coefficient(lambda, coeff);
	if (lambda == basis_.last_wavelet(jmax, pmax))
	  break;
        if (lambda == basis_.last_wavelet(jmax, p)){
            ++p;
            lambda = basis_.first_generator(j0,p);
        }
        else
            ++lambda;
            
      }
    fnorm_sqr = l2_norm_sqr(fhelp);
    
    // sort the coefficients into fcoeffs
    fcoeffs.resize(fhelp.size());
    unsigned int id(0);
    for (typename InfiniteVector<double,Index>::const_iterator it(fhelp.begin()), itend(fhelp.end());
	 it != itend; ++it, ++id)
      fcoeffs[id] = std::pair<Index,double>(it.index(), *it);
    sort(fcoeffs.begin(), fcoeffs.end(), typename InfiniteVector<double,Index>::decreasing_order());
    
    rhs_precomputed = true;
  }
  
  template <class RFRAME>
  double
  PeriodicFrameIntervalGramian<RFRAME>::f(const Index& lambda) const
  {
    // f(v) = \int_0^1 g(t)v(t) dt

    double r = 0;
    
    // first we compute the support of psi_lambda
    const int j = lambda.j()+lambda.e();
    int k1, k2;
    basis_.support(lambda, k1, k2); // note: k2 may be less or equal to k1 in the case of overlap
    const int length = (k2 > k1 ? k2-k1 : k2+(1<<j)-k1); // number of subintervals
    
    // setup Gauss points and weights for a composite quadrature formula:
    const unsigned int N_Gauss = 5;
    const double h = 1.0/(1<<j);

    Array1D<double> gauss_points (N_Gauss*length), vvalues;
    int k = k1;
    for (int patch = 0; patch < length; patch++, k = dyadic_modulo(++k,j)) // work on 2^{-j}[k,k+1]
      for (unsigned int n = 0; n < N_Gauss; n++)
 	gauss_points[patch*N_Gauss+n] = h*(2*k+1+GaussPoints[N_Gauss-1][n])/2;
    //cout << h*k1 << ", " << h*k2 << endl;
    //cout << gauss_points << endl;
    
    
    // - compute point values of the integrand
    basis_.evaluate(0, lambda, gauss_points, vvalues);
    //cout << vvalues << endl;
    // - add all integral shares
    for (int patch = k1, id = 0; patch < k1+length; patch++)
      for (unsigned int n = 0; n < N_Gauss; n++, id++) {
	const double t = gauss_points[id];
	const double gauss_weight = GaussWeights[N_Gauss-1][n] * h;
	    
	const double gt = plp_.g(t);
	if (gt != 0)
	  r += gt
	    * vvalues[id]
	    * gauss_weight;
      }
    //cout << r << endl;
    return r;
  }
}
