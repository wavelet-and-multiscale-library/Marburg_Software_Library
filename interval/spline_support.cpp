// implementation for p_support.h

#include <Rd/cdf_utils.h>

namespace WaveletTL
{
  template <int d, int dT, SplineBasisFlavor flavor, int s0, int s1, int sT0, int sT1>
  inline
  void support(const SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>& basis,
	       const typename SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::Index& lambda,
	       int& k1, int& k2)
  {
    basis.support(lambda, k1, k2);
  }
  
  template <int d, int dT, SplineBasisFlavor flavor, int s0, int s1, int sT0, int sT1>
  bool intersect_supports(const SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>& basis,
			  const typename SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::Index& lambda,
			  const int m, const int a, const int b,
			  int& j, int& k1, int& k2)
  {
    // universal routine for any value of "flavor"
    const int j_lambda = lambda.j() + lambda.e();
    j = std::max(j_lambda, m);
    int k1_lambda, k2_lambda;
    basis.support(lambda, k1_lambda, k2_lambda);
    const int jmb = b<<(j-m);
    const int jjk1 = k1_lambda<<(j-j_lambda);
    if (jjk1 >= jmb)
      return false;
    else {
      const int jma = a<<(j-m);
      const int jjk2 = k2_lambda<<(j-j_lambda);
      if (jma >= jjk2)
	return false;
      else {
	k1 = std::max(jjk1, jma);
	k2 = std::min(jjk2, jmb);
      }
    }
    return true;
  }

  template <int d, int dT, SplineBasisFlavor flavor, int s0, int s1, int sT0, int sT1>
  inline
  bool intersect_supports(const SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>& basis,
			  const typename SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::Index& lambda,
			  const typename SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::Index& nu,
			  typename SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::Support& supp)
  {
    // universal routine for any value of "flavor"
    int k1_nu, k2_nu;
    basis.support(nu, k1_nu, k2_nu);
    return intersect_supports(basis, lambda, nu.j()+nu.e(), k1_nu, k2_nu,
			      supp.j, supp.k1, supp.k2);
  }

  // special variant of intersecting_wavelets() for flavor==P_construction:
  template <int d, int dT, int s0, int s1, int sT0, int sT1>
  void intersecting_wavelets(const SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1>& basis,
			     const typename SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1>::Index& lambda,
			     const int j, const bool generators,
			     std::list<std::pair<typename SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1>::Index,
			     typename SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1>::Support> >& intersecting)
  {
    typedef typename SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1>::Index Index;
    typedef typename SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1>::Support Support;

    intersecting.clear();

    // compute support of \psi_\lambda
    const int j_lambda = lambda.j() + lambda.e();
    int k1_lambda, k2_lambda;
    support(basis, lambda, k1_lambda, k2_lambda);
    
#if 1
    // new code
    Support supp;
    if (generators) {
      // the leftmost generator on the level j, s.th. its support intersects [a,b], fulfills
      //   2^{-j}(k+ell2) > a  but  2^{-j}(k-1+ell2) <= a,
      // so that ...
      const int firstk = std::max(basis.DeltaLmin(), (int) floor(ldexp(1.0,j-j_lambda)*k1_lambda-ell2<d>())+1);
      
      // the rightmost generator on the level j, s.th. its support intersects [a,b], fulfills
      //   2^{-j}(k+ell1) < b  but  2^{-j}(k+1+ell1) >= b,
      // so that ...
      const int lastk  = std::min(basis.DeltaRmax(j), (int) ceil(ldexp(1.0,j-j_lambda)*k2_lambda-ell1<d>())-1);
      
      Index nu(j, 0, firstk);
      for (int k = firstk; k <= lastk; k++, ++nu) {
	intersect_supports(basis, nu, j_lambda, k1_lambda, k2_lambda, supp.j, supp.k1, supp.k2);
	intersecting.push_back(std::make_pair(nu, supp));
      }
    } else {
      // if the left boundary wavelets are not involved, then
      // the leftmost wavelet on the level j, s.th. its support intersects [a,b], fulfills
      //   2^{-j}(k+(d+dT)/2) > a  but  2^{-j}(k-1+(d+dT)/2) <= a,
      // so that ...
      const int firstk = (ldexp(1.0,j-j_lambda)*k1_lambda < d+dT-1 // overestimate, TODO!
			  ? 0
			  : (int) floor(ldexp(1.0,j-j_lambda)*k1_lambda-(d+dT)/2)+1);

      // if the right boundary wavelets are not involved, then
      // the rightmost wavelet on the level j, s.th. its support intersects [a,b], fulfills
      //   2^{-j}(k-(d+dT)/2+1) < b  but  2^{-j}(k-(d+dT)/2+2) >= a,
      // so that ...
      const int lastk = (ldexp(1.0,-j_lambda)*k2_lambda > 1 - ldexp(1.0,-j)*(d+dT-1) // overestimate, TODO!
			 ? (1<<j)-1
			 : (int) ceil(ldexp(1.0,j-j_lambda)*k2_lambda+(d+dT)/2)-2);

      Index nu(j, 1, firstk);
      for (int k = firstk; k <= lastk; k++, ++nu) {
	intersect_supports(basis, nu, j_lambda, k1_lambda, k2_lambda, supp.j, supp.k1, supp.k2);
	intersecting.push_back(std::make_pair(nu, supp));
      }
    }
#else
    // old code, a brute force solution
    Support supp;
    if (generators) {
      for (Index nu = basis.first_generator(j);; ++nu) {
	if (intersect_supports(basis, nu, j_lambda, k1_lambda, k2_lambda, supp.j, supp.k1, supp.k2))
	  intersecting.push_back(std::make_pair(nu, supp));
	if (nu == basis.last_generator(j)) break;
      }
    } else {
      for (Index nu = basis.first_wavelet(j);; ++nu) {
	if (intersect_supports(basis, nu, j_lambda, k1_lambda, k2_lambda, supp.j, supp.k1, supp.k2))
	  intersecting.push_back(std::make_pair(nu, supp));
	if (nu == basis.last_wavelet(j)) break;
      }
    }
#endif
  }

  // generic variant of intersecting_wavelets() for DS type wavelets
  template <int d, int dT, SplineBasisFlavor flavor, int s0, int s1, int sT0, int sT1>
  void intersecting_wavelets(const SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>& basis,
			     const typename SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::Index& lambda,
			     const int j, const bool generators,
			     std::list<std::pair<typename SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::Index,
			     typename SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::Support> >& intersecting)
  {
    typedef typename SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::Index Index;
    typedef typename SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::Support Support;

    intersecting.clear();

    // compute support of \psi_\lambda
    const int j_lambda = lambda.j() + lambda.e();
    int k1_lambda, k2_lambda;
    support(basis, lambda, k1_lambda, k2_lambda);
    
    // a brute force solution
    Support supp;
    if (generators) {
      for (Index nu = basis.first_generator(j);; ++nu) {
	if (intersect_supports(basis, nu, j_lambda, k1_lambda, k2_lambda, supp.j, supp.k1, supp.k2))
	  intersecting.push_back(std::make_pair(nu, supp));
	if (nu == basis.last_generator(j)) break;
      }
    } else {
      for (Index nu = basis.first_wavelet(j);; ++nu) {
	if (intersect_supports(basis, nu, j_lambda, k1_lambda, k2_lambda, supp.j, supp.k1, supp.k2))
	  intersecting.push_back(std::make_pair(nu, supp));
	if (nu == basis.last_wavelet(j)) break;
      }
    }
  }

  // special variant of intersecting_wavelets() for flavor==P_construction
  template <int d, int dT, int s0, int s1, int sT0, int sT1>
  void intersecting_wavelets(const SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1>& basis,
			     const typename SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1>::Index& lambda,
			     const int j, const bool generators,
			     std::list<typename SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1>::Index>& intersecting)
  {
    typedef typename SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1>::Index Index;
    typedef typename SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1>::Support Support;

    intersecting.clear();

    // compute support of \psi_\lambda
    const int j_lambda = lambda.j() + lambda.e();
    int k1_lambda, k2_lambda;
    support(basis, lambda, k1_lambda, k2_lambda);
    
#if 1
    // new code
    if (generators) {
      // the leftmost generator on the level j, s.th. its support intersects [a,b], fulfills
      //   2^{-j}(k+ell2) > a  but  2^{-j}(k-1+ell2) <= a,
      // so that ...
      const int firstk = std::max(basis.DeltaLmin(), (int) floor(ldexp(1.0,j-j_lambda)*k1_lambda-ell2<d>())+1);
      
      // the rightmost generator on the level j, s.th. its support intersects [a,b], fulfills
      //   2^{-j}(k+ell1) < b  but  2^{-j}(k+1+ell1) >= b,
      // so that ...
      const int lastk  = std::min(basis.DeltaRmax(j), (int) ceil(ldexp(1.0,j-j_lambda)*k2_lambda-ell1<d>())-1);

      for (int k = firstk; k <= lastk; k++)
	intersecting.push_back(Index(j, 0, k));
      
    } else {
      // if the left boundary wavelets are not involved, then
      // the rightmost wavelet on the level j, s.th. its support intersects [a,b], fulfills
      //   2^{-j}(k+(d+dT)/2) > a  but  2^{-j}(k-1+(d+dT)/2) <= a,
      // so that ...
      const int firstk = (ldexp(1.0,j-j_lambda)*k1_lambda < d+dT-1 // overestimate, TODO!
			  ? 0
			  : (int) floor(ldexp(1.0,j-j_lambda)*k1_lambda-(d+dT)/2)+1);

      // if the right boundary wavelets are not involved, then
      // the leftmost wavelet on the level j, s.th. its support intersects [a,b], fulfills
      //   2^{-j}(k-(d+dT)/2+1) < b  but  2^{-j}(k-(d+dT)/2+2) >= a,
      // so that ...
      const int lastk = (ldexp(1.0,-j_lambda)*k2_lambda > 1 - ldexp(1.0,-j)*(d+dT-1) // overestimate, TODO!
			 ? (1<<j)-1
			 : (int) ceil(ldexp(1.0,j-j_lambda)*k2_lambda+(d+dT)/2)-2);

      for (int k = firstk; k <= lastk; k++)
	intersecting.push_back(Index(j, 1, k));
    }
#else
    // old code, a brute force solution
    Support supp;
    if (generators) {
      for (Index nu = basis.first_generator(j);; ++nu) {
	if (intersect_supports(basis, nu, j_lambda, k1_lambda, k2_lambda, supp.j, supp.k1, supp.k2))
	  intersecting.push_back(nu);
	if (nu == basis.last_generator(j)) break;
      }
    } else {
      for (Index nu = basis.first_wavelet(j);; ++nu) {
	if (intersect_supports(basis, nu, j_lambda, k1_lambda, k2_lambda, supp.j, supp.k1, supp.k2))
	  intersecting.push_back(nu);
	if (nu == basis.last_wavelet(j)) break;
      }
    }
#endif
  }

  // generic variant of intersecting_wavelets() for DS type wavelets
  template <int d, int dT, SplineBasisFlavor flavor, int s0, int s1, int sT0, int sT1>
  void intersecting_wavelets(const SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>& basis,
			     const typename SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::Index& lambda,
			     const int j, const bool generators,
			     std::list<typename SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::Index>& intersecting)
  {
    typedef typename SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::Index Index;
    typedef typename SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::Support Support;

    intersecting.clear();
    
    // compute support of \psi_\lambda
    const int j_lambda = lambda.j() + lambda.e();
    int k1_lambda, k2_lambda;
    support(basis, lambda, k1_lambda, k2_lambda);

    // a brute force solution
    Support supp;
    if (generators) {
      for (Index nu = basis.first_generator(j);; ++nu) {
	if (intersect_supports(basis, nu, j_lambda, k1_lambda, k2_lambda, supp.j, supp.k1, supp.k2))
	  intersecting.push_back(nu);
	if (nu == basis.last_generator(j)) break;
      }
    } else {
      for (Index nu = basis.first_wavelet(j);; ++nu) {
	if (intersect_supports(basis, nu, j_lambda, k1_lambda, k2_lambda, supp.j, supp.k1, supp.k2))
	  intersecting.push_back(nu);
	if (nu == basis.last_wavelet(j)) break;
      }
    }
  }



  template <int d, int dT, SplineBasisFlavor flavor, int s0, int s1, int sT0, int sT1>
  bool intersect_singular_support(const SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>& basis,
				  const typename SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::Index& lambda,
				  const int m, const int a, const int b,
				  int& j, int& k1, int& k2)
  {
    // universal routine for any value of "flavor"

    const int j_lambda = lambda.j() + lambda.e();
    
    if (j_lambda < m) {
      // granularity of support intersection is m
      return intersect_supports(basis, lambda, m, a, b, j, k1, k2);
    } else {
      // Here j_lambda >= m, i.e., the support of psi_lambda determines the
      // granularity of the support intersection
      
      int k1_lambda, k2_lambda;
      basis.support(lambda, k1_lambda, k2_lambda);
      
      // first exclude the case that supp(psi_lambda) and 2^{-m}[a,b] do not intersect at all
      const int b_help = (1<<(j_lambda-m)) * b;
      if (k1_lambda >= b_help)
	return false; // supp(psi_lambda) is completely right of 2^{-m}[a,b]
      else {
	const int a_help = (1<<(j_lambda-m)) * a;
	if (k2_lambda <= a_help)
	  return false; // supp(psi_lambda) is completely left of 2^{-m}[a,b]
	else {
	  if (k1_lambda >= a_help && k2_lambda <= b_help) {
	    // subset criterion, necessary for supp(psi_lambda) not hitting singsupp(2^{-m}[a,b])
	    
	    // determine largest integer n, such that 2^{-j}k1_lambda >= 2^{-m}n
	    // (we know that n < b)
	    const int n = (int) floor(ldexp(1.0, m-j_lambda) * k1_lambda);
	    if (k2_lambda <= (1<<(j_lambda-m)) * (n+1))
	      return false; // 2^{-m}k2_lambda <= 2^{-j}(n+1), i.e., no singuar support intersection
	  }
	  j = j_lambda;
	  k1 = std::max(k1_lambda, a_help);
	  k2 = std::min(k2_lambda, b_help);
	}
      }
    }
    return true;
  }

  template <int d, int dT, SplineBasisFlavor flavor, int s0, int s1, int sT0, int sT1>
  inline
  bool intersect_singular_support(const SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>& basis,
				  const typename SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::Index& lambda,
				  const typename SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::Index& nu,
				  int& j, int& k1, int& k2)
  {
    // universal routine for any value of "flavor"

    const int j_lambda = lambda.j() + lambda.e();
    const int j_nu = nu.j() + nu.e();

    if (j_lambda < j_nu)
      return intersect_singular_support(basis, nu, lambda, j, k1, k2);
    
    int k1_nu, k2_nu;
    basis.support(nu, k1_nu, k2_nu);
    return intersect_singular_support(basis, lambda, j_nu, k1_nu, k2_nu, j, k1, k2);
  }

  template <int d, int dT, SplineBasisFlavor flavor, int s0, int s1, int sT0, int sT1>
  inline
  bool intersect_singular_support(const SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>& basis,
				  const typename SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::Index& lambda,
				  const typename SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::Index& nu)
  {
    // universal routine for any value of "flavor"

    int j, k1, k2;
    return intersect_singular_support(basis, lambda, nu, j, k1, k2);
  }
}
