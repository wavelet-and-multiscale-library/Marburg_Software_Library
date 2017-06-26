// implementation for periodic_frame_support.h

#include <algebra/matrix.h>
#include <Rd/cdf_utils.h>

#include "periodic_frame.h"

namespace WaveletTL
{
    
  template <class RFRAME>
  void intersecting_quarklets(const PeriodicFrame<RFRAME>& basis,
                             const typename PeriodicFrame<RFRAME>::Index& lambda,
			     const int j, const bool generators, 
                             std::list<typename PeriodicFrame<RFRAME>::Index>& intersecting,
                             const int p)
  {
      typedef PeriodicFrame<RFRAME> Basis;
      typedef typename Basis::Index Index;
    intersecting.clear();
      
    // a brute force solution
    if (generators) {
      for (Index nu = first_q_generator<Basis>(j,p);; ++nu) {
	if (basis.intersect_supports(lambda, nu))
	  intersecting.push_back(nu);
	if (nu == last_q_generator<Basis>(j,p)) break;
      }
    } else {
      for (Index nu = first_quarklet<Basis>(j,p);; ++nu) {
	if (basis.intersect_supports(lambda,nu))
	  intersecting.push_back(nu);
	if (nu == last_quarklet<Basis>(j,p)) break;
      }
    }
  }
  
#if 0  
//work in progress NOT NECESSARY FOR OUR PURPOSES
  template <class RBasis>
  bool intersect_supports(const PeriodicBasis<RBasis>& basis,
                          const typename PeriodicBasis<RBasis>::Index& lambda,
			  const int m, const int a, const int b,
			  int& j, int& k1, int& k2)
  {
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
  
  //computing the support doesn't work until now because of difficulties with periodicity
  template <class RBasis>
  bool intersect_supports(const PeriodicBasis<RBasis>& basis,
                          const typename PeriodicBasis<RBasis>::Index& lambda,
			  const typename PeriodicBasis<RBasis>::Index& nu,
			  typename PeriodicBasis<RBasis>::Support& supp)
  {
    int k1_nu, k2_nu;
    support(basis, nu, k1_nu, k2_nu);
    return intersect_supports(basis, lambda, nu.j()+nu.e(), k1_nu, k2_nu,
			      supp.j, supp.k1, supp.k2);
  }
  //computing the support doesn't work until now because of difficulties with periodicity 

  template <class RBasis>
  void intersecting_quarklets(const PeriodicBasis<RBasis>& basis,
                             const typename PeriodicBasis<RBasis>::Index& lambda,
			     const int j, const bool generators, std::list<std::pair<typename PeriodicBasis<RBasis>::Index,
			     typename PeriodicBasis<RBasis>::Support> >& intersecting)
  {
    typedef typename PeriodicBasis<RBasis>::Index Index;
    typedef typename PeriodicBasis<RBasis>::Support Support;

    intersecting.clear();

    // compute support of \psi_\lambda
    const int j_lambda = lambda.j() + lambda.e();
    int k1_lambda, k2_lambda;
    support(basis, lambda, k1_lambda, k2_lambda);
    
    // a brute force solution
    Support supp;
    if (generators) {
      for (Index nu = first_generator(&basis, j);; ++nu) {
	if (intersect_supports(basis, nu, j_lambda, k1_lambda, k2_lambda, supp.j, supp.k1, supp.k2))
	  intersecting.push_back(std::make_pair(nu, supp));
	if (nu == last_generator(&basis, j)) break;
      }
    } else {
      for (Index nu = first_wavelet(&basis, j);; ++nu) {
	if (intersect_supports(basis, nu, j_lambda, k1_lambda, k2_lambda, supp.j, supp.k1, supp.k2))
	  intersecting.push_back(std::make_pair(nu, supp));
	if (nu == last_wavelet(&basis, j)) break;
      }
    }
  }
  */
  
  /*
  template <int d, int dT, DSBiorthogonalizationMethod BIO>
  bool intersect_singular_support(const DSBasis<d,dT,BIO>& basis,
				  const typename DSBasis<d,dT,BIO>::Index& lambda,
				  const int m, const int a, const int b,
				  int& j, int& k1, int& k2)
  {
    const int j_lambda = lambda.j() + lambda.e();
    
    if (j_lambda < m) {
      // granularity of support intersection is m
      return intersect_supports(basis, lambda, m, a, b, j, k1, k2);
    } else {
      // Here j_lambda >= m, i.e., the support of psi_lambda determines the
      // granularity of the support intersection
      
      int k1_lambda, k2_lambda;
      support(basis, lambda, k1_lambda, k2_lambda);
      
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

  template <int d, int dT, DSBiorthogonalizationMethod BIO>
  bool intersect_singular_support(const DSBasis<d,dT,BIO>& basis,
				  const typename DSBasis<d,dT,BIO>::Index& lambda,
				  const typename DSBasis<d,dT,BIO>::Index& nu,
				  int& j, int& k1, int& k2)
  {
    const int j_lambda = lambda.j() + lambda.e();
    const int j_nu = nu.j() + nu.e();

    if (j_lambda < j_nu)
      return intersect_singular_support(basis, nu, lambda, j, k1, k2);

    int k1_nu, k2_nu;
    support(basis, nu, k1_nu, k2_nu);
    return intersect_singular_support(basis, lambda, j_nu, k1_nu, k2_nu, j, k1, k2);
  }

  

  template <int d, int dT, DSBiorthogonalizationMethod BIO>
  void relevant_wavelets(const DSBasis<d,dT,BIO>& basis,
			 const typename DSBasis<d,dT,BIO>::Index& lambda,
			 const int j, const bool generators,
			 std::list<typename DSBasis<d,dT,BIO>::Index>& relevant)
  {
    typedef typename DSBasis<d,dT,BIO>::Index Index;
    typedef typename DSBasis<d,dT,BIO>::Support Support;
    
    relevant.clear();

    // compute support of \psi_\lambda
    const int j_lambda = lambda.j() + lambda.e();
    int k1_lambda, k2_lambda;
    support(basis, lambda, k1_lambda, k2_lambda);
    
    // a brute force solution
    if (generators) {
      Support supp;
      for (Index nu = first_generator(&basis, j);; ++nu) {
	if (intersect_singular_support(basis, nu, j_lambda, k1_lambda, k2_lambda, supp.j, supp.k1, supp.k2))
	  relevant.push_back(nu);
	if (nu == last_generator(&basis, j)) break;
      }
    } else {
      Support supp;
      for (Index nu = first_wavelet(&basis, j);; ++nu) {
	if (intersect_singular_support(basis, nu, j_lambda, k1_lambda, k2_lambda, supp.j, supp.k1, supp.k2))
	  relevant.push_back(nu);
	if (nu == last_wavelet(&basis, j)) break;
      }
    }
  }

  template <int d, int dT, DSBiorthogonalizationMethod BIO>
  void relevant_wavelets(const DSBasis<d,dT,BIO>& basis,
			 const typename DSBasis<d,dT,BIO>::Index& lambda,
			 const int j, const bool generators,
			 std::list<std::pair<typename DSBasis<d,dT,BIO>::Index,
			 typename DSBasis<d,dT,BIO>::Support> >& relevant)
  {
    typedef typename DSBasis<d,dT,BIO>::Index Index;
    typedef typename DSBasis<d,dT,BIO>::Support Support;
    
    relevant.clear();

    // compute support of \psi_\lambda
    const int j_lambda = lambda.j() + lambda.e();
    int k1_lambda, k2_lambda;
    support(basis, lambda, k1_lambda, k2_lambda);
    
    // a brute force solution
    if (generators) {
      Support supp;
      for (Index nu = first_generator(&basis, j);; ++nu) {
	if (intersect_singular_support(basis, nu, j_lambda, k1_lambda, k2_lambda, supp.j, supp.k1, supp.k2))
	  relevant.push_back(std::make_pair(nu, supp));
	if (nu == last_generator(&basis, j)) break;
      }
    } else {
      Support supp;
      for (Index nu = first_wavelet(&basis, j);; ++nu) {
	if (intersect_singular_support(basis, nu, j_lambda, k1_lambda, k2_lambda, supp.j, supp.k1, supp.k2))
	  relevant.push_back(std::make_pair(nu, supp));
	if (nu == last_wavelet(&basis, j)) break;
      }
    }
  }
   */#endif
}
