// implementation for jl_support.h

namespace WaveletTL
{
  void support(const JLBasis& basis,
 	       const Index& lambda,
 	       int& k1, int& k2)
  {
    // we use supp(phi_i)=supp(psi_i)=[-1,1]
    if (lambda.e() == 0) // generator
      {
	k1 = std::max(0, lambda.k()-1);
	k2 = std::min(1<<lambda.j(), lambda.k()+1);
      }
    else // wavelet
      {
	k1 = std::max(0, lambda.k()-1)*2;
	k2 = std::min(1<<lambda.j(), lambda.k()+1)*2;
      }
  }
  
  bool intersect_supports(const JLBasis& basis,
 			  const Index& lambda,
 			  const int m, const int a, const int b,
 			  int& j, int& k1, int& k2)
  {
    const int j_lambda = lambda.j() + lambda.e();
    j = std::max(j_lambda, m);
    int k1_lambda, k2_lambda;
    support(basis, lambda, k1_lambda, k2_lambda);
    const int jmb = (1<<(j-m)) * b;
    const int jjk1 = (1<<(j-j_lambda)) * k1_lambda;
    if (jjk1 >= jmb)
      return false;
    else {
      const int jma = (1<<(j-m)) * a;
      const int jjk2 = (1<<(j-j_lambda)) * k2_lambda;
      if (jma >= jjk2)
 	return false;
      else {
 	k1 = std::max(jjk1, jma);
 	k2 = std::min(jjk2, jmb);
      }
    }
    return true;
  }
  
  inline
  bool intersect_supports(const JLBasis& basis,
			  const Index& lambda,
			  const Index& nu,
			  JLBasis::Support& supp)
  {
    int k1_nu, k2_nu;
    support(basis, nu, k1_nu, k2_nu);
    return intersect_supports(basis, lambda, nu.j()+nu.e(), k1_nu, k2_nu,
			      supp.j, supp.k1, supp.k2);
  }
  
  void intersecting_wavelets(const JLBasis& basis,
			     const Index& lambda,
			     const int j, const bool generators,
			     std::list<std::pair<Index, JLBasis::Support> >& intersecting)
  {
    typedef JLBasis::Support Support;

    intersecting.clear();

    // compute support of \psi_\lambda
    const int j_lambda = lambda.j() + lambda.e();
    int k1_lambda, k2_lambda;
    support(basis, lambda, k1_lambda, k2_lambda);
    
#if 1
    // new code
    Support supp;
    supp.k1 = supp.k2 = 0;
    if (generators) {
      // the leftmost generator (of type phi_0) on the level j, s.th. its support intersects [a,b], fulfills
      //   2^{-j}(k+1) > a  but  2^{-j}k <= a,
      // so that ...
      const int firstk0 = std::max(1, (int) floor(ldexp(1.0,j-j_lambda)*k1_lambda));
      
      // the rightmost generator (of type phi_0) on the level j, s.th. its support intersects [a,b], fulfills
      //   2^{-j}(k-1) < b  but  2^{-j}k >= b,
      // so that ...
      const int lastk0  = std::min((1<<j)-1, (int) ceil(ldexp(1.0,j-j_lambda)*k2_lambda));
      
      // the leftmost generator (of type phi_1) on the level j, s.th. its support intersects [a,b], fulfills
      //   2^{-j}(k+1) > a  but  2^{-j}k <= a,
      // so that ...
      const int firstk1 = std::max(0, (int) floor(ldexp(1.0,j-j_lambda)*k1_lambda));
      
      // the rightmost generator (of type phi_1) on the level j, s.th. its support intersects [a,b], fulfills
      //   2^{-j}(k-1) < b  but  2^{-j}k >= b,
      // so that ...
      const int lastk1  = std::min(1<<j, (int) ceil(ldexp(1.0,j-j_lambda)*k2_lambda));

      Index nu(j, 0, 0, firstk0);
      for (int k = firstk0; k <= lastk0; k++, ++nu) {
  	intersect_supports(basis, nu, j_lambda, k1_lambda, k2_lambda, supp.j, supp.k1, supp.k2);
  	intersecting.push_back(std::make_pair(nu, supp));
      }
      nu = Index(j, 0, 1, firstk1);
      for (int k = firstk1; k <= lastk1; k++, ++nu) {
  	intersect_supports(basis, nu, j_lambda, k1_lambda, k2_lambda, supp.j, supp.k1, supp.k2);
  	intersecting.push_back(std::make_pair(nu, supp));
      }
    } else {
      // the leftmost wavelet (of type psi_0) on the level j, s.th. its support intersects [a,b], fulfills
      //   2^{-j}(k+1) > a  but  2^{-j}k <= a,
      // so that ...
      const int firstk0 = std::max(1, (int) floor(ldexp(1.0,j-j_lambda)*k1_lambda));
      
      // the rightmost wavelet (of type psi_0) on the level j, s.th. its support intersects [a,b], fulfills
      //   2^{-j}(k-1) < b  but  2^{-j}k >= b,
      // so that ...
      const int lastk0  = std::min((1<<j)-1, (int) ceil(ldexp(1.0,j-j_lambda)*k2_lambda));
      
      // the leftmost wavelet (of type psi_1) on the level j, s.th. its support intersects [a,b], fulfills
      //   2^{-j}(k+1) > a  but  2^{-j}k <= a,
      // so that ...
      const int firstk1 = std::max(0, (int) floor(ldexp(1.0,j-j_lambda)*k1_lambda));
      
      // the rightmost wavelet (of type psi_1) on the level j, s.th. its support intersects [a,b], fulfills
      //   2^{-j}(k-1) < b  but  2^{-j}k >= b,
      // so that ...
      const int lastk1  = std::min(1<<j, (int) ceil(ldexp(1.0,j-j_lambda)*k2_lambda));

      Index nu(j, 1, 0, firstk0);
      for (int k = firstk0; k <= lastk0; k++, ++nu) {
  	intersect_supports(basis, nu, j_lambda, k1_lambda, k2_lambda, supp.j, supp.k1, supp.k2);
  	intersecting.push_back(std::make_pair(nu, supp));
      }
      nu = Index(j, 1, 1, firstk1);
      for (int k = firstk1; k <= lastk1; k++, ++nu) {
  	intersect_supports(basis, nu, j_lambda, k1_lambda, k2_lambda, supp.j, supp.k1, supp.k2);
  	intersecting.push_back(std::make_pair(nu, supp));
      }
    }
#else
    // old code, a brute force solution
    Support supp;
    if (generators) {
      for (Index nu = first_generator(j);; ++nu) {
 	if (intersect_supports(basis, nu, j_lambda, k1_lambda, k2_lambda, supp.j, supp.k1, supp.k2))
 	  intersecting.push_back(std::make_pair(nu, supp));
 	if (nu == last_generator(j)) break;
      }
    } else {
      for (Index nu = first_wavelet(j);; ++nu) {
 	if (intersect_supports(basis, nu, j_lambda, k1_lambda, k2_lambda, supp.j, supp.k1, supp.k2))
 	  intersecting.push_back(std::make_pair(nu, supp));
 	if (nu == last_wavelet(j)) break;
      }
    }
#endif
  }

  void intersecting_wavelets(const JLBasis& basis,
  			     const Index& lambda,
  			     const int j, const bool generators,
  			     std::list<Index>& intersecting)
  {
    typedef JLBasis::Support Support;

    intersecting.clear();

    // compute support of \psi_\lambda
    const int j_lambda = lambda.j() + lambda.e();
    int k1_lambda, k2_lambda;
    support(basis, lambda, k1_lambda, k2_lambda);
        
#if 1
    // new code
    if (generators) {
      // the leftmost generator (of type phi_0) on the level j, s.th. its support intersects [a,b], fulfills
      //   2^{-j}(k+1) > a  but  2^{-j}k <= a,
      // so that ...
      const int firstk0 = std::max(1, (int) floor(ldexp(1.0,j-j_lambda)*k1_lambda));
      
      // the rightmost generator (of type phi_0) on the level j, s.th. its support intersects [a,b], fulfills
      //   2^{-j}(k-1) < b  but  2^{-j}k >= b,
      // so that ...
      const int lastk0  = std::min((1<<j)-1, (int) ceil(ldexp(1.0,j-j_lambda)*k2_lambda));
      
      // the leftmost generator (of type phi_1) on the level j, s.th. its support intersects [a,b], fulfills
      //   2^{-j}(k+1) > a  but  2^{-j}k <= a,
      // so that ...
      const int firstk1 = std::max(0, (int) floor(ldexp(1.0,j-j_lambda)*k1_lambda));
      
      // the rightmost generator (of type phi_1) on the level j, s.th. its support intersects [a,b], fulfills
      //   2^{-j}(k-1) < b  but  2^{-j}k >= b,
      // so that ...
      const int lastk1  = std::min(1<<j, (int) ceil(ldexp(1.0,j-j_lambda)*k2_lambda));

      for (int k = firstk0; k <= lastk0; k++)
  	intersecting.push_back(Index(j, 0, 0, k));
      for (int k = firstk1; k <= lastk1; k++)
  	intersecting.push_back(Index(j, 0, 1, k));
    } else {
      // the leftmost wavelet (of type psi_0) on the level j, s.th. its support intersects [a,b], fulfills
      //   2^{-j}(k+1) > a  but  2^{-j}k <= a,
      // so that ...
      const int firstk0 = std::max(1, (int) floor(ldexp(1.0,j-j_lambda)*k1_lambda));
      
      // the rightmost wavelet (of type psi_0) on the level j, s.th. its support intersects [a,b], fulfills
      //   2^{-j}(k-1) < b  but  2^{-j}k >= b,
      // so that ...
      const int lastk0  = std::min((1<<j)-1, (int) ceil(ldexp(1.0,j-j_lambda)*k2_lambda));
      
      // the leftmost wavelet (of type psi_1) on the level j, s.th. its support intersects [a,b], fulfills
      //   2^{-j}(k+1) > a  but  2^{-j}k <= a,
      // so that ...
      const int firstk1 = std::max(0, (int) floor(ldexp(1.0,j-j_lambda)*k1_lambda));
      
      // the rightmost wavelet (of type psi_1) on the level j, s.th. its support intersects [a,b], fulfills
      //   2^{-j}(k-1) < b  but  2^{-j}k >= b,
      // so that ...
      const int lastk1  = std::min(1<<j, (int) ceil(ldexp(1.0,j-j_lambda)*k2_lambda));

      for (int k = firstk0; k <= lastk0; k++)
  	intersecting.push_back(Index(j, 1, 0, k));
      for (int k = firstk1; k <= lastk1; k++)
  	intersecting.push_back(Index(j, 1, 1, k));
    }
#else
    // old code, a brute force solution
    Support supp;
    if (generators) {
      for (Index nu = first_generator(j);; ++nu) {
 	if (intersect_supports(basis, nu, j_lambda, k1_lambda, k2_lambda, supp.j, supp.k1, supp.k2))
 	  intersecting.push_back(nu);
 	if (nu == last_generator(j)) break;
      }
    } else {
      for (Index nu = first_wavelet(j);; ++nu) {
 	if (intersect_supports(basis, nu, j_lambda, k1_lambda, k2_lambda, supp.j, supp.k1, supp.k2))
 	  intersecting.push_back(nu);
 	if (nu == last_wavelet(j)) break;
      }
    }
#endif
  }
  
  bool intersect_singular_support(const JLBasis& basis,
				  const JLBasis::Index& lambda,
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
	      return false; // 2^{-m}k2_lambda <= 2^{-j}(n+1), i.e., no singular support intersection
	  }
	  j = j_lambda;
	  k1 = std::max(k1_lambda, a_help);
	  k2 = std::min(k2_lambda, b_help);
	}
      }
    }
    return true;
  }

  inline
  bool intersect_singular_support(const JLBasis& basis,
				  const JLBasis::Index& lambda,
				  const JLBasis::Index& nu,
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

  inline
  bool intersect_singular_support(const JLBasis& basis,
 				  const JLBasis::Index& lambda,
 				  const JLBasis::Index& nu)
  {
    int j, k1, k2;
    return intersect_singular_support(basis, lambda, nu, j, k1, k2);
  }

}
