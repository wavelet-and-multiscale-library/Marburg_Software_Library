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
	// phi_{j,s0},...,phi_{j,2^j-s1}             <-> 2^{j/2}phi_0(2^j*x-k), k=s0,...,s^j-s1
	// phi_{j,2^j-s1+1},...,phi_{j,2^{j+1}-s1+1} <-> 2^{j/2}phi_1(2^j*x-k), k=0,...,2^j

	const int first_half = (1<<lambda.j())-basis.get_s1();
	if (lambda.k() <= first_half) {
	  // use phi_0
	  k1 = std::max(0, lambda.k()-1);
	  k2 = std::min(1<<lambda.j(), lambda.k()+1);
	} else {
	  // use phi_1
	  k1 = std::max(0, lambda.k()-first_half-2);
	  k2 = std::min(1<<lambda.j(), lambda.k()-first_half);
	}
      }
    else // wavelet
      {
	// psi_{j,1},...,psi_{j,2^j-1}     <-> 2^{j/2}psi_0(2^j*x-k), k=1,...,2^j-1
	// psi_{j,2^j},...,psi_{j,2^{j+1}} <-> 2^{j/2}psi_1(2^j*x-k), k=0,...,2^j

	const int first_half = (1<<lambda.j())-1;
	if (lambda.k() <= first_half) {
	  // use psi_0
	  k1 = std::max(0, lambda.k()-1)*2;
	  k2 = std::min(1<<lambda.j(), lambda.k()+1)*2;
	} else {
	  // use psi_1
	  k1 = std::max(0, lambda.k()-first_half-2)*2;
	  k2 = std::min(1<<lambda.j(), lambda.k()-first_half)*2;
	}
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
    if (generators) {
      // last index of phi_0-like generators
      const int first_half = (1<<lambda.j())-basis.get_s1();

      // the leftmost generator (of type phi_0) on the level j, s.th. its support intersects [a,b], fulfills
      //   2^{-j}(k+1) > a  but  2^{-j}k <= a,
      // so that ...
      const int firstk0 = std::max(basis.DeltaLmin(), (int) floor(ldexp(1.0,j-j_lambda)*k1_lambda-1)+1);
      
      // the rightmost generator (of type phi_0) on the level j, s.th. its support intersects [a,b], fulfills
      //   2^{-j}(k-1) < b  but  2^{-j}k >= b,
      // so that ...
      const int lastk0  = std::min(first_half, (int) ceil(ldexp(1.0,j-j_lambda)*k2_lambda+1)-1);
      
      // the leftmost generator (of type phi_1) on the level j, s.th. its support intersects [a,b], fulfills
      //   2^{-j}(k-first_half) > a  but  2^{-j}(k-first_half-1) <= a,
      // so that ...
      const int firstk1 = std::max(first_half+1, (int) floor(ldexp(1.0,j-j_lambda)*k1_lambda)+first_half+1);
      
      // the rightmost generator (of type phi_1) on the level j, s.th. its support intersects [a,b], fulfills
      //   2^{-j}(k-first_half-2) < b  but  2^{-j}(k-first_half-3) >= b,
      // so that ...
      const int lastk1  = std::min(basis.DeltaRmax(j), (int) ceil(ldexp(1.0,j-j_lambda)*k2_lambda)+first_half+1);

      Index nu(j, 0, firstk0, &basis);
      for (int k = firstk0; k <= lastk0; k++, ++nu) {
 	intersect_supports(basis, nu, j_lambda, k1_lambda, k2_lambda, supp.j, supp.k1, supp.k2);
 	intersecting.push_back(std::make_pair(nu, supp));
      }
      nu = Index(j, 0, firstk1, &basis);
      for (int k = firstk1; k <= lastk1; k++, ++nu) {
 	intersect_supports(basis, nu, j_lambda, k1_lambda, k2_lambda, supp.j, supp.k1, supp.k2);
 	intersecting.push_back(std::make_pair(nu, supp));
      }
    } else {
      // last index of psi_0-like wavelets
      const int first_half = (1<<lambda.j())-1;

      // the leftmost wavelet (of type psi_0) on the level j, s.th. its support intersects [a,b], fulfills
      //   2^{-j}(k+1) > a  but  2^{-j}k <= a,
      // so that ...
      const int firstk0 = std::max(1, (int) floor(ldexp(1.0,j-j_lambda)*k1_lambda-1)+1);
      
      // the rightmost wavelet (of type psi_0) on the level j, s.th. its support intersects [a,b], fulfills
      //   2^{-j}(k-1) < b  but  2^{-j}k >= b,
      // so that ...
      const int lastk0  = std::min(first_half, (int) ceil(ldexp(1.0,j-j_lambda)*k2_lambda+1)-1);
      
      // the leftmost wavelet (of type psi_1) on the level j, s.th. its support intersects [a,b], fulfills
      //   2^{-j}(k-first_half) > a  but  2^{-j}(k-first_half-1) <= a,
      // so that ...
      const int firstk1 = std::max(first_half+1, (int) floor(ldexp(1.0,j-j_lambda)*k1_lambda)+first_half+1);
      
      // the rightmost generator (of type phi_1) on the level j, s.th. its support intersects [a,b], fulfills
      //   2^{-j}(k-first_half-2) < b  but  2^{-j}(k-first_half-3) >= b,
      // so that ...
      const int lastk1  = std::min(basis.Nablamax(j), (int) ceil(ldexp(1.0,j-j_lambda)*k2_lambda)+first_half+1);

      Index nu(j, 1, firstk0, &basis);
      for (int k = firstk0; k <= lastk0; k++, ++nu) {
 	intersect_supports(basis, nu, j_lambda, k1_lambda, k2_lambda, supp.j, supp.k1, supp.k2);
 	intersecting.push_back(std::make_pair(nu, supp));
      }
      nu = Index(j, 1, firstk1, &basis);
      for (int k = firstk1; k <= lastk1; k++, ++nu) {
 	intersect_supports(basis, nu, j_lambda, k1_lambda, k2_lambda, supp.j, supp.k1, supp.k2);
 	intersecting.push_back(std::make_pair(nu, supp));
      }
    }
#else
    // old code, a brute force solution
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
      const int first_half = (1<<lambda.j())-basis.get_s1();

      // the leftmost generator (of type phi_0) on the level j, s.th. its support intersects [a,b], fulfills
      //   2^{-j}(k+1) > a  but  2^{-j}k <= a,
      // so that ...
      const int firstk0 = std::max(basis.DeltaLmin(), (int) floor(ldexp(1.0,j-j_lambda)*k1_lambda-1)+1);
      
      // the rightmost generator (of type phi_0) on the level j, s.th. its support intersects [a,b], fulfills
      //   2^{-j}(k-1) < b  but  2^{-j}k >= b,
      // so that ...
      const int lastk0  = std::min(first_half, (int) ceil(ldexp(1.0,j-j_lambda)*k2_lambda+1)-1);

      // the leftmost generator (of type phi_1) on the level j, s.th. its support intersects [a,b], fulfills
      //   2^{-j}(k-first_half) > a  but  2^{-j}(k-first_half-1) <= a,
      // so that ...
      const int firstk1 = std::max(first_half+1, (int) floor(ldexp(1.0,j-j_lambda)*k1_lambda)+first_half+1);
      
      // the rightmost generator (of type phi_1) on the level j, s.th. its support intersects [a,b], fulfills
      //   2^{-j}(k-first_half-2) < b  but  2^{-j}(k-first_half-3) >= b,
      // so that ...
      const int lastk1  = std::min(basis.DeltaRmax(j), (int) ceil(ldexp(1.0,j-j_lambda)*k2_lambda)+first_half+1);

      for (int k = firstk0; k <= lastk0; k++)
 	intersecting.push_back(Index(j, 0, k, &basis));
      for (int k = firstk1; k <= lastk1; k++)
 	intersecting.push_back(Index(j, 0, k, &basis));
    } else {
      // last index of psi_0-like wavelets
      const int first_half = (1<<lambda.j())-1;
      
      // the leftmost wavelet (of type psi_0) on the level j, s.th. its support intersects [a,b], fulfills
      //   2^{-j}(k+1) > a  but  2^{-j}k <= a,
      // so that ...
      const int firstk0 = std::max(1, (int) floor(ldexp(1.0,j-j_lambda)*k1_lambda-1)+1);
      
      // the rightmost wavelet (of type psi_0) on the level j, s.th. its support intersects [a,b], fulfills
      //   2^{-j}(k-1) < b  but  2^{-j}k >= b,
      // so that ...
      const int lastk0  = std::min(first_half, (int) ceil(ldexp(1.0,j-j_lambda)*k2_lambda+1)-1);
      
      // the leftmost wavelet (of type psi_1) on the level j, s.th. its support intersects [a,b], fulfills
      //   2^{-j}(k-first_half) > a  but  2^{-j}(k-first_half-1) <= a,
      // so that ...
      const int firstk1 = std::max(first_half+1, (int) floor(ldexp(1.0,j-j_lambda)*k1_lambda)+first_half+1);
      
      // the rightmost generator (of type phi_1) on the level j, s.th. its support intersects [a,b], fulfills
      //   2^{-j}(k-first_half-2) < b  but  2^{-j}(k-first_half-3) >= b,
      // so that ...
      const int lastk1  = std::min(basis.Nablamax(j), (int) ceil(ldexp(1.0,j-j_lambda)*k2_lambda)+first_half+1);

      for (int k = firstk0; k <= lastk0; k++)
 	intersecting.push_back(Index(j, 1, k, &basis));
      for (int k = firstk1; k <= lastk1; k++)
 	intersecting.push_back(Index(j, 1, k, &basis));
    }
#else
    // old code, a brute force solution
    Support supp;
    if (generators) {
      for (Index nu = first_generator(&basis, j);; ++nu) {
 	if (intersect_supports(basis, nu, j_lambda, k1_lambda, k2_lambda, supp.j, supp.k1, supp.k2))
 	  intersecting.push_back(nu);
 	if (nu == last_generator(&basis, j)) break;
      }
    } else {
      for (Index nu = first_wavelet(&basis, j);; ++nu) {
 	if (intersect_supports(basis, nu, j_lambda, k1_lambda, k2_lambda, supp.j, supp.k1, supp.k2))
 	  intersecting.push_back(nu);
 	if (nu == last_wavelet(&basis, j)) break;
      }
    }
#endif
  }

//   template <int d, int dT>
//   bool intersect_singular_support(const PBasis<d,dT>& basis,
// 				  const typename PBasis<d,dT>::Index& lambda,
// 				  const int m, const int a, const int b,
// 				  int& j, int& k1, int& k2)
//   {
//     const int j_lambda = lambda.j() + lambda.e();
    
//     if (j_lambda < m) {
//       // granularity of support intersection is m
//       return intersect_supports(basis, lambda, m, a, b, j, k1, k2);
//     } else {
//       // Here j_lambda >= m, i.e., the support of psi_lambda determines the
//       // granularity of the support intersection
      
//       int k1_lambda, k2_lambda;
//       support(basis, lambda, k1_lambda, k2_lambda);
      
//       // first exclude the case that supp(psi_lambda) and 2^{-m}[a,b] do not intersect at all
//       const int b_help = (1<<(j_lambda-m)) * b;
//       if (k1_lambda >= b_help)
// 	return false; // supp(psi_lambda) is completely right of 2^{-m}[a,b]
//       else {
// 	const int a_help = (1<<(j_lambda-m)) * a;
// 	if (k2_lambda <= a_help)
// 	  return false; // supp(psi_lambda) is completely left of 2^{-m}[a,b]
// 	else {
// 	  if (k1_lambda >= a_help && k2_lambda <= b_help) {
// 	    // subset criterion, necessary for supp(psi_lambda) not hitting singsupp(2^{-m}[a,b])
	    
// 	    // determine largest integer n, such that 2^{-j}k1_lambda >= 2^{-m}n
// 	    // (we know that n < b)
// 	    const int n = (int) floor(ldexp(1.0, m-j_lambda) * k1_lambda);
// 	    if (k2_lambda <= (1<<(j_lambda-m)) * (n+1))
// 	      return false; // 2^{-m}k2_lambda <= 2^{-j}(n+1), i.e., no singuar support intersection
// 	  }
// 	  j = j_lambda;
// 	  k1 = std::max(k1_lambda, a_help);
// 	  k2 = std::min(k2_lambda, b_help);
// 	}
//       }
//     }
//     return true;
//   }

//   template <int d, int dT>
//   inline
//   bool intersect_singular_support(const PBasis<d,dT>& basis,
// 				  const typename PBasis<d,dT>::Index& lambda,
// 				  const typename PBasis<d,dT>::Index& nu,
// 				  int& j, int& k1, int& k2)
//   {
//     const int j_lambda = lambda.j() + lambda.e();
//     const int j_nu = nu.j() + nu.e();

//     if (j_lambda < j_nu)
//       return intersect_singular_support(basis, nu, lambda, j, k1, k2);
    
//     int k1_nu, k2_nu;
//     support(basis, nu, k1_nu, k2_nu);
//     return intersect_singular_support(basis, lambda, j_nu, k1_nu, k2_nu, j, k1, k2);
//   }

//   template <int d, int dT>
//   inline
//   bool intersect_singular_support(const PBasis<d,dT>& basis,
// 				  const typename PBasis<d,dT>::Index& lambda,
// 				  const typename PBasis<d,dT>::Index& nu)
//   {
//     int j, k1, k2;
//     return intersect_singular_support(basis, lambda, nu, j, k1, k2);
//   }

}
