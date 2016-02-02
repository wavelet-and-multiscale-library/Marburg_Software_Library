// implementation for p_support.h

#include <Rd/cdf_utils.h>

#include "i_index.h"
#include "p_basis.h"

namespace WaveletTL
{
  template <int d, int dT>
  inline
  void support(const PBasis<d,dT>& basis,
	       const typename PBasis<d,dT>::Index& lambda,
	       int& k1, int& k2)
  {
#if 1
    basis.support(lambda, k1, k2);
#else
    // the following code has moved to PBasis:

    if (lambda.e() == 0) // generator
      {
	// phi_{j,k}(x) = 2^{j/2} B_{j,k-ell_1}
	k1 = std::max(0, lambda.k() + ell1<d>());
	k2 = std::min(1<<lambda.j(), lambda.k() + ell2<d>());
      }
    else // wavelet
      {
	// cf. [P, p. 125]
	if (lambda.k() < (d+dT)/2-1) {
	  // left boundary wavelet
	  k1 = 0;

 	  // hard code special case
//  	  if ( (d==3 && dT==3) &&  ((basis.get_s0()==1) && (basis.get_s1()==1)) ) {
// 	    k2 = 8;
//  	  }
//  	  else

	  k2 = 2*(d+dT)-2; // overestimate, TODO

	} else {
	  if ((1<<lambda.j())-lambda.k() <= (d+dT)/2-1) {
	    // right boundary wavelet

	    // hard code special case
// 	    if ( (d==3 && dT==3) &&  ((basis.get_s0()==1) && (basis.get_s1()==1)) ) {
// 	      k1 = (1<<(lambda.j()+1)) - 8;
//  	    }
//  	    else

	    k1 = (1<<(lambda.j()+1))-(2*(d+dT)-2); // overestimate, TODO

	    k2 = 1<<(lambda.j()+1);
	  } else {
	    // interior wavelet (CDF)
	    k1 = 2*(lambda.k()-(d+dT)/2+1);
	    k2 = k1+2*(d+dT)-2;
	  }
	}
      }
#endif
  }
  
    template <int d, int dT>
    inline
    void support(const PBasis<d,dT>& basis,
                 const int j_, const int e_, const int k_,
                 int& k1, int& k2)
    {
        basis.support(j_, e_, k_, k1, k2);
    }
    
  template <int d, int dT>
  bool intersect_supports(const PBasis<d,dT>& basis,
			  const typename PBasis<d,dT>::Index& lambda,
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
  
  template <int d, int dT>
  inline
  bool intersect_supports(const PBasis<d,dT>& basis,
			  const typename PBasis<d,dT>::Index& lambda,
			  const typename PBasis<d,dT>::Index& nu,
			  typename PBasis<d,dT>::Support& supp)
  {
    int k1_nu, k2_nu;
    support(basis, nu, k1_nu, k2_nu);
    return intersect_supports(basis, lambda, nu.j()+nu.e(), k1_nu, k2_nu,
			      supp.j, supp.k1, supp.k2);
  }

  
    template <int d, int dT>
    void intersecting_wavelets(const PBasis<d,dT>& basis,
                               const typename PBasis<d,dT>::Index& lambda,
                               const int j, const bool generators,
                               std::list<typename PBasis<d,dT>::Index>& intersecting)
    {
        typedef typename PBasis<d,dT>::Index Index;
        typedef typename PBasis<d,dT>::Support Support;
        intersecting.clear();
        // compute support of \psi_\lambda
        const int j_lambda = lambda.j() + lambda.e();
        int k1_lambda, k2_lambda;
        support(basis, lambda, k1_lambda, k2_lambda);
#if 1
        if (generators) 
        {
            // Generator on level j has support 2^{-j}[k+ell1 ; k+ell2]
            // the leftmost generator on the level j, s.th. its support intersects [2^{-j_lambda}k1; 2^{-j_lambda}k2], fulfills
            //   2^{-j}(k+ell2-1) <= 2^{-j_lambda}k1 < 2^{-j}(k+ell2)
            const int firstk = std::max(basis.DeltaLmin(),
                                        (j < j_lambda) ? (k1_lambda / (1<<(j_lambda-j))-ell2<d>()+1 )
                                                       : ((k1_lambda << (j-j_lambda)) -ell2<d>()+1) );
            //assert (firstk == std::max(basis.DeltaLmin(), (int) floor(ldexp(1.0,j-j_lambda)*k1_lambda-ell2<d>())+1) );
            //const int firstk = std::max(basis.DeltaLmin(), (int) floor(ldexp(1.0,j-j_lambda)*k1_lambda-ell2<d>())+1);
            // the rightmost generator on the level j, s.th. its support intersects [2^{-j_lambda}k1; 2^{-j_lambda}k2], fulfills
            //   2^{-j}(k+ell1) < 2^{-j_lambda}k2 <= 2^{-j}(k+1+ell1)
            const int lastk = std::min(basis.DeltaRmax(j), (j < j_lambda) ? (k2_lambda / (1<<(j_lambda-j)) + ( ((k2_lambda % (1<<(j_lambda-j))) == 0)? 0 : 1 ) -1-ell1<d>() )
                                                                          : ((k2_lambda << (j-j_lambda)) -1-ell1<d>() )   );
            //assert (lastk == std::min(basis.DeltaRmax(j), (int) ceil(ldexp(1.0,j-j_lambda)*k2_lambda-ell1<d>())-1));
            //const int lastk  = std::min(basis.DeltaRmax(j), (int) ceil(ldexp(1.0,j-j_lambda)*k2_lambda-ell1<d>())-1);
            for (int k = firstk; k <= lastk; k++)
                intersecting.push_back(Index(j, 0, k, &basis));
        } 
        else 
        {
            // Inner wavelet on level j has support 2^{-j}[k-(d+dt)/2+1; k +(d+dt)/2]
            // the leftmost inner wavelet on the level j, s.th. its support intersects [2^{-j_lambda}k1; 2^{-j_lambda}k2], fulfills
            //   2^{-j}(k+(d+dt)/2-1) <= 2^{-j_lambda}k1 < 2^{-j}(k+(d+dt)/2)
            // compute:   2^{-j}n <= 2^{-j_lambda}k1 < 2^{-j}(n+1)
            int n = (j < j_lambda) ? (k1_lambda / (1<<(j_lambda-j)))
                                   : (k1_lambda << (j-j_lambda));
            int mink = d+dT-2; //used as temporary storage
            if (n > mink)
            {
                // no left boundary wavelets intersect
                mink = n - mink/2; // = n+1-(d+dt)/2;
            }
            else
            {
                if (n < mink)
                {
                    // all left boundary wavelets intersect
                    mink = 0;
                }
                else
                {
                    // left boundary wavelets may intersect 
                    if (basis.get_s0() == d-1)
                    {
                        // left boundary wavelets have the same support as the leftmost inner wavelet
                        mink = 0;
                    }
                    else
                    {
                        // left boundary wavelets do not intersect
                        mink /= 2;
                        // newmink = n+1-(d+dt)/2;  // it holds that n == d+dt-2
                    }
                }
            }
            /*
            assert (mink >= (ldexp(1.0,j-j_lambda)*k1_lambda < d+dT-1 // overestimate, TODO!
			  ? 0
			  : (int) floor(ldexp(1.0,j-j_lambda)*k1_lambda-(d+dT)/2)+1));
            */
            
            // Inner wavelet on level j has support 2^{-j}[k-(d+dt)/2+1; k +(d+dt)/2]
            // the rightmost inner wavelet on the level j, s.th. its support intersects [2^{-j_lambda}k1; 2^{-j_lambda}k2], fulfills
            //        2^{-j}(k-(d+dt)/2+1) < 2^{-j_lambda}k2 <= 2^{-j}(k-(d+dt)/2+2)
            // biggest n such that 2^{-j}n < 2^{-j_lambda}k2_lambda <= 2^{-j}(n+1)
            n = (j < j_lambda) ? (k2_lambda / (1<<(j_lambda-j)) + (((k2_lambda % (1<<j_lambda-j)) == 0)?0:1) ) -1 
                               : (k2_lambda << (j-j_lambda)) -1;
            int maxk = (1<<j) -(d+dT)+1; // used as temporary value;
            if (n < maxk)
            {
                // no right boundary wavelets intersect
                maxk = n-1+(d+dT)/2;
            }
            else
            {
                if (n > maxk)
                {
                    // all right boundary wavelets intersect
                    maxk = (1<<j)-1;
                }
                else
                {
                    // right boundary wavelets may intersect 
                    if (basis.get_s1() == d-1)
                    {
                        // right boundary wavelets have the same support as the rightmost inner wavelet
                        maxk = (1<<j)-1;
                    }
                    else
                    {
                        // right boundary wavelets do not intersect
                        maxk = n-1+(d+dT)/2; // with n = newmaxk, but formula doesnt get cheaper without introducing another temp variable
                    }
                }
            }
            /*
            assert (maxk <= (ldexp(1.0,-j_lambda)*k2_lambda > 1 - ldexp(1.0,-j)*(d+dT-1) // overestimate, TODO!
			 ? (1<<j)-1
			 : (int) ceil(ldexp(1.0,j-j_lambda)*k2_lambda+(d+dT)/2)-2) );
             */
            for (int k = mink; k <= maxk; k++)
                intersecting.push_back(Index(j, 1, k, &basis));
        }
#else
        // old code, a brute force solution
        Support supp;
        if (generators) 
        {
            for (Index nu = first_generator(&basis, j);; ++nu) 
            {
                if (intersect_supports(basis, nu, j_lambda, k1_lambda, k2_lambda, supp.j, supp.k1, supp.k2))
                    intersecting.push_back(nu);
                if (nu == last_generator(&basis, j)) break;
            }
        } 
        else 
        {
            for (Index nu = first_wavelet(&basis, j);; ++nu) 
            {
                if (intersect_supports(basis, nu, j_lambda, k1_lambda, k2_lambda, supp.j, supp.k1, supp.k2))
                    intersecting.push_back(nu);
                if (nu == last_wavelet(&basis, j)) break;
            }
        }
#endif
    }
    
    template <int d, int dT>
    void intersecting_wavelets(const PBasis<d,dT>& basis,
                               const typename PBasis<d,dT>::Index& lambda,
                               const int j, const bool generators,
                               std::list<std::pair<typename PBasis<d,dT>::Index, typename PBasis<d,dT>::Support> >& intersecting)
    {

        typedef typename PBasis<d,dT>::Index Index;
        typedef typename PBasis<d,dT>::Support Support;
        intersecting.clear();
        // compute support of \psi_\lambda
        const int j_lambda = lambda.j() + lambda.e();
        int k1_lambda, k2_lambda;
        support(basis, lambda, k1_lambda, k2_lambda);
#if 1
        Support supp;
        if (generators) 
        {
            // Generator on level j has support 2^{-j}[k+ell1 ; k+ell2]
            // the leftmost generator on the level j, s.th. its support intersects [2^{-j_lambda}k1; 2^{-j_lambda}k2], fulfills
            //   2^{-j}(k+ell2-1) <= 2^{-j_lambda}k1 < 2^{-j}(k+ell2)
            const int firstk = std::max(basis.DeltaLmin(),
                                        (j < j_lambda) ? (k1_lambda / (1<<(j_lambda-j))-ell2<d>()+1 )
                                                       : ((k1_lambda << (j-j_lambda)) -ell2<d>()+1) );
            //assert (firstk == std::max(basis.DeltaLmin(), (int) floor(ldexp(1.0,j-j_lambda)*k1_lambda-ell2<d>())+1) );
            //const int firstk = std::max(basis.DeltaLmin(), (int) floor(ldexp(1.0,j-j_lambda)*k1_lambda-ell2<d>())+1);
            // the rightmost generator on the level j, s.th. its support intersects [2^{-j_lambda}k1; 2^{-j_lambda}k2], fulfills
            //   2^{-j}(k+ell1) < 2^{-j_lambda}k2 <= 2^{-j}(k+1+ell1)
            const int lastk = std::min(basis.DeltaRmax(j), (j < j_lambda) ? (k2_lambda / (1<<(j_lambda-j)) + ( ((k2_lambda % (1<<(j_lambda-j))) == 0)? 0 : 1 ) -1-ell1<d>() )
                                                                          : ((k2_lambda << (j-j_lambda)) -1-ell1<d>() )   );
            //assert (lastk == std::min(basis.DeltaRmax(j), (int) ceil(ldexp(1.0,j-j_lambda)*k2_lambda-ell1<d>())-1));
            //const int lastk  = std::min(basis.DeltaRmax(j), (int) ceil(ldexp(1.0,j-j_lambda)*k2_lambda-ell1<d>())-1);
// INEFFICIENT:
            Index nu(j, 0, firstk, &basis);
            for (int k = firstk; k <= lastk; k++, ++nu) 
            {
                intersect_supports(basis, nu, j_lambda, k1_lambda, k2_lambda, supp.j, supp.k1, supp.k2);
                intersecting.push_back(std::make_pair(nu, supp));
            }
        }
        else 
        {
            // Inner wavelet on level j has support 2^{-j}[k-(d+dt)/2+1; k +(d+dt)/2]
            // the leftmost inner wavelet on the level j, s.th. its support intersects [2^{-j_lambda}k1; 2^{-j_lambda}k2], fulfills
            //   2^{-j}(k+(d+dt)/2-1) <= 2^{-j_lambda}k1 < 2^{-j}(k+(d+dt)/2)
            // compute:   2^{-j}n <= 2^{-j_lambda}k1 < 2^{-j}(n+1)
            int n = (j < j_lambda) ? (k1_lambda / (1<<(j_lambda-j)))
                                   : (k1_lambda << (j-j_lambda));
            int mink = d+dT-2; //used as temporary storage
            if (n > mink)
            {
                // no left boundary wavelets intersect
                mink = n - mink/2; // = n+1-(d+dt)/2;
            }
            else
            {
                if (n < mink)
                {
                    // all left boundary wavelets intersect
                    mink = 0;
                }
                else
                {
                    // left boundary wavelets may intersect 
                    if (basis.get_s0() == d-1)
                    {
                        // left boundary wavelets have the same support as the leftmost inner wavelet
                        mink = 0;
                    }
                    else
                    {
                        // left boundary wavelets do not intersect
                        mink /= 2;
                        // newmink = n+1-(d+dt)/2;  // it holds that n == d+dt-2
                    }
                }
            }
            /*
            assert (mink >= (ldexp(1.0,j-j_lambda)*k1_lambda < d+dT-1 // overestimate, TODO!
			  ? 0
			  : (int) floor(ldexp(1.0,j-j_lambda)*k1_lambda-(d+dT)/2)+1));
            */
            
            // Inner wavelet on level j has support 2^{-j}[k-(d+dt)/2+1; k +(d+dt)/2]
            // the rightmost inner wavelet on the level j, s.th. its support intersects [2^{-j_lambda}k1; 2^{-j_lambda}k2], fulfills
            //        2^{-j}(k-(d+dt)/2+1) < 2^{-j_lambda}k2 <= 2^{-j}(k-(d+dt)/2+2)
            // biggest n such that 2^{-j}n < 2^{-j_lambda}k2_lambda <= 2^{-j}(n+1)
            n = (j < j_lambda) ? (k2_lambda / (1<<(j_lambda-j)) + (((k2_lambda % (1<<j_lambda-j)) == 0)?0:1) ) -1 
                               : (k2_lambda << (j-j_lambda)) -1;
            int maxk = (1<<j) -(d+dT)+1; // used as temporary value;
            if (n < maxk)
            {
                // no right boundary wavelets intersect
                maxk = n-1+(d+dT)/2;
            }
            else
            {
                if (n > maxk)
                {
                    // all right boundary wavelets intersect
                    maxk = (1<<j)-1;
                }
                else
                {
                    // right boundary wavelets may intersect 
                    if (basis.get_s1() == d-1)
                    {
                        // right boundary wavelets have the same support as the rightmost inner wavelet
                        maxk = (1<<j)-1;
                    }
                    else
                    {
                        // right boundary wavelets do not intersect
                        maxk = n-1+(d+dT)/2; // with n = newmaxk, but formula doesnt get cheaper without introducing another temp variable
                    }
                }
            }
            /*
            assert (maxk <= (ldexp(1.0,-j_lambda)*k2_lambda > 1 - ldexp(1.0,-j)*(d+dT-1) // overestimate, TODO!
			 ? (1<<j)-1
			 : (int) ceil(ldexp(1.0,j-j_lambda)*k2_lambda+(d+dT)/2)-2) );
            */
// INEFFICIENT            
            Index nu(j, 1, mink, &basis);
            for (int k = mink; k <= maxk; k++, ++nu) 
            {
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


    

    template <int d, int dT>
    bool intersect_singular_support(const PBasis<d,dT>& basis,
                                    const typename PBasis<d,dT>::Index& lambda,
                                    const int m, const int a, const int b,
                                    int& j, int& k1, int& k2)
    {
        const int j_lambda = lambda.j() + lambda.e();
        if (j_lambda < m) 
        {
            // granularity of support intersection is m
            return intersect_supports(basis, lambda, m, a, b, j, k1, k2);
        } 
        else 
        {
            // Here j_lambda >= m, i.e., the support of psi_lambda determines the
            // granularity of the support intersection
            int k1_lambda, k2_lambda;
            support(basis, lambda, k1_lambda, k2_lambda);
            // first exclude the case that supp(psi_lambda) and 2^{-m}[a,b] do not intersect at all
            const int b_help = (1<<(j_lambda-m)) * b;
            if (k1_lambda >= b_help)
                return false; // supp(psi_lambda) is completely right of 2^{-m}[a,b]
            else 
            {
                const int a_help = (1<<(j_lambda-m)) * a;
                if (k2_lambda <= a_help)
                    return false; // supp(psi_lambda) is completely left of 2^{-m}[a,b]
                else 
                {
                    if (k1_lambda >= a_help && k2_lambda <= b_help) 
                    {
                        // subset criterion, necessary for supp(psi_lambda) not hitting singsupp(2^{-m}[a,b])
                        // smallest n such that 2^{-m}n <= 2^{-j}k1_lambda < 2^{-m}(n+1)
                        const int n = (m < j_lambda) ? (k1_lambda / (1<<(j_lambda-m)))
                                                     : (k1_lambda << (m-j_lambda));
                        // (we know that n < b)
                        // assert (n == (int) floor(ldexp(1.0, m-j_lambda) * k1_lambda));
                        // determine largest integer n, such that 2^{-j}k1_lambda >= 2^{-m}n
                        //const int n = (int) floor(ldexp(1.0, m-j_lambda) * k1_lambda);
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

    template <int d, int dT>
    inline
    bool intersect_singular_support(const PBasis<d,dT>& basis,
                                    const typename PBasis<d,dT>::Index& lambda,
                                    const typename PBasis<d,dT>::Index& nu,
                                    int& j, int& k1, int& k2)
    {
#if 1
        return intersect_singular_support(basis, lambda.j(), lambda.e(), lambda.k(), nu.j(), nu.e(), nu.k(), j, k1, k2);
#else             
        const int j_lambda = lambda.j() + lambda.e();
        const int j_nu = nu.j() + nu.e();
        if (j_lambda < j_nu)
            return intersect_singular_support(basis, nu, lambda, j, k1, k2);
        int k1_nu, k2_nu;
        support(basis, nu, k1_nu, k2_nu);
        return intersect_singular_support(basis, lambda, j_nu, k1_nu, k2_nu, j, k1, k2);
#endif
    }

    template <int d, int dT>
    inline
    bool intersect_singular_support(const PBasis<d,dT>& basis,
                                    const typename PBasis<d,dT>::Index& lambda,
                                    const typename PBasis<d,dT>::Index& nu)
    {
        int j, k1, k2;
        return intersect_singular_support(basis, lambda.j(), lambda.e(), lambda.k(), nu.j(), nu.e(), nu.k(), j, k1, k2);
    }
  
    template <int d, int dT>
    bool intersect_singular_support(const PBasis<d,dT>& basis,
                                    const int lamj, const int lame, const int lamk,
                                    const int muj, const int mue, const int muk,
                                    int& j, int& k1, int& k2)
    {
        const int j_lambda = lamj+lame;
        const int j_mu = muj+mue;
        if (j_lambda < j_mu) 
        {
            return intersect_singular_support(basis, muj, mue, muk, lamj,lame,lamk, j, k1, k2);
        } 
        else 
        {
            // Here j_lambda >= j_mu, i.e., the support of psi_lambda determines the
            // granularity of the support intersection
            int k1_lambda, k2_lambda, k1_mu, k2_mu;
            support(basis, lamj,lame,lamk, k1_lambda, k2_lambda);
            support(basis, muj,mue,muk, k1_mu, k2_mu);
            // first exclude the case that supp(psi_lambda) and 2^{-m}[a,b] do not intersect at all
            const int r_help = (k2_mu << (j_lambda-j_mu));
            if (k1_lambda >= r_help)
                return false; // supp(psi_lambda) is completely right of 2^{-j_mu}[k1_mu,k2_mu]
            else 
            {
                const int l_help = (k1_mu << (j_lambda-j_mu));
                if (k2_lambda <= l_help)
                    return false; // supp(psi_lambda) is completely left of 2^{-j_mu}[k1_mu,k2_mu]
                else 
                {
                    if (k1_lambda >= l_help && k2_lambda <= r_help) 
                    {
                        // subset criterion, necessary for supp(psi_lambda) not hitting singsupp(2^{-j_mu}[k1_mu,k2_mu])
                        // determine largest integer n, such that 2^{-j_lambda}k1_lambda >= 2^{-j_mu}n
                        // (we know that k1_mu <= n < k2_mu)
                        const int n = (k1_lambda / (1<< (j_lambda-j_mu)));
                        if (k2_lambda <= ((n+1)<<(j_lambda-j_mu)))
                        {
                            return false; // 2^{-j_mu}k2_lambda <= 2^{-j_lambda}(n+1), i.e., no singuar support intersection
                        }
                    }
                    j = j_lambda;
                    k1 = std::max(k1_lambda, l_help);
                    k2 = std::min(k2_lambda, r_help);
                }
            }
        }
        return true;
    }
    
    template <int d, int dT>
    void get_intersecting_wavelets_on_level(const PBasis<d,dT>& basis,
                                            const typename PBasis<d,dT>::Index& lambda,
                                            const int j, 
                                            const bool generators,
                                            int& mink, 
                                            int& maxk)
    { 
        //  the following code is equivalent to
        //get_intersecting_wavelets_on_level(basis,lambda.j(), lambda.e(), lambda.k(), j, generators, mink, maxk);
        // however, this call is slower!
        typedef typename PBasis<d,dT>::Index Index;
        typedef typename PBasis<d,dT>::Support Support;
        // compute support of \psi_\lambda
        const int j_lambda = lambda.j() + lambda.e();
        int k1_lambda, k2_lambda;
        support(basis, lambda, k1_lambda, k2_lambda);
        if (generators) 
        {
            // smallest n such that 2^{-j}n <= 2^{-j_lambda}k1_lambda < 2^{-j}(n+1)
            mink = (j < j_lambda) ? (k1_lambda / (1<<(j_lambda-j)))
                                  : (k1_lambda << (j-j_lambda));
            // for generators: n = k-1+ell2<d>()
            mink += 1 - (d - d/2); 
            mink = std::max(basis.DeltaLmin(), mink);
            // old code
            // mink = std::max(basis.DeltaLmin(), (int) floor(ldexp(1.0,j-j_lambda)*k1_lambda-(d-d/2))+1);
            
            // biggest n such that 2^{-j}n < 2^{-j_lambda}k2_lambda <= 2^{-j}(n+1)
            maxk = (j < j_lambda) ? (k2_lambda / (1<<(j_lambda-j)) + ( ((k2_lambda % (1<<(j_lambda-j))) == 0)? 0 : 1) )
                                   : (k2_lambda << (j-j_lambda));
            maxk += d/2 -1;
            maxk = std::min(basis.DeltaRmax(j), maxk);
            // old code
            //maxk = std::min(basis.DeltaRmax(j), (int) ceil(ldexp(1.0,j-j_lambda)*k2_lambda-(-d/2))-1);
        } 
        else 
        {
            //old code relies on an overestimate of the suppor size of the boundary wavelets
            /*
            mink = (ldexp(1.0,j-j_lambda)*k1_lambda < d+dT-1 // overestimate, TODO!
                    ? 0
                    : (int) floor(ldexp(1.0,j-j_lambda)*k1_lambda-(d+dT)/2)+1);
            maxk = (ldexp(1.0,-j_lambda)*k2_lambda > 1 - ldexp(1.0,-j)*(d+dT-1) // overestimate, TODO!
                    ? (1<<j)-1
                    : (int) ceil(ldexp(1.0,j-j_lambda)*k2_lambda+(d+dT)/2)-2);
             */
            // new code:
            // find first inner wavelet that intersects
            // smallest n such that 2^{-j}n <= 2^{-j_lambda}k1_lambda < 2^{-j}(n+1)
            int n = (j < j_lambda) ? (k1_lambda / (1<<(j_lambda-j)))
                                   : (k1_lambda << (j-j_lambda));
            /*
            // the following code is equivalent, looks better but is a bit slower (left here for understanding)
            if (basis.get_s0() == d-1)
            {
                mink = d+dT-2;
                mink = (n > mink)? (n - mink/2) : 0;
            }
            else
            {
                mink = d+dT-2;
                mink = (n < mink)? 0 : (n - mink/2);
            }
            */
            mink = d+dT-2; //used as temporary storage
            if (n > mink)
            {
                // no left boundary wavelets intersect
                mink = n - mink/2; // = n+1-(d+dt)/2;
            }
            else
            {
                if (n < mink)
                {
                    // all left boundary wavelets intersect
                    mink = 0;
                }
                else
                {
                    // left boundary wavelets may intersect 
                    if (basis.get_s0() == d-1)
                    {
                        // left boundary wavelets have the same support as the leftmost inner wavelet
                        mink = 0;
                    }
                    else
                    {
                        // left boundary wavelets do not intersect
                        mink /= 2;
                        // newmink = n+1-(d+dt)/2;  // with n = newmink
                    }
                }
            }
            // find last inner wavelet that intersects
            // biggest n such that 2^{-j}n < 2^{-j_lambda}k2_lambda <= 2^{-j}(n+1)
            n = (j < j_lambda) ? (k2_lambda / (1<<(j_lambda-j)) + (((k2_lambda % (1<<j_lambda-j)) == 0)?0:1) ) -1 
                               : (k2_lambda << (j-j_lambda)) -1;
            /*
            // again some clear, but slightly slower code:
            maxk = (1<<j) -(d+dT)+1; // used as temporary value;
            if (basis.get_s1() == d-1)
            {
                maxk = (n < maxk) ?  (n-1+(d+dT)/2) :( (1<<j) -1);
            }
            else
            {
                maxk = (n > maxk) ?  ( (1<<j) -1) : (n-1+(d+dT)/2);
            }
             
            */
            
            maxk = (1<<j) -(d+dT)+1; // used as temporary value;
            if (n < maxk)
            {
                // no right boundary wavelets intersect
                maxk = n-1+(d+dT)/2;
            }
            else
            {
                if (n > maxk)
                {
                    // all right boundary wavelets intersect
                    maxk = (1<<j)-1;
                }
                else
                {
                    // right boundary wavelets may intersect 
                    if (basis.get_s1() == d-1)
                    {
                        // right boundary wavelets have the same support as the rightmost inner wavelet
                        maxk = (1<<j)-1;
                    }
                    else
                    {
                        // right boundary wavelets do not intersect
                        maxk = n-1+(d+dT)/2; // with n = newmaxk, but formula doesnt get cheaper without introducing another temp variable
                    }
                }
            }
        }
    }
    
    template <int d, int dT>
    void get_intersecting_wavelets_on_level(const PBasis<d,dT>& basis,
                                            const int lamj, const int lame, const int lamk,
                                            const int j, 
                                            const bool generators,
                                            int& mink, 
                                            int& maxk)
    {
        typedef typename PBasis<d,dT>::Index Index;
        typedef typename PBasis<d,dT>::Support Support;
        // compute support of \psi_\lambda
        const int j_lambda = lamj+lame;
        int k1_lambda, k2_lambda;
        support(basis, lamj, lame, lamk, k1_lambda, k2_lambda);
        if (generators) 
        {
            // the leftmost generator on the level j, s.th. its support intersects [a,b], fulfills
            //   2^{-j}(k+ell2) > a  but  2^{-j}(k-1+ell2) <= a,
            // so that ...
            
            // smallest n such that 2^{-j}n <= 2^{-j_lambda}k1_lambda < 2^{-j}(n+1)
            mink = (j < j_lambda) ? (k1_lambda / (1<<(j_lambda-j)))
                                  : (k1_lambda << (j-j_lambda));
            // for generators: n = k-1+ell2<d>()
            mink += 1 - (d - d/2); 
            //assert(mink == (int) floor(ldexp(1.0,j-j_lambda)*k1_lambda-(d-d/2))+1);
                
            mink = std::max(basis.DeltaLmin(), mink);
            //old code
            //mink = std::max(basis.DeltaLmin(), (int) floor(ldexp(1.0,j-j_lambda)*k1_lambda-(d-d/2))+1);
            
            // the rightmost generator on the level j, s.th. its support intersects [a,b], fulfills
            //   2^{-j}(k+ell1) < b  but  2^{-j}(k+1+ell1) >= b,
            // so that ...
            maxk = (j < j_lambda) ? (k2_lambda / (1<<(j_lambda-j)) + ( ((k2_lambda % (1<<(j_lambda-j))) == 0)? 0 : 1) )
                                  : (k2_lambda << (j-j_lambda));
            maxk += d/2 -1;
            //assert (maxk == (int) ceil(ldexp(1.0,j-j_lambda)*k2_lambda-(-d/2))-1);
            maxk = std::min(basis.DeltaRmax(j), maxk);
            //old code
            //maxk = std::min(basis.DeltaRmax(j), (int) ceil(ldexp(1.0,j-j_lambda)*k2_lambda-(-d/2))-1);
        }
        else 
        {
            // find first inner wavelet that intersects
            // smallest n such that 2^{-j}n <= 2^{-j_lambda}k1_lambda < 2^{-j}(n+1)
            int n = (j < j_lambda) ? (k1_lambda / (1<<(j_lambda-j)))
                                   : (k1_lambda << (j-j_lambda));
            /*
            // the following code is equivalent, looks better but is a bit slower (left here for understanding)
            if (basis.get_s0() == d-1)
            {
                mink = d+dT-2;
                mink = (n > mink)? (n - mink/2) : 0;
            }
            else
            {
                mink = d+dT-2;
                mink = (n < mink)? 0 : (n - mink/2);
            }
            */
            mink = d+dT-2; //used as temporary storage
            if (n > mink)
            {
                // no left boundary wavelets intersect
                mink = n - mink/2; // = n+1-(d+dt)/2;
            }
            else
            {
                if (n < mink)
                {
                    // all left boundary wavelets intersect
                    mink = 0;
                }
                else
                {
                    // left boundary wavelets may intersect 
                    if (basis.get_s0() == d-1)
                    {
                        // left boundary wavelets have the same support as the leftmost inner wavelet
                        mink = 0;
                    }
                    else
                    {
                        // left boundary wavelets do not intersect
                        mink /= 2;
                        // newmink = n+1-(d+dt)/2;  // with n = newmink
                    }
                }
            }
            
            
            // find last inner wavelet that intersects
            // biggest n such that 2^{-j}n < 2^{-j_lambda}k2_lambda <= 2^{-j}(n+1)
            n = (j < j_lambda) ? (k2_lambda / (1<<(j_lambda-j)) + (((k2_lambda % (1<<j_lambda-j)) == 0)?0:1) ) -1 
                               : (k2_lambda << (j-j_lambda)) -1;
            /*
            // again some clear, but slightly slower code:
            maxk = (1<<j) -(d+dT)+1; // used as temporary value;
            if (basis.get_s1() == d-1)
            {
                maxk = (n < maxk) ?  (n-1+(d+dT)/2) :( (1<<j) -1);
            }
            else
            {
                maxk = (n > maxk) ?  ( (1<<j) -1) : (n-1+(d+dT)/2);
            }
             
            */
            
            maxk = (1<<j) -(d+dT)+1; // used as temporary value;
            if (n < maxk)
            {
                // no right boundary wavelets intersect
                maxk = n-1+(d+dT)/2;
            }
            else
            {
                if (n > maxk)
                {
                    // all right boundary wavelets intersect
                    maxk = (1<<j)-1;
                }
                else
                {
                    // right boundary wavelets may intersect 
                    if (basis.get_s1() == d-1)
                    {
                        // right boundary wavelets have the same support as the rightmost inner wavelet
                        maxk = (1<<j)-1;
                    }
                    else
                    {
                        // right boundary wavelets do not intersect
                        maxk = n-1+(d+dT)/2; // with n = newmaxk, but formula doesnt get cheaper without introducing another temp variable
                    }
                }
            }
            
            // old code is much shorter, but relies on the fact that the outmost inner wavelets have the same support than the respective boundary wavelets
            /*
            // if the left boundary wavelets are not involved, then
            // the rightmost wavelet on the level j, s.th. its support intersects [a,b], fulfills
            //   2^{-j}(k+(d+dT)/2) > a  but  2^{-j}(k-1+(d+dT)/2) <= a,
            // so that ...
            mink = (ldexp(1.0,j-j_lambda)*k1_lambda < d+dT-1 // overestimate, TODO!
                    ? 0
                    : (int) floor(ldexp(1.0,j-j_lambda)*k1_lambda-(d+dT)/2)+1);
            // if the right boundary wavelets are not involved, then
            // the leftmost wavelet on the level j, s.th. its support intersects [a,b], fulfills
            //   2^{-j}(k-(d+dT)/2+1) < b  but  2^{-j}(k-(d+dT)/2+2) >= a,
            // so that ...
            maxk = (ldexp(1.0,-j_lambda)*k2_lambda > 1 - ldexp(1.0,-j)*(d+dT-1) // overestimate, TODO!
                    ? (1<<j)-1
                    : (int) ceil(ldexp(1.0,j-j_lambda)*k2_lambda+(d+dT)/2)-2);
             */
        }
    }
    
    template <int d, int dT>
    void get_intersecting_wavelets_on_level(const PBasis<d,dT>& basis,
                                            const int lamj, const int lame, const int k1_lambda, const int k2_lambda,
                                            const int j, 
                                            const bool generators,
                                            int& mink, 
                                            int& maxk)
    {
        typedef typename PBasis<d,dT>::Index Index;
        typedef typename PBasis<d,dT>::Support Support;
        // compute support of \psi_\lambda
        const int j_lambda = lamj+lame;
        //int k1_lambda, k2_lambda;
        //support(basis, lamj, lame, lamk, k1_lambda, k2_lambda);
        if (generators) 
        {
            // the leftmost generator on the level j, s.th. its support intersects [a,b], fulfills
            //   2^{-j}(k+ell2) > a  but  2^{-j}(k-1+ell2) <= a,
            // so that ...
            
            // smallest n such that 2^{-j}n <= 2^{-j_lambda}k1_lambda < 2^{-j}(n+1)
            mink = (j < j_lambda) ? (k1_lambda / (1<<(j_lambda-j)))
                                  : (k1_lambda << (j-j_lambda));
            // for generators: n = k-1+ell2<d>()
            mink += 1 - (d - d/2); 
            mink = std::max(basis.DeltaLmin(), mink);
            // the rightmost generator on the level j, s.th. its support intersects [a,b], fulfills
            //   2^{-j}(k+ell1) < b  but  2^{-j}(k+1+ell1) >= b,
            // so that ...
            maxk = (j < j_lambda) ? (k2_lambda / (1<<(j_lambda-j)) + ( ((k2_lambda % (1<<(j_lambda-j))) == 0)? 0 : 1) )
                                  : (k2_lambda << (j-j_lambda));
            maxk += d/2 -1;
            //assert (maxk == (int) ceil(ldexp(1.0,j-j_lambda)*k2_lambda-(-d/2))-1);
            maxk = std::min(basis.DeltaRmax(j), maxk);
            //old code
            //maxk = std::min(basis.DeltaRmax(j), (int) ceil(ldexp(1.0,j-j_lambda)*k2_lambda-(-d/2))-1);
        }
        else 
        {
            // find first inner wavelet that intersects
            // smallest n such that 2^{-j}n <= 2^{-j_lambda}k1_lambda < 2^{-j}(n+1)
            int n = (j < j_lambda) ? (k1_lambda / (1<<(j_lambda-j)))
                                   : (k1_lambda << (j-j_lambda));
            mink = d+dT-2; //used as temporary storage
            if (n > mink)
            {
                // no left boundary wavelets intersect
                mink = n - mink/2; // = n+1-(d+dt)/2;
            }
            else
            {
                if (n < mink)
                {
                    // all left boundary wavelets intersect
                    mink = 0;
                }
                else
                {
                    // left boundary wavelets may intersect 
                    if (basis.get_s0() == d-1)
                    {
                        // left boundary wavelets have the same support as the leftmost inner wavelet
                        mink = 0;
                    }
                    else
                    {
                        // left boundary wavelets do not intersect
                        mink /= 2;
                        // newmink = n+1-(d+dt)/2;  // with n = newmink
                    }
                }
            }
            // find last inner wavelet that intersects
            // biggest n such that 2^{-j}n < 2^{-j_lambda}k2_lambda <= 2^{-j}(n+1)
            n = (j < j_lambda) ? (k2_lambda / (1<<(j_lambda-j)) + (((k2_lambda % (1<<j_lambda-j)) == 0)?0:1) ) -1 
                               : (k2_lambda << (j-j_lambda)) -1;
            maxk = (1<<j) -(d+dT)+1; // used as temporary value;
            if (n < maxk)
            {
                // no right boundary wavelets intersect
                maxk = n-1+(d+dT)/2;
            }
            else
            {
                if (n > maxk)
                {
                    // all right boundary wavelets intersect
                    maxk = (1<<j)-1;
                }
                else
                {
                    // right boundary wavelets may intersect 
                    if (basis.get_s1() == d-1)
                    {
                        // right boundary wavelets have the same support as the rightmost inner wavelet
                        maxk = (1<<j)-1;
                    }
                    else
                    {
                        // right boundary wavelets do not intersect
                        maxk = n-1+(d+dT)/2; // with n = newmaxk, but formula doesnt get cheaper without introducing another temp variable
                    }
                }
            }
        }
    }
    
}
