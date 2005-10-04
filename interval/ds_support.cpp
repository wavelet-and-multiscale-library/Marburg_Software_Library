// implementation for dku_support.h

#include <algebra/matrix.h>
#include <Rd/cdf_utils.h>

namespace WaveletTL
{
  template <int d, int dT>
  void support(const DSBasis<d,dT>& basis,
	       const typename DSBasis<d,dT>::Index& lambda,
	       int& k1, int& k2)
  {
    if (lambda.e() == 0) // generator
      {
	const Matrix<double>& CLA = basis.get_CLA();
	if (lambda.k() < basis.DeltaLmin()+(int)CLA.column_dimension())
	  {
 	    // left boundary generator
	    k1 = 0;
	    k2 = CLA.row_dimension();
	  }
	else
	  {
	    const Matrix<double>& CRA = basis.get_CRA();
	    if (lambda.k() > basis.DeltaRmax(lambda.j())-(int)CRA.column_dimension())
	      {
		// right boundary generator
		k1 = (1<<lambda.j()) - CRA.row_dimension();
		k2 = 1<<lambda.j();
	      }
	    else
	      {
		k1 = lambda.k() + ell1<d>();
		k2 = lambda.k() + ell2<d>();
	      }
	  }
      }
    else // wavelet
      {
	// To determine which generators would be necessary to create the
	// wavelet in question, we mimic a reconstruct_1() call:
	
	typedef Vector<double>::size_type size_type;

	const size_type row = lambda.k();
	size_type row_j0 = row;
	size_type offset = 0;

	const int j0 = basis.j0();

	if (lambda.j() > j0) {
	  const size_type rows_top = 1<<(j0-1);
	  if (row >= rows_top) {
	    const size_type bottom = (1<<lambda.j())-(1<<(j0-1));
	    if (row >= bottom) {
	      row_j0 = row+rows_top-bottom;
	      offset = basis.Deltasize(lambda.j()+1)-basis.Deltasize(j0+1);
	    } else {
	      if ((int)row < (1<<(lambda.j()-1))) {
		row_j0 = rows_top-1;
		offset = 2*(row-rows_top)+2;
	      } else {
		row_j0 = 1<<(j0-1);
		offset = basis.Deltasize(lambda.j()+1)-basis.Deltasize(j0+1)+2*((int)row-bottom);
	      }
	    }
	  }
	}

	const SparseMatrix<double>& Mj1_t = basis.get_Mj1_t();
	const size_type kleft  = basis.DeltaLmin() + Mj1_t.get_nth_index(row_j0,0) + offset;
	const size_type kright = basis.DeltaLmin() + Mj1_t.get_nth_index(row_j0,Mj1_t.entries_in_row(row_j0)-1) + offset;
	
	int dummy;
	support(basis, typename DSBasis<d,dT>::Index(lambda.j()+1, 0, kleft, &basis), k1, dummy);
	support(basis, typename DSBasis<d,dT>::Index(lambda.j()+1, 0, kright, &basis), dummy, k2);
      }
  }

  template <int d, int dT>
  bool intersect_supports(const DSBasis<d,dT>& basis,
			  const typename DSBasis<d,dT>::Index& lambda,
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
  
  template <int d, int dT>
  bool intersect_supports(const DSBasis<d,dT>& basis,
			  const typename DSBasis<d,dT>::Index& lambda,
			  const typename DSBasis<d,dT>::Index& nu,
			  typename DSBasis<d,dT>::Support& supp)
  {
    int k1_nu, k2_nu;
    support(basis, nu, k1_nu, k2_nu);
    return intersect_supports(basis, lambda, nu.j()+nu.e(), k1_nu, k2_nu,
			      supp.j, supp.k1, supp.k2);
  }

  template <int d, int dT>
  void intersecting_wavelets(const DSBasis<d,dT>& basis,
			     const typename DSBasis<d,dT>::Index& lambda,
			     const int j, const bool generators,
			     std::list<std::pair<typename DSBasis<d,dT>::Index,
			     typename DSBasis<d,dT>::Support> >& intersecting)
  {
    typedef typename DSBasis<d,dT>::Index Index;
    typedef typename DSBasis<d,dT>::Support Support;

    intersecting.clear();

    // compute support of \psi_\lambda
    const int j_lambda = lambda.j() + lambda.e();
    int k1_lambda, k2_lambda;
    support(basis, lambda, k1_lambda, k2_lambda);
    
    // a brute force solution
    if (generators) {
      Support supp;
      for (Index nu = first_generator(&basis, j);; ++nu) {
	if (intersect_supports(basis, nu, j_lambda, k1_lambda, k2_lambda, supp.j, supp.k1, supp.k2))
	  intersecting.push_back(std::make_pair(nu, supp));
	if (nu == last_generator(&basis, j)) break;
      }
    } else {
      Support supp;
      for (Index nu = first_wavelet(&basis, j);; ++nu) {
	if (intersect_supports(basis, nu, j_lambda, k1_lambda, k2_lambda, supp.j, supp.k1, supp.k2))
	  intersecting.push_back(std::make_pair(nu, supp));
	if (nu == last_wavelet(&basis, j)) break;
      }
    }
  }

  template <int d, int dT>
  void intersecting_wavelets(const DSBasis<d,dT>& basis,
			     const typename DSBasis<d,dT>::Index& lambda,
			     const int j, const bool generators,
			     std::list<typename DSBasis<d,dT>::Index>& intersecting)
  {
    typedef typename DSBasis<d,dT>::Index Index;
    typedef typename DSBasis<d,dT>::Support Support;

    intersecting.clear();

    // compute support of \psi_\lambda
    const int j_lambda = lambda.j() + lambda.e();
    int k1_lambda, k2_lambda;
    support(basis, lambda, k1_lambda, k2_lambda);
    
    // a brute force solution
    if (generators) {
      Support supp;
      for (Index nu = first_generator(&basis, j);; ++nu) {
	if (intersect_supports(basis, nu, j_lambda, k1_lambda, k2_lambda, supp.j, supp.k1, supp.k2))
	  intersecting.push_back(nu);
	if (nu == last_generator(&basis, j)) break;
      }
    } else {
      Support supp;
      for (Index nu = first_wavelet(&basis, j);; ++nu) {
	if (intersect_supports(basis, nu, j_lambda, k1_lambda, k2_lambda, supp.j, supp.k1, supp.k2))
	  intersecting.push_back(nu);
	if (nu == last_wavelet(&basis, j)) break;
      }
    }
  }

  template <int d, int dT>
  bool intersect_singular_support(const DSBasis<d,dT>& basis,
				  const typename DSBasis<d,dT>::Index& lambda,
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

  template <int d, int dT>
  bool intersect_singular_support(const DSBasis<d,dT>& basis,
				  const typename DSBasis<d,dT>::Index& lambda,
				  const typename DSBasis<d,dT>::Index& nu,
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

  template <int d, int dT>
  void relevant_wavelets(const DSBasis<d,dT>& basis,
			 const typename DSBasis<d,dT>::Index& lambda,
			 const int j, const bool generators,
			 std::list<std::pair<typename DSBasis<d,dT>::Index, typename DSBasis<d,dT>::Support> >& relevant)
  {
    typedef typename DSBasis<d,dT>::Index Index;
    typedef typename DSBasis<d,dT>::Support Support;
    
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
}
