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
    if (jjk1 < jmb) {
      const int jma = (1<<(j-m)) * a;
      const int jjk2 = (1<<(j-j_lambda)) * k2_lambda;
      if (jma < jjk2) {
	k1 = std::max(jjk1, jma);
	k2 = std::min(jjk2, jmb);
	return true;
      }
    }
    return false;
  }
  
  template <int d, int dT>
  bool intersect_supports(const DSBasis<d,dT>& basis,
			  const typename DSBasis<d,dT>::Index& lambda,
			  const typename DSBasis<d,dT>::Index& nu,
			  int& j, int& k1, int& k2)
  {
    int k1_nu, k2_nu;
    support(basis, nu, k1_nu, k2_nu);
    return intersect_supports(basis, lambda, nu.j()+nu.e(), k1_nu, k2_nu, j, k1, k2);
  }

  template <int d, int dT>
  void intersecting_wavelets(const DSBasis<d,dT>& basis,
			     const typename DSBasis<d,dT>::Index& lambda,
			     const int j, const bool generators,
			     std::list<std::pair<typename DSBasis<d,dT>::Index, Support1D> >& intersecting)
  {
    typedef typename DSBasis<d,dT>::Index Index;

    intersecting.clear();

    // compute support of \psi_\lambda
    const int j_lambda = lambda.j() + lambda.e();
    int k1_lambda, k2_lambda;
    support(basis, lambda, k1_lambda, k2_lambda);
    
    // a brute force solution
    if (generators) {
      for (Index nu = first_generator(&basis, j);; ++nu) {
	Support1D supp;
	if (intersect_supports(basis, nu, j_lambda, k1_lambda, k2_lambda, supp.j, supp.k1, supp.k2))
	  intersecting.push_back(std::make_pair(nu, supp));
	if (nu == last_generator(&basis, j)) break;
      }
    } else {
      for (Index nu = first_wavelet(&basis, j);; ++nu) {
	Support1D supp;
	if (intersect_supports(basis, nu, j_lambda, k1_lambda, k2_lambda, supp.j, supp.k1, supp.k2))
	  intersecting.push_back(std::make_pair(nu, supp));
	if (nu == last_wavelet(&basis, j)) break;
      }
    }
  }
}
