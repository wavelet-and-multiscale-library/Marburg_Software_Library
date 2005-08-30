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
			  const typename DSBasis<d,dT>::Index& nu,
			  int& j, int& k1, int& k2)
  {
    const int j_lambda = lambda.j() + lambda.e();
    const int j_nu     = nu.j() + nu.e();
    j = std::max(j_lambda, j_nu);
    int k1_lambda, k2_lambda, k1_nu, k2_nu;
    support(basis, lambda, k1_lambda, k2_lambda);
    support(basis, nu    , k1_nu    , k2_nu    );
    if ((1<<(j-j_lambda))*k1_lambda < (1<<(j-j_nu))*k2_nu) {
      if ((1<<(j-j_nu))*k1_nu < (1<<(j-j_lambda))*k2_lambda) {
	k1 = std::max((1<<(j-j_lambda))*k1_lambda, (1<<(j-j_nu))*k1_nu);
	k2 = std::min((1<<(j-j_lambda))*k2_lambda, (1<<(j-j_nu))*k2_nu);
	return true;
      }
    }
    return false;
  }
}
