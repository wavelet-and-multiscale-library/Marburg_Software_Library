// implementation for p_support.h

namespace WaveletTL
{
  template <int d, int dT>
  void support(const PBasis<d,dT>& basis,
	       const typename PBasis<d,dT>::Index& lambda,
	       int& k1, int& k2)
  {
    if (lambda.e() == 0) // generator
      {
	// phi_{j,k}(x) = 2^{j/2} B_{j,k-ell_1}
	k1 = std::max(0, lambda.k() + ell1<d>());
	k2 = std::min(1<<lambda.j(), lambda.k() + ell2<d>());
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
	support(basis, typename PBasis<d,dT>::Index(lambda.j()+1, 0, kleft, &basis), k1, dummy);
	support(basis, typename PBasis<d,dT>::Index(lambda.j()+1, 0, kright, &basis), dummy, k2);
      }
  }
}
