// implementation for full_gramian.h

namespace WaveletTL
{
  template <int d, int dT>
  FullGramian<d,dT>::FullGramian(const SplineBasis<d,dT>& sb)
    : sb_(sb), j_(sb.j0())
  {
  }

  template <int d, int dT>
  inline
  const typename FullGramian<d,dT>::size_type
  FullGramian<d,dT>::row_dimension() const
  {
    return sb_.Deltasize(j_);
  }  
  
  template <int d, int dT>
  inline
  const typename FullGramian<d,dT>::size_type
  FullGramian<d,dT>::column_dimension() const
  {
    return row_dimension(); // square
  }  

  template <int d, int dT>
  void FullGramian<d,dT>::set_level(const int j) const
  {
    assert(j >= sb_.j0());
    j_ = j;
  }

  template <int d, int dT>
  const double FullGramian<d,dT>::get_entry(const size_type row, const size_type column) const
  {
    assert(row < row_dimension() && column < column_dimension());
    
    Vector<double> ecol(column_dimension()), col(row_dimension());
    ecol[column] = 1.0;
    apply(ecol, col);

    return col[row];
  }

  template <int d, int dT>
  template <class VECTOR>
  void FullGramian<d,dT>::apply(const VECTOR& x, VECTOR& Mx) const
  {
    assert(Mx.size() == row_dimension());

    VECTOR y(x);

    // apply wavelet transformation T_{j-1}
    // (does nothing if j==j0)
    if (j_ > sb_.j0())
      sb_.apply_Tj(j_-1, y, Mx);
    else
      Mx.swap(y);

    // apply Gramian w.r.t the B-Splines in V_j
    if (d == 2) {
      // apply tridiag(1/6,2/3,1/6)
      y[0] = (4*Mx[0] + Mx[1])/6;
      y[row_dimension()-1] = (4*Mx[row_dimension()-1] + Mx[row_dimension()-2])/6;
      for (size_type i(1); i < row_dimension()-1; i++)
	y[i] = (Mx[i-1] + 4*Mx[i] + Mx[i+1])/6;
    } else {
      if (d == 3) {
	// not yet implemented
      }
    }
    
    // apply transposed wavelet transformation T_{j-1}^T
    // (does nothing if j==j0)
    if (j_ > sb_.j0())
      sb_.apply_Tj_transposed(j_-1, y, Mx);
    else
      Mx.swap(y);
  }

  template <int d, int dT>
  void FullGramian<d,dT>::print(std::ostream &os,
				  const unsigned int tabwidth,
				  const unsigned int precision) const
  {
    if (row_dimension() == 0)
      os << "[]" << std::endl; // Matlab style
    else
      {
	unsigned int old_precision = os.precision(precision);
	for (typename FullGramian<d,dT>::size_type i(0); i < row_dimension(); ++i)
	  {
	    for (typename FullGramian<d,dT>::size_type j(0); j < column_dimension(); ++j)
	      os << std::setw(tabwidth) << std::setprecision(precision)
		 << get_entry(i, j);
	    os << std::endl;
	  }
	os.precision(old_precision);
      }
  }

  template <int d, int dT>
  inline
  std::ostream& operator << (std::ostream& os, const FullGramian<d,dT>& M)
  {
    M.print(os);
    return os;
  }

}
