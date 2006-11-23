// implementation for full_laplacian.h

namespace WaveletTL
{
  template <int d, int dT>
  FullLaplacian<d,dT>::FullLaplacian(const SplineBasis<d,dT>& sb)
    : sb_(sb), j_(sb.j0())
  {
  }

  template <int d, int dT>
  inline
  const typename FullLaplacian<d,dT>::size_type
  FullLaplacian<d,dT>::row_dimension() const
  {
    return sb_.Deltasize(j_);
  }  
  
  template <int d, int dT>
  inline
  const typename FullLaplacian<d,dT>::size_type
  FullLaplacian<d,dT>::column_dimension() const
  {
    return row_dimension(); // square
  }  

  template <int d, int dT>
  void FullLaplacian<d,dT>::set_level(const int j) const
  {
    assert(j >= sb_.j0());
    j_ = j;
  }

  template <int d, int dT>
  const double FullLaplacian<d,dT>::get_entry(const size_type row, const size_type column) const
  {
    assert(row < row_dimension() && column < column_dimension());
    
    Vector<double> ecol(column_dimension()), col(row_dimension());
    ecol[column] = 1.0;
    apply(ecol, col);

    return col[row];
  }

  template <int d, int dT>
  template <class VECTOR>
  void FullLaplacian<d,dT>::apply(const VECTOR& x, VECTOR& Mx) const
  {
    assert(Mx.size() == row_dimension());

    VECTOR y(x);

    // apply diagonal preconditioner D^{-1}
    // (does nothing if j==j0)
    for (int j = sb_.j0(); j < j_; j++) {
      for (int k(sb_.Deltasize(j)); k < sb_.Deltasize(j+1); k++)
	y[k] /= (1<<j);
    }

    // apply wavelet transformation T_{j-1}
    // (does nothing if j==j0)
    if (j_ > sb_.j0())
      sb_.apply_Tj(j_-1, y, Mx);
    else
      Mx.swap(y);

    // apply Laplacian w.r.t the B-Splines in V_j
    if (d == 2) {
      // apply 2^{2j}*tridiag(-1,2,-1)
      y[0] = ldexp(1.0, 2*j_) * (2*Mx[0] - Mx[1]);
      y[row_dimension()-1] = ldexp(1.0, 2*j_) * (2*Mx[row_dimension()-1]-Mx[row_dimension()-2]);
      for (size_type i(1); i < row_dimension()-1; i++)
	y[i] = ldexp(1.0, 2*j_) * (2*Mx[i] - Mx[i-1] - Mx[i+1]);
    } else {
      if (d == 3) {
	// cf. [P, Bsp. 3.26]
	y[0] = ldexp(1.0, 2*j_) * (4*Mx[0]/3 - Mx[1]/6 - Mx[2]/6);
	y[1] = ldexp(1.0, 2*j_) * (-Mx[0]/6 + Mx[1] - Mx[2]/3 - Mx[3]/6);
	const size_type m = row_dimension();
	y[m-1] = ldexp(1.0, 2*j_) * (-Mx[m-3]/6 - Mx[m-2]/6 + 4*Mx[m-1]/3);
	y[m-2] = ldexp(1.0, 2*j_) * (-Mx[m-4]/6 - Mx[m-3]/3 + Mx[m-2] - Mx[m-1]/6);
	for (size_type i(2); i < m-2; i++)
	  y[i] = ldexp(1.0, 2*j_) * (-Mx[i-2]/6 - Mx[i-1]/3 + Mx[i] - Mx[i+1]/3 - Mx[i+2]/6);
      }
    }
    
    // apply transposed wavelet transformation T_{j-1}^T
    // (does nothing if j==j0)
    if (j_ > sb_.j0())
      sb_.apply_Tj_transposed(j_-1, y, Mx);
    else
      Mx.swap(y);
    
    // apply diagonal preconditioner D^{-1}
    // (does nothing if j==j0)
    for (int j = sb_.j0(); j < j_; j++) {
      for (int k(sb_.Deltasize(j)); k < sb_.Deltasize(j+1); k++)
	Mx[k] /= (1<<j);
    }
  }

  template <int d, int dT>
  void FullLaplacian<d,dT>::print(std::ostream &os,
				  const unsigned int tabwidth,
				  const unsigned int precision) const
  {
    if (row_dimension() == 0)
      os << "[]" << std::endl; // Matlab style
    else
      {
	unsigned int old_precision = os.precision(precision);
	for (typename FullLaplacian<d,dT>::size_type i(0); i < row_dimension(); ++i)
	  {
	    for (typename FullLaplacian<d,dT>::size_type j(0); j < column_dimension(); ++j)
	      os << std::setw(tabwidth) << std::setprecision(precision)
		 << get_entry(i, j);
	    os << std::endl;
	  }
	os.precision(old_precision);
      }
  }

  template <int d, int dT>
  inline
  std::ostream& operator << (std::ostream& os, const FullLaplacian<d,dT>& M)
  {
    M.print(os);
    return os;
  }

}
