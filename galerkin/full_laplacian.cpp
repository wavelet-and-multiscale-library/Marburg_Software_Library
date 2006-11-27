// implementation for full_laplacian.h

#include <map>

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

#if 0
    Vector<double> ecol(column_dimension()), col(row_dimension());
    ecol[column] = 1.0;
    apply(ecol, col);

    return col[row];
#else
    std::map<size_type,double> ecol, col;
    ecol[column] = 1.0;
    apply(ecol, col);

    return col[row];
#endif
  }

  template <int d, int dT>
  template <class VECTOR>
  void FullLaplacian<d,dT>::apply(const VECTOR& x, VECTOR& Mx) const
  {
    assert(Mx.size() == row_dimension());

    VECTOR y(x);

    // apply diagonal preconditioner D^{-1}
    for (int k(0); k < sb_.Deltasize(sb_.j0()); k++)
      y[k] /= (1<<sb_.j0());
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
    const double factor = ldexp(1.0, 2*j_); // not "1<<(2*j_)" !
    if (d == 2) {
      // apply 2^{2j}*tridiag(-1,2,-1)
      y[0] = factor * (2*Mx[0] - Mx[1]);
      const size_type m = row_dimension();
      y[m-1] = factor * (2*Mx[m-1]-Mx[m-2]);
      for (size_type i(1); i < m-1; i++)
	y[i] = factor * (2*Mx[i] - Mx[i-1] - Mx[i+1]);
    } else {
      if (d == 3) {
	// cf. [P, Bsp. 3.26]
	y[0] = factor * (4*Mx[0]/3 - Mx[1]/6 - Mx[2]/6);
	y[1] = factor * (-Mx[0]/6 + Mx[1] - Mx[2]/3 - Mx[3]/6);
	const size_type m = row_dimension();
	y[m-1] = factor * (-Mx[m-3]/6 - Mx[m-2]/6 + 4*Mx[m-1]/3);
	y[m-2] = factor * (-Mx[m-4]/6 - Mx[m-3]/3 + Mx[m-2] - Mx[m-1]/6);
	for (size_type i(2); i < m-2; i++)
	  y[i] = factor * (-Mx[i-2]/6 - Mx[i-1]/3 + Mx[i] - Mx[i+1]/3 - Mx[i+2]/6);
      }
    }
    
    // apply transposed wavelet transformation T_{j-1}^T
    // (does nothing if j==j0)
    if (j_ > sb_.j0())
      sb_.apply_Tj_transposed(j_-1, y, Mx);
    else
      Mx.swap(y);
    
    // apply diagonal preconditioner D^{-1}
    for (int k(0); k < sb_.Deltasize(sb_.j0()); k++)
      Mx[k] /= (1<<sb_.j0());
    for (int j = sb_.j0(); j < j_; j++) {
      for (int k(sb_.Deltasize(j)); k < sb_.Deltasize(j+1); k++)
	Mx[k] /= (1<<j);
    }
  }

  template <int d, int dT>
  void FullLaplacian<d,dT>::apply(const std::map<size_type,double>& x,
				  std::map<size_type,double>& Mx) const
  {
    std::map<size_type,double> y(x);

    // apply diagonal preconditioner D^{-1}
    for (int k(0); k < sb_.Deltasize(sb_.j0()); k++)
      y[k] /= (1<<sb_.j0());
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
    const double factor = ldexp(1.0, 2*j_); // not "1<<(2*j_)" !
    y.clear();
    if (d == 2) {
      // apply 2^{2j}*tridiag(-1,2,-1)
      for (std::map<size_type,double>::const_iterator it(Mx.begin());
	   it != Mx.end(); ++it) {
	y[it->first] += factor * 2*it->second;
	if (it->first > 0)
	  y[it->first-1] -= factor * it->second;
	if (it->first < column_dimension()-1)
	  y[it->first+1] -= factor * it->second;
      }
    } else {
      if (d == 3) {
	// cf. [P, Bsp. 3.26]
	for (std::map<size_type,double>::const_iterator it(Mx.begin());
	     it != Mx.end(); ++it) {
	  const size_type m = row_dimension();
	  switch(it->first) {
	  case 0:
	    y[0] += factor * 4*it->second/3;
	    y[1] -= factor * it->second/6;
	    y[2] -= factor * it->second/6;
	    break;
	  case 1:
	    y[0] -= factor * it->second/6;
	    y[1] += factor * it->second;
	    y[2] -= factor * it->second/3;
	    y[3] -= factor * it->second/6;
	    break;
	  default: // >= 2
	    switch(m-1-it->first) {
	    case 0: // m-1
	      y[m-1] += factor * 4*it->second/3;
	      y[m-2] -= factor * it->second/6;
	      y[m-3] -= factor * it->second/6;
	      break;
	    case 1: // m-2
	      y[m-1] -= factor * it->second/6;
	      y[m-2] += factor * it->second;
	      y[m-3] -= factor * it->second/3;
	      y[m-4] -= factor * it->second/6;
	      break;
	    default: // < m-2
	      y[it->first-2] -= factor * it->second/6;
	      y[it->first-1] -= factor * it->second/3;
	      y[it->first]   += factor * it->second;
	      y[it->first+1] -= factor * it->second/3;
	      y[it->first+2] -= factor * it->second/6;
	      break;
	    }
	    break;
	  }
	}
      }
    }
    
    // apply transposed wavelet transformation T_{j-1}^T
    // (does nothing if j==j0)
    if (j_ > sb_.j0())
      sb_.apply_Tj_transposed(j_-1, y, Mx);
    else
      Mx.swap(y);
    
    // apply diagonal preconditioner D^{-1}
    for (int k(0); k < sb_.Deltasize(sb_.j0()); k++)
      Mx[k] /= (1<<sb_.j0());
    for (int j = sb_.j0(); j < j_; j++) {
      for (int k(sb_.Deltasize(j)); k < sb_.Deltasize(j+1); k++)
	Mx[k] /= (1<<j);
    }

    // remove unnecessary zeros
    for (typename std::map<size_type,double>::iterator it(Mx.begin()); it != Mx.end();) {
      if (it->second == 0)
	Mx.erase(it++);
      else
	++it;
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
