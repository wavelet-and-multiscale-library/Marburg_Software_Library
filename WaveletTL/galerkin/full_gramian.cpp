// implementation for full_gramian.h

#include <list>
#include <map>

namespace WaveletTL
{
  template <int d, int dT, int s0, int s1, int sT0, int sT1, int J0>
  FullGramian<d,dT,s0,s1,sT0,sT1,J0>::FullGramian(const SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1,J0>& sb)
    : sb_(sb), j_(sb.j0())
  {
  }

  template <int d, int dT, int s0, int s1, int sT0, int sT1, int J0>
  inline
  const typename FullGramian<d,dT,s0,s1,sT0,sT1,J0>::size_type
  FullGramian<d,dT,s0,s1,sT0,sT1,J0>::row_dimension() const
  {
    return sb_.Deltasize(j_);
  }  
  
  template <int d, int dT, int s0, int s1, int sT0, int sT1, int J0>
  inline
  const typename FullGramian<d,dT,s0,s1,sT0,sT1,J0>::size_type
  FullGramian<d,dT,s0,s1,sT0,sT1,J0>::column_dimension() const
  {
    return row_dimension(); // square
  }  

  template <int d, int dT, int s0, int s1, int sT0, int sT1, int J0>
  void FullGramian<d,dT,s0,s1,sT0,sT1,J0>::set_level(const int j) const
  {
    assert(j >= sb_.j0());
    j_ = j;
  }

  template <int d, int dT, int s0, int s1, int sT0, int sT1, int J0>
  const double FullGramian<d,dT,s0,s1,sT0,sT1,J0>::get_entry(const size_type row, const size_type col) const
  {
    assert(row < row_dimension() && col < column_dimension());
    
#if 0
    Vector<double> ecol(column_dimension()), vcol(row_dimension(), false);
    ecol[col] = 1.0;
    apply(ecol, vcol);

    return vcol[row];
#else
    const int j0 = sb_.j0();

    // determine level of "row" and "col"
    int jrow = j0;
    if (row >= (size_type) sb_.Deltasize(j0)) {
      jrow += 1+(int)floor(log(((double)(row-sb_.Deltasize(j0)))/(1<<j0)+1)/M_LN2);
    }
    int jcol = j0;
    if (col >= (size_type) sb_.Deltasize(j0)) {
      jcol += 1+(int)floor(log(((double)(col-sb_.Deltasize(j0)))/(1<<j0)+1)/M_LN2);
    }

    // determine generator coeffs
    std::map<size_type,double> e_k_row; e_k_row[row] = 1.0;
    std::map<size_type,double> e_k_col; e_k_col[col] = 1.0;
    std::map<size_type,double> y_row, y_col;

    int j = std::max(jrow,jcol);
    if (j > j0) {
      sb_.apply_Tj(j-1, e_k_row, y_row);
      sb_.apply_Tj(j-1, e_k_col, y_col);
    } else {
      y_row.swap(e_k_row);
      y_col.swap(e_k_col);
    }

    // compute Ay_row (see "apply")
    std::map<size_type,double> Ay;
    if (d == 2) {
      if (s0==1 && s1==1) {
	// homogeneous b.c., apply tridiag(1/6,2/3,1/6)
	for (std::map<size_type,double>::const_iterator it(y_row.begin());
	     it != y_row.end(); ++it) {
	  Ay[it->first] += 2*it->second/3;
	  if (it->first > 0)
	    Ay[it->first-1] += it->second/6;
	  if (it->first < column_dimension()-1)
	    Ay[it->first+1] += it->second/6;
	}
      } else {
	if (s0==1) {
	  // homogeneous b.c. only at x=0, make a modification at right boundary
	  for (std::map<size_type,double>::const_iterator it(y_row.begin());
	       it != y_row.end(); ++it) {
	    Ay[it->first] +=
	      (it->first == column_dimension()-1 ? 1 : 2)*it->second/3;
	    if (it->first > 0)
	      Ay[it->first-1] += it->second/6;
	    if (it->first < column_dimension()-1)
	      Ay[it->first+1] += it->second/6;
	  }
	} else {
	  // homogeneous b.c. only at x=1, make a modification at left boundary
	  for (std::map<size_type,double>::const_iterator it(y_row.begin());
	       it != y_row.end(); ++it) {
	    Ay[it->first] +=
	      (it->first == 0 ? 1 : 2)*it->second/3;
	    if (it->first > 0)
	      Ay[it->first-1] += it->second/6;
	    if (it->first < column_dimension()-1)
	      Ay[it->first+1] += it->second/6;
	  }
	}
      }
    } else {
      if (d == 3) {
	// cf. [P, Bsp. 3.23]
	if (s0==1 && s1==1) {
	  // homogeneous b.c.
	  for (std::map<size_type,double>::const_iterator it(y_row.begin());
	       it != y_row.end(); ++it) {
	    const size_type m = sb_.Deltasize(j);
	    switch(it->first) {
	    case 0:
	      Ay[0] += it->second/3;
	      Ay[1] += 5*it->second/24;
	      Ay[2] += it->second/120;
	      break;
	    case 1:
	      Ay[0] += 5*it->second/24;
	      Ay[1] += 11*it->second/20;
	      Ay[2] += 13*it->second/60;
	      Ay[3] += it->second/120;
	      break;
	    default: // >= 2
	      switch(m-1-it->first) {
	      case 0: // m-1
		Ay[m-1] += it->second/3;
		Ay[m-2] += 5*it->second/24;
		Ay[m-3] += it->second/120;
		break;
	      case 1: // m-2
		Ay[m-1] += 5*it->second/24;
		Ay[m-2] += 11*it->second/20;
		Ay[m-3] += 13*it->second/60;
		Ay[m-4] += it->second/120;
		break;
	      default: // < m-2
		Ay[it->first-2] += it->second/120;
		Ay[it->first-1] += 13*it->second/60;
		Ay[it->first]   += 11*it->second/20;
		Ay[it->first+1] += 13*it->second/60;
		Ay[it->first+2] += it->second/120;
		break;
	      }
	      break;
	    }
	  }
	} else {
	  if (s0==1) {
	    // homogeneous b.c. at x=0
	    for (std::map<size_type,double>::const_iterator it(y_row.begin());
		 it != y_row.end(); ++it) {
	      const size_type m = sb_.Deltasize(j);
	      switch(it->first) {
	      case 0:
		Ay[0] += it->second/3;
		Ay[1] += 5*it->second/24;
		Ay[2] += it->second/120;
		break;
	      case 1:
		Ay[0] += 5*it->second/24;
		Ay[1] += 11*it->second/20;
		Ay[2] += 13*it->second/60;
		Ay[3] += it->second/120;
		break;
	      default: // >= 2
		switch(m-1-it->first) {
		case 0: // m-1
		  Ay[m-1] += it->second/5;
		  Ay[m-2] += 7*it->second/60;
		  Ay[m-3] += it->second/60;
		  break;
		case 1: // m-2
		  Ay[m-1] += 7*it->second/60;
		  Ay[m-2] += it->second/3;
		  Ay[m-3] += 5*it->second/24;
		  Ay[m-4] += it->second/120;
		  break;
		case 2: // m-3
		  Ay[m-1] += it->second/60;
		  Ay[m-2] += 5*it->second/24;
		  Ay[m-3] += 11*it->second/20;
		  Ay[m-4] += 13*it->second/60;
		  Ay[m-5] += it->second/120;
		  break;
		default: // < m-3
		  Ay[it->first-2] += it->second/120;
		  Ay[it->first-1] += 13*it->second/60;
		  Ay[it->first]   += 11*it->second/20;
		  Ay[it->first+1] += 13*it->second/60;
		  Ay[it->first+2] += it->second/120;
		  break;
		}
		break;
	      }
	    }
	  } else {
	    // homogeneous b.c. at x=1
	    for (std::map<size_type,double>::const_iterator it(y_row.begin());
		 it != y_row.end(); ++it) {
	      const size_type m = sb_.Deltasize(j);
	      switch(it->first) {
	      case 0:
		Ay[0] += it->second/5;
		Ay[1] += 7*it->second/60;
		Ay[2] += it->second/60;
		break;
	      case 1:
		Ay[0] += 7*it->second/60;
		Ay[1] += it->second/3;
		Ay[2] += 5*it->second/24;
		Ay[3] += it->second/120;
		break;
	      case 2:
		Ay[0] += it->second/60;
		Ay[1] += 5*it->second/24;
		Ay[2] += 11*it->second/20;
		Ay[3] += 13*it->second/60;
		Ay[4] += it->second/120;
		break;
	    default: // >= 2
	      switch(m-1-it->first) {
	      case 0: // m-1
		Ay[m-1] += it->second/3;
		Ay[m-2] += 5*it->second/24;
		Ay[m-3] += it->second/120;
		break;
	      case 1: // m-2
		Ay[m-1] += 5*it->second/24;
		Ay[m-2] += 11*it->second/20;
		Ay[m-3] += 13*it->second/60;
		Ay[m-4] += it->second/120;
		break;
	      default: // < m-2
		Ay[it->first-2] += it->second/120;
		Ay[it->first-1] += 13*it->second/60;
		Ay[it->first]   += 11*it->second/20;
		Ay[it->first+1] += 13*it->second/60;
		Ay[it->first+2] += it->second/120;
		break;
	      }
	      break;
	      }
	    }
	  }
	}
      }
    }
 
    // compute inner product <y_col,Ay_row>
    double r = 0;
    for (typename std::map<size_type,double>::const_iterator
	   ity(y_col.begin()),
	   ityend(y_col.end()),
	   itAy(Ay.begin()),
	   itAyend(Ay.end());
	 ity != ityend && itAy != itAyend; ++itAy)
      {
 	while (ity != ityend && ity->first < itAy->first) ++ity;
	if (ity != ityend)
	  if (itAy->first == ity->first)
	    r += ity->second * itAy->second;
      }
    return r;
 
//     std::map<size_type,double> ecol, vcol;
//     ecol[col] = 1.0;
//     apply(ecol, vcol);

//     return vcol[row];
#endif
  }

  template <int d, int dT, int s0, int s1, int sT0, int sT1, int J0>
  const double
  FullGramian<d,dT,s0,s1,sT0,sT1,J0>::diagonal(const size_type row) const
  {
    std::map<size_type,double> e_k; e_k[row] = 1.0;
    std::map<size_type,double> y;
    
    // determine corresponding level of "row"
    int jrow = sb_.j0();
    if (row >= (size_type) sb_.Deltasize(sb_.j0())) {
      jrow += 1+(int)floor(log(((double)(row-sb_.Deltasize(sb_.j0())))/(1<<sb_.j0())+1)/M_LN2);
    }
    
    // apply wavelet transformation y=T_{jrow-1}e_k
    // (does nothing if jrow==j0)
    if (jrow > sb_.j0())
      sb_.apply_Tj(jrow-1, e_k, y);
    else
      y.swap(e_k);
    
    // compute Ay (see "apply")
    std::map<size_type,double> Ay;
    if (d == 2) {
      if (s0==1 && s1==1) {
	// homogeneous b.c., apply tridiag(1/6,2/3,1/6)
	for (std::map<size_type,double>::const_iterator it(y.begin());
	     it != y.end(); ++it) {
	  Ay[it->first] += 2*it->second/3;
	  if (it->first > 0)
	    Ay[it->first-1] += it->second/6;
	  if (it->first < column_dimension()-1)
	    Ay[it->first+1] += it->second/6;
	}
      } else {
	if (s0==1) {
	  // homogeneous b.c. only at x=0, make a modification at right boundary
	  for (std::map<size_type,double>::const_iterator it(y.begin());
	       it != y.end(); ++it) {
	    Ay[it->first] +=
	      (it->first == column_dimension()-1 ? 1 : 2)*it->second/3;
	    if (it->first > 0)
	      Ay[it->first-1] += it->second/6;
	    if (it->first < column_dimension()-1)
	      Ay[it->first+1] += it->second/6;
	  }
	} else {
	  // homogeneous b.c. only at x=1, make a modification at left boundary
	  for (std::map<size_type,double>::const_iterator it(y.begin());
	       it != y.end(); ++it) {
	    Ay[it->first] +=
	      (it->first == 0 ? 1 : 2)*it->second/3;
	    if (it->first > 0)
	      Ay[it->first-1] += it->second/6;
	    if (it->first < column_dimension()-1)
	      Ay[it->first+1] += it->second/6;
	  }
	}
      }
    } else {
      if (d == 3) {
	// cf. [P, Bsp. 3.23]
	if (s0==1 && s1==1) {
	  // homogeneous b.c.
	  for (std::map<size_type,double>::const_iterator it(y.begin());
	       it != y.end(); ++it) {
	    const size_type m = sb_.Deltasize(jrow);
	    switch(it->first) {
	    case 0:
	      Ay[0] += it->second/3;
	      Ay[1] += 5*it->second/24;
	      Ay[2] += it->second/120;
	      break;
	    case 1:
	      Ay[0] += 5*it->second/24;
	      Ay[1] += 11*it->second/20;
	      Ay[2] += 13*it->second/60;
	      Ay[3] += it->second/120;
	      break;
	    default: // >= 2
	      switch(m-1-it->first) {
	      case 0: // m-1
		Ay[m-1] += it->second/3;
		Ay[m-2] += 5*it->second/24;
		Ay[m-3] += it->second/120;
		break;
	      case 1: // m-2
		Ay[m-1] += 5*it->second/24;
		Ay[m-2] += 11*it->second/20;
		Ay[m-3] += 13*it->second/60;
		Ay[m-4] += it->second/120;
		break;
	      default: // < m-2
		Ay[it->first-2] += it->second/120;
		Ay[it->first-1] += 13*it->second/60;
		Ay[it->first]   += 11*it->second/20;
		Ay[it->first+1] += 13*it->second/60;
		Ay[it->first+2] += it->second/120;
		break;
	      }
	      break;
	    }
	  }
	} else {
	  if (s0==1) {
	    // homogeneous b.c. at x=0
	    for (std::map<size_type,double>::const_iterator it(y.begin());
		 it != y.end(); ++it) {
	      const size_type m = sb_.Deltasize(jrow);
	      switch(it->first) {
	      case 0:
		Ay[0] += it->second/3;
		Ay[1] += 5*it->second/24;
		Ay[2] += it->second/120;
		break;
	      case 1:
		Ay[0] += 5*it->second/24;
		Ay[1] += 11*it->second/20;
		Ay[2] += 13*it->second/60;
		Ay[3] += it->second/120;
		break;
	      default: // >= 2
		switch(m-1-it->first) {
		case 0: // m-1
		  Ay[m-1] += it->second/5;
		  Ay[m-2] += 7*it->second/60;
		  Ay[m-3] += it->second/60;
		  break;
		case 1: // m-2
		  Ay[m-1] += 7*it->second/60;
		  Ay[m-2] += it->second/3;
		  Ay[m-3] += 5*it->second/24;
		  Ay[m-4] += it->second/120;
		  break;
		case 2: // m-3
		  Ay[m-1] += it->second/60;
		  Ay[m-2] += 5*it->second/24;
		  Ay[m-3] += 11*it->second/20;
		  Ay[m-4] += 13*it->second/60;
		  Ay[m-5] += it->second/120;
		  break;
		default: // < m-3
		  Ay[it->first-2] += it->second/120;
		  Ay[it->first-1] += 13*it->second/60;
		  Ay[it->first]   += 11*it->second/20;
		  Ay[it->first+1] += 13*it->second/60;
		  Ay[it->first+2] += it->second/120;
		  break;
		}
		break;
	      }
	    }
	  } else {
	    // homogeneous b.c. at x=1
	    for (std::map<size_type,double>::const_iterator it(y.begin());
		 it != y.end(); ++it) {
	      const size_type m = sb_.Deltasize(jrow);
	      switch(it->first) {
	      case 0:
		Ay[0] += it->second/5;
		Ay[1] += 7*it->second/60;
		Ay[2] += it->second/60;
		break;
	      case 1:
		Ay[0] += 7*it->second/60;
		Ay[1] += it->second/3;
		Ay[2] += 5*it->second/24;
		Ay[3] += it->second/120;
		break;
	      case 2:
		Ay[0] += it->second/60;
		Ay[1] += 5*it->second/24;
		Ay[2] += 11*it->second/20;
		Ay[3] += 13*it->second/60;
		Ay[4] += it->second/120;
		break;
	    default: // >= 2
	      switch(m-1-it->first) {
	      case 0: // m-1
		Ay[m-1] += it->second/3;
		Ay[m-2] += 5*it->second/24;
		Ay[m-3] += it->second/120;
		break;
	      case 1: // m-2
		Ay[m-1] += 5*it->second/24;
		Ay[m-2] += 11*it->second/20;
		Ay[m-3] += 13*it->second/60;
		Ay[m-4] += it->second/120;
		break;
	      default: // < m-2
		Ay[it->first-2] += it->second/120;
		Ay[it->first-1] += 13*it->second/60;
		Ay[it->first]   += 11*it->second/20;
		Ay[it->first+1] += 13*it->second/60;
		Ay[it->first+2] += it->second/120;
		break;
	      }
	      break;
	      }
	    }
	  }
	}
      }
    }
    
    // compute inner product <y,Ay>
    double r = 0;
    for (typename std::map<size_type,double>::const_iterator
	   ity(y.begin()),
	   ityend(y.end()),
	   itAy(Ay.begin()),
	   itAyend(Ay.end());
	 ity != ityend && itAy != itAyend; ++itAy)
      {
 	while (ity != ityend && ity->first < itAy->first) ++ity;
	if (ity != ityend)
	  if (itAy->first == ity->first)
	    r += ity->second * itAy->second;
      }
    return r;
  }
  
  template <int d, int dT, int s0, int s1, int sT0, int sT1, int J0>
  template <class VECTOR>
  void FullGramian<d,dT,s0,s1,sT0,sT1,J0>::apply(const VECTOR& x, VECTOR& Mx) const
  {
    assert(Mx.size() == row_dimension());

    VECTOR y(x.size(), false); // no initialization necessary
    
    // apply wavelet transformation T_{j-1}
    // (does nothing if j==j0)
    if (j_ > sb_.j0())
      sb_.apply_Tj(j_-1, x, Mx);
    else
      Mx = x;

    // apply Gramian w.r.t the B-Splines in V_j
    if (d == 2) {
      if (s0==1 && s1==1) {
	// homogeneous b.c., apply tridiag(1/6,2/3,1/6)
	y[0] = (4*Mx[0] + Mx[1])/6;
	const size_type m = row_dimension();
	y[m-1] = (4*Mx[m-1] + Mx[m-2])/6;
	for (size_type i(1); i < m-1; i++)
	  y[i] = (Mx[i-1] + 4*Mx[i] + Mx[i+1])/6;
      } else {
	if (s0==1) {
	  // homogeneous b.c. at x=0
	  y[0] = (4*Mx[0] + Mx[1])/6;
	  const size_type m = row_dimension();
	  y[m-1] = (2*Mx[m-1] + Mx[m-2])/6;
	  for (size_type i(1); i < m-1; i++)
	    y[i] = (Mx[i-1] + 4*Mx[i] + Mx[i+1])/6;
	} else {
	  // homogeneous b.c. at x=1
	  y[0] = (2*Mx[0] + Mx[1])/6;
	  const size_type m = row_dimension();
	  y[m-1] = (4*Mx[m-1] + Mx[m-2])/6;
	  for (size_type i(1); i < m-1; i++)
	    y[i] = (Mx[i-1] + 4*Mx[i] + Mx[i+1])/6;
	}
      }
    } else {
      if (d == 3) {
	// cf. [P, Bsp. 3.23]
	if (s0==1 && s1==1) {
	  // homogeneous b.c.
	  y[0] = Mx[0]/3 + 5*Mx[1]/24 + Mx[2]/120;
	  y[1] = 5*Mx[0]/24 + 11*Mx[1]/20 + 13*Mx[2]/60 + Mx[3]/120;
	  const size_type m = row_dimension();
	  y[m-1] = Mx[m-3]/120 + 5*Mx[m-2]/24 + Mx[m-1]/3;
	  y[m-2] = Mx[m-4]/120 + 13*Mx[m-3]/60 + 11*Mx[m-2]/20 + 5*Mx[m-1]/24;
	  for (size_type i(2); i < m-2; i++)
	    y[i] = Mx[i-2]/120 + 13*Mx[i-1]/60 + 11*Mx[i]/20 + 13*Mx[i+1]/60 + Mx[i+2]/120;
	} else {
	  if (s0==1) {
	    // homogeneous b.c. at x=0	    
	    y[0] = Mx[0]/3 + 5*Mx[1]/24 + Mx[2]/120;
	    y[1] = 5*Mx[0]/24 + 11*Mx[1]/20 + 13*Mx[2]/60 + Mx[3]/120;
	    const size_type m = row_dimension();
	    y[m-1] = Mx[m-3]/60 + 7*Mx[m-2]/60 + Mx[m-1]/5;
	    y[m-2] = Mx[m-4]/120 + 5*Mx[m-3]/24 + Mx[m-2]/3 + 7*Mx[m-1]/60;
	    y[m-3] = Mx[m-5]/120 + 13*Mx[m-4]/60 + 11*Mx[m-3]/20 + 5*Mx[m-2]/24 + Mx[m-1]/60;
	    for (size_type i(2); i < m-3; i++)
	      y[i] = Mx[i-2]/120 + 13*Mx[i-1]/60 + 11*Mx[i]/20 + 13*Mx[i+1]/60 + Mx[i+2]/120;
	  } else {
	    // homogeneous b.c. at x=1
	    y[0] = Mx[0]/5 + 7*Mx[1]/60 + Mx[2]/60;
	    y[1] = 7*Mx[0]/60 + Mx[1]/3 + 5*Mx[2]/24 + Mx[3]/120;
	    y[2] = Mx[0]/60 + 5*Mx[1]/24 + 11*Mx[2]/20 + 13*Mx[3]/60 + Mx[4]/120;
	    const size_type m = row_dimension();
	    y[m-1] = Mx[m-3]/120 + 5*Mx[m-2]/24 + Mx[m-1]/3;
	    y[m-2] = Mx[m-4]/120 + 13*Mx[m-3]/60 + 11*Mx[m-2]/20 + 5*Mx[m-1]/24;
	    for (size_type i(3); i < m-2; i++)
	      y[i] = Mx[i-2]/120 + 13*Mx[i-1]/60 + 11*Mx[i]/20 + 13*Mx[i+1]/60 + Mx[i+2]/120;
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
  }

  template <int d, int dT, int s0, int s1, int sT0, int sT1, int J0>
  void FullGramian<d,dT,s0,s1,sT0,sT1,J0>::apply(const std::map<size_type,double>& x,
				std::map<size_type,double>& Mx) const
  {
    std::map<size_type,double> y; // no initialization necessary
    
    // apply wavelet transformation T_{j-1}
    // (does nothing if j==j0)
    if (j_ > sb_.j0())
      sb_.apply_Tj(j_-1, x, Mx);
    else
      Mx = x;

    // apply Gramian w.r.t the B-Splines in V_j
    if (d == 2) {
      if (s0==1 && s1==1) {
	// homogeneous b.c., apply tridiag(1/6,2/3,1/6)
	for (std::map<size_type,double>::const_iterator it(Mx.begin());
	     it != Mx.end(); ++it) {
	  y[it->first] += 2*it->second/3;
	  if (it->first > 0)
	    y[it->first-1] += it->second/6;
	  if (it->first < column_dimension()-1)
	    y[it->first+1] += it->second/6;
	}
      } else {
	if (s0==1) {
	  // homogeneous b.c. only at x=0, make a modification at right boundary
	  for (std::map<size_type,double>::const_iterator it(Mx.begin());
	       it != Mx.end(); ++it) {
	    y[it->first] +=
	      (it->first == column_dimension()-1 ? 1 : 2)*it->second/3;
	    if (it->first > 0)
	      y[it->first-1] += it->second/6;
	    if (it->first < column_dimension()-1)
	      y[it->first+1] += it->second/6;
	  }
	} else {
	  // homogeneous b.c. only at x=1, make a modification at left boundary
	  for (std::map<size_type,double>::const_iterator it(Mx.begin());
	       it != Mx.end(); ++it) {
	    y[it->first] +=
	      (it->first == 0 ? 1 : 2)*it->second/3;
	    if (it->first > 0)
	      y[it->first-1] += it->second/6;
	    if (it->first < column_dimension()-1)
	      y[it->first+1] += it->second/6;
	  }
	}
      }
    } else {
      if (d == 3) {
	// cf. [P, Bsp. 3.23]
	if (s0==1 && s1==1) {
	  // homogeneous b.c.
	  for (std::map<size_type,double>::const_iterator it(Mx.begin());
	       it != Mx.end(); ++it) {
	    const size_type m = row_dimension();
	    switch(it->first) {
	    case 0:
	      y[0] += it->second/3;
	      y[1] += 5*it->second/24;
	      y[2] += it->second/120;
	      break;
	    case 1:
	      y[0] += 5*it->second/24;
	      y[1] += 11*it->second/20;
	      y[2] += 13*it->second/60;
	      y[3] += it->second/120;
	      break;
	    default: // >= 2
	      switch(m-1-it->first) {
	      case 0: // m-1
		y[m-1] += it->second/3;
		y[m-2] += 5*it->second/24;
		y[m-3] += it->second/120;
		break;
	      case 1: // m-2
		y[m-1] += 5*it->second/24;
		y[m-2] += 11*it->second/20;
		y[m-3] += 13*it->second/60;
		y[m-4] += it->second/120;
		break;
	      default: // < m-2
		y[it->first-2] += it->second/120;
		y[it->first-1] += 13*it->second/60;
		y[it->first]   += 11*it->second/20;
		y[it->first+1] += 13*it->second/60;
		y[it->first+2] += it->second/120;
		break;
	      }
	      break;
	    }
	  }
	} else {
	  if (s0==1) {
	    // homogeneous b.c. at x=0
	    for (std::map<size_type,double>::const_iterator it(Mx.begin());
		 it != Mx.end(); ++it) {
	      const size_type m = row_dimension();
	      switch(it->first) {
	      case 0:
		y[0] += it->second/3;
		y[1] += 5*it->second/24;
		y[2] += it->second/120;
		break;
	      case 1:
		y[0] += 5*it->second/24;
		y[1] += 11*it->second/20;
		y[2] += 13*it->second/60;
		y[3] += it->second/120;
		break;
	      default: // >= 2
		switch(m-1-it->first) {
		case 0: // m-1
		  y[m-1] += it->second/5;
		  y[m-2] += 7*it->second/60;
		  y[m-3] += it->second/60;
		  break;
		case 1: // m-2
		  y[m-1] += 7*it->second/60;
		  y[m-2] += it->second/3;
		  y[m-3] += 5*it->second/24;
		  y[m-4] += it->second/120;
		  break;
		case 2: // m-3
		  y[m-1] += it->second/60;
		  y[m-2] += 5*it->second/24;
		  y[m-3] += 11*it->second/20;
		  y[m-4] += 13*it->second/60;
		  y[m-5] += it->second/120;
		  break;
		default: // < m-3
		  y[it->first-2] += it->second/120;
		  y[it->first-1] += 13*it->second/60;
		  y[it->first]   += 11*it->second/20;
		  y[it->first+1] += 13*it->second/60;
		  y[it->first+2] += it->second/120;
		  break;
		}
		break;
	      }
	    }
	  } else {
	    // homogeneous b.c. at x=1
	    for (std::map<size_type,double>::const_iterator it(Mx.begin());
		 it != Mx.end(); ++it) {
	      const size_type m = row_dimension();
	      switch(it->first) {
	      case 0:
		y[0] += it->second/5;
		y[1] += 7*it->second/60;
		y[2] += it->second/60;
		break;
	      case 1:
		y[0] += 7*it->second/60;
		y[1] += it->second/3;
		y[2] += 5*it->second/24;
		y[3] += it->second/120;
		break;
	      case 2:
		y[0] += it->second/60;
		y[1] += 5*it->second/24;
		y[2] += 11*it->second/20;
		y[3] += 13*it->second/60;
		y[4] += it->second/120;
		break;
	    default: // >= 2
	      switch(m-1-it->first) {
	      case 0: // m-1
		y[m-1] += it->second/3;
		y[m-2] += 5*it->second/24;
		y[m-3] += it->second/120;
		break;
	      case 1: // m-2
		y[m-1] += 5*it->second/24;
		y[m-2] += 11*it->second/20;
		y[m-3] += 13*it->second/60;
		y[m-4] += it->second/120;
		break;
	      default: // < m-2
		y[it->first-2] += it->second/120;
		y[it->first-1] += 13*it->second/60;
		y[it->first]   += 11*it->second/20;
		y[it->first+1] += 13*it->second/60;
		y[it->first+2] += it->second/120;
		break;
	      }
	      break;
	      }
	    }
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

    // remove unnecessary zeros
    for (typename std::map<size_type,double>::iterator it(Mx.begin()); it != Mx.end();) {
      if (it->second == 0)
	Mx.erase(it++);
      else
	++it;
    }  
  }
  
  template <int d, int dT, int s0, int s1, int sT0, int sT1, int J0>
  void
  FullGramian<d,dT,s0,s1,sT0,sT1,J0>::to_sparse(SparseMatrix<double>& S) const
  {
    // we utilize that the Gramian is symmetric
    const size_type m = row_dimension();
    S.resize(m, m);
    for (size_type row(0); row < m; row++) {
      std::map<size_type,double> x, Mx;
      x[row] = 1.0;
      apply(x, Mx);
      typename std::list<size_type> indices;
      std::list<double> entries;
      for (typename std::map<size_type,double>::const_iterator it(Mx.begin()); it != Mx.end(); ++it) {
	indices.push_back(it->first);
	entries.push_back(it->second);
      }
      S.set_row(row, indices, entries);
    }
  }
  
  template <int d, int dT, int s0, int s1, int sT0, int sT1, int J0>
  void FullGramian<d,dT,s0,s1,sT0,sT1,J0>::print(std::ostream &os,
				  const unsigned int tabwidth,
				  const unsigned int precision) const
  {
    if (row_dimension() == 0)
      os << "[]" << std::endl; // Matlab style
    else
      {
	unsigned int old_precision = os.precision(precision);
	for (size_type i(0); i < row_dimension(); ++i)
	  {
	    for (size_type j(0); j < column_dimension(); ++j)
	      os << std::setw(tabwidth) << std::setprecision(precision)
		 << get_entry(i, j);
	    os << std::endl;
	  }
	os.precision(old_precision);
      }
  }

  template <int d, int dT, int s0, int s1, int sT0, int sT1, int J0>
  inline
  std::ostream& operator << (std::ostream& os, const FullGramian<d,dT,s0,s1,sT0,sT1,J0>& M)
  {
    M.print(os);
    return os;
  }

}
