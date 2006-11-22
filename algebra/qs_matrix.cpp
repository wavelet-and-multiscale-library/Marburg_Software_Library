// implementation for qs_matrix.h

namespace MathTL
{
  template <class C>
  QuasiStationaryMatrix<C>
  ::QuasiStationaryMatrix(const int j0,
			  const unsigned int mj0, const unsigned int nj0,
			  const Matrix<C>& ML, const Matrix<C>& MR,
			  const Vector<C>& bandL, const Vector<C>& bandR,
			  const int offsetL, const int offsetR)
    : j0_(j0), mj0_(mj0), nj0_(nj0),
      ML_(ML), MR_(MR),
      bandL_(bandL), bandR_(bandR),
      offsetL_(offsetL), offsetR_(offsetR)
  {
    j_ = j0_;
  }

  template <class C>
  inline
  const typename QuasiStationaryMatrix<C>::size_type
  QuasiStationaryMatrix<C>::row_dimension() const
  {
    return mj0_-(1<<j0_)+(1<<j_);
  }  

  template <class C>
  inline
  const typename QuasiStationaryMatrix<C>::size_type
  QuasiStationaryMatrix<C>::column_dimension() const
  {
    return nj0_-(1<<j0_)+(1<<j_);
  }  

  template <class C>
  const C QuasiStationaryMatrix<C>::get_entry(const size_type row, const size_type column) const
  {
    assert(row < row_dimension() && column < column_dimension());
    
    // for readability:
    const unsigned int ml = ML_.row_dimension();
    const unsigned int nl = ML_.column_dimension();
    const unsigned int mr = MR_.row_dimension();
    const unsigned int nr = MR_.column_dimension();
    const unsigned int m = row_dimension();
    const unsigned int n = column_dimension();

    if (column < nl) { // left block column region
      if (row < ml)
	return ML_.get_entry(row, column);
    } else {
      if (column < n/2) { // left band column region
	const unsigned int colbegin = offsetL_+2*(column-nl);
	if (row >= colbegin && row < colbegin+bandL_.size())
	  return bandL_[row-colbegin];
      } else {
	if (column >= n-nr) { // right block column region
	  if (row >= m-mr)
	    return MR_.get_entry(row-(m-mr),
				 column-(n-nr));
	} else { // right band column region
	  const unsigned int colendplus1 = m-offsetR_-2*(n-nr-1-column);
	  if (row < colendplus1 && row >= colendplus1-bandR_.size())
	    return bandR_[(row+bandR_.size())-colendplus1];
	}
      }
    }

    return 0;
  }

  template <class C>
  void QuasiStationaryMatrix<C>::set_level(const int j)
  {
    assert(j >= j0_);
    j_ = j;
  }

  template <class C>
  template <class VECTOR>
  void QuasiStationaryMatrix<C>::apply(const VECTOR& x, VECTOR& Mx) const
  {
    assert(Mx.size() == row_dimension());
    
    // for readability:
    const unsigned int ml = ML_.row_dimension();
    const unsigned int nl = ML_.column_dimension();
    const unsigned int mr = MR_.row_dimension();
    const unsigned int nr = MR_.column_dimension();
    const unsigned int m = row_dimension();
    const unsigned int n = column_dimension();

    for (size_type i(0); i < row_dimension(); i++) {
      C help(0);

      // contribution from upper left corner block
      if (i < ml)
 	for (size_type j(0); j < nl; j++)
 	  help += ML_.get_entry(i, j) * x[j];
      
      // contribution from left band
      const unsigned int jbeginlefthelp = (i+2*nl+1)-bandL_.size()-offsetL_;
      for (size_type j = std::max(nl, jbeginlefthelp-jbeginlefthelp/2);
	   j <= std::min(n/2-1, ((i+2*nl)-offsetL_)/2); j++)
	help += bandL_[i-offsetL_-2*(j-nl)] * x[j];
      
      // contribution from right band
      const unsigned int jbeginrighthelp = (i+offsetR_+2*(n-nr-1)+1)-m;
      for (size_type j = std::max(n/2, jbeginrighthelp-jbeginrighthelp/2);
	   j <= std::min(n-nr-1, ((i+offsetR_+2*(n-nr-1)+bandR_.size())-m)/2); j++)
	help += bandR_[i-(m-offsetR_-2*(n-nr-1-j)-bandR_.size())] * x[j];

      // contribution from lower right corner block
      if (i >= m-mr)
	for (size_type j(0); j < nr; j++)
	  help += MR_.get_entry(i-(m-mr), j) * x[n-nr+j];
      
      Mx[i] = help;
    }
  }

  template <class C>
  template <class VECTOR>
  void QuasiStationaryMatrix<C>::apply_transposed(const VECTOR& x, VECTOR& Mtx) const
  {
    assert(Mtx.size() == column_dimension());
    
    // for readability:
    const unsigned int ml = ML_.row_dimension();
    const unsigned int nl = ML_.column_dimension();
    const unsigned int mr = MR_.row_dimension();
    const unsigned int nr = MR_.column_dimension();
    const unsigned int m = row_dimension();
    const unsigned int n = column_dimension();

    // contribution from upper left corner block
    for (size_type i(0); i < nl; i++) {
      C help(0);
      for (size_type j(0); j < ml; j++) {
	help += ML_.get_entry(j, i) * x[j];
      }
      Mtx[i] = help;
    }

    // contribution from left band
    for (size_type i(nl); i < n/2; i++) {
      C help(0);
      const unsigned int jbegin = offsetL_+2*(i-nl);
      for (size_type j(jbegin); j < jbegin+bandL_.size(); j++)
	help += bandL_[j-jbegin] * x[j];
      Mtx[i] = help;
    }
    
    // contribution from right band
    for (size_type i(n/2); i < n-nr; i++) {
      C help(0);
      const unsigned int jendplus1 = m-offsetR_-2*(n-nr-1-i);
      for (size_type j(jendplus1-bandR_.size()); j < jendplus1; j++)
	help += bandR_[(j+bandR_.size())-jendplus1] * x[j];
      Mtx[i] = help;
    }

    // contribution from lower right corner block
    for (size_type i(n-nr); i < n; i++) {
      C help(0);
      for (size_type j(m-mr); j < m; j++) {
	help += MR_.get_entry(j-(m-mr), i-(n-nr)) * x[j];
      }
      Mtx[i] = help;
    }
    
  }
  
  template <class C>
  void
  QuasiStationaryMatrix<C>::to_sparse(SparseMatrix<C>& S) const
  {
    // for readability:
    const unsigned int nl = ML_.column_dimension();
    const unsigned int mr = MR_.row_dimension();
    const unsigned int nr = MR_.column_dimension();
    const unsigned int m = row_dimension();
    const unsigned int n = column_dimension();

    S.resize(m, n);

    // corner blocks
    S.set_block(0, 0, ML_);
    S.set_block(m-mr, n-nr, MR_);
    
    for (unsigned int i(0); i < m; i++) {
      // left band
      const unsigned int jbeginlefthelp = (i+2*nl+1)-bandL_.size()-offsetL_;
      for (size_type j = std::max(nl, jbeginlefthelp-jbeginlefthelp/2);
	   j <= std::min(n/2-1, ((i+2*nl)-offsetL_)/2); j++)
	S.set_entry(i, j, bandL_[i-offsetL_-2*(j-nl)]);
      
      // right band
      const unsigned int jbeginrighthelp = (i+offsetR_+2*(n-nr-1)+1)-m;
      for (size_type j = std::max(n/2, jbeginrighthelp-jbeginrighthelp/2);
	   j <= std::min(n-nr-1, ((i+offsetR_+2*(n-nr-1)+bandR_.size())-m)/2); j++)
	S.set_entry(i, j, bandR_[i-(m-offsetR_-2*(n-nr-1-j)-bandR_.size())]);
    }
  }
    
  template <class C>
  void QuasiStationaryMatrix<C>::print(std::ostream &os,
				       const unsigned int tabwidth,
				       const unsigned int precision) const
  {
    if (row_dimension() == 0)
      os << "[]" << std::endl; // Matlab style
    else
      {
	unsigned int old_precision = os.precision(precision);
	for (typename QuasiStationaryMatrix<C>::size_type i(0); i < row_dimension(); ++i)
	  {
	    for (typename QuasiStationaryMatrix<C>::size_type j(0); j < column_dimension(); ++j)
	      os << std::setw(tabwidth) << std::setprecision(precision)
		 << get_entry(i, j);
	    os << std::endl;
	  }
	os.precision(old_precision);
      }
  }

  template <class C>
  inline
  std::ostream& operator << (std::ostream& os, const QuasiStationaryMatrix<C>& M)
  {
    M.print(os);
    return os;
  }
}
