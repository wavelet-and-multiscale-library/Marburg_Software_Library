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
	  return bandL_[row-(offsetL_+2*(column-nl))];
      } else {
	if (column >= n-nr) { // right block column region
	  if (row >= m-mr)
	    return MR_.get_entry(row-(m-mr),
				 column-(n-nr));
	} else { // right band column region
	  const unsigned int colendplus1 = m-offsetR_-2*(n-nr-1-column);
	  if (row < colendplus1 && row >= colendplus1-bandR_.size())
	    return bandR_[0];
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
