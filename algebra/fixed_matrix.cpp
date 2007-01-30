// implementation of MathTL::FixedMatrix inline functions

#include <iomanip>
#include <sstream>

namespace MathTL
{
  template <class C, unsigned int ROW_DIM, unsigned int COL_DIM>
  inline
  FixedMatrix<C, ROW_DIM, COL_DIM>::FixedMatrix(const C value)
    : entries_(value)
  {
  }

  template <class C, unsigned int ROW_DIM, unsigned int COL_DIM>
  inline
  FixedMatrix<C, ROW_DIM, COL_DIM>::FixedMatrix(const FixedMatrix<C, ROW_DIM, COL_DIM>& M)
    : entries_(M.entries_)
  {
  }

  template <class C, unsigned int ROW_DIM, unsigned int COL_DIM>
  FixedMatrix<C, ROW_DIM, COL_DIM>::FixedMatrix(const char* str, const bool byrow)
    : entries_()
  {
    setup_matrix(str, byrow);
  }

  template <class C, unsigned int ROW_DIM, unsigned int COL_DIM>
  FixedMatrix<C, ROW_DIM, COL_DIM>::FixedMatrix(const char* str, const C denominator, const bool byrow)
    : entries_()
  {
    setup_matrix(str, byrow);
    scale(C(1)/denominator);
  }

  template <class C, unsigned int ROW_DIM, unsigned int COL_DIM>
  void FixedMatrix<C, ROW_DIM, COL_DIM>::setup_matrix(const char* str, const bool byrow)
  {
    std::istringstream ins(str);
    if (byrow)
      {
        for (size_type i(0); i < ROW_DIM && ins.good(); i++)
          for (size_type j(0); j < COL_DIM && ins.good(); j++)
            ins >> (*this).operator () (i,j);
      }
    else
      {
        for (size_type j(0); j < COL_DIM && ins.good(); j++)
          for (size_type i(0); i < ROW_DIM && ins.good(); i++)
            ins >> (*this).operator () (i,j);
      }
  }

  template <class C, unsigned int ROW_DIM, unsigned int COL_DIM>
  inline
  FixedMatrix<C, ROW_DIM, COL_DIM>::~FixedMatrix()
  {
  }

  template <class C, unsigned int ROW_DIM, unsigned int COL_DIM>
  inline
  const typename FixedMatrix<C, ROW_DIM, COL_DIM>::size_type
  FixedMatrix<C, ROW_DIM, COL_DIM>::row_dimension() const
  {
    return ROW_DIM;
  }

  template <class C, unsigned int ROW_DIM, unsigned int COL_DIM>
  inline
  const typename FixedMatrix<C, ROW_DIM, COL_DIM>::size_type
  FixedMatrix<C, ROW_DIM, COL_DIM>::column_dimension() const
  {
    return COL_DIM;
  }

  template <class C, unsigned int ROW_DIM, unsigned int COL_DIM>
  inline
  const typename FixedMatrix<C, ROW_DIM, COL_DIM>::size_type
  FixedMatrix<C, ROW_DIM, COL_DIM>::size() const
  {
    return row_dimension()*column_dimension();
  }

  template <class C, unsigned int ROW_DIM, unsigned int COL_DIM>
  inline
  const typename FixedMatrix<C, ROW_DIM, COL_DIM>::size_type
  FixedMatrix<C, ROW_DIM, COL_DIM>::memory_consumption() const
  {
    return sizeof(*this) + row_dimension()*column_dimension()*sizeof(C);
  }

  template <class C, unsigned int ROW_DIM, unsigned int COL_DIM>
  inline
  bool FixedMatrix<C, ROW_DIM, COL_DIM>::empty() const
  {
    return size()==0;
  }

  template <class C, unsigned int ROW_DIM, unsigned int COL_DIM>
  inline
  const C FixedMatrix<C, ROW_DIM, COL_DIM>::operator () (const size_type row, const size_type column) const
  {
    return entries_[row+column*ROW_DIM];
  }

  template <class C, unsigned int ROW_DIM, unsigned int COL_DIM>
  inline
  const C FixedMatrix<C, ROW_DIM, COL_DIM>::get_entry(const size_type row, const size_type column) const
  {
    return this->operator () (row, column);
  }

  template <class C, unsigned int ROW_DIM, unsigned int COL_DIM>
  inline
  C& FixedMatrix<C, ROW_DIM, COL_DIM>::operator () (const size_type row, const size_type column)
  {
    return entries_[row+column*ROW_DIM];
  }

  template <class C, unsigned int ROW_DIM, unsigned int COL_DIM>
  inline
  void FixedMatrix<C, ROW_DIM, COL_DIM>::set_entry(const size_type row, const size_type column, const C value)
  {
    this->operator () (row, column) = value;
  }

  template <class C, unsigned int ROW_DIM, unsigned int COL_DIM>
  template <class C2>
  bool FixedMatrix<C, ROW_DIM, COL_DIM>::operator == (const FixedMatrix<C2, ROW_DIM, COL_DIM>& M) const
  {
    if (ROW_DIM != M.row_dimension() || COL_DIM != M.column_dimension()) return false;
    return std::equal(entries_.begin(), entries_.end(), M.entries_.begin());
  }
  
  template <class C, unsigned int ROW_DIM, unsigned int COL_DIM>
  template <class C2>
  inline
  bool FixedMatrix<C, ROW_DIM, COL_DIM>::operator != (const FixedMatrix<C2, ROW_DIM, COL_DIM>& M) const
  {
    return !((*this) == M);
  }

  template <class C, unsigned int ROW_DIM, unsigned int COL_DIM>
  FixedMatrix<C, ROW_DIM, COL_DIM>& FixedMatrix<C, ROW_DIM, COL_DIM>::operator = (const FixedMatrix<C, ROW_DIM, COL_DIM>& M)
  {
    entries_ = M.entries_;

    return *this;
  }

  template <class C, unsigned int ROW_DIM, unsigned int COL_DIM>
  inline
  void FixedMatrix<C, ROW_DIM, COL_DIM>::scale(const C s)
  {
    entries_.scale(s);
  }

  template <class C, unsigned int ROW_DIM, unsigned int COL_DIM>
  template <class VECTOR>
  void FixedMatrix<C, ROW_DIM, COL_DIM>::apply(const VECTOR& x, VECTOR& Mx) const
  {
    assert(Mx.size() == ROW_DIM);
    
    for (typename FixedMatrix<C, ROW_DIM, COL_DIM>::size_type i(0); i < ROW_DIM; i++)
      {
        Mx[i] = 0;
        for (typename FixedMatrix<C, ROW_DIM, COL_DIM>::size_type j(0);
             j < COL_DIM; j++)
          Mx[i] += this->operator () (i, j) * x[j];
      }
  }

  template <class C, unsigned int ROW_DIM, unsigned int COL_DIM>
  void FixedMatrix<C, ROW_DIM, COL_DIM>::apply(const FixedVector<C, COL_DIM>& x, FixedVector<C, ROW_DIM>& Mx) const
  {
    assert(Mx.size() == ROW_DIM);
    
    for (typename FixedMatrix<C, ROW_DIM, COL_DIM>::size_type i(0); i < ROW_DIM; i++)
      {
        Mx[i] = 0;
        for (typename FixedMatrix<C, ROW_DIM, COL_DIM>::size_type j(0);
             j < COL_DIM; j++)
          Mx[i] += this->operator () (i, j) * x[j];
      }
  }

  template <class C, unsigned int ROW_DIM, unsigned int COL_DIM>
  template <class VECTOR>
  void FixedMatrix<C, ROW_DIM, COL_DIM>::apply_transposed(const VECTOR& x, VECTOR& Mtx) const
  {
    assert(Mtx.size() == COL_DIM);
    
    for (typename FixedMatrix<C, ROW_DIM, COL_DIM>::size_type i(0); i < COL_DIM; i++)
      {
        Mtx[i] = 0;
        for (typename FixedMatrix<C, ROW_DIM, COL_DIM>::size_type j(0);
             j < ROW_DIM; j++)
          Mtx[i] += this->operator () (j, i) * x[j];
      }
  }

  template <class C, unsigned int ROW_DIM, unsigned int COL_DIM>
  void FixedMatrix<C, ROW_DIM, COL_DIM>::apply_transposed(const FixedVector<C, ROW_DIM>& x, FixedVector<C, COL_DIM>& Mtx) const
  {
    assert(Mtx.size() == COL_DIM);
    
    for (typename FixedMatrix<C, ROW_DIM, COL_DIM>::size_type i(0); i < COL_DIM; i++)
      {
        Mtx[i] = 0;
        for (typename FixedMatrix<C, ROW_DIM, COL_DIM>::size_type j(0);
             j < ROW_DIM; j++)
          Mtx[i] += this->operator () (j, i) * x[j];
      }
  }
  
  template <class C, unsigned int ROW_DIM, unsigned int COL_DIM>
  void FixedMatrix<C, ROW_DIM, COL_DIM>::print(std::ostream &os,
                        const unsigned int tabwidth,
                        const unsigned int precision) const
  {
    if (empty())
      os << "[]" << std::endl; // Matlab style
    else
      {
        unsigned int old_precision = os.precision(precision);
        for (typename FixedMatrix<C, ROW_DIM, COL_DIM>::size_type i(0); i < row_dimension(); ++i)
          {
            for (typename FixedMatrix<C, ROW_DIM, COL_DIM>::size_type j(0); j < column_dimension(); ++j)
              os << std::setw(tabwidth) << std::setprecision(precision)
                 << this->operator () (i, j);
            os << std::endl;
          }
        os.precision(old_precision);
      }
  }
  
  template <class C, unsigned int ROW_DIM, unsigned int COL_DIM>
  FixedMatrix<C, ROW_DIM, COL_DIM> operator - (const FixedMatrix<C, ROW_DIM, COL_DIM>& M, const FixedMatrix<C, ROW_DIM, COL_DIM>& N)
  {
    assert(M.column_dimension() == N.column_dimension());
    assert(M.row_dimension() == N.row_dimension());
    typedef typename FixedMatrix<C, ROW_DIM, COL_DIM>::size_type size_type;

    FixedMatrix<C, ROW_DIM, COL_DIM> R;
    for (size_type i(0); i < M.row_dimension(); i++)
      for (size_type j(0); j < M.column_dimension(); j++)
        {
          R(i, j) = M(i, j) - N(i, j);
        }

    return R;
  }

  template <class C, unsigned int ROW_DIM, unsigned int COL_DIM, unsigned int COL_DIM2>
  FixedMatrix<C, ROW_DIM, COL_DIM> operator * (const FixedMatrix<C, ROW_DIM, COL_DIM>& M, const FixedMatrix<C, COL_DIM, COL_DIM2>& N)
  {
    assert(M.column_dimension() == N.row_dimension());
    typedef typename FixedMatrix<C, ROW_DIM, COL_DIM>::size_type size_type;

    FixedMatrix<C, ROW_DIM, COL_DIM2> R;
    for (size_type i(0); i < M.row_dimension(); i++)
      for (size_type j(0); j < N.column_dimension(); j++)
        {
          double help(0);
          for (size_type k(0); k < N.row_dimension(); k++)
            help += M(i, k) * N(k, j);
          R(i, j) = help;
        }

    return R;
  }

  template <class C, unsigned int ROW_DIM, unsigned int COL_DIM>
  FixedMatrix<C, COL_DIM, ROW_DIM> transpose(const FixedMatrix<C, ROW_DIM, COL_DIM>& M)
  {
    typedef typename FixedMatrix<C, ROW_DIM, COL_DIM>::size_type size_type;

    FixedMatrix<C, COL_DIM, ROW_DIM> R;
    for (size_type i(0); i < M.row_dimension(); i++)
      for (size_type j(0); j < M.column_dimension(); j++)
        R(j, i) = M(i, j);

    return R;
  }

  template <class C, unsigned int ROW_DIM, unsigned int COL_DIM>
  inline
  std::ostream& operator << (std::ostream& os, const FixedMatrix<C, ROW_DIM, COL_DIM>& M)
  {
    M.print(os);
    return os;
  }
}
