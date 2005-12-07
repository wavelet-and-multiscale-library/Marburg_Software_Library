// implementation of MathTL::SparseMatrix inline functions

#include <cassert>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <algebra/vector.h>
#include <algebra/matrix.h>

using namespace std;

namespace MathTL
{
  template <class C>
  SparseMatrix<C>::SparseMatrix(const size_type m)
    : rowdim_(m), coldim_(m)
  {
    assert(m >= 1);

    entries_ = new C*[m];
    indices_ = new size_type*[m];
    
    assert(entries_ != NULL);
    assert(indices_ != NULL);

    for (size_type row(0); row < m; row++) {
      entries_[row] = NULL;
      indices_[row] = NULL;
    }
  }

  template <class C>
  SparseMatrix<C>::SparseMatrix(const size_type m, const size_type n)
    : rowdim_(m), coldim_(n)
  {
    assert(m >= 1 && n >= 1);

    entries_ = new C*[m];
    indices_ = new size_type*[m];
    
    assert(entries_ != NULL);
    assert(indices_ != NULL);

    for (size_type row(0); row < m; row++) {
      entries_[row] = NULL;
      indices_[row] = NULL;
    }
  }

  template <class C>
  SparseMatrix<C>::SparseMatrix(const SparseMatrix<C>& M)
    : rowdim_ (M.row_dimension()), coldim_(M.column_dimension())
  {
    if (this != &M) {
      entries_ = new C*[rowdim_];
      indices_ = new size_type*[rowdim_];
      
      assert(entries_ != NULL);
      assert(indices_ != NULL);
      
      for (size_type row(0); row < rowdim_; row++) {
	if (M.indices_[row]) {
	  indices_[row] = new size_type[M.indices_[row][0]+1];
	  assert(indices_[row] != NULL);
	  for (size_type j(0); j <= M.indices_[row][0]; j++)
	    indices_[row][j] = M.indices_[row][j];
	  
	  entries_[row] = new C[indices_[row][0]];
	  assert(entries_[row] != NULL);
	  for (size_type j(0); j < indices_[row][0]; j++)
	    entries_[row][j] = M.entries_[row][j];
	} else {
	  indices_[row] = NULL;
	  entries_[row] = NULL;
	}
      }
    }
  }

  template <class C>
  SparseMatrix<C>::~SparseMatrix()
  {
//     cout << "~SparseMatrix() called for matrix" << endl;
//     cout << *this << endl;
//     cout << "indices_[0]=" << indices_[0] << endl;
    kill();
//     cout << "... done!" << endl;
  }

  template <class C>
  void SparseMatrix<C>::kill()
  {
    for (size_type row(0); row < rowdim_; row++) {
      if (indices_[row]) {
	delete[] indices_[row];
	delete[] entries_[row];
      }
    }
    delete[] indices_;
    delete[] entries_;
  }
  
  template <class C>
  inline
  const typename SparseMatrix<C>::size_type
  SparseMatrix<C>::row_dimension() const
  {
    return rowdim_;
  }
  
  template <class C>
  inline
  const typename SparseMatrix<C>::size_type
  SparseMatrix<C>::column_dimension() const
  {
    return coldim_;
  }

  template <class C>
  const typename SparseMatrix<C>::size_type
  SparseMatrix<C>::size() const
  {
    size_type nz(0);
    for (size_type row(0); row < rowdim_; row++)
      if (indices_[row])
	nz += indices_[row][0];

    return nz;
  }

  template <class C>
  void SparseMatrix<C>::resize(const size_type m, const size_type n)
  {
    assert(m >= 1 && n >= 1);

    kill();

    rowdim_ = m;
    coldim_ = n;

    indices_ = new size_type*[m];
    entries_ = new C*[m];
    
    assert(indices_ != NULL);
    assert(entries_ != NULL);

    for (size_type row(0); row < m; row++) {
      indices_[row] = NULL;
      entries_[row] = NULL;
    }
  }

  template <class C>
  inline
  const typename SparseMatrix<C>::size_type
  SparseMatrix<C>::entries_in_row(const size_type row) const
  {
    assert(row < rowdim_);
    return indices_[row] ? indices_[row][0] : 0;
  }
  
  template <class C>
  inline
  const typename SparseMatrix<C>::size_type
  SparseMatrix<C>::get_nth_index(const size_type row, const size_type n) const
  {
    return indices_[row][n+1];
  }
  
  template <class C>
  inline
  const C SparseMatrix<C>::get_nth_entry(const size_type row, const size_type n) const
  {
    return entries_[row][n];
  }
  
  template <class C>
  const C SparseMatrix<C>::get_entry(const size_type row, const size_type column) const
  {
    assert(row < rowdim_ && column < coldim_);

    if (indices_[row])
      {
	for (size_type k(1); k <= indices_[row][0]; k++)
	  {
	    if (indices_[row][k] == column)
	      return entries_[row][k-1];
	  }
      }

    return 0;
  }

  template <class C>
  template <class MATRIX>
  void SparseMatrix<C>::get_block(const size_type firstrow, const size_type firstcolumn,
				  const size_type rows, const size_type columns,
				  MATRIX& M) const
  {
    assert(firstrow+rows <= rowdim_ && firstcolumn+columns <= coldim_);

    M.resize(rows, columns);
    for (size_type i(0); i < rows; i++) // brute force
      for (size_type j(0); j < columns; j++)
	{
	  M.set_entry(i, j, get_entry(firstrow+i, firstcolumn+j));
	}
  }
  
  template <class C>
  void SparseMatrix<C>::get_row(const size_type row, InfiniteVector<C, size_type>& v, const size_type offset) const
  {
    assert(row < rowdim_);

    v.clear();
    if (indices_[row])
      {
	for (size_type k(1); k <= indices_[row][0]; k++)
	  v.set_coefficient(indices_[row][k]+offset, entries_[row][k-1]);
      }
  }

  template <class C>
  void SparseMatrix<C>::set_entry(const size_type row, const size_type column,
 				  const C value)
  {
    assert(row < rowdim_ && column < coldim_);

    if (indices_[row]) // there are entries in the desired row
      {
	size_type ind(0);
	for (size_type k(1); k <= indices_[row][0]; k++)
	  if (indices_[row][k] == column)
	    ind = k; // can be sped up a bit (by immediate break after a positive search result)
	
	if (ind)
	  entries_[row][ind-1] = value;
	else {
	  ind = 0; // redundant
	  
	  size_type* helpi = new size_type[indices_[row][0]+2];
	  C*         helpd = new C        [indices_[row][0]+1];
	  
	  assert(helpi != NULL);
	  assert(helpd != NULL);

	  for (size_type k(1); k <= indices_[row][0]; k++)
	    if (indices_[row][k] < column)
	      ind = k;
	  for (size_type k(1); k <= ind; k++) {
	    helpi[k]   = indices_[row][k];
	    helpd[k-1] = entries_[row][k-1];
	  }
	  for (size_type k(indices_[row][0]); k >= ind+1; k--) {
	    helpi[k+1] = indices_[row][k];
	    helpd[k]   = entries_[row][k-1];
	  }
	  helpi[ind+1] = column;
	  helpi[0] = indices_[row][0]+1;
	  
	  delete[] indices_[row];
	  delete[] entries_[row];

	  indices_[row] = helpi;
	  entries_[row] = helpd;

	  entries_[row][ind] = value;
	}
      }
    else
      {
	indices_[row] = new size_type[2];
	entries_[row] = new C[1]; // TODO: check pointers

	assert(indices_[row] != NULL);
	assert(entries_[row] != NULL);

	indices_[row][0] = 1;
	indices_[row][1] = column;
	entries_[row][0] = value;
      }
  }

  template <class C>
  void
  SparseMatrix<C>::set_row(const size_type row,
			   const std::list<size_type>& indices,
			   const std::list<C>& entries)
  {
    assert(row < rowdim_);
    assert(indices.size() == entries.size());
    
    if (indices_[row]) {
      delete[] indices_[row];
      delete[] entries_[row];
    }

    if (indices.size() > 0) {
      indices_[row] = new size_type[indices.size()+1];
      assert(indices_[row] != NULL);
      indices_[row][0] = indices.size();
      size_type id = 1;
      for (typename std::list<size_type>::const_iterator it(indices.begin()), itend(indices.end());
	   it != itend; ++it, ++id)
	indices_[row][id] = *it;
      
      entries_[row] = new C[indices_[row][0]];
      assert(entries_[row] != NULL);
      id = 0;
      for (typename std::list<C>::const_iterator it(entries.begin()), itend(entries.end());
	   it != itend; ++it, ++id)
	entries_[row][id] = *it;
    } else {
      entries_[row] = NULL;
      indices_[row] = NULL;
    }
  }

  template <class C>
  template <class MATRIX>
  void SparseMatrix<C>::set_block(const size_type firstrow, const size_type firstcolumn,
				  const MATRIX& M, const bool mirror)
  {
    assert(firstrow+M.row_dimension() <= row_dimension());
    assert(firstcolumn+M.column_dimension() <= column_dimension());

    // the following code can be optimized
    for (size_type row(0); row < M.row_dimension(); row++)
      for (size_type column(0); column < M.column_dimension(); column++)
	{
	  if (mirror)
	    set_entry(firstrow+M.row_dimension()-1-row,
		      firstcolumn+M.column_dimension()-1-column,
		      M.get_entry(row, column));
	  else
	    set_entry(row+firstrow, column+firstcolumn, M.get_entry(row, column));
	}
  }

  template <class C>
  SparseMatrix<C>&
  SparseMatrix<C>::operator = (const SparseMatrix<C>& M)
  {
    if (this != &M) {
      resize(M.row_dimension(), M.column_dimension());
      
      for (size_type row(0); row < rowdim_; row++) {
	if (M.indices_[row]) {
	  indices_[row] = new size_type[M.indices_[row][0]+1];
	  assert(indices_[row] != NULL);
	  for (size_type j(0); j <= M.indices_[row][0]; j++)
	    indices_[row][j] = M.indices_[row][j];
	  
	  entries_[row] = new C[indices_[row][0]];
	  assert(entries_[row] != NULL);
	  for (size_type j(0); j < indices_[row][0]; j++)
	    entries_[row][j] = M.entries_[row][j];
	} else {
	  indices_[row] = NULL;
	  entries_[row] = NULL;
	}
      }
    }

    return *this;
  }

  template <class C>
  void SparseMatrix<C>::diagonal(const size_type n, const C diag)
  {
    resize(n, n);

    for (size_type row(0); row < n; row++)
      set_entry(row, row, diag);
  }

  template <class C>
  template <class VECTOR>
  void SparseMatrix<C>::apply(const VECTOR& x, VECTOR& Mx) const
  {
    assert(Mx.size() == rowdim_);
    
    for (size_type i(0); i < rowdim_; i++) {
      C help(0);
      if (indices_[i]) {
	for (size_type j(1); j <= indices_[i][0]; j++)
	  help += entries_[i][j-1] * x[indices_[i][j]];
      }
      Mx[i] = help;
    }
  }

  template <class C>
  template <class VECTOR>
  void SparseMatrix<C>::apply_transposed(const VECTOR& x, VECTOR& Mtx) const
  {
    assert(Mtx.size() == coldim_);
    
    for (size_type i(0); i < coldim_; i++)
      Mtx[i] = C(0);

    for (size_type i(0); i < rowdim_; i++) {
      if (indices_[i]) {
	for (size_type j(1); j <= indices_[i][0]; j++)
	  Mtx[indices_[i][j]] += entries_[i][j-1] * x[i];
      }
    }
  }

  template <class C>
  void SparseMatrix<C>::compress(const double eta)
  {
    SparseMatrix<C> M(rowdim_, coldim_);

    for (size_type i(0); i < rowdim_; i++) {
      if (indices_[i]) {
	for (size_type j(1); j <= indices_[i][0]; j++)
	  if (fabs(entries_[i][j-1])>=eta)
	    M.set_entry(i, indices_[i][j], entries_[i][j-1]);
      }
    }

    this->operator = (M);
  }

  template <class C>
  void SparseMatrix<C>::print(std::ostream &os,
			      const unsigned int tabwidth,
			      const unsigned int precision) const
  {
    if (row_dimension() == 0)
      os << "[]" << std::endl; // Matlab style
    else
      {
	unsigned int old_precision = os.precision(precision);
	for (typename SparseMatrix<C>::size_type i(0); i < row_dimension(); ++i)
	  {
	    for (typename SparseMatrix<C>::size_type j(0); j < column_dimension(); ++j)
	      os << std::setw(tabwidth) << std::setprecision(precision)
		 << get_entry(i, j);
	    os << std::endl;
	  }
	os.precision(old_precision);
      }
  }

  template <class C>
  void SparseMatrix<C>::matlab_output (char *file, char *Matrixname, int binary) const
  {
    if (binary)
      {
	unsigned int i,j;

	char Filename[200];

	Filename[0] = '\x0';

	strcat(Filename, file);
	strcat(Filename, ".m");

	std::ofstream m_file(Filename);


	Filename[0] = '\x0';

	strcat(Filename, file);
	strcat(Filename, ".bin");

	std::ofstream bin_file(Filename);

  
	m_file << "fid_vector=fopen('" << file 
	       << ".bin', 'r');" << endl
	       << "n=fread(fid_vector, 1, 'int');" << endl 
	       << "m=fread(fid_vector, 1, 'int');" << endl
	       << "l=fread(fid_vector, 1, 'int');" << endl
	       << "i=fread(fid_vector, l, 'int');" << endl
	       << "j=fread(fid_vector, l, 'int');" << endl
	       << "s=fread(fid_vector, l, 'double');" << endl
	       << Matrixname << "=sparse(i,j,s,n,m);" << endl
	       << "fclose(fid_vector);" << endl;

	int r=row_dimension(), c=column_dimension(), l=size();

	bin_file.write((char*)(&r), sizeof(int));
	bin_file.write((char*)(&c), sizeof(int));
	bin_file.write((char*)(&l), sizeof(int));

	for (i=0; i<row_dimension(); i++)
	  {
	    if (indices_[i])
	      {
		int ii=i+1;
	       
		for (j=1; j<=indices_[i][0]; j++)
		  bin_file.write((char*)(&ii), sizeof(int));
	      }
	  }   

	for (i=0; i<row_dimension(); i++)
	  {
	    if (indices_[i])
	      {               
		for (j=1; j<=indices_[i][0]; j++)
		  {
                    int jj=indices_[i][j]+1;
                    bin_file.write((char*)(&jj), sizeof(int));
		  }
	      }
	  }   

        for (i=0; i<row_dimension(); i++)
	  {
            if (indices_[i])
	      bin_file.write((char*)(entries_[i]), indices_[i][0]*sizeof(C));
	  }   
      }
  else
  {
    unsigned int i;
	 
    char *Filename = new char[200];
    Filename[0] = '\x0';
	 
    strcat(Filename, file);
    strcat(Filename, ".m");
	 
    std::ofstream s;
    s.open(Filename);
	 
    delete[] Filename;
	 
	 
    s.setf(std::ios_base::scientific, std::ios_base::fixed);
    s.precision(12);
	 
    s << Matrixname << "=sparse(" << row_dimension() << "," << column_dimension() << ");" << endl;
	
    for (i=0; i<row_dimension(); i++) {
      if (indices_[i]) {
	for (size_type k(1); k <= indices_[i][0]; k++) {
	  s << Matrixname << "(" << i+1 << "," << indices_[i][k]+1 << ")="
	    << entries_[i][k-1] << ";" << endl;
	}
      }
    }   
  }
}
  

  template <class C>
  SparseMatrix<C> operator - (const SparseMatrix<C>& M, const SparseMatrix<C>& N)
  {
    assert(M.column_dimension() == N.column_dimension());
    assert(M.row_dimension() == N.row_dimension());
    typedef typename Matrix<C>::size_type size_type;

    SparseMatrix<C> R(M.row_dimension(), M.column_dimension());
    for (size_type i(0); i < M.row_dimension(); i++)
      for (size_type j(0); j < M.column_dimension(); j++)
	{
	  const double help = M.get_entry(i, j) - N.get_entry(i, j);
	  if (help != 0)
	    R.set_entry(i, j, help);
	}

    return R;
  }
  
  template <class C>
  SparseMatrix<C> operator * (const SparseMatrix<C>& M, const SparseMatrix<C>& N)
  {
    assert(M.column_dimension() == N.row_dimension());
    typedef typename SparseMatrix<C>::size_type size_type;

    SparseMatrix<C> R(M.row_dimension(), N.column_dimension());
    for (size_type i(0); i < M.row_dimension(); i++)
      for (size_type j(0); j < N.column_dimension(); j++)
	{
	  double help(0);
	  for (size_type k(0); k < N.row_dimension(); k++)
	    help += M.get_entry(i, k) * N.get_entry(k, j);

	  if (help != 0)
	    R.set_entry(i, j, help);
	}

    return R;
  }

  template <class C>
  SparseMatrix<C> transpose(const SparseMatrix<C>& M)
  {
    typedef typename SparseMatrix<C>::size_type size_type;

    SparseMatrix<C> R(M.column_dimension(), M.row_dimension());
    for (size_type i(0); i < M.row_dimension(); i++) // TODO: speedup this brute force hack
      for (size_type j(0); j < M.column_dimension(); j++)
	{
	  double help(M.get_entry(i, j));
	  if (help != 0)
	    R.set_entry(j, i, help);
	}

    return R;
  }

  template <class C>
  inline
  std::ostream& operator << (std::ostream& os, const SparseMatrix<C>& M)
  {
    M.print(os);
    return os;
  }
}
