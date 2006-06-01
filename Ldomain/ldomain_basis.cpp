// implementation for ldomain_basis.h

#include <cmath>
#include <iostream>

using std::cout;
using std::endl;

namespace WaveletTL
{
  template <class IBASIS>
  LDomainBasis<IBASIS>::LDomainBasis()
    : basis1d_(false, false)
  {
  }

  template <class IBASIS>
  const int
  LDomainBasis<IBASIS>::Deltasize(const int j) const {
    const unsigned int Deltaj = basis1d().Deltasize(j);
    return 3*(Deltaj-2)*(Deltaj-2)+2*(Deltaj-2);
  }

  template <class IBASIS>
  const int
  LDomainBasis<IBASIS>::Nabla01size(const int j) const {
    const unsigned int Deltaj = basis1d().Deltasize(j);
    return 3*(Deltaj-2)*(1<<j)+(1<<j);
  }
  
  template <class IBASIS>
  const int
  LDomainBasis<IBASIS>::Nabla10size(const int j) const {
    const unsigned int Deltaj = basis1d().Deltasize(j);
    return 3*(Deltaj-2)*(1<<j)+(1<<j);
  }
  
  template <class IBASIS>
  const int
  LDomainBasis<IBASIS>::Nabla11size(const int j) const {
    return 3*(1<<(2*j));
  }
  
  template <class IBASIS>
  typename LDomainBasis<IBASIS>::Index
  LDomainBasis<IBASIS>::first_generator(const int j) const
  {
    assert(j >= j0());

    typename Index::type_type e;

    // setup lowest translation index for e=(0,0), p=0
    typename Index::translation_type k;
    k[0] = k[1] = basis1d().DeltaLmin()+1;
    
    return Index(j, e, 0, k, this);
  }

  template <class IBASIS>
  typename LDomainBasis<IBASIS>::Index
  LDomainBasis<IBASIS>::last_generator(const int j) const
  {
    assert(j >= j0());

    typename Index::type_type e;

    // setup highest translation index for e=(0,0), p=4
    typename Index::translation_type k; // k[0]=0 by convention
    k[1] = basis1d().DeltaRmax(j)-1;
    
    return Index(j, e, 4, k, this);
  }

  template <class IBASIS>
  typename LDomainBasis<IBASIS>::Index
  LDomainBasis<IBASIS>::first_wavelet(const int j) const
  {
    assert(j >= j0());

    typename Index::type_type e;
    e[1] = 1;

    // setup lowest translation index for e=(0,1), p=0
    typename Index::translation_type k;
    k[0] = basis1d().DeltaLmin()+1;
    k[1] = basis1d().Nablamin();
    
    return Index(j, e, 0, k, this);
  }

  template <class IBASIS>
  typename LDomainBasis<IBASIS>::Index
  LDomainBasis<IBASIS>::first_wavelet(const int j, const typename Index::type_type& ewish) const
  {
    assert(j >= j0());
    
    typename Index::type_type e(ewish);
    
    // setup lowest translation index appropriately
    typename Index::translation_type k;
    const int ecode(e[0]+2*e[1]);
    if (ecode == 0) {
      // e = (0,0)
      k[0] = k[1] = basis1d().DeltaLmin()+1;
    } else {
      if (ecode == 1) {
	// e = (1,0)
	k[0] = basis1d().Nablamin();
	k[1] = basis1d().DeltaLmin()+1;
      } else {
	if (ecode == 2) {
	  // e = (0,1)
	  k[0] = basis1d().DeltaLmin()+1;
	  k[1] = basis1d().Nablamin();
	} else {
	  // e = (1,1)
	  k[0] = k[1] = basis1d().Nablamin();
	}
      }
    }
    
    return Index(j, e, 0, k, this);

  }

  template <class IBASIS>
  typename LDomainBasis<IBASIS>::Index
  LDomainBasis<IBASIS>::last_wavelet(const int j) const
  {
    assert(j >= j0());
    
    typename Index::type_type e;
    e[0] = e[1] = 1;

    // setup highest translation index for e=(1,1), p=2
    typename Index::translation_type k;
    k[0] = k[1] = basis1d().Nablamax(j);
    
    return Index(j, e, 2, k, this);
  }

  template <class IBASIS>
  const BlockMatrix<double>&
  LDomainBasis<IBASIS>::get_Mj0 (const int j) const {
    // check whether the j-th matrix already exists in the cache
    typename MatrixCache::iterator matrix_lb(Mj0_cache.lower_bound(j));
    typename MatrixCache::iterator matrix_it(matrix_lb);
    if (matrix_lb == Mj0_cache.end() ||
	Mj0_cache.key_comp()(j, matrix_lb->first))
      {
// 	cout << "LDomainBasis::get_Mj0() cache miss" << endl;

	// compute Mj0 and insert it into the cache
	BlockMatrix<double> Mj0(5, 5);
	typedef typename MatrixCache::value_type value_type;
 	matrix_it = Mj0_cache.insert(matrix_lb, value_type(j, Mj0));
	
	const unsigned int Deltaj   = basis1d().Deltasize(j);
	const unsigned int Deltajp1 = basis1d().Deltasize(j+1);

	// row/column 0,1,2 <-> patches 0,1,2
	for (unsigned int patch = 0; patch <= 2; patch++) {
	  matrix_it->second.resize_block_row   (patch, (Deltajp1-2)*(Deltajp1-2));
	  matrix_it->second.resize_block_column(patch, (Deltaj-2)*(Deltaj-2));
	}
	
 	// row/column 3,4 <-> interface 3,4
	for (unsigned int patch = 3; patch <= 4; patch++) {
	  matrix_it->second.resize_block_row   (patch, Deltajp1-2);
	  matrix_it->second.resize_block_column(patch, Deltaj-2);
	}

 	// prepare 1d matrices
 	SparseMatrix<double> Mj0_1d; basis1d().assemble_Mj0(j, Mj0_1d);
   	SparseMatrix<double> Mj0_1d_interior(Deltajp1-2, Deltaj-2);
 	for (unsigned int row = 0; row < Deltajp1-2; row++)
 	  for (unsigned int column = 0; column < Deltaj-2; column++)
 	    Mj0_1d_interior.set_entry(row, column, Mj0_1d.get_entry(row+1, column+1));
	SparseMatrix<double> Mj0_1d_left(Deltajp1-2, 1);
 	for (unsigned int row = 0; row < Deltajp1-2; row++)
	  Mj0_1d_left.set_entry(row, 0, Mj0_1d.get_entry(row+1, 0));
	SparseMatrix<double> Mj0_1d_right(Deltajp1-2, 1);
 	for (unsigned int row = 0; row < Deltajp1-2; row++)
	  Mj0_1d_right.set_entry(row, 0, Mj0_1d.get_entry(row+1, Deltaj-1));
	SparseMatrix<double> Mj0_1d_left_top(1, 1); Mj0_1d_left_top.set_entry(0, 0, Mj0_1d.get_entry(0, 0));
	
 	// patch generators decompose only into themselves
	for (int patch = 0; patch <= 2; patch++)
	  matrix_it->second.set_block(patch, patch,
				      new KroneckerMatrix<double,SparseMatrix<double>,SparseMatrix<double> >
				      (Mj0_1d_interior, Mj0_1d_interior));
	
	// interface generators decompose into themselves and patch generators from the neighboring patches
 	matrix_it->second.set_block(0, 3,
 				    new KroneckerMatrix<double,SparseMatrix<double>,SparseMatrix<double> >
 				    (Mj0_1d_interior, Mj0_1d_left, M_SQRT1_2));
 	matrix_it->second.set_block(1, 3,
 				    new KroneckerMatrix<double,SparseMatrix<double>,SparseMatrix<double> >
 				    (Mj0_1d_interior, Mj0_1d_right, M_SQRT1_2));
 	matrix_it->second.set_block(3, 3,
 				    new KroneckerMatrix<double,SparseMatrix<double>,SparseMatrix<double> >
 				    (Mj0_1d_left_top, Mj0_1d_interior));

 	matrix_it->second.set_block(1, 4,
 				    new KroneckerMatrix<double,SparseMatrix<double>,SparseMatrix<double> >
 				    (Mj0_1d_right, Mj0_1d_interior, M_SQRT1_2));
 	matrix_it->second.set_block(2, 4,
 				    new KroneckerMatrix<double,SparseMatrix<double>,SparseMatrix<double> >
 				    (Mj0_1d_left, Mj0_1d_interior, M_SQRT1_2));
 	matrix_it->second.set_block(4, 4,
 				    new KroneckerMatrix<double,SparseMatrix<double>,SparseMatrix<double> >
 				    (Mj0_1d_left_top, Mj0_1d_interior));
      }
    else
      {
// 	cout << "LDomainBasis::get_Mj0() cache hit" << endl;
      }
    
    return matrix_it->second;
  }
  
  template <class IBASIS>
  const BlockMatrix<double>&
  LDomainBasis<IBASIS>::get_Mj0T (const int j) const {
    // check whether the j-th matrix already exists in the cache
    typename MatrixCache::iterator matrix_lb(Mj0T_cache.lower_bound(j));
    typename MatrixCache::iterator matrix_it(matrix_lb);
    if (matrix_lb == Mj0T_cache.end() ||
	Mj0T_cache.key_comp()(j, matrix_lb->first))
      {
// 	cout << "LDomainBasis::get_Mj0T() cache miss" << endl;

	// compute Mj0T and insert it into the cache
	BlockMatrix<double> Mj0T(5, 5);
	typedef typename MatrixCache::value_type value_type;
 	matrix_it = Mj0T_cache.insert(matrix_lb, value_type(j, Mj0T));
	
	const unsigned int Deltaj   = basis1d().Deltasize(j);
	const unsigned int Deltajp1 = basis1d().Deltasize(j+1);

	// row/column 0,1,2 <-> patches 0,1,2
	for (unsigned int patch = 0; patch <= 2; patch++) {
	  matrix_it->second.resize_block_row   (patch, (Deltajp1-2)*(Deltajp1-2));
	  matrix_it->second.resize_block_column(patch, (Deltaj-2)*(Deltaj-2));
	}
	
 	// row/column 3,4 <-> interface 3,4
	for (unsigned int patch = 3; patch <= 4; patch++) {
	  matrix_it->second.resize_block_row   (patch, Deltajp1-2);
	  matrix_it->second.resize_block_column(patch, Deltaj-2);
	}

 	// prepare 1d matrices
 	SparseMatrix<double> Mj0T_1d; basis1d().assemble_Mj0T(j, Mj0T_1d);
   	SparseMatrix<double> Mj0T_1d_interior(Deltajp1-2, Deltaj-2);
 	for (unsigned int row = 0; row < Deltajp1-2; row++)
 	  for (unsigned int column = 0; column < Deltaj-2; column++)
 	    Mj0T_1d_interior.set_entry(row, column, Mj0T_1d.get_entry(row+1, column+1));
	SparseMatrix<double> Mj0T_1d_left(Deltajp1-2, 1);
 	for (unsigned int row = 0; row < Deltajp1-2; row++)
	  Mj0T_1d_left.set_entry(row, 0, Mj0T_1d.get_entry(row+1, 0));
	SparseMatrix<double> Mj0T_1d_right(Deltajp1-2, 1);
 	for (unsigned int row = 0; row < Deltajp1-2; row++)
	  Mj0T_1d_right.set_entry(row, 0, Mj0T_1d.get_entry(row+1, Deltaj-1));
	SparseMatrix<double> Mj0T_1d_left_top(1, 1); Mj0T_1d_left_top.set_entry(0, 0, Mj0T_1d.get_entry(0, 0));
	
 	// patch generators decompose only into themselves
	for (int patch = 0; patch <= 2; patch++)
	  matrix_it->second.set_block(patch, patch,
				      new KroneckerMatrix<double,SparseMatrix<double>,SparseMatrix<double> >
				      (Mj0T_1d_interior, Mj0T_1d_interior));
	
	// interface generators decompose into themselves and patch generators from the neighboring patches
 	matrix_it->second.set_block(0, 3,
 				    new KroneckerMatrix<double,SparseMatrix<double>,SparseMatrix<double> >
 				    (Mj0T_1d_interior, Mj0T_1d_left, M_SQRT1_2));
 	matrix_it->second.set_block(1, 3,
 				    new KroneckerMatrix<double,SparseMatrix<double>,SparseMatrix<double> >
 				    (Mj0T_1d_interior, Mj0T_1d_right, M_SQRT1_2));
 	matrix_it->second.set_block(3, 3,
 				    new KroneckerMatrix<double,SparseMatrix<double>,SparseMatrix<double> >
 				    (Mj0T_1d_left_top, Mj0T_1d_interior));

 	matrix_it->second.set_block(1, 4,
 				    new KroneckerMatrix<double,SparseMatrix<double>,SparseMatrix<double> >
 				    (Mj0T_1d_right, Mj0T_1d_interior, M_SQRT1_2));
 	matrix_it->second.set_block(2, 4,
 				    new KroneckerMatrix<double,SparseMatrix<double>,SparseMatrix<double> >
 				    (Mj0T_1d_left, Mj0T_1d_interior, M_SQRT1_2));
 	matrix_it->second.set_block(4, 4,
 				    new KroneckerMatrix<double,SparseMatrix<double>,SparseMatrix<double> >
 				    (Mj0T_1d_left_top, Mj0T_1d_interior));
      }
    else
      {
// 	cout << "LDomainBasis::get_Mj0T() cache hit" << endl;
      }
    
    return matrix_it->second;
  }

  template <class IBASIS>
  const SparseMatrix<double>&
  LDomainBasis<IBASIS>::get_Mj1c_1d (const int j) const {
    // check whether the j-th matrix already exists in the cache
    typename Matrix1DCache::iterator matrix_lb(Mj1c_1d_cache.lower_bound(j));
    typename Matrix1DCache::iterator matrix_it(matrix_lb);
    if (matrix_lb == Mj1c_1d_cache.end() ||
	Mj1c_1d_cache.key_comp()(j, matrix_lb->first))
      {
// 	cout << "LDomainBasis::get_Mj1c_1d() cache miss" << endl;

	// compute Mj1c_1d and insert it into the cache, it is
	//   Mj1c = (I-Mj0*<Theta_{j+1},Lambdatilde>^T)Mj1

 	const unsigned int Deltaj   = basis1d().Deltasize(j);
	const unsigned int Deltajp1 = basis1d().Deltasize(j+1);

	SparseMatrix<double> Mj0;  basis1d().assemble_Mj0 (j, Mj0 );
	SparseMatrix<double> Mj1;  basis1d().assemble_Mj1 (j, Mj1 );
	SparseMatrix<double> Mj0T; basis1d().assemble_Mj0T(j, Mj0T);
	SparseMatrix<double> Mj0T_modified(Deltajp1, Deltaj);
	Mj0T_modified.set_entry(0, 0, M_SQRT2);
	Mj0T_modified.set_entry(Deltajp1-1, Deltaj-1, M_SQRT2);
	for (unsigned int row = 0; row < Deltajp1; row++)
	  for (unsigned int column = 1; column < Deltaj-1; column++)
	    Mj0T_modified.set_entry(row, column, Mj0T.get_entry(row, column));
	SparseMatrix<double> Mj1c = Mj1 - (Mj0*transpose(Mj0T_modified)*Mj1);
	Mj1c.compress(1e-12);

	typedef typename Matrix1DCache::value_type value_type;
 	matrix_it = Mj1c_1d_cache.insert(matrix_lb, value_type(j, Mj1c));
      }
    else
      {
// 	cout << "LDomainBasis::get_Mj1c_1d() cache hit" << endl;
      }

    return matrix_it->second;
  }

  template <class IBASIS>
  const BlockMatrix<double>&
  LDomainBasis<IBASIS>::get_Mj1c_01 (const int j) const {
    // check whether the j-th matrix already exists in the cache
    typename MatrixCache::iterator matrix_lb(Mj1c_01_cache.lower_bound(j));
    typename MatrixCache::iterator matrix_it(matrix_lb);
    if (matrix_lb == Mj1c_01_cache.end() ||
	Mj1c_01_cache.key_comp()(j, matrix_lb->first))
      {
// 	cout << "LDomainBasis::get_Mj1c_01() cache miss" << endl;

 	// compute Mj1c_01 and insert it into the cache
 	BlockMatrix<double> Mj1c_01(5, 4);
 	typedef typename MatrixCache::value_type value_type;
  	matrix_it = Mj1c_01_cache.insert(matrix_lb, value_type(j, Mj1c_01));
	
 	const unsigned int Deltaj   = basis1d().Deltasize(j);
 	const unsigned int Deltajp1 = basis1d().Deltasize(j+1);

 	// row/column 0,1,2 <-> patches 0,1,2
 	for (unsigned int patch = 0; patch <= 2; patch++) {
 	  matrix_it->second.resize_block_row   (patch, (Deltajp1-2)*(Deltajp1-2));
 	  matrix_it->second.resize_block_column(patch, (Deltaj-2)*(1<<j));
 	}
	
  	// row 3/4 <-> interface 3,4
 	for (unsigned int patch = 3; patch <= 4; patch++) {
 	  matrix_it->second.resize_block_row   (patch, Deltajp1-2);
 	}

	// column 3 <-> interface 4
	matrix_it->second.resize_block_column(3, 1<<j);

 	// prepare 1d matrices
 	SparseMatrix<double> Mj0_1d; basis1d().assemble_Mj0(j, Mj0_1d);
   	SparseMatrix<double> Mj0_1d_interior(Deltajp1-2, Deltaj-2);
 	for (unsigned int row = 0; row < Deltajp1-2; row++)
 	  for (unsigned int column = 0; column < Deltaj-2; column++)
 	    Mj0_1d_interior.set_entry(row, column, Mj0_1d.get_entry(row+1, column+1));
	SparseMatrix<double> Mj0_1d_left(Deltajp1-2, 1);
 	for (unsigned int row = 0; row < Deltajp1-2; row++)
	  Mj0_1d_left.set_entry(row, 0, Mj0_1d.get_entry(row+1, 0));
	SparseMatrix<double> Mj0_1d_right(Deltajp1-2, 1);
 	for (unsigned int row = 0; row < Deltajp1-2; row++)
	  Mj0_1d_right.set_entry(row, 0, Mj0_1d.get_entry(row+1, Deltaj-1));
	SparseMatrix<double> Mj0_1d_left_top(1, 1); Mj0_1d_left_top.set_entry(0, 0, Mj0_1d.get_entry(0, 0));
	
	const SparseMatrix<double>& Mj1c_1d = get_Mj1c_1d(j);
	SparseMatrix<double> Mj1c_1d_interior(Deltajp1-2, 1<<j);
 	for (unsigned int row = 0; row < Deltajp1-2; row++)
 	  for (int column = 0; column < 1<<j; column++)
 	    Mj1c_1d_interior.set_entry(row, column, Mj1c_1d.get_entry(row+1, column));	
	
   	// patch wavelets decompose only into patch generators
  	for (int patch = 0; patch <= 2; patch++)
  	  matrix_it->second.set_block(patch, patch,
  				      new KroneckerMatrix<double,SparseMatrix<double>,SparseMatrix<double> >
  				      (Mj0_1d_interior, Mj1c_1d_interior));
	
  	// interface generators decompose into themselves and patch generators from the neighboring patches
  	matrix_it->second.set_block(1, 3,
  				    new KroneckerMatrix<double,SparseMatrix<double>,SparseMatrix<double> >
  				    (Mj0_1d_right, Mj1c_1d_interior, M_SQRT1_2));
  	matrix_it->second.set_block(2, 3,
  				    new KroneckerMatrix<double,SparseMatrix<double>,SparseMatrix<double> >
  				    (Mj0_1d_left, Mj1c_1d_interior, M_SQRT1_2));
  	matrix_it->second.set_block(4, 3,
  				    new KroneckerMatrix<double,SparseMatrix<double>,SparseMatrix<double> >
  				    (Mj0_1d_left_top, Mj1c_1d_interior));
      }
    else
      {
// 	cout << "LDomainBasis::get_Mj1c_01() cache hit" << endl;
      }
    
    return matrix_it->second;
  }

  template <class IBASIS>
  const BlockMatrix<double>&
  LDomainBasis<IBASIS>::get_Mj1c_10 (const int j) const {
    // check whether the j-th matrix already exists in the cache
    typename MatrixCache::iterator matrix_lb(Mj1c_10_cache.lower_bound(j));
    typename MatrixCache::iterator matrix_it(matrix_lb);
    if (matrix_lb == Mj1c_10_cache.end() ||
	Mj1c_10_cache.key_comp()(j, matrix_lb->first))
      {
// 	cout << "LDomainBasis::get_Mj1c_10() cache miss" << endl;

 	// compute Mj1c_10 and insert it into the cache
 	BlockMatrix<double> Mj1c_10(5, 4);
 	typedef typename MatrixCache::value_type value_type;
  	matrix_it = Mj1c_10_cache.insert(matrix_lb, value_type(j, Mj1c_10));
	
 	const unsigned int Deltaj   = basis1d().Deltasize(j);
 	const unsigned int Deltajp1 = basis1d().Deltasize(j+1);

 	// row/column 0,1,2 <-> patches 0,1,2
 	for (unsigned int patch = 0; patch <= 2; patch++) {
 	  matrix_it->second.resize_block_row   (patch, (Deltajp1-2)*(Deltajp1-2));
 	  matrix_it->second.resize_block_column(patch, (Deltaj-2)*(1<<j));
 	}
	
  	// row 3/4 <-> interface 3,4
 	for (unsigned int patch = 3; patch <= 4; patch++) {
 	  matrix_it->second.resize_block_row   (patch, Deltajp1-2);
 	}

	// column 3 <-> interface 3
	matrix_it->second.resize_block_column(3, 1<<j);

 	// prepare 1d matrices
 	SparseMatrix<double> Mj0_1d; basis1d().assemble_Mj0(j, Mj0_1d);
   	SparseMatrix<double> Mj0_1d_interior(Deltajp1-2, Deltaj-2);
 	for (unsigned int row = 0; row < Deltajp1-2; row++)
 	  for (unsigned int column = 0; column < Deltaj-2; column++)
 	    Mj0_1d_interior.set_entry(row, column, Mj0_1d.get_entry(row+1, column+1));
	SparseMatrix<double> Mj0_1d_left(Deltajp1-2, 1);
 	for (unsigned int row = 0; row < Deltajp1-2; row++)
	  Mj0_1d_left.set_entry(row, 0, Mj0_1d.get_entry(row+1, 0));
	SparseMatrix<double> Mj0_1d_right(Deltajp1-2, 1);
 	for (unsigned int row = 0; row < Deltajp1-2; row++)
	  Mj0_1d_right.set_entry(row, 0, Mj0_1d.get_entry(row+1, Deltaj-1));
	SparseMatrix<double> Mj0_1d_left_top(1, 1); Mj0_1d_left_top.set_entry(0, 0, Mj0_1d.get_entry(0, 0));
	
	const SparseMatrix<double>& Mj1c_1d = get_Mj1c_1d(j);
	SparseMatrix<double> Mj1c_1d_interior(Deltajp1-2, 1<<j);
 	for (unsigned int row = 0; row < Deltajp1-2; row++)
 	  for (int column = 0; column < 1<<j; column++)
 	    Mj1c_1d_interior.set_entry(row, column, Mj1c_1d.get_entry(row+1, column));	
	
    	// patch wavelets decompose only into themselves
   	for (int patch = 0; patch <= 2; patch++)
   	  matrix_it->second.set_block(patch, patch,
   				      new KroneckerMatrix<double,SparseMatrix<double>,SparseMatrix<double> >
   				      (Mj1c_1d_interior, Mj0_1d_interior));
	
  	// interface generators decompose into themselves and patch generators from the neighboring patches
   	matrix_it->second.set_block(0, 3,
   				    new KroneckerMatrix<double,SparseMatrix<double>,SparseMatrix<double> >
   				    (Mj1c_1d_interior, Mj0_1d_right, M_SQRT1_2));
   	matrix_it->second.set_block(1, 3,
   				    new KroneckerMatrix<double,SparseMatrix<double>,SparseMatrix<double> >
   				    (Mj1c_1d_interior, Mj0_1d_left, M_SQRT1_2));
   	matrix_it->second.set_block(3, 3,
   				    new KroneckerMatrix<double,SparseMatrix<double>,SparseMatrix<double> >
   				    (Mj0_1d_left_top, Mj1c_1d_interior));
      }
    else
      {
// 	cout << "LDomainBasis::get_Mj1c_10() cache hit" << endl;
      }
    
    return matrix_it->second;
  }

  template <class IBASIS>
  const BlockMatrix<double>&
  LDomainBasis<IBASIS>::get_Mj1c_11 (const int j) const {
    // check whether the j-th matrix already exists in the cache
    typename MatrixCache::iterator matrix_lb(Mj1c_11_cache.lower_bound(j));
    typename MatrixCache::iterator matrix_it(matrix_lb);
    if (matrix_lb == Mj1c_11_cache.end() ||
	Mj1c_11_cache.key_comp()(j, matrix_lb->first))
      {
// 	cout << "LDomainBasis::get_Mj1c_11() cache miss" << endl;

 	// compute Mj1c_11 and insert it into the cache
 	BlockMatrix<double> Mj1c_11(5, 3);
 	typedef typename MatrixCache::value_type value_type;
  	matrix_it = Mj1c_11_cache.insert(matrix_lb, value_type(j, Mj1c_11));
	
 	const unsigned int Deltajp1 = basis1d().Deltasize(j+1);

 	// row/column 0,1,2 <-> patches 0,1,2
 	for (unsigned int patch = 0; patch <= 2; patch++) {
 	  matrix_it->second.resize_block_row   (patch, (Deltajp1-2)*(Deltajp1-2));
 	  matrix_it->second.resize_block_column(patch, 1<<(2*j));
 	}
	
  	// row 3/4 <-> interface 3,4
 	for (unsigned int patch = 3; patch <= 4; patch++) {
 	  matrix_it->second.resize_block_row   (patch, Deltajp1-2);
 	}

 	// prepare 1d matrices
	const SparseMatrix<double>& Mj1c_1d = get_Mj1c_1d(j);
	SparseMatrix<double> Mj1c_1d_interior(Deltajp1-2, 1<<j);
 	for (unsigned int row = 0; row < Deltajp1-2; row++)
 	  for (int column = 0; column < 1<<j; column++)
 	    Mj1c_1d_interior.set_entry(row, column, Mj1c_1d.get_entry(row+1, column));	
	
   	// patch wavelets decompose only into themselves
  	for (int patch = 0; patch <= 2; patch++)
  	  matrix_it->second.set_block(patch, patch,
  				      new KroneckerMatrix<double,SparseMatrix<double>,SparseMatrix<double> >
  				      (Mj1c_1d_interior, Mj1c_1d_interior));
      }
    else
      {
// 	cout << "LDomainBasis::get_Mj1c_11() cache hit" << endl;
      }
    
    return matrix_it->second;
  }
  
  template <class IBASIS>
  void
  LDomainBasis<IBASIS>::reconstruct(const InfiniteVector<double, Index>& c,
				    const int j,
				    InfiniteVector<double, Index>& d) const {
    for (typename InfiniteVector<double, Index>::const_iterator it(c.begin()), itend(c.end());
	 it != itend; ++it) {
      InfiniteVector<double, Index> help;
      reconstruct_1(it.index(), j, help);
      d.add(*it, help);
    }
  }
  
  template <class IBASIS>
  void
  LDomainBasis<IBASIS>::reconstruct_t(const InfiniteVector<double, Index>& c,
				      const int j,
				      InfiniteVector<double, Index>& d) const {
    for (typename InfiniteVector<double, Index>::const_iterator it(c.begin()), itend(c.end());
	 it != itend; ++it) {
      InfiniteVector<double, Index> help;
      reconstruct_t_1(it.index(), j, help);
      d.add(*it, help);
    }
  }

  template <class IBASIS>
  void
  LDomainBasis<IBASIS>::reconstruct_1(const Index& lambda,
				      const int j,
				      InfiniteVector<double, Index>& c) const {
    typedef typename IBASIS::Index IIndex;
    
    if (lambda.j() >= j) {
      // then we can just copy psi_lambda
      c.add_coefficient(lambda, 1.0);
    } else {
      // For the reconstruction of psi_lambda, we have to compute
      // the corresponding column of the transformation matrix Mj=(Mj0, Mj1).
      // For the left half (generators), this is comparatively easy. The wavelet case
      // will cause much more effort due to the biorthogonalization equation
      //   Mj1 = (I-Mj0*(Mj0T^T))*Mj1c
      //       = Mj1c - Mj0*(Mj0T^T)*Mj1c
      // and the specific structure of the initial stable completion
      //   Mj1c = tensor product of factors (I-Mj0*(sqrt(2)*e_1 cut(Mj0T^T) sqrt(2)*e_n))*Mj1
      
      // storage for the corresponding column of Mj1
      Vector<double> generators(Deltasize(lambda.j()+1));

      const int ecode(lambda.e()[0]+2*lambda.e()[1]);
      if (ecode == 0) {
	// generator

	// compute the corresponding column of Mj0
	const BlockMatrix<double>& Mj0 = get_Mj0(lambda.j());
	Vector<double> unitvector(Deltasize(lambda.j()));
	unsigned int id = 0;
	for (Index mu = first_generator(lambda.j()); mu != lambda; ++mu) id++;
 	unitvector[id] = 1.0;
  	Mj0.apply(unitvector, generators);	
      } else {
	if (ecode == 1) {
	  // (1,0)-wavelet
	  
	  // compute the corresponding column of Mj1c_10
	  const BlockMatrix<double>& Mj1c_10 = get_Mj1c_10(lambda.j());
	  Vector<double> unitvector(Nabla10size(lambda.j()));
	  unsigned int id = 0;
	  typename Index::type_type e;
	  e[0] = 1; // e = (1,0)
	  for (Index mu = first_wavelet(lambda.j(), e); mu != lambda; ++mu) id++;
 	  unitvector[id] = 1.0;
 	  Mj1c_10.apply(unitvector, generators);
	} else {
 	  if (ecode == 2) {
 	    // (0,1)-wavelet
	    
	    // compute the corresponding column of Mj1c_01
 	    const BlockMatrix<double>& Mj1c_01 = get_Mj1c_01(lambda.j());
 	    Vector<double> unitvector(Nabla01size(lambda.j()));
 	    unsigned int id = 0;
 	    for (Index mu = first_wavelet(lambda.j()); mu != lambda; ++mu) id++;
 	    unitvector[id] = 1.0;
 	    Mj1c_01.apply(unitvector, generators);
	  } else {
	    // (1,1)-wavelet

	    // compute the corresponding column of Mj1c_11
 	    const BlockMatrix<double>& Mj1c_11 = get_Mj1c_11(lambda.j());
	    
 	    Vector<double> unitvector(Nabla11size(lambda.j()));
 	    unsigned int id = 0;
	    typename Index::type_type e;
	    e[0] = e[1] = 1; // e = (1,1)
	    for (Index mu = first_wavelet(lambda.j(), e); mu != lambda; ++mu) id++;
 	    unitvector[id] = 1.0;
 	    Mj1c_11.apply(unitvector, generators);
	  }
	}
      }

      // Now that the corresponding column of Mj1 has been computed, we collect
      // all relevant generators from the scale |lambda|+1
      // (this is the identity part of the biorthogonalization equation)
      unsigned int id = 0;
      for (Index mu = first_generator(lambda.j()+1); id < generators.size(); id++, ++mu) {
	if (fabs(generators[id]) > 1e-12) {
	  InfiniteVector<double, Index> d;
	  reconstruct_1(mu, j, d);
	  c.add(generators[id], d);
	}
      }

      if (ecode > 0) {
 	// second part of the biorthogonalization equation,
 	// compute the corresponding column of Mj0*Mj0T^T*Mj1c
	Vector<double> help(Deltasize(lambda.j()));
 	get_Mj0T(lambda.j()).apply_transposed(generators, help);
 	get_Mj0 (lambda.j()).apply(help, generators);

 	// collect all relevant generators from the scale |lambda+1|
 	unsigned int id = 0;
 	for (Index mu = first_generator(lambda.j()+1); id < generators.size(); id++, ++mu) {
  	  if (fabs(generators[id]) > 1e-12) {
 	    InfiniteVector<double, Index> d;
 	    reconstruct_1(mu, j, d);
 	    c.add(-generators[id], d);
  	  }
 	}
      }
    }
  }
  
  template <class IBASIS>
  void
  LDomainBasis<IBASIS>::reconstruct_t_1(const Index& lambda,
					const int j,
					InfiniteVector<double, Index>& c) const {
    // not implemented (yet)
  }

}
