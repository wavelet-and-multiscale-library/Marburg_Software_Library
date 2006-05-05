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
  LDomainBasis<IBASIS>::Deltasize(const int j) const
  {
    
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
	cout << "LDomainBasis::get_Mj0() cache miss" << endl;

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
 	SparseMatrix<double> Mj0_1d_interior(Mj0_1d.row_dimension()-2, Mj0_1d.column_dimension()-2);
 	for (unsigned int row = 0; row < Mj0_1d_interior.row_dimension(); row++)
 	  for (unsigned int column = 0; column < Mj0_1d_interior.column_dimension(); column++)
 	    Mj0_1d_interior.set_entry(row, column, Mj0_1d.get_entry(row+1, column+1));
	SparseMatrix<double> Mj0_1d_left(Mj0_1d.row_dimension()-2, 1);
 	for (unsigned int row = 0; row < Mj0_1d_left.row_dimension(); row++)
	  Mj0_1d_left.set_entry(row, 0, Mj0_1d.get_entry(row+1, 0));
	SparseMatrix<double> Mj0_1d_right(Mj0_1d.row_dimension()-2, 1);
 	for (unsigned int row = 0; row < Mj0_1d_right.row_dimension(); row++)
	  Mj0_1d_right.set_entry(row, 0, Mj0_1d.get_entry(row+1, Mj0_1d.column_dimension()-1));
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
	cout << "LDomainBasis::get_Mj0() cache hit" << endl;
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
	cout << "LDomainBasis::get_Mj0T() cache miss" << endl;

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
 	SparseMatrix<double> Mj0T_1d_interior(Mj0T_1d.row_dimension()-2, Mj0T_1d.column_dimension()-2);
 	for (unsigned int row = 0; row < Mj0T_1d_interior.row_dimension(); row++)
 	  for (unsigned int column = 0; column < Mj0T_1d_interior.column_dimension(); column++)
 	    Mj0T_1d_interior.set_entry(row, column, Mj0T_1d.get_entry(row+1, column+1));
	SparseMatrix<double> Mj0T_1d_left(Mj0T_1d.row_dimension()-2, 1);
 	for (unsigned int row = 0; row < Mj0T_1d_left.row_dimension(); row++)
	  Mj0T_1d_left.set_entry(row, 0, Mj0T_1d.get_entry(row+1, 0));
	SparseMatrix<double> Mj0T_1d_right(Mj0T_1d.row_dimension()-2, 1);
 	for (unsigned int row = 0; row < Mj0T_1d_right.row_dimension(); row++)
	  Mj0T_1d_right.set_entry(row, 0, Mj0T_1d.get_entry(row+1, Mj0T_1d.column_dimension()-1));
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
	cout << "LDomainBasis::get_Mj0T() cache hit" << endl;
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
	cout << "LDomainBasis::get_Mj1c_1d() cache miss" << endl;

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
	cout << "LDomainBasis::get_Mj1c_1d() cache hit" << endl;
      }

    return matrix_it->second;
  }

  template <class IBASIS>
  const BlockMatrix<double>&
  LDomainBasis<IBASIS>::get_Mj1c_01 (const int j) const {
    // check whether the j-th matrix already exists in the cache
    typename MatrixCache::iterator matrix_lb(Mj1c_cache.lower_bound(j));
    typename MatrixCache::iterator matrix_it(matrix_lb);
    if (matrix_lb == Mj1c_cache.end() ||
	Mj1c_cache.key_comp()(j, matrix_lb->first))
      {
	cout << "LDomainBasis::get_Mj1c() cache miss" << endl;

// 	// compute Mj0T and insert it into the cache
// 	BlockMatrix<double> Mj0T(5, 5);
// 	typedef typename MatrixCache::value_type value_type;
//  	matrix_it = Mj0T_cache.insert(matrix_lb, value_type(j, Mj0T));
	
// 	const unsigned int Deltaj   = basis1d().Deltasize(j);
// 	const unsigned int Deltajp1 = basis1d().Deltasize(j+1);

// 	// row/column 0,1,2 <-> patches 0,1,2
// 	for (unsigned int patch = 0; patch <= 2; patch++) {
// 	  matrix_it->second.resize_block_row   (patch, (Deltajp1-2)*(Deltajp1-2));
// 	  matrix_it->second.resize_block_column(patch, (Deltaj-2)*(Deltaj-2));
// 	}
	
//  	// row/column 3,4 <-> interface 3,4
// 	for (unsigned int patch = 3; patch <= 4; patch++) {
// 	  matrix_it->second.resize_block_row   (patch, Deltajp1-2);
// 	  matrix_it->second.resize_block_column(patch, Deltaj-2);
// 	}

//  	// prepare 1d matrices
//  	SparseMatrix<double> Mj0T_1d; basis1d().assemble_Mj0T(j, Mj0T_1d);
//  	SparseMatrix<double> Mj0T_1d_interior(Mj0T_1d.row_dimension()-2, Mj0T_1d.column_dimension()-2);
//  	for (unsigned int row = 0; row < Mj0T_1d_interior.row_dimension(); row++)
//  	  for (unsigned int column = 0; column < Mj0T_1d_interior.column_dimension(); column++)
//  	    Mj0T_1d_interior.set_entry(row, column, Mj0T_1d.get_entry(row+1, column+1));
// 	SparseMatrix<double> Mj0T_1d_left(Mj0T_1d.row_dimension()-2, 1);
//  	for (unsigned int row = 0; row < Mj0T_1d_left.row_dimension(); row++)
// 	  Mj0T_1d_left.set_entry(row, 0, Mj0T_1d.get_entry(row+1, 0));
// 	SparseMatrix<double> Mj0T_1d_right(Mj0T_1d.row_dimension()-2, 1);
//  	for (unsigned int row = 0; row < Mj0T_1d_right.row_dimension(); row++)
// 	  Mj0T_1d_right.set_entry(row, 0, Mj0T_1d.get_entry(row+1, Mj0T_1d.column_dimension()-1));
// 	SparseMatrix<double> Mj0T_1d_left_top(1, 1); Mj0T_1d_left_top.set_entry(0, 0, Mj0T_1d.get_entry(0, 0));
	
//  	// patch generators decompose only into themselves
// 	for (int patch = 0; patch <= 2; patch++)
// 	  matrix_it->second.set_block(patch, patch,
// 				      new KroneckerMatrix<double,SparseMatrix<double>,SparseMatrix<double> >
// 				      (Mj0T_1d_interior, Mj0T_1d_interior));
	
// 	// interface generators decompose into themselves and patch generators from the neighboring patches
//  	matrix_it->second.set_block(0, 3,
//  				    new KroneckerMatrix<double,SparseMatrix<double>,SparseMatrix<double> >
//  				    (Mj0T_1d_interior, Mj0T_1d_left, M_SQRT1_2));
//  	matrix_it->second.set_block(1, 3,
//  				    new KroneckerMatrix<double,SparseMatrix<double>,SparseMatrix<double> >
//  				    (Mj0T_1d_interior, Mj0T_1d_right, M_SQRT1_2));
//  	matrix_it->second.set_block(3, 3,
//  				    new KroneckerMatrix<double,SparseMatrix<double>,SparseMatrix<double> >
//  				    (Mj0T_1d_left_top, Mj0T_1d_interior));

//  	matrix_it->second.set_block(1, 4,
//  				    new KroneckerMatrix<double,SparseMatrix<double>,SparseMatrix<double> >
//  				    (Mj0T_1d_right, Mj0T_1d_interior, M_SQRT1_2));
//  	matrix_it->second.set_block(2, 4,
//  				    new KroneckerMatrix<double,SparseMatrix<double>,SparseMatrix<double> >
//  				    (Mj0T_1d_left, Mj0T_1d_interior, M_SQRT1_2));
//  	matrix_it->second.set_block(4, 4,
//  				    new KroneckerMatrix<double,SparseMatrix<double>,SparseMatrix<double> >
//  				    (Mj0T_1d_left_top, Mj0T_1d_interior));
      }
    else
      {
	cout << "LDomainBasis::get_Mj1c() cache hit" << endl;
      }
    
    return matrix_it->second;
  }

//   template <class IBASIS>
//   const BlockMatrix<double>&
//   LDomainBasis<IBASIS>::get_Mj1 (const int j) const {
//   }

//   template <class IBASIS>
//   const BlockMatrix<double>&
//   LDomainBasis<IBASIS>::get_Mj1T  (const int j) const {
//   }

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
      
      const int ecode(lambda.e()[0]+2*lambda.e()[1]);

      if (ecode == 0) {
	// generator

	// compute the corresponding column of Mj0
	const BlockMatrix<double>& Mj0 = get_Mj0(lambda.j());
	const Vector<double> unitvector, generators;

	// collect all relevant generators from the scale |lambda|+1
	


//  	    for (int l(2*lambda.k()+abegin); l <= 2*lambda.k()+aend; l++)
//  	      {
//  		InfiniteVector<double, Index> d;
//  		reconstruct_1(Index(lambda.j()+1, 0, l), j, d);
//  		c.add(M_SQRT1_2 * a().get_coefficient(MultiIndex<int, 1>(l-2*lambda.k())), d);
//  	      }

      } else {

      
//       InfiniteVector<double,IIndex> gcoeffs0, gcoeffs1, gcoeffs2;
//       InfiniteVector<double,Index> psic; // psi_lambda^check

      
//       // compute c from psic
//       c.swap(psic); // TODO: implement biorthogonalization

      }









      
//       switch(ecode) {
//       case 0:
//  	// generator
	
// 	// collect all relevant generators from the scale |lambda|+1
	

//  	break;
//       case 1:
//  	// (1,0)-wavelet
//  	break;
//       case 2:
//  	// (0,1)-wavelet
	
//  	// First treat the identity part of the biorthogonalization equation, i.e.,
//  	// compute the generator coefficients of the initial stable completion
// 	//   Mj1c = (I-Mj0*<Theta_{j+1},Lambda_j^tilde>^T)*Mj1
// 	// So, the wavelets which use one of the nonvaninshing boundary generators
// 	// will be modified to be zero at the endpoints.
	
//  	switch(lambda.p()) {
// 	case 0:
// 	  // psic_lambda is a tensor product of a generator and a wavelet on patch 0
// 	  basis00().reconstruct_1(IIndex(lambda.j(), 0, lambda.k()[0], &basis00()),
// 				  lambda.j()+1, gcoeffs0);
// 	  basis10().reconstruct_1(IIndex(lambda.j(), 1, lambda.k()[1], &basis10()),
// 				  lambda.j()+1, gcoeffs1);

// 	  // add the tensor product generators needed to construct psic_lambda
// 	  for (typename InfiniteVector<double,IIndex>::const_iterator it1(gcoeffs0.begin()),
// 		 it1end(gcoeffs0.end());
// 	       it1 != it1end; ++it1)
// 	    for (typename InfiniteVector<double,IIndex>::const_iterator it2(gcoeffs1.begin()),
// 		   it2end(gcoeffs1.end());
// 		 it2 != it2end; ++it2) {
// 	      if (it2.index().k() == basis10().DeltaLmin()) {
// 		// interface generator
// 		psic.add_coefficient(Index(lambda.j()+1,
// 					   typename Index::type_type(),
// 					3,
// 					   typename Index::translation_type(it1.index().k(), 0),
// 					   this),
// 				     *it1 * *it2); // no factor!
// 	      } else {
// 		// patch generator
// 		psic.add_coefficient(Index(lambda.j()+1,
// 					   typename Index::type_type(),
// 					   0,
// 					   typename Index::translation_type(it1.index().k(), it2.index().k()),
// 					   this),
// 				     *it1 * *it2);
// 	      }
// 	    }
// 	  break;
// 	case 1:
// 	  // psic_lambda is a tensor product of a generator and a wavelet on patch 1
// 	  basis01().reconstruct_1(IIndex(lambda.j(), 0, lambda.k()[0], &basis01()),
// 				  lambda.j()+1, gcoeffs0);
// 	  basis01().reconstruct_1(IIndex(lambda.j(), 1, lambda.k()[1], &basis01()),
// 				  lambda.j()+1, gcoeffs1);

// 	  // add the tensor product generators needed to construct psic_lambda
// 	  for (typename InfiniteVector<double,IIndex>::const_iterator it1(gcoeffs0.begin()),
// 		 it1end(gcoeffs0.end());
// 	       it1 != it1end; ++it1)
// 	    for (typename InfiniteVector<double,IIndex>::const_iterator it2(gcoeffs1.begin()),
// 		   it2end(gcoeffs1.end());
// 		 it2 != it2end; ++it2) {
// 	      if (it2.index().k() == basis01().DeltaRmax(lambda.j()+1)) {
// 		// interface generator
// 		psic.add_coefficient(Index(lambda.j()+1,
// 					   typename Index::type_type(),
// 					   3,
// 					   typename Index::translation_type(it1.index().k(), 0),
// 					   this),
// 				     *it1 * *it2); // no factor!
// 	      } else {
// 		// patch generator
// 		psic.add_coefficient(Index(lambda.j()+1,
// 					   typename Index::type_type(),
// 					   1,
// 					   typename Index::translation_type(it1.index().k(), it2.index().k()),
// 					   this),
// 				     *it1 * *it2);
// 	      }
// 	    }
// 	  break;
// 	case 2:
// 	  // psic_lambda is a tensor product of a generator and a wavelet on patch 2
// 	  basis10().reconstruct_1(IIndex(lambda.j(), 0, lambda.k()[0], &basis10()),
// 				  lambda.j()+1, gcoeffs0);
// 	  basis00().reconstruct_1(IIndex(lambda.j(), 1, lambda.k()[1], &basis00()),
// 				  lambda.j()+1, gcoeffs1);

// 	  // add the tensor product generators needed to construct psic_lambda
// 	  for (typename InfiniteVector<double,IIndex>::const_iterator it1(gcoeffs0.begin()),
// 		 it1end(gcoeffs0.end());
// 	       it1 != it1end; ++it1)
// 	    for (typename InfiniteVector<double,IIndex>::const_iterator it2(gcoeffs1.begin()),
// 		   it2end(gcoeffs1.end());
// 		 it2 != it2end; ++it2) {
// 	      // (always!) a patch generator
// 	      psic.add_coefficient(Index(lambda.j()+1,
// 					 typename Index::type_type(),
// 					 2,
// 					 typename Index::translation_type(it1.index().k(), it2.index().k()),
// 					 this),
// 				   *it1 * *it2);
// 	    }
// 	  break;
// 	case 4:
// 	  // psic_lambda decomposes into two tensor products of a generator and a wavelet,
// 	  // on patches 1 and 2
// 	  basis01().reconstruct_1(IIndex(lambda.j(), 0, basis01().DeltaRmax(lambda.j()), &basis01()),
// 				  lambda.j()+1, gcoeffs0);
// 	  basis10().reconstruct_1(IIndex(lambda.j(), 0, basis10().DeltaLmin(), &basis10()),
// 				  lambda.j()+1, gcoeffs1);
// 	  basis01().reconstruct_1(IIndex(lambda.j(), 1, lambda.k()[1], &basis01()),
// 				  lambda.j()+1, gcoeffs2);
	  
// 	  // add the tensor product generators needed to construct psic_lambda
// 	  for (typename InfiniteVector<double,IIndex>::const_iterator it1(gcoeffs0.begin()),
// 		 it1end(gcoeffs0.end());
// 	       it1 != it1end; ++it1)
// 	    for (typename InfiniteVector<double,IIndex>::const_iterator it2(gcoeffs2.begin()),
// 		   it2end(gcoeffs2.end());
// 		 it2 != it2end; ++it2) {
// 	      if (it1.index().k() == basis01().DeltaRmax(lambda.j()+1)) {
// 		// interface generator
// 		psic.add_coefficient(Index(lambda.j()+1,
// 					   typename Index::type_type(),
// 					   4,
// 					   typename Index::translation_type(0, it2.index().k()),
// 					   this),
// 				     *it1 * *it2); // no factor!
// 	      } else {
// 		// patch generator
// 		psic.add_coefficient(Index(lambda.j()+1,
// 					   typename Index::type_type(),
// 					   1,
// 					   typename Index::translation_type(it1.index().k(), it2.index().k()),
// 					   this),
// 				     *it1 * *it2);
// 	      }
// 	    }
//  	  for (typename InfiniteVector<double,IIndex>::const_iterator it1(gcoeffs1.begin()),
//  		 it1end(gcoeffs1.end());
//  	       it1 != it1end; ++it1)
//  	    for (typename InfiniteVector<double,IIndex>::const_iterator it2(gcoeffs2.begin()),
//  		   it2end(gcoeffs2.end());
//  		 it2 != it2end; ++it2) {
//  	      if (it1.index().k() > basis10().DeltaLmin()) {
//  		// patch generator (interface generators already processed above)
//  		psic.add_coefficient(Index(lambda.j()+1,
// 					   typename Index::type_type(),
// 					   2,
// 					   typename Index::translation_type(it1.index().k(), it2.index().k()),
// 					   this),
// 				     *it1 * *it2);
//  	      }
//  	    }
// 	  break;
// 	} // end switch(lambda.p())
	
//  	break;
//       case 3:
//  	// (1,1)-wavelet
//  	break;
//       }


    }
  }
  
  template <class IBASIS>
  void
  LDomainBasis<IBASIS>::reconstruct_t_1(const Index& lambda,
					const int j,
					InfiniteVector<double, Index>& c) const {
  }

}
