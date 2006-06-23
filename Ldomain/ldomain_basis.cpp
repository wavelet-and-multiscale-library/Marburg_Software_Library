// implementation for ldomain_basis.h

#include <cmath>
#include <time.h>
#include <iostream>

using std::cout;
using std::endl;

namespace WaveletTL
{
  template <class IBASIS>
  LDomainBasis<IBASIS>::LDomainBasis()
    : basis1d_(false, false)
  {
#if _WAVELETTL_LDOMAINBASIS_VERBOSITY >= 1
    supp_hits = 0;
    supp_misses = 0;
#endif
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
    return (3*basis1d().Deltasize(j)-5)*(1<<j);
  }
  
  template <class IBASIS>
  const int
  LDomainBasis<IBASIS>::Nabla10size(const int j) const {
    return (3*basis1d().Deltasize(j)-5)*(1<<j);
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
    typename Index::translation_type k(basis1d().DeltaLmin()+1, basis1d().DeltaLmin()+1);
    
    return Index(j, e, 0, k, this);
  }

  template <class IBASIS>
  typename LDomainBasis<IBASIS>::Index
  LDomainBasis<IBASIS>::last_generator(const int j) const
  {
    assert(j >= j0());

    typename Index::type_type e;

    // setup highest translation index for e=(0,0), p=4
    typename Index::translation_type k(0, basis1d().DeltaRmax(j)-1);
    
    return Index(j, e, 4, k, this);
  }

  template <class IBASIS>
  typename LDomainBasis<IBASIS>::Index
  LDomainBasis<IBASIS>::first_wavelet(const int j) const
  {
    assert(j >= j0());

    typename Index::type_type e(0, 1);

    // setup lowest translation index for e=(0,1), p=0
    typename Index::translation_type k(basis1d().DeltaLmin()+1, basis1d().Nablamin());
    
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
    
    typename Index::type_type e(1, 1);

    // setup highest translation index for e=(1,1), p=2
    typename Index::translation_type k(basis1d().Nablamax(j), basis1d().Nablamax(j));
    
    return Index(j, e, 2, k, this);
  }

  template <class IBASIS>
  void
  LDomainBasis<IBASIS>::support(const Index& lambda, Support& supp) const
  {
    // check whether the supp(psi_lambda) already exists in the cache
    typename SupportCache::iterator supp_lb(supp_cache.lower_bound(lambda));
    typename SupportCache::iterator supp_it(supp_lb);
    if (supp_lb == supp_cache.end() ||
	supp_cache.key_comp()(lambda, supp_lb->first))
      {
	// compute supp(psi_lambda) and insert it into the cache
	typedef typename SupportCache::value_type value_type;

	const int ecode = lambda.e()[0]+2*lambda.e()[1];
	const int lambdaj = lambda.j();
	
	if (ecode == 0) {
	  // psi_lambda is a generator. Here we know by construction of the
	  // composite basis that per patch, psi_lambda looks like a single
	  // tensor product of 1D generators (possibly weighted by a factor).
	  
	  supp.j = lambdaj;
	  
	  switch (lambda.p()) {
	  case 0:
	    // psi_lambda completely lives on patch 0
	    
	    basis1d().support(typename IBASIS::Index(lambdaj,
						     0,
						     lambda.k()[0],
						     &basis1d()),
			      supp.xmin[0],
			      supp.xmax[0]);
	    
	    basis1d().support(typename IBASIS::Index(lambdaj,
						     0,
						     lambda.k()[1],
						     &basis1d()),
			      supp.ymin[0],
			      supp.ymax[0]);
	    
	    supp.xmin[1] = supp.xmin[2] = -1;
	    
	    break;
	  case 1:
	    // psi_lambda completely lives on patch 1
	    
	    basis1d().support(typename IBASIS::Index(lambdaj,
						     0,
						     lambda.k()[0],
						     &basis1d()),
			      supp.xmin[1],
			      supp.xmax[1]);
	    
	    basis1d().support(typename IBASIS::Index(lambdaj,
						     0,
						     lambda.k()[1],
						     &basis1d()),
			      supp.ymin[1],
			      supp.ymax[1]);
	    
	    supp.xmin[0] = supp.xmin[2] = -1;
	    
	    break;
	  case 2:
	    // psi_lambda completely lives on patch 2
	    
	    basis1d().support(typename IBASIS::Index(lambdaj,
						     0,
						     lambda.k()[0],
						     &basis1d()),
			      supp.xmin[2],
			      supp.xmax[2]);
	    
	    basis1d().support(typename IBASIS::Index(lambdaj,
						     0,
						     lambda.k()[1],
						     &basis1d()),
			      supp.ymin[2],
			      supp.ymax[2]);
	    
	    supp.xmin[0] = supp.xmin[1] = -1;
	    
	    break;
	  case 3:
	    // psi_lambda lives on patches 0 and 1
	    
	    basis1d().support(typename IBASIS::Index(lambdaj,
						     0,
						     lambda.k()[0],
						     &basis1d()),
			      supp.xmin[0],
			      supp.xmax[0]);
	    
	    basis1d().support(typename IBASIS::Index(lambdaj,
						     0,
						     basis1d().DeltaLmin(),
						     &basis1d()),
			      supp.ymin[0],
			      supp.ymax[0]);
	    
	    basis1d().support(typename IBASIS::Index(lambdaj,
						     0,
						     lambda.k()[0],
						     &basis1d()),
			      supp.xmin[1],
			      supp.xmax[1]);
	    
	    basis1d().support(typename IBASIS::Index(lambdaj,
						     0,
						     basis1d().DeltaRmax(lambdaj),
						     &basis1d()),
			      supp.ymin[1],
			      supp.ymax[1]);
	    
	    supp.xmin[2] = -1;
	    
	    break;
	  case 4:
	    // psi_lambda lives on patches 1 and 2
	    
	    basis1d().support(typename IBASIS::Index(lambdaj,
						     0,
						     basis1d().DeltaRmax(lambdaj),
						     &basis1d()),
			      supp.xmin[1],
			      supp.xmax[1]);
	    
	    basis1d().support(typename IBASIS::Index(lambdaj,
						     0,
						     lambda.k()[1],
						     &basis1d()),
			      supp.ymin[1],
			      supp.ymax[1]);
	    
	    basis1d().support(typename IBASIS::Index(lambdaj,
						     0,
						     basis1d().DeltaLmin(),
						     &basis1d()),
			      supp.xmin[2],
			      supp.xmax[2]);
	    
	    basis1d().support(typename IBASIS::Index(lambdaj,
						     0,
						     lambda.k()[1],
						     &basis1d()),
			      supp.ymin[2],
			      supp.ymax[2]);
	    
	    supp.xmin[0] = -1;
	    
	    break;
	  }
	} else {
	  // wavelet
	  
	  supp.j = lambdaj+1;
	  
	  // compute the expansion coefficients of psi_lambda w.r.t. the
	  // generators of the next higher scale, then aggregating all the supports
	  // (of course, this is a brute force solution...)
	  InfiniteVector<double, Index> gcoeffs;
	  reconstruct_1(lambda, lambdaj+1, gcoeffs);
	  
	  Support tempsupp;
	  
	  // initialize the support with an "empty" set
	  for (int p = 0; p <= 2; p++) {
	    supp.xmin[p] = -1;
	  }
	  
	  for (typename InfiniteVector<double,Index>::const_iterator it(gcoeffs.begin()),
		 itend(gcoeffs.end()); it != itend; ++it)
	    {
	      // compute supp(psi_mu)
	      support(it.index(), tempsupp);
	      
	      // for each patch p, update the corresponding support estimate
	      for (int p = 0; p <= 2; p++) {
		if (tempsupp.xmin[p] != -1) {
		  // a nontrivial new support share, we have to do something
		  if (supp.xmin[p] == -1) {
		    // previous support estimate was "empty", we have to insert a nontrivial new one
		    supp.xmin[p] = tempsupp.xmin[p];
		    supp.xmax[p] = tempsupp.xmax[p];
		    supp.ymin[p] = tempsupp.ymin[p];
		    supp.ymax[p] = tempsupp.ymax[p];
		  } else {
		    // previous support estimate was nontrivial, we have to compute a new one
		    supp.xmin[p] = std::min(supp.xmin[p], tempsupp.xmin[p]);
		    supp.xmax[p] = std::max(supp.xmax[p], tempsupp.xmax[p]);
		    supp.ymin[p] = std::min(supp.ymin[p], tempsupp.ymin[p]);
		    supp.ymax[p] = std::max(supp.ymax[p], tempsupp.ymax[p]);
		  }
		}
	      }
	    }
	}
	
	supp_it = supp_cache.insert(supp_lb, value_type(lambda, supp));

#if _WAVELETTL_LDOMAINBASIS_VERBOSITY >= 1
	supp_misses++;
	if ((supp_hits+supp_misses)%1000000 == 0)
	  {
	    cout << "[LdomainBasis support cache (hits/misses/total/hit ratio): ("
		 << supp_hits << "/"
		 << supp_misses << "/"
		 << supp_hits+supp_misses << "/"
		 << (double) supp_hits/(supp_hits+supp_misses)
		 << ")"
		 << "]"
		 << endl;
	  }
#endif
      }
    else
      {
	// cache hit, copy the precomputed support
  	supp.j = supp_it->second.j;
  	for (unsigned int i = 0; i < 3; i++) {
  	  supp.xmin[i] = supp_it->second.xmin[i];
  	  supp.xmax[i] = supp_it->second.xmax[i];
  	  supp.ymin[i] = supp_it->second.ymin[i];
  	  supp.ymax[i] = supp_it->second.ymax[i];
  	}
	
#if _WAVELETTL_LDOMAINBASIS_VERBOSITY >= 1
	supp_hits++;
	if ((supp_hits+supp_misses)%1000000 == 0)
	  {
	    cout << "[LdomainBasis support cache (hits/misses/total/hit ratio): ("
		 << supp_hits << "/"
		 << supp_misses << "/"
		 << supp_hits+supp_misses << "/"
		 << (double) supp_hits/(supp_hits+supp_misses)
		 << ")"
		 << "]"
		 << endl;
	  }
#endif
      }  
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

#if _WAVELETTL_LDOMAINBASIS_VERBOSITY >= 2
    clock_t tstart, tend, tmiddle1, tmiddle2;
#endif

    const int lambdaj = lambda.j();

    if (lambdaj >= j) {
      // then we can just copy psi_lambda
      c.add_coefficient(lambda, 1.0);
    } else {

#if _WAVELETTL_LDOMAINBASIS_VERBOSITY >= 2
      cout << "LDomainBasis::reconstruct_1() nontrivially called with lambda=" << lambda << endl;
#endif

      // For the reconstruction of psi_lambda, we have to compute
      // the corresponding column of the transformation matrix Mj=(Mj0, Mj1).
      // For the left half (generators), this is comparatively easy. The wavelet case
      // will cause much more effort due to the biorthogonalization equation
      //   Mj1 = (I-Mj0*(Mj0T^T))*Mj1c
      //       = Mj1c - Mj0*(Mj0T^T)*Mj1c
      // and the specific structure of the initial stable completion
      //   Mj1c = tensor product of factors (I-Mj0*(sqrt(2)*e_1 cut(Mj0T^T) sqrt(2)*e_n))*Mj1
      
      // storage for the corresponding column of Mj1
      Vector<double> generators(Deltasize(lambdaj+1));

#if _WAVELETTL_LDOMAINBASIS_VERBOSITY >= 2
      tstart = clock();
#endif

      const int ecode(lambda.e()[0]+2*lambda.e()[1]);
      
      if (ecode == 0) {
	// generator

	// compute the corresponding column of Mj0
	const BlockMatrix<double>& Mj0 = get_Mj0(lambdaj);
	Vector<double> unitvector(Deltasize(lambdaj));
#if 0
	// check the number
	int id = lambda.number();
	int idcheck = 0;
	for (Index mu = first_generator(lambdaj); mu != lambda; ++mu) idcheck++;
	if (id != idcheck) {
	  cout << "in LDomainBasis::reconstruct_1(), id != idcheck!" << endl;
	  abort();
	}
#endif
 	unitvector[lambda.number()] = 1.0;
  	Mj0.apply(unitvector, generators);	
      } else {
	if (ecode == 1) {
	  // (1,0)-wavelet
	  
	  // compute the corresponding column of Mj1c_10
	  const BlockMatrix<double>& Mj1c_10 = get_Mj1c_10(lambdaj);
	  Vector<double> unitvector(Nabla10size(lambdaj));
	  const typename Index::type_type e(1, 0);
#if 0
	  // check the number
	  int id = lambda.number()-first_wavelet(lambdaj,e).number();
	  int idcheck = 0;
	  for (Index mu = first_wavelet(lambdaj, e); mu != lambda; ++mu) idcheck++;
	  if (id != idcheck) {
	    cout << "in LDomainBasis::reconstruct_1(), id != idcheck!" << endl;
	    abort();
	  }
#endif
 	  unitvector[lambda.number()-first_wavelet(lambdaj,e).number()] = 1.0;
 	  Mj1c_10.apply(unitvector, generators);	  
	} else {
 	  if (ecode == 2) {
 	    // (0,1)-wavelet
	    
	    // compute the corresponding column of Mj1c_01
 	    const BlockMatrix<double>& Mj1c_01 = get_Mj1c_01(lambdaj);
 	    Vector<double> unitvector(Nabla01size(lambdaj));
#if 0
	    // check the number
	    int id = lambda.number()-first_wavelet(lambdaj).number();
	    int idcheck = 0;
	    for (Index mu = first_wavelet(lambdaj); mu != lambda; ++mu) idcheck++;
	    if (id != idcheck) {
	      cout << "in LDomainBasis::reconstruct_1(), id != idcheck!" << endl;
	      abort();
	    }
#endif
	    unitvector[lambda.number()-first_wavelet(lambdaj).number()] = 1.0;
 	    Mj1c_01.apply(unitvector, generators);
	  } else {
	    // (1,1)-wavelet

	    // compute the corresponding column of Mj1c_11
 	    const BlockMatrix<double>& Mj1c_11 = get_Mj1c_11(lambdaj);
	    
 	    Vector<double> unitvector(Nabla11size(lambdaj));
	    
	    const typename Index::type_type e(1, 1);
#if 0
	  // check the number
	    int id = lambda.number()-first_wavelet(lambdaj,e).number();
	    int idcheck = 0;
	    for (Index mu = first_wavelet(lambdaj, e); mu != lambda; ++mu) idcheck++;
	    if (id != idcheck) {
	      cout << "in LDomainBasis::reconstruct_1(), id != idcheck!" << endl;
	      abort();
	    }
#endif
 	    unitvector[lambda.number()-first_wavelet(lambdaj,e).number()] = 1.0;
 	    Mj1c_11.apply(unitvector, generators);
	  }
	}
      }

#if _WAVELETTL_LDOMAINBASIS_VERBOSITY >= 2
      tend = clock();
      cout << "* in reconstruct_1(), time needed for the column of Mj1: "
	   << (double)(tend-tstart)/CLOCKS_PER_SEC
	   << "s" << endl;

      tstart = clock();
#endif

      // Now that the corresponding column of Mj1 has been computed, we collect
      // all relevant generators from the scale |lambda|+1
      // (this is the identity part of the biorthogonalization equation)
      unsigned int id = 0;
      if (lambdaj+1 >= j) {
	for (Index mu(first_generator(lambdaj+1)); id < generators.size(); id++, ++mu) {
	  c.add_coefficient(mu, generators[id]);
	}
      } else {
	for (Index mu(first_generator(lambdaj+1)); id < generators.size(); id++, ++mu) {
	  InfiniteVector<double, Index> d;
	  reconstruct_1(mu, j, d);
	  c.add(generators[id], d);
	}
      }
      
#if _WAVELETTL_LDOMAINBASIS_VERBOSITY >= 2
      tend = clock();
      cout << "* in reconstruct_1(), time needed for the I part       : "
	   << (double)(tend-tstart)/CLOCKS_PER_SEC
	   << "s" << endl;

      tstart = clock();
#endif

      if (ecode > 0) {
 	// second part of the biorthogonalization equation,
 	// compute the corresponding column of -Mj0*Mj0T^T*Mj1c
	Vector<double> help(Deltasize(lambdaj));
 	get_Mj0T(lambdaj).apply_transposed(generators, help);
#if _WAVELETTL_LDOMAINBASIS_VERBOSITY >= 2
	tmiddle1 = clock();
#endif
	
 	get_Mj0 (lambdaj).apply(help, generators);
#if _WAVELETTL_LDOMAINBASIS_VERBOSITY >= 2
	tmiddle2 = clock();
#endif

 	// collect all relevant generators from the scale |lambda+1|
 	unsigned int id = 0;
	if (lambdaj+1 >= j) {
	  for (Index mu(first_generator(lambdaj+1)); id < generators.size(); id++, ++mu) {
	    c.add_coefficient(mu, -generators[id]);
	  }
	} else {
	  for (Index mu(first_generator(lambdaj+1)); id < generators.size(); id++, ++mu) {
	    InfiniteVector<double, Index> d;
	    reconstruct_1(mu, j, d);
	    c.add(-generators[id], d);
	  }
 	}
      }

#if _WAVELETTL_LDOMAINBASIS_VERBOSITY >= 2
      tend = clock();
      cout << "* in reconstruct_1(), time needed to apply Mj0T^T      : "
	   << (double)(tmiddle1-tstart)/CLOCKS_PER_SEC
	   << "s" << endl;
      cout << "* in reconstruct_1(), time needed to apply Mj0         : "
	   << (double)(tmiddle2-tmiddle1)/CLOCKS_PER_SEC
	   << "s" << endl;
      cout << "* in reconstruct_1(), time needed for the bio. part    : "
	   << (double)(tend-tmiddle2)/CLOCKS_PER_SEC
	   << "s" << endl;
#endif
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
