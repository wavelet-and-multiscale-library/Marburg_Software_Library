// implementation for cube_basis.h

#include <map>

namespace WaveletTL
{
  template <class IBASIS, unsigned int DIM>
  CubeBasis<IBASIS,DIM>::CubeBasis()
    : bases_()
  {
    // we only need one instance of IBASIS, without b.c.
    IBASIS* b = new IBASIS();
    bases_infact.push_back(b);
    for (unsigned int i = 0; i < DIM; i++)
      bases_[i] = b;
    j0_ = bases_[0]->j0();
  }

  template <class IBASIS, unsigned int DIM>
  CubeBasis<IBASIS,DIM>::CubeBasis(const FixedArray1D<int,2*DIM>& s,
				   const FixedArray1D<int,2*DIM>& sT) {
    for (unsigned int i = 0; i < DIM; i++) {
      // check whether the corresponding 1d basis already exists
      IBASIS* b = 0;
      for (typename list<IBASIS*>::const_iterator it(bases_infact.begin());
	   it != bases_infact.end(); ++it) {
	if ((*it)->get_s0() == s[2*i]
	    && (*it)->get_s1() == s[2*i+1]
	    && (*it)->get_sT0() == sT[2*i]
	    && (*it)->get_sT1() == sT[2*i+1]) {
	  b = *it;
	  break;
	}
      }
      if (b == 0) {
	b = new IBASIS(s[2*i], s[2*i+1], sT[2*i], sT[2*i+1]);
	bases_infact.push_back(b);
      }
      bases_[i] = b;
    }
    
    j0_ = bases_[0]->j0();
  }
  
  template <class IBASIS, unsigned int DIM>
  CubeBasis<IBASIS,DIM>::CubeBasis(const FixedArray1D<bool,2*DIM>& bc) {
    for (unsigned int i = 0; i < DIM; i++) {
      // check whether the corresponding 1d basis already exists
      IBASIS* b = 0;
      for (typename list<IBASIS*>::const_iterator it(bases_infact.begin());
	   it != bases_infact.end(); ++it) {
	if (((*it)->get_s0()==1) == bc[2*i]
	    && ((*it)->get_s1()==1) == bc[2*i+1]) {
	  b = *it;
	  break;
	}
      }
      if (b == 0) {
	b = new IBASIS(bc[2*i], bc[2*i+1]);
	bases_infact.push_back(b);
      }
      bases_[i] = b;
    }
    
    j0_ = bases_[0]->j0();
  }
  
  template <class IBASIS, unsigned int DIM>
  CubeBasis<IBASIS,DIM>::~CubeBasis()
  {
    for (typename list<IBASIS*>::const_iterator it(bases_infact.begin());
	 it != bases_infact.end(); ++it)
      delete *it;
  }

  template <class IBASIS, unsigned int DIM>
  typename CubeBasis<IBASIS,DIM>::Index
  CubeBasis<IBASIS,DIM>::first_generator(const int j) const
  {
    assert(j >= j0());

    typename Index::type_type e;
    typename Index::translation_type k;
    for (unsigned int i = 0; i < DIM; i++) {
      typename IBASIS::Index lambda(bases()[i]->first_generator(j));
      k[i] = lambda.k();
    }
    
    return Index(j, e, k, this);
  }

  template <class IBASIS, unsigned int DIM>
  typename CubeBasis<IBASIS,DIM>::Index
  CubeBasis<IBASIS,DIM>::last_generator(const int j) const
  {
    assert(j >= j0());

    typename Index::type_type e;
    typename Index::translation_type k;
    for (unsigned int i = 0; i < DIM; i++) {
      typename IBASIS::Index lambda(bases()[i]->last_generator(j));
      k[i] = lambda.k();
    }

    return Index(j, e, k, this);
  }

  template <class IBASIS, unsigned int DIM>
  typename CubeBasis<IBASIS,DIM>::Index
  CubeBasis<IBASIS,DIM>::first_wavelet(const int j) const
  {
    assert(j >= j0());
    
    typename Index::type_type e;
    typename Index::translation_type k;
    for (unsigned int i = 0; i < DIM-1; i++) {
      typename IBASIS::Index lambda(bases()[i]->first_generator(j));
      k[i] = lambda.k();
    }
    typename IBASIS::Index lambda(bases()[DIM-1]->first_wavelet(j));
    k[DIM-1] = lambda.k();
    e[DIM-1] = 1;

    return Index(j, e, k, this);
  }

  template <class IBASIS, unsigned int DIM>
  typename CubeBasis<IBASIS,DIM>::Index
  CubeBasis<IBASIS,DIM>::last_wavelet(const int j) const
  {
    assert(j >= j0());
    
    typename Index::type_type e;
    typename Index::translation_type k;
    for (unsigned int i = 0; i < DIM; i++) {
      typename IBASIS::Index lambda(bases()[i]->last_wavelet(j));
      k[i] = lambda.k();
      e[i] = 1;
    }

    return Index(j, e, k, this);
  }

  template <class IBASIS, unsigned int DIM>
  void
  CubeBasis<IBASIS,DIM>::decompose(const InfiniteVector<double, Index>& c,
				   const int j0,
				   InfiniteVector<double, Index>& d) const {
    InfiniteVector<double, Index> help;
    for (typename InfiniteVector<double, Index>::const_iterator it(c.begin()), itend(c.end());
	 it != itend; ++it) {
      decompose_1(it.index(), j0, help); // calls help.clear() first
      d.add(*it, help);
    }
  }
  
  template <class IBASIS, unsigned int DIM>
  void
  CubeBasis<IBASIS,DIM>::decompose_t(const InfiniteVector<double, Index>& c,
				     const int j0,
				     InfiniteVector<double, Index>& d) const {
    InfiniteVector<double, Index> help;
    for (typename InfiniteVector<double, Index>::const_iterator it(c.begin()), itend(c.end());
	 it != itend; ++it) {
      decompose_t_1(it.index(), j0, help); // calls help.clear() first
      d.add(*it, help);
    }
  }

  template <class IBASIS, unsigned int DIM>
  void
  CubeBasis<IBASIS,DIM>::reconstruct(const InfiniteVector<double, Index>& c,
				     const int j,
				     InfiniteVector<double, Index>& d) const {
    for (typename InfiniteVector<double, Index>::const_iterator it(c.begin()), itend(c.end());
	 it != itend; ++it) {
      InfiniteVector<double, Index> help;
      reconstruct_1(it.index(), j, help);
      d.add(*it, help);
    }
  }
  
  template <class IBASIS, unsigned int DIM>
  void
  CubeBasis<IBASIS,DIM>::reconstruct_t(const InfiniteVector<double, Index>& c,
				       const int j,
				       InfiniteVector<double, Index>& d) const {
    for (typename InfiniteVector<double, Index>::const_iterator it(c.begin()), itend(c.end());
	 it != itend; ++it) {
      InfiniteVector<double, Index> help;
      reconstruct_t_1(it.index(), j, help);
      d.add(*it, help);
    }
  }

  template <class IBASIS, unsigned int DIM>
  void
  CubeBasis<IBASIS,DIM>::decompose_1(const Index& lambda,
				     const int j0,
				     InfiniteVector<double,Index>& c) const {
    assert(lambda.j() >= j0);
    c.clear();

    static MultiIndex<unsigned int,DIM> zero;
    if (lambda.e() == zero) {
      // a generator on a (possibly) fine level
      if (lambda.j() == j0) {
	// generators from the coarsest level can be copied
	c.set_coefficient(lambda, 1.0);
      } else {
	// j>j0, perform multiscale decomposition
	
	typedef typename IBASIS::Index IIndex;
	FixedArray1D<InfiniteVector<double,IIndex>,DIM> coeffs;
	for (unsigned int i = 0; i < DIM; i++)
	  bases_[i]->decompose_1(IIndex(lambda.j(), lambda.e()[i], lambda.k()[i], bases_[i]), lambda.j()-1, coeffs[i]);
	
	if (DIM == 2) {
	  // directly add all tensor product wavelets needed
	  for (typename InfiniteVector<double,IIndex>::const_iterator it1(coeffs[0].begin()), it1end(coeffs[0].end());
	       it1 != it1end; ++it1)
	    for (typename InfiniteVector<double,IIndex>::const_iterator it2(coeffs[1].begin()), it2end(coeffs[1].end());
		 it2 != it2end; ++it2) {
	      InfiniteVector<double,Index> d;
	      decompose_1(Index(it1.index().j(),
				typename Index::type_type(it1.index().e(), it2.index().e()),
				typename Index::translation_type(it1.index().k(), it2.index().k()),
				this),
			  j0, d);
	      c.add(*it1 * *it2, d);
	    }
	} else {
	  // prepare all tensor product wavelet indices needed (+values)
	  typedef std::list<std::pair<FixedArray1D<IIndex,DIM>,double> > list_type;
	  list_type indices;
	  FixedArray1D<IIndex,DIM> helpindex;
	  for (typename InfiniteVector<double,IIndex>::const_iterator itv(coeffs[0].begin()), itvend(coeffs[0].end());
	       itv != itvend; ++itv) {
	    helpindex[0] = itv.index();
	    indices.push_back(std::make_pair(helpindex, *itv));
	  }
	  for (unsigned int i(1); i < DIM; i++) {
	    list_type sofar;
	    sofar.swap(indices);
	    for (typename list_type::const_iterator it(sofar.begin()), itend(sofar.end());
		 it != itend; ++it) {
	      helpindex = it->first;
	      for (typename InfiniteVector<double,IIndex>::const_iterator itv(coeffs[i].begin()), itvend(coeffs[i].end());
		   itv != itvend; ++itv) {
		helpindex[i] = itv.index();
		indices.push_back(std::make_pair(helpindex, *itv * it->second));
	      }
	    }
	  }

	  // write results into c
	  typename Index::type_type help_e;
	  typename Index::translation_type help_k;
	  for (typename list_type::const_iterator it(indices.begin()), itend(indices.end());
	       it != itend; ++it) {
	    InfiniteVector<double,Index> d;
	    for (unsigned int i = 0; i < DIM; i++) {
	      help_e[i] = it->first[i].e();
	      help_k[i] = it->first[i].k();
	    }
	    decompose_1(Index(it->first[0].j(), help_e, help_k, this), j0, d);
	    c.add(it->second, d);
	  }
	}
      }
    } else {
      // the true wavelet coefficients don't have to be modified
      c.set_coefficient(lambda, 1.0);
    }
  }
  
  template <class IBASIS, unsigned int DIM>
  void
  CubeBasis<IBASIS,DIM>::decompose_t_1(const Index& lambda,
				       const int j0,
				       InfiniteVector<double, Index>& c) const {
    assert(lambda.j() >= j0);
    c.clear();

    static MultiIndex<unsigned int,DIM> zero;
    if (lambda.e() == zero) {
      // a generator on a (possibly) fine level
      if (lambda.j() == j0) {
	// generators from the coarsest level can be copied
	c.set_coefficient(lambda, 1.0);
      } else {
	// j>j0, perform multiscale decomposition
	
	typedef typename IBASIS::Index IIndex;
	FixedArray1D<InfiniteVector<double,IIndex>,DIM> coeffs;
	for (unsigned int i = 0; i < DIM; i++)
	  bases_[i]->decompose_t_1(IIndex(lambda.j(), lambda.e()[i], lambda.k()[i], bases_[i]), lambda.j()-1, coeffs[i]);
	
	if (DIM == 2) {
	  // directly add all tensor product wavelets needed
	  for (typename InfiniteVector<double,IIndex>::const_iterator it1(coeffs[0].begin()), it1end(coeffs[0].end());
	       it1 != it1end; ++it1)
	    for (typename InfiniteVector<double,IIndex>::const_iterator it2(coeffs[1].begin()), it2end(coeffs[1].end());
		 it2 != it2end; ++it2) {
	      InfiniteVector<double,Index> d;
	      decompose_t_1(Index(it1.index().j(),
				  typename Index::type_type(it1.index().e(), it2.index().e()),
				  typename Index::translation_type(it1.index().k(), it2.index().k()),
				  this),
			    j0, d);
	      c.add(*it1 * *it2, d);
	    }
	} else {
	  // prepare all tensor product wavelet indices needed (+values)
	  typedef std::list<std::pair<FixedArray1D<IIndex,DIM>,double> > list_type;
	  list_type indices;
	  FixedArray1D<IIndex,DIM> helpindex;
	  for (typename InfiniteVector<double,IIndex>::const_iterator itv(coeffs[0].begin()), itvend(coeffs[0].end());
	       itv != itvend; ++itv) {
	    helpindex[0] = itv.index();
	    indices.push_back(std::make_pair(helpindex, *itv));
	  }
	  for (unsigned int i(1); i < DIM; i++) {
	    list_type sofar;
	    sofar.swap(indices);
	    for (typename list_type::const_iterator it(sofar.begin()), itend(sofar.end());
		 it != itend; ++it) {
	      helpindex = it->first;
	      for (typename InfiniteVector<double,IIndex>::const_iterator itv(coeffs[i].begin()), itvend(coeffs[i].end());
		   itv != itvend; ++itv) {
		helpindex[i] = itv.index();
		indices.push_back(std::make_pair(helpindex, *itv * it->second));
	      }
	    }
	  }

	  // write results into c
	  typename Index::type_type help_e;
	  typename Index::translation_type help_k;
	  for (typename list_type::const_iterator it(indices.begin()), itend(indices.end());
	       it != itend; ++it) {
	    InfiniteVector<double,Index> d;
	    for (unsigned int i = 0; i < DIM; i++) {
	      help_e[i] = it->first[i].e();
	      help_k[i] = it->first[i].k();
	    }
	    decompose_t_1(Index(it->first[0].j(), help_e, help_k, this), j0, d);
	    c.add(it->second, d);
	  }
	}
      }
    } else {
      // the true wavelet coefficients don't have to be modified
      c.set_coefficient(lambda, 1.0);
    }
  }

  template <class IBASIS, unsigned int DIM>
  void
  CubeBasis<IBASIS,DIM>::reconstruct_1(const Index& lambda,
				       const int j,
				       InfiniteVector<double, Index>& c) const {
    if (lambda.j() >= j) {
      // then we can just copy \psi_\lambda
      c.add_coefficient(lambda, 1.0);
    } else {
      // reconstruct by recursion

      typedef typename IBASIS::Index IIndex;
      FixedArray1D<InfiniteVector<double,IIndex>,DIM> coeffs;
      for (unsigned int i = 0; i < DIM; i++)
	bases_[i]->reconstruct_1(IIndex(lambda.j(), lambda.e()[i], lambda.k()[i], bases_[i]), lambda.j()+1, coeffs[i]);
      
      if (DIM == 2) {
	// directly add all tensor product wavelets needed
	for (typename InfiniteVector<double,IIndex>::const_iterator it1(coeffs[0].begin()), it1end(coeffs[0].end());
	     it1 != it1end; ++it1)
	  for (typename InfiniteVector<double,IIndex>::const_iterator it2(coeffs[1].begin()), it2end(coeffs[1].end());
	       it2 != it2end; ++it2) {
	    InfiniteVector<double,Index> d;
	    reconstruct_1(Index(it1.index().j(),
				typename Index::type_type(it1.index().e(), it2.index().e()),
				typename Index::translation_type(it1.index().k(), it2.index().k()),
				this),
			  j, d);
	    c.add(*it1 * *it2, d);
	  }
      } else {
	// prepare all tensor product wavelet indices needed (+values)
	typedef std::list<std::pair<FixedArray1D<IIndex,DIM>,double> > list_type;
	list_type indices;
	FixedArray1D<IIndex,DIM> helpindex;
	for (typename InfiniteVector<double,IIndex>::const_iterator itv(coeffs[0].begin()), itvend(coeffs[0].end());
	     itv != itvend; ++itv) {
	  helpindex[0] = itv.index();
	  indices.push_back(std::make_pair(helpindex, *itv));
	}
	for (unsigned int i(1); i < DIM; i++) {
	  list_type sofar;
	  sofar.swap(indices);
	  for (typename list_type::const_iterator it(sofar.begin()), itend(sofar.end());
	       it != itend; ++it) {
	    helpindex = it->first;
	    for (typename InfiniteVector<double,IIndex>::const_iterator itv(coeffs[i].begin()), itvend(coeffs[i].end());
		 itv != itvend; ++itv) {
	      helpindex[i] = itv.index();
	      indices.push_back(std::make_pair(helpindex, *itv * it->second));
	    }
	  }
	}
	
	// write results into c
	typename Index::type_type help_e;
	typename Index::translation_type help_k;
	for (typename list_type::const_iterator it(indices.begin()), itend(indices.end());
	     it != itend; ++it) {
	  InfiniteVector<double,Index> d;
	  for (unsigned int i = 0; i < DIM; i++) {
	    help_e[i] = it->first[i].e();
	    help_k[i] = it->first[i].k();
	  }
	  reconstruct_1(Index(it->first[0].j(), help_e, help_k, this), j, d);
	  c.add(it->second, d);
	}
      }
    }
  }
  
  template <class IBASIS, unsigned int DIM>
  void
  CubeBasis<IBASIS,DIM>::reconstruct_t_1(const Index& lambda,
					 const int j,
					 InfiniteVector<double, Index>& c) const {
    if (lambda.j() >= j) {
      // then we can just copy \psi_\lambda
      c.add_coefficient(lambda, 1.0);
    } else {
      // reconstruct by recursion

      typedef typename IBASIS::Index IIndex;
      FixedArray1D<InfiniteVector<double,IIndex>,DIM> coeffs;
      for (unsigned int i = 0; i < DIM; i++)
	bases_[i]->reconstruct_t_1(IIndex(lambda.j(), lambda.e()[i], lambda.k()[i], bases_[i]), lambda.j()+1, coeffs[i]);
      
      if (DIM == 2) {
	// directly add all tensor product wavelets needed
	for (typename InfiniteVector<double,IIndex>::const_iterator it1(coeffs[0].begin()), it1end(coeffs[0].end());
	     it1 != it1end; ++it1)
	  for (typename InfiniteVector<double,IIndex>::const_iterator it2(coeffs[1].begin()), it2end(coeffs[1].end());
	       it2 != it2end; ++it2) {
	    InfiniteVector<double,Index> d;
	    reconstruct_t_1(Index(it1.index().j(),
				  typename Index::type_type(it1.index().e(), it2.index().e()),
				  typename Index::translation_type(it1.index().k(), it2.index().k()),
				  this),
			    j, d);
	    c.add(*it1 * *it2, d);
	  }
      } else {
	// prepare all tensor product wavelet indices needed (+values)
	typedef std::list<std::pair<FixedArray1D<IIndex,DIM>,double> > list_type;
	list_type indices;
	FixedArray1D<IIndex,DIM> helpindex;
	for (typename InfiniteVector<double,IIndex>::const_iterator itv(coeffs[0].begin()), itvend(coeffs[0].end());
	     itv != itvend; ++itv) {
	  helpindex[0] = itv.index();
	  indices.push_back(std::make_pair(helpindex, *itv));
	}
	for (unsigned int i(1); i < DIM; i++) {
	  list_type sofar;
	  sofar.swap(indices);
	  for (typename list_type::const_iterator it(sofar.begin()), itend(sofar.end());
	       it != itend; ++it) {
	    helpindex = it->first;
	    for (typename InfiniteVector<double,IIndex>::const_iterator itv(coeffs[i].begin()), itvend(coeffs[i].end());
		 itv != itvend; ++itv) {
	      helpindex[i] = itv.index();
	      indices.push_back(std::make_pair(helpindex, *itv * it->second));
	    }
	  }
	}
	
	// write results into c
	typename Index::type_type help_e;
	typename Index::translation_type help_k;
	for (typename list_type::const_iterator it(indices.begin()), itend(indices.end());
	     it != itend; ++it) {
	  InfiniteVector<double,Index> d;
	  for (unsigned int i = 0; i < DIM; i++) {
	    help_e[i] = it->first[i].e();
	    help_k[i] = it->first[i].k();
	  }
	  reconstruct_t_1(Index(it->first[0].j(), help_e, help_k, this), j, d);
	  c.add(it->second, d);
	}
      }
    }
  }
}
