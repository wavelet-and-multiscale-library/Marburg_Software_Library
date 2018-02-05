// implementation for cube_basis.h

#include <map>
#include <numerics/gauss_data.h>

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
    delete_pointers = true;
  }

  template <class IBASIS, unsigned int DIM>
  CubeBasis<IBASIS,DIM>::CubeBasis(const CubeBasis<IBASIS, DIM>& other)
    : full_collection(other.full_collection), j0_(other.j0_), jmax_(other.jmax_),
      bases_infact(), bases_()
  {
    for (typename list<IBASIS*>::const_iterator it(other.bases_infact.begin());
         it != other.bases_infact.end(); ++it)
    {
        // make deep copy of the 1D bases:
        IBASIS* ibasis = new IBASIS(*(*it));
        bases_infact.push_back(ibasis);
        for (int i = 0; i < DIM; ++i)
        {
			if (*it == other.bases_[i])
				bases_[i] = ibasis;
		}
    }
    delete_pointers = true;
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
	//b = new IBASIS("",s[2*i], s[2*i+1], sT[2*i], sT[2*i+1]);
	b = new IBASIS(s[2*i], s[2*i+1], sT[2*i], sT[2*i+1]);
	bases_infact.push_back(b);
      }
      bases_[i] = b;
    }
    
    j0_ = bases_[0]->j0();
    delete_pointers = true;
  }

  template <class IBASIS, unsigned int DIM>
  CubeBasis<IBASIS,DIM>::CubeBasis(const FixedArray1D<int,2*DIM>& s) {
    for (unsigned int i = 0; i < DIM; i++) {
      // check whether the corresponding 1d basis already exists
      IBASIS* b = 0;
      for (typename list<IBASIS*>::const_iterator it(bases_infact.begin());
	   it != bases_infact.end(); ++it) {
	if ((*it)->get_s0() == s[2*i]
	    && (*it)->get_s1() == s[2*i+1]) {
	  b = *it;
	  break;
	}
      }
      if (b == 0) {
	b = new IBASIS(s[2*i], s[2*i+1]);
	bases_infact.push_back(b);
      }
      bases_[i] = b;
    }
    
    j0_ = bases_[0]->j0();
    delete_pointers = true;
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
    delete_pointers = true;
  }
  
  template <class IBASIS, unsigned int DIM>
  CubeBasis<IBASIS,DIM>::CubeBasis(const FixedArray1D<IBASIS*,DIM> bases) {
    for (unsigned int i = 0; i < DIM; i++) {
      bases_infact.push_back(bases[i]);
      bases_[i] = bases[i];
    }
    
    j0_ = bases_[0]->j0();
    delete_pointers = false;
  }

  template <class IBASIS, unsigned int DIM>
  CubeBasis<IBASIS,DIM>::~CubeBasis()
  {
    if (delete_pointers) {
      for (typename list<IBASIS*>::const_iterator it(bases_infact.begin());
	   it != bases_infact.end(); ++it)
	delete *it;
    }
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

    static typename Index::type_type zero;
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

    static typename Index::type_type zero;
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
      
      if (DIM == 2) { // TODO: template specialization here!!!
	// directly add all tensor product wavelets needed
 	if (lambda.j()+1 >= j) {
	  for (typename InfiniteVector<double,IIndex>::const_iterator it1(coeffs[0].begin()), it1end(coeffs[0].end());
	       it1 != it1end; ++it1)
	    for (typename InfiniteVector<double,IIndex>::const_iterator it2(coeffs[1].begin()), it2end(coeffs[1].end());
		 it2 != it2end; ++it2) {
	      c.add_coefficient(Index(lambda.j()+1,
				      typename Index::type_type(it1.index().e(), it2.index().e()),
				      typename Index::translation_type(it1.index().k(), it2.index().k()),
				      this),
				*it1 * *it2);
	    }
 	} else {
	  for (typename InfiniteVector<double,IIndex>::const_iterator it1(coeffs[0].begin()), it1end(coeffs[0].end());
	       it1 != it1end; ++it1)
	    for (typename InfiniteVector<double,IIndex>::const_iterator it2(coeffs[1].begin()), it2end(coeffs[1].end());
		 it2 != it2end; ++it2) {
	      InfiniteVector<double,Index> d;
	      reconstruct_1(Index(lambda.j()+1,
				  typename Index::type_type(it1.index().e(), it2.index().e()),
				  typename Index::translation_type(it1.index().k(), it2.index().k()),
				  this),
			    j, d);
	      c.add(*it1 * *it2, d);
	    }
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
 	if (lambda.j()+1 >= j) {
	  for (typename list_type::const_iterator it(indices.begin()), itend(indices.end());
	       it != itend; ++it) {
	    for (unsigned int i = 0; i < DIM; i++) {
	      help_e[i] = it->first[i].e();
	      help_k[i] = it->first[i].k();
	    }
	    c.add_coefficient(Index(lambda.j()+1, help_e, help_k, this),
			      it->second);
	  }
	} else {
	  for (typename list_type::const_iterator it(indices.begin()), itend(indices.end());
	       it != itend; ++it) {
	    InfiniteVector<double,Index> d;
	    for (unsigned int i = 0; i < DIM; i++) {
	      help_e[i] = it->first[i].e();
	      help_k[i] = it->first[i].k();
	    }
	    reconstruct_1(Index(lambda.j()+1, help_e, help_k, this), j, d);
	    c.add(it->second, d);
	  }
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
 	if (lambda.j()+1 >= j) {
	  for (typename InfiniteVector<double,IIndex>::const_iterator it1(coeffs[0].begin()), it1end(coeffs[0].end());
	       it1 != it1end; ++it1)
	    for (typename InfiniteVector<double,IIndex>::const_iterator it2(coeffs[1].begin()), it2end(coeffs[1].end());
		 it2 != it2end; ++it2) {
	      c.add_coefficient(Index(lambda.j()+1,
				      typename Index::type_type(it1.index().e(), it2.index().e()),
				      typename Index::translation_type(it1.index().k(), it2.index().k()),
				      this),
				*it1 * *it2);
	    }
	} else {
	  for (typename InfiniteVector<double,IIndex>::const_iterator it1(coeffs[0].begin()), it1end(coeffs[0].end());
	       it1 != it1end; ++it1)
	    for (typename InfiniteVector<double,IIndex>::const_iterator it2(coeffs[1].begin()), it2end(coeffs[1].end());
		 it2 != it2end; ++it2) {
	      InfiniteVector<double,Index> d;
	      reconstruct_t_1(Index(lambda.j()+1,
				    typename Index::type_type(it1.index().e(), it2.index().e()),
				    typename Index::translation_type(it1.index().k(), it2.index().k()),
				    this),
			      j, d);
	      c.add(*it1 * *it2, d);
	    }
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
 	if (lambda.j()+1 >= j) {
	  for (typename list_type::const_iterator it(indices.begin()), itend(indices.end());
	       it != itend; ++it) {
	    for (unsigned int i = 0; i < DIM; i++) {
	      help_e[i] = it->first[i].e();
	      help_k[i] = it->first[i].k();
	    }
	    c.add_coefficient(Index(lambda.j()+1, help_e, help_k, this),
			      it->second);
	  }
	} else {
	  for (typename list_type::const_iterator it(indices.begin()), itend(indices.end());
	       it != itend; ++it) {
	    InfiniteVector<double,Index> d;
	    for (unsigned int i = 0; i < DIM; i++) {
	      help_e[i] = it->first[i].e();
	      help_k[i] = it->first[i].k();
	    }
	    reconstruct_t_1(Index(lambda.j()+1, help_e, help_k, this), j, d);
	    c.add(it->second, d);
	  }
	}
      }
    }
  }

  template <class IBASIS, unsigned int DIM>
  void
  CubeBasis<IBASIS,DIM>::support(const Index& lambda, Support& supp) const
  {
    const unsigned int jplus = multi_degree(lambda.e()) > 0 ? 1 : 0;
    supp.j = lambda.j() + jplus;
    for (unsigned int i(0); i < DIM; i++) {
      bases()[i]->support(typename IBASIS::Index(lambda.j(),
						 lambda.e()[i],
						 lambda.k()[i],
						 bases()[i]),
			  supp.a[i], supp.b[i]);
      if (lambda.e()[i] == 0 && jplus > 0) {
	supp.a[i] *= 2;
	supp.b[i] *= 2;
      }
    }
  }
  
  template <class IBASIS, unsigned int DIM>
  void
  CubeBasis<IBASIS,DIM>::expand
  (const Function<DIM>* f,
   const bool primal,
   const int jmax,
   InfiniteVector<double,Index>& coeffs) const
  {
    if (primal == false)        //! integrate against primal wavelets and generators
    {
      for (Index lambda = first_generator(j0());;++lambda)
      {
	    const double coeff = integrate(f, lambda);
	    if (fabs(coeff)>1e-15)
	    {
	      coeffs.set_coefficient(lambda, coeff);
 	    }
 	    if (lambda == last_wavelet(jmax))
 	    {
 	      break;
        }
      }
    }    
    else                        //! integrate against dual wavelets and generators (we integrate against the primal wavelets
                                //! and generators and multiply the resulting coefficients with the inverse of the primal gramian)

    {
#ifdef P_POISSON
      //! integrate against the primal wavelets and generators
      InfiniteVector<double,Index> dual_coeffs;
      for (Index lambda = first_generator(j0());;++lambda)
      {
	    const double coeff = integrate(f, lambda);
        dual_coeffs.set_coefficient(lambda, coeff);
 	    if (lambda == last_wavelet(jmax))
 	    {
 	      break;
        }
      }

      //! setup the primal gramian
      SparseMatrix<double> A_Lambda;
      //set<typename PROBLEM::WaveletBasis::Index>& Lambda;
      std::set<Index> Lambda;

      dual_coeffs.support(Lambda);
      A_Lambda.resize(Lambda.size(), Lambda.size());

      typedef typename SparseMatrix<double>::size_type size_type;

      size_type row = 0;
      //typedef typename PROBLEM::Index Index;

      for (typename std::set<Index>::const_iterator it1(Lambda.begin()), itend(Lambda.end()); it1 != itend; ++it1, ++row)
      {
	    // double d1 = preconditioned ? P.D(*it1) : 1.0;
	    std::list<size_type> indices;
	    std::list<double> entries;

	    size_type column = 0;

	    for (typename std::set<Index>::const_iterator it2(Lambda.begin()); it2 != itend; ++it2, ++column)
	    {
          double entry = integrate(*it2, *it1);
	      if (fabs(entry) > 1e-15)
	      {
		    indices.push_back(column);
		    //entries.push_back(entry / (preconditioned ? d1 * P.D(*it2) : 1.0));
		    entries.push_back(entry);
	      }
	    }

	    A_Lambda.set_row(row, indices, entries);
      }


      //! solve system of linear equations
      Vector<double> F_Lambda(Lambda.size());

      // copy dual_coeffs into F_Lambda
      unsigned int id = 0;
      for (typename std::set<Index>::const_iterator it = Lambda.begin(), itend = Lambda.end(); it != itend; ++it, ++id)
      {
	    F_Lambda[id] = dual_coeffs.get_coefficient(*it);
      }

      // setup initial approximation xk
      Vector<double> xk(Lambda.size());
      id = 0;
      for (typename std::set<Index>::const_iterator it = Lambda.begin(), itend = Lambda.end(); it != itend; ++it, ++id)
      {
	    //xk[id] = dual_coeffs.get_coefficient(*it);
	    xk[id] = 0.0;
      }

      unsigned int iterations = 0;
      CG(A_Lambda, F_Lambda, xk, 1e-15, 150, iterations);

      id = 0;
      for (typename std::set<Index>::const_iterator it = Lambda.begin(), itend = Lambda.end(); it != itend; ++it, ++id)
      {
	    coeffs.set_coefficient(*it, xk[id]);
      }

#endif
    }  // else END

  
  }
  
  template <class IBASIS, unsigned int DIM>
  double
  CubeBasis<IBASIS,DIM>::integrate
  (const Function<DIM>* f,
   const Index& lambda) const
  {
    // f(v) = \int_0^1 g(t)v(t) dt
    
    double r = 0;
    
    // first compute supp(psi_lambda)
    Support supp;
    support(lambda, supp);
    
    // setup Gauss points and weights for a composite quadrature formula:
    const int N_Gauss = 5;
    const double h = ldexp(1.0, -supp.j); // granularity for the quadrature
    FixedArray1D<Array1D<double>,DIM> gauss_points, gauss_weights, v_values;
    for (unsigned int i = 0; i < DIM; i++) {
      gauss_points[i].resize(N_Gauss*(supp.b[i]-supp.a[i]));
      gauss_weights[i].resize(N_Gauss*(supp.b[i]-supp.a[i]));
      for (int patch = supp.a[i]; patch < supp.b[i]; patch++)
	for (int n = 0; n < N_Gauss; n++) {
	  gauss_points[i][(patch-supp.a[i])*N_Gauss+n]
	    = h*(2*patch+1+GaussPoints[N_Gauss-1][n])/2.;
	  gauss_weights[i][(patch-supp.a[i])*N_Gauss+n]
	    = h*GaussWeights[N_Gauss-1][n];
	}
    }

    // compute the point values of the integrand (where we use that it is a tensor product)
    for (unsigned int i = 0; i < DIM; i++)
      bases()[i]->evaluate(0,
			   typename IBASIS::Index(lambda.j(),
						  lambda.e()[i],
						  lambda.k()[i],
						  bases()[i]),
			   gauss_points[i], v_values[i]);
    
    // iterate over all points and sum up the integral shares
    int index[DIM]; // current multiindex for the point values
    for (unsigned int i = 0; i < DIM; i++)
      index[i] = 0;
    
    Point<DIM> x;
    while (true) {
      for (unsigned int i = 0; i < DIM; i++)
	x[i] = gauss_points[i][index[i]];
      double share = f->value(x);
      for (unsigned int i = 0; i < DIM; i++)
	share *= gauss_weights[i][index[i]] * v_values[i][index[i]];
      r += share;

      // "++index"
      bool exit = false;
      for (unsigned int i = 0; i < DIM; i++) {
	if (index[i] == N_Gauss*(supp.b[i]-supp.a[i])-1) {
	  index[i] = 0;
	  exit = (i == DIM-1);
	} else {
	  index[i]++;
	  break;
	}
      }
      if (exit) break;
    }
    
    return r;
  }

  template <class IBASIS, unsigned int DIM>
  double
  CubeBasis<IBASIS,DIM>::evaluate(const unsigned int derivative, const Index& lambda, const Point<DIM> x) const
  {
    double value = 1.0;
    
    if (derivative > 0)
    {
      cout << "CubeBasis.evaluate: evaluation of derivative not implemented correctly!!" << endl;
      exit(1);
    }

    for (unsigned int i = 0; i < DIM; i++) // loop through components of the tensor product
      value *= bases_[i]->evaluate(derivative,
                                   typename IBASIS::Index(lambda.j(), lambda.e()[i], lambda.k()[i], bases_[i]),
                                   x[i]);

    return value;
  }
  
  template <class IBASIS, unsigned int DIM>
  void
  CubeBasis<IBASIS,DIM>::setup_full_collection()
  {
    if (jmax_ == -1 || jmax_ < j0_) {
      cout << "CubeBasis<IBASIS,DIM>::setup_full_collection(): specify a maximal level of resolution first!" << endl;
      abort();
    }   
    
    int degrees_of_freedom = 1;

    for (unsigned int i = 0; i < DIM; i++) {
      degrees_of_freedom *= bases_[i]->Deltasize(jmax_+1);
    }
    cout << "total degrees of freedom between j0_ and jmax_ is " << degrees_of_freedom << endl;

    cout << "setting up collection of wavelet indices..." << endl;
    full_collection.resize(degrees_of_freedom);
    int k = 0;
    for (Index ind = first_generator(j0_); ind <= last_wavelet(jmax_); ++ind) {
      //cout << ind << " " << ind.number() << endl;
      full_collection[k] = ind;
      k++;
    }
    cout << "done setting up collection of wavelet indices..." << endl;

  }
  
  template <class IBASIS, unsigned int DIM>
  double
  CubeBasis<IBASIS,DIM>::integrate(const Index& lambda, const Index& mu) const
  {

    double r = 0.0;

    // first decide whether the supports of psi_lambda and psi_mu intersect
    Support supp;

    if (intersect_supports(*this, lambda, mu, supp))
    {

	  // setup Gauss points and weights for a composite quadrature formula:
	  const int doe = DIM * primal_polynomial_degree();
	  const int N_Gauss = (doe+1)/2;
	  const double h = ldexp(1.0, -supp.j); // granularity for the quadrature (h=2^-j)
	  FixedArray1D<Array1D<double>,DIM> gauss_points, gauss_weights;
	  for (unsigned int i = 0; i < DIM; i++)
	  {
	    gauss_points[i].resize(N_Gauss*(supp.b[i]-supp.a[i]));
	    gauss_weights[i].resize(N_Gauss*(supp.b[i]-supp.a[i]));
	    for (int patch = supp.a[i]; patch < supp.b[i]; patch++)
	    {
	      for (int n = 0; n < N_Gauss; n++)
	      {
	        gauss_points[i][(patch-supp.a[i])*N_Gauss+n] = h*(2*patch+1+GaussPoints[N_Gauss-1][n])/2.;
	        gauss_weights[i][(patch-supp.a[i])*N_Gauss+n] = h*GaussWeights[N_Gauss-1][n];
	      }
        }
  	  }


	  // compute point values of the integrand (where we use that it is a tensor product)
	  FixedArray1D<Array1D<double>,DIM>
	    psi_lambda_values,     // values of the components of psi_lambda at gauss_points[i]
	    psi_mu_values;         // -"-, for psi_mu

	  for (unsigned int i = 0; i < DIM; i++)
      {
	    bases()[i]->evaluate(0,
                             typename IBASIS::Index(lambda.j(), lambda.e()[i], lambda.k()[i], bases()[i]),
                             gauss_points[i], psi_lambda_values[i]);

	    bases()[i]->evaluate(0,
                             typename IBASIS::Index(mu.j(), mu.e()[i], mu.k()[i], bases()[i]),
                             gauss_points[i], psi_mu_values[i]);
	  }

	  // iterate over all points and sum up the integral shares
	  int index[DIM]; // current multiindex for the point values
	  for (unsigned int i = 0; i < DIM; i++)
	  {
	    index[i] = 0;
      }

	  double weights, share;

	  while (true)
	  {

	    // product of current Gauss weights
	    weights = 1.0;
	    for (unsigned int i = 0; i < DIM; i++)
	    {
	      weights *= gauss_weights[i][index[i]];
        }

	    // compute the share psi_lambda(x)psi_mu(x)
	    share = weights;
	    for (unsigned int i = 0; i < DIM; i++)
	    {
	      share *= psi_lambda_values[i][index[i]] * psi_mu_values[i][index[i]];
	    }
	    r += share;

	    // "++index"
	    bool exit = false;
	    for (unsigned int i = 0; i < DIM; i++)
	    {
	      if (index[i] == N_Gauss*(supp.b[i]-supp.a[i])-1)
	      {
		    index[i] = 0;
		    exit = (i == DIM-1);
	      }
	      else
	      {
		    index[i]++;
		    break;
	      }
	    }
	    if (exit) break;
	  }

    }     // if (intersect_supports(basis_, lambda, mu, supp)) ENDE

    return r;
  }

}
