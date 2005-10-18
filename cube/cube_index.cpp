// implementation for cube_index.h

namespace WaveletTL
{
  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  CubeIndex<IBASIS,DIM,CUBEBASIS>::CubeIndex(const CUBEBASIS* basis)
    : basis_(basis)
  {
    if (basis_ == 0) {
      j_ = 0; // invalid (e and k are initialized by zero automatically)
    } else {
      j_ = basis_->j0(); // coarsest level;
      // e_ is zero by default: generator
      for (unsigned int i = 0; i < DIM; i++)
	k_[i] = basis_->bases()[i]->DeltaLmin();
    }
  }

  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  CubeIndex<IBASIS,DIM,CUBEBASIS>::CubeIndex(const int j,
				   const type_type& e,
				   const translation_type& k,
				   const CUBEBASIS* basis)
    : basis_(basis), j_(j), e_(e), k_(k)
  {
  }

  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  CubeIndex<IBASIS,DIM,CUBEBASIS>::CubeIndex(const CubeIndex& lambda)
    : basis_(lambda.basis_), j_(lambda.j_), e_(lambda.e_), k_(lambda.k_)
  {
  }

  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  bool
  CubeIndex<IBASIS,DIM,CUBEBASIS>::operator == (const CubeIndex& lambda) const
  {
    return (j_ == lambda.j() &&
	    e_ == lambda.e() &&
	    k_ == lambda.k());
  }

  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  bool
  CubeIndex<IBASIS,DIM,CUBEBASIS>::operator < (const CubeIndex& lambda) const
  {
    // standard lexicographic order on (j,e,k),
    // we assume that e and k are already lexicographically ordered (cf. MultiIndex)
    return (j_ < lambda.j() ||
	    (j_ == lambda.j() && (e_ < lambda.e() ||
				  (e_ == lambda.e() && k_ < lambda.k()))));
  }

  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  CubeIndex<IBASIS,DIM,CUBEBASIS>&
  CubeIndex<IBASIS,DIM,CUBEBASIS>::operator ++ ()
  {
    bool eplusplus = false;
    for (int i = DIM-1; i >= 0; i--) {
      const int last_index = (e_[i] == 0
			      ? basis_->bases()[i]->DeltaRmax(j_)
			      : basis_->bases()[i]->Nablamax(j_));
      if (k_[i] == last_index) {
	k_[i] = (e_[i] == 0
		 ? basis_->bases()[i]->DeltaLmin()
		 : basis_->bases()[i]->Nablamin());
	eplusplus = (i == 0);
      } else {
	++k_[i];
	break;
      }
    }

    bool jplusplus = false;
    if (eplusplus) {
      for (int i = DIM-1; i >= 0; i--) {
	if (e_[i] == 1) {
 	  e_[i] = 0;
 	  jplusplus = (i == 0);
	} else {
	  ++e_[i];
	  break;
	}
      }
     
      if (!jplusplus) // then choose lowest translation index k=k(j,e)
	for (unsigned int i = 0; i < DIM; i++)
	  k_[i] = (e_[i] == 0
		   ? basis_->bases()[i]->DeltaLmin()
		   : basis_->bases()[i]->Nablamin());
    }

    if (jplusplus) {
      ++j_;
      // choose lowest type e=(0,...,0,1) and lowest translation index k=k(j,e)
      for (unsigned int i = 0; i < DIM-1; i++) {
	e_[i] = 0;
	k_[i] = basis_->bases()[i]->DeltaLmin();
      }
      e_[DIM-1] = 1;
      k_[DIM-1] = basis_->bases()[DIM-1]->Nablamin();
    }
    
    return *this;
  }

  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  CubeIndex<IBASIS,DIM,CUBEBASIS>
  first_generator(const CUBEBASIS* basis, const int j)
  {
    assert(j >= basis->j0());

    typename CubeIndex<IBASIS,DIM,CUBEBASIS>::type_type e;
    typename CubeIndex<IBASIS,DIM,CUBEBASIS>::translation_type k;
    for (unsigned int i = 0; i < DIM; i++) {
      typename IBASIS::Index lambda(first_generator<IBASIS>(basis->bases()[i], j));
      k[i] = lambda.k();
    }

    return CubeIndex<IBASIS,DIM,CUBEBASIS>(j, e, k, basis);
  }

  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  CubeIndex<IBASIS,DIM,CUBEBASIS>
  last_generator(const CUBEBASIS* basis, const int j)
  {
    assert(j >= basis->j0());

    typename CubeIndex<IBASIS,DIM,CUBEBASIS>::type_type e;
    typename CubeIndex<IBASIS,DIM,CUBEBASIS>::translation_type k;
    for (unsigned int i = 0; i < DIM; i++) {
      typename IBASIS::Index lambda(last_generator<IBASIS>(basis->bases()[i], j));
      k[i] = lambda.k();
    }

    return CubeIndex<IBASIS,DIM,CUBEBASIS>(j, e, k, basis);
  }

  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  CubeIndex<IBASIS,DIM,CUBEBASIS>
  first_wavelet(const CUBEBASIS* basis, const int j)
  {
    assert(j >= basis->j0());
    
    typename CubeIndex<IBASIS,DIM,CUBEBASIS>::type_type e;
    typename CubeIndex<IBASIS,DIM,CUBEBASIS>::translation_type k;
    for (unsigned int i = 0; i < DIM-1; i++) {
      typename IBASIS::Index lambda(first_generator<IBASIS>(basis->bases()[i], j));
      k[i] = lambda.k();
    }
    typename IBASIS::Index lambda(first_wavelet<IBASIS>(basis->bases()[DIM-1], j));
    k[DIM-1] = lambda.k();
    e[DIM-1] = 1;

    return CubeIndex<IBASIS,DIM,CUBEBASIS>(j, e, k, basis);
  }

  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  CubeIndex<IBASIS,DIM,CUBEBASIS>
  last_wavelet(const CUBEBASIS* basis, const int j)
  {
    assert(j >= basis->j0());
    
    typename CubeIndex<IBASIS,DIM,CUBEBASIS>::type_type e;
    typename CubeIndex<IBASIS,DIM,CUBEBASIS>::translation_type k;
    for (unsigned int i = 0; i < DIM; i++) {
      typename IBASIS::Index lambda(last_wavelet<IBASIS>(basis->bases()[i], j));
      k[i] = lambda.k();
      e[i] = 1;
    }

    return CubeIndex<IBASIS,DIM,CUBEBASIS>(j, e, k, basis);
  }
}
