// implementation for ldomain_index.h

namespace WaveletTL
{
  template <class IBASIS>
  LDomainIndex<IBASIS>::LDomainIndex(const LDomainBasis<IBASIS>* basis)
    : basis_(basis), p_(0)
  {
    if (basis_ == 0) {
      j_ = 0; // invalid (e and k are initialized by zero automatically)
    } else {
      j_ = basis_->j0(); // coarsest level;
      // e_ is zero by default: generator
      k_[0] = basis->basis00().DeltaLmin();
      k_[1] = basis->basis10().DeltaLmin()+1;
    }
  }

  template <class IBASIS>
  LDomainIndex<IBASIS>::LDomainIndex(const int j,
				     const type_type& e,
				     const int p,
				     const translation_type& k,
				     const LDomainBasis<IBASIS>* basis)
    : basis_(basis), j_(j), e_(e), p_(p), k_(k)
  {
  }

  template <class IBASIS>
  LDomainIndex<IBASIS>::LDomainIndex(const LDomainIndex& lambda)
    : basis_(lambda.basis_), j_(lambda.j_), e_(lambda.e_), p_(lambda.p_), k_(lambda.k_)
  {
  }

  template <class IBASIS>
  bool
  LDomainIndex<IBASIS>::operator == (const LDomainIndex& lambda) const
  {
    return (j_ == lambda.j() &&
	    e_ == lambda.e() &&
	    p_ == lambda.p() &&
	    k_ == lambda.k());
  }

  template <class IBASIS>
  LDomainIndex<IBASIS>&
  LDomainIndex<IBASIS>::operator ++ ()
  {
    // decide whether the patch number has to be increased
    bool pplusplus = false;
    for (int i = 1; i >= 0; i--) {
      // determine last possible translation index into the i-th direction,
      // this will in general depend on the current patch number
      int last_index = 0;
      if (e_[i] == 0) { // generator
	switch(p_) {
	case 0:
	  last_index = (i == 0
			? basis_->basis00().DeltaRmax(j_)
			: basis_->basis10().DeltaRmax(j_));
	  break;
	case 1:
	  last_index = basis_->basis01().DeltaRmax(j_)-1; // independent from i
	  break;
	case 2:
	  last_index = (i == 0
			? basis_->basis10().DeltaRmax(j_)
			: basis_->basis00().DeltaRmax(j_));
	  break;
	case 3:
	  last_index = (i == 0
			? basis_->basis00().DeltaRmax(j_)
			: 0); // by convention
	  break;
	case 4:
	  last_index = (i == 0
			? 0 // by convention
			: basis_->basis00().DeltaRmax(j_));
	  break;
	}
      } else { // wavelet, maximal translation index is independent from the patch number
	last_index = basis_->basis00().Nablamax(j_);
      }

      if (k_[i] == last_index) {
	// set k_[i] to the lowest possible translation index
	if (e_[i] == 0) { // generator, 
	  switch(p_) {
	  case 0:
	    k_[i] = (i == 0
		     ? basis_->basis00().DeltaLmin()
		     : basis_->basis10().DeltaLmin()+1);
	    break;
	  case 1:
	    k_[i] = basis_->basis01().DeltaLmin(); // independent from i
	    break;
	  case 2:
	    k_[i] = (i == 0
		     ? basis_->basis10().DeltaLmin()+1
		     : basis_->basis00().DeltaLmin());
	    break;
	  case 3:
	    k_[i] = (i == 0
		     ? basis_->basis00().DeltaLmin()
		     : 0); // by convention
	    break;
	  case 4:
	    k_[i] = (i == 0
		     ? 0 // by convention
		     : basis_->basis00().DeltaLmin());
	    break;
	  }
	} else { // wavelet, minimal translation index is independent from the patch number
	  k_[i] = basis_->basis00().Nablamin();
	}
	pplusplus = (i == 0);
      } else {
	++k_[i];
	break;
      }
    }

    bool eplusplus = false;
    if (pplusplus) {
      switch (p_) {
      case 0:
      case 1:
	p_++;
	break;
      case 2:
	if (e_[1] == 1) {
	  if (e_[0] == 0)
	    p_ = 4;
	  else
	    eplusplus = true; // no (1,1) wavelets on the interfaces
	} else p_ = 3;
	break;
      case 3:
	if (e_[0] == 1)
	  eplusplus = true;
	else
	  p_ = 4;
	break;
      case 4:
	eplusplus = true;
	break;
      }

      if (!eplusplus) { // then choose lowest translation index k=k(j,e,p)
	switch(p_) { // we know that p_>0
	case 1:
	  k_[0] = (e_[0] == 0
		   ? basis_->basis01().DeltaLmin()
		   : basis_->basis01().Nablamin());
	  k_[1] = (e_[1] == 0
		   ? basis_->basis01().DeltaLmin()
		   : basis_->basis01().Nablamin());
	  break;
	case 2:
	  k_[0] = (e_[0] == 0
		   ? basis_->basis10().DeltaLmin()+1
		   : basis_->basis10().Nablamin());
	  k_[1] = (e_[1] == 0
		   ? basis_->basis00().DeltaLmin()
		   : basis_->basis00().Nablamin());
	  break;
	case 3:
	  k_[0] = (e_[0] == 0
		   ? basis_->basis00().DeltaLmin()
		   : basis_->basis00().Nablamin());
	  k_[1] = 0; // by convention;
	  break;
	case 4:
	  k_[0] = 0; // by convention
	  k_[1] = (e_[1] == 0
		   ? basis_->basis00().DeltaLmin()
		   : basis_->basis00().Nablamin());
	}
      }
    } else return *this;

    bool jplusplus = false;
    if (eplusplus) {
      // advance e
      if (e_[0] == 1 && e_[1] == 1)
	jplusplus = true;
      else {
	if (e_[0] == 1)
	  e_[1] = 1;
	else {
	  if (e_[1] == 0)
	    e_[1] = 1;
	  else {
	    e_[0] = 1;
	    e_[1] = 0;
	  }
	}

	// choose lowest patch number ...
	p_ = 0;

	// ... and lowest translation index k = k(j,e,0)
	k_[0] = (e_[0] == 0
		 ? basis_->basis00().DeltaLmin()
		 : basis_->basis00().Nablamin());
	k_[1] = (e_[1] == 0
		 ? basis_->basis10().DeltaLmin()+1
		 : basis_->basis10().Nablamin());
      }

    } else return *this;

    if (jplusplus) {
      ++j_;

      // choose lowest type (we have to advance to a wavelet) ...
      e_[0] = 0;
      e_[1] = 1;
      
      // ... lowest patch number ...
      p_ = 0;
      
      // ... and lowest translation index k = k(j,(0,1),0)
      k_[0] = basis_->basis00().DeltaLmin();
      k_[1] = basis_->basis10().Nablamin();
    }
    
    return *this;
  }

  template <class IBASIS>
  bool
  LDomainIndex<IBASIS>::operator < (const LDomainIndex& lambda) const
  {
    // standard lexicographic order on (j,e,p,k),
    // we assume that e and k are already lexicographically ordered (cf. MultiIndex)
    return (j_ < lambda.j() ||
	    (j_ == lambda.j() &&
	     (e_ < lambda.e() ||
	      (e_ == lambda.e() &&
	       (p_ < lambda.p() ||
		(p_ == lambda.p() && k_ < lambda.k()
		 )
		)
	       )
	      )
	     )
	    );
  }
  
  template <class IBASIS>
  LDomainIndex<IBASIS>
  first_generator(const LDomainBasis<IBASIS>* basis, const int j)
  {
    assert(j >= basis->j0());

    typename LDomainIndex<IBASIS>::type_type e;

    // setup lowest translation index for e=(0,0), p=0
    typename LDomainIndex<IBASIS>::translation_type k;
    k[0] = basis->basis00().DeltaLmin();
    k[1] = basis->basis10().DeltaLmin()+1;
    
    return LDomainIndex<IBASIS>(j, e, 0, k, basis);
  }

  template <class IBASIS>
  LDomainIndex<IBASIS>
  last_generator(const LDomainBasis<IBASIS>* basis, const int j)
  {
    assert(j >= basis->j0());

    typename LDomainIndex<IBASIS>::type_type e;

    // setup highest translation index for e=(0,0), p=4
    typename LDomainIndex<IBASIS>::translation_type k; // k[0]=0 by convention
    k[1] = basis->basis00().DeltaRmax(j);
    
    return LDomainIndex<IBASIS>(j, e, 4, k, basis);
  }

  template <class IBASIS>
  LDomainIndex<IBASIS>
  first_wavelet(const LDomainBasis<IBASIS>* basis, const int j)
  {
    assert(j >= basis->j0());

    typename LDomainIndex<IBASIS>::type_type e;
    e[1] = 1;

    // setup lowest translation index for e=(0,1), p=0
    typename LDomainIndex<IBASIS>::translation_type k;
    k[0] = basis->basis00().DeltaLmin();
    k[1] = basis->basis10().Nablamin();
    
    return LDomainIndex<IBASIS>(j, e, 0, k, basis);
  }

  template <class IBASIS>
  LDomainIndex<IBASIS>
  last_wavelet(const LDomainBasis<IBASIS>* basis, const int j)
  {
    assert(j >= basis->j0());
    
    typename LDomainIndex<IBASIS>::type_type e;
    e[0] = e[1] = 1;

    // setup highest translation index for e=(1,1), p=2
    typename LDomainIndex<IBASIS>::translation_type k;
    k[0] = basis->basis10().DeltaRmax(j);
    k[1] = basis->basis00().DeltaRmax(j);
    
    return LDomainIndex<IBASIS>(j, e, 2, k, basis);
  }
}
