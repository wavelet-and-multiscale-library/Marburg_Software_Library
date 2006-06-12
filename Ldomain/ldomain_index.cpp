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
      k_[0] = basis->basis1d().DeltaLmin();
      k_[1] = k_[0]+1;
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
      // determine the highest possible translation index into the i-th direction,
      // this will in general depend on the current patch number
      int last_index = 0;
      if (e_[i] == 0) { // generator in the i-th direction
	switch(p_) {
	case 0:
	case 1:
	case 2:
	  last_index = basis_->basis1d().DeltaRmax(j_)-1;
	  break;
	case 3:
	  last_index = (i == 0
			? basis_->basis1d().DeltaRmax(j_)-1
			: 0); // by convention
	  break;
	case 4:
	  last_index = (i == 0
			? 0 // by convention
			: basis_->basis1d().DeltaRmax(j_)-1);
	  break;
	}
      } else {
	// wavelet in the i-th direction,
	// the maximal translation index is independent from the patch number
	last_index = basis_->basis1d().Nablamax(j_); // should be (1<<j)-1
      }

      if (k_[i] == last_index) {
	// reset k_[i] to the lowest possible translation index
	if (e_[i] == 0) { // generator, 
	  switch(p_) {
	  case 0:
	  case 1:
	  case 2:
	    k_[i] = basis_->basis1d().DeltaLmin()+1;
	    break;
	  case 3:
	    k_[i] = (i == 0
		     ? basis_->basis1d().DeltaLmin()+1
		     : 0); // by convention
	    break;
	  case 4:
	    k_[i] = (i == 0
		     ? 0 // by convention
		     : basis_->basis1d().DeltaLmin()+1);
	    break;
	  }
	} else { // wavelet, minimal translation index is independent from the patch number
	  k_[i] = basis_->basis1d().Nablamin(); // should be 0
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
	    p_ = 4; // there are no (0,1) wavelets on the interface 3
	  else
	    eplusplus = true; // there are no (1,1) wavelets on the interfaces
	} else p_ = 3;
	break;
      case 3:
	if (e_[0] == 1)
	  eplusplus = true; // there are no (1,*) wavelets on the interface 4
	else
	  p_ = 4;
	break;
      case 4:
	eplusplus = true; // highest patch number reached
	break;
      }

      if (!eplusplus) { // then choose lowest translation index k=k(j,e,p)
	switch(p_) { // we know that p_>0
	case 1:
	case 2:
	  k_[0] = (e_[0] == 0
		   ? basis_->basis1d().DeltaLmin()+1
		   : basis_->basis1d().Nablamin());
	  k_[1] = (e_[1] == 0
		   ? basis_->basis1d().DeltaLmin()+1
		   : basis_->basis1d().Nablamin());
	  break;
	case 3:
	  k_[0] = (e_[0] == 0
		   ? basis_->basis1d().DeltaLmin()+1
		   : basis_->basis1d().Nablamin());
	  k_[1] = 0; // by convention;
	  break;
	case 4:
	  k_[0] = 0; // by convention
	  k_[1] = (e_[1] == 0
		   ? basis_->basis1d().DeltaLmin()+1
		   : basis_->basis1d().Nablamin());
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
		 ? basis_->basis1d().DeltaLmin()+1
		 : basis_->basis1d().Nablamin());
	k_[1] = (e_[1] == 0
		 ? basis_->basis1d().DeltaLmin()+1
		 : basis_->basis1d().Nablamin());
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
      k_[0] = basis_->basis1d().DeltaLmin()+1;
      k_[1] = basis_->basis1d().Nablamin();
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
  const int
  LDomainIndex<IBASIS>::number() const
  {
    const int ecode(e()[0]+2*e()[1]);
    
    const int DeltaLmin = basis_->basis1d().DeltaLmin();
    const int Nablamin  = basis_->basis1d().Nablamin();
    const int Deltasize = basis_->basis1d().Deltasize(j());
    const int Nablasize = basis_->basis1d().Nablasize(j());

    switch(ecode) {
    case 0: // e=(0,0)
      switch(p()) {
      case 0:
	return
	  (k()[0]-DeltaLmin-1)*(Deltasize-2)
	  + (k()[1]-DeltaLmin-1);
	break;
      case 1:
	return
	  (Deltasize-2)*(Deltasize-2)
	  + (k()[0]-DeltaLmin-1)*(Deltasize-2)
	  + (k()[1]-DeltaLmin-1);
	break;
      case 2:
	return
	  2*(Deltasize-2)*(Deltasize-2)
	  + (k()[0]-DeltaLmin-1)*(Deltasize-2)
	  + (k()[1]-DeltaLmin-1);
	break;
      case 3:
	return
	  3*(Deltasize-2)*(Deltasize-2)
	  + (k()[0]-DeltaLmin-1);
	break;
      case 4:
	return
	  3*(Deltasize-2)*(Deltasize-2)
	  + (Deltasize-2)
	  + (k()[1]-DeltaLmin-1);
	break;
      default:
	return 42; // for the compiler, this will never happen
      }
      break;
    case 1: // e=(1,0)
      switch(p()) {
      case 0:
	return
	  basis_->Deltasize(j())
	  + basis_->Nabla01size(j())
	  + (k()[0]-Nablamin)*(Deltasize-2)
	  + (k()[1]-DeltaLmin-1);
	break;
      case 1:
	return
	  basis_->Deltasize(j())
	  + basis_->Nabla01size(j())
	  + (Deltasize-2)*Nablasize
	  + (k()[0]-Nablamin)*(Deltasize-2)
	  + (k()[1]-DeltaLmin-1);
	break;
      case 2:
	return
	  basis_->Deltasize(j())
	  + basis_->Nabla01size(j())
	  + 2*(Deltasize-2)*Nablasize
	  + (k()[0]-Nablamin)*(Deltasize-2)
	  + (k()[1]-DeltaLmin-1);
	break;
      case 3:
	return
	  basis_->Deltasize(j())
	  + basis_->Nabla01size(j())
	  + 3*(Deltasize-2)*Nablasize
	  + (k()[0]-Nablamin);
	break;
      default:
	return 42; // for the compiler, this will never happen
      }
      break;
    case 2: // e=(0,1)
      switch(p()) {
      case 0:
	return
	  basis_->Deltasize(j())
	  + (k()[0]-DeltaLmin-1)*Nablasize
	  + (k()[1]-Nablamin);
	break;
      case 1:
	return
	  basis_->Deltasize(j())
	  + (Deltasize-2)*Nablasize
	  + (k()[0]-DeltaLmin-1)*Nablasize
	  + (k()[1]-Nablamin);
	break;
      case 2:
	return
	  basis_->Deltasize(j())
	  + 2*(Deltasize-2)*Nablasize
	  + (k()[0]-DeltaLmin-1)*Nablasize
	  + (k()[1]-Nablamin);
	break;
      case 4:
	return
	  basis_->Deltasize(j())
	  + 3*(Deltasize-2)*Nablasize
	  + (k()[1]-Nablamin);
	break;
      default:
	return 42; // for the compiler, this will never happen
      }
    case 3: // e=(1,1)
      switch(p()) {
      case 0:
	return
	  basis_->Deltasize(j())
	  + basis_->Nabla01size(j())
	  + basis_->Nabla10size(j())
	  + (k()[0]-Nablamin)*Nablasize
	  + (k()[1]-Nablamin);
	break;
      case 1:
	return
	  basis_->Deltasize(j())
	  + basis_->Nabla01size(j())
	  + basis_->Nabla10size(j())
	  + Nablasize*Nablasize
	  + (k()[0]-Nablamin)*Nablasize
	  + (k()[1]-Nablamin);
	break;
      case 2:
	return
	  basis_->Deltasize(j())
	  + basis_->Nabla01size(j())
	  + basis_->Nabla10size(j())
	  + 2*Nablasize*Nablasize
	  + (k()[0]-Nablamin)*Nablasize
	  + (k()[1]-Nablamin);
	break;
      default:
	return 42; // for the compiler, this will never happen
      }
      break;
    default:
      return 42; // for the compiler, this will never happen
    }
  }
  
  template <class IBASIS>
  LDomainIndex<IBASIS>
  first_generator(const LDomainBasis<IBASIS>* basis, const int j)
  {
    assert(j >= basis->j0());

    typename LDomainIndex<IBASIS>::type_type e;

    // setup lowest translation index for e=(0,0), p=0
    typename LDomainIndex<IBASIS>::translation_type k(basis->basis1d().DeltaLmin()+1,
						      basis->basis1d().DeltaLmin()+1);
    
    return LDomainIndex<IBASIS>(j, e, 0, k, basis);
  }

  template <class IBASIS>
  LDomainIndex<IBASIS>
  last_generator(const LDomainBasis<IBASIS>* basis, const int j)
  {
    assert(j >= basis->j0());

    typename LDomainIndex<IBASIS>::type_type e;

    // setup highest translation index for e=(0,0), p=4
    typename LDomainIndex<IBASIS>::translation_type k(0, basis->basis1d().DeltaRmax(j)-1);
    
    return LDomainIndex<IBASIS>(j, e, 4, k, basis);
  }

  template <class IBASIS>
  LDomainIndex<IBASIS>
  first_wavelet(const LDomainBasis<IBASIS>* basis, const int j)
  {
    assert(j >= basis->j0());

    typename LDomainIndex<IBASIS>::type_type e(0, 1);

    // setup lowest translation index for e=(0,1), p=0
    typename LDomainIndex<IBASIS>::translation_type k(basis->basis1d().DeltaLmin()+1,
						      basis->basis1d().Nablamin());
    
    return LDomainIndex<IBASIS>(j, e, 0, k, basis);
  }

  template <class IBASIS>
  LDomainIndex<IBASIS>
  first_wavelet(const LDomainBasis<IBASIS>* basis,
		const int j,
		const typename LDomainIndex<IBASIS>::type_type& ewish)
  {
    assert(j >= basis->j0());

    typename LDomainIndex<IBASIS>::type_type e(ewish);
    
    // setup lowest translation index appropriately
    typename LDomainIndex<IBASIS>::translation_type k;
    const int ecode(e[0]+2*e[1]);
    if (ecode == 0) {
      // e = (0,0)
      k[0] = k[1] = basis->basis1d().DeltaLmin()+1;
    } else {
      if (ecode == 1) {
	// e = (1,0)
	k[0] = basis->basis1d().Nablamin();
	k[1] = basis->basis1d().DeltaLmin()+1;
      } else {
	if (ecode == 2) {
	  // e = (0,1)
	  k[0] = basis->basis1d().DeltaLmin()+1;
	  k[1] = basis->basis1d().Nablamin();
	} else {
	  // e = (1,1)
	  k[0] = k[1] = basis->basis1d().Nablamin();
	}
      }
    }
    
    return LDomainIndex<IBASIS>(j, e, 0, k, basis);
  }
  
  template <class IBASIS>
  LDomainIndex<IBASIS>
  last_wavelet(const LDomainBasis<IBASIS>* basis, const int j)
  {
    assert(j >= basis->j0());
    
    typename LDomainIndex<IBASIS>::type_type e(1, 1);

    // setup highest translation index for e=(1,1), p=2
    typename LDomainIndex<IBASIS>::translation_type k(basis->basis1d().Nablamax(j),
						      basis->basis1d().Nablamax(j));
    
    return LDomainIndex<IBASIS>(j, e, 2, k, basis);
  }
}
