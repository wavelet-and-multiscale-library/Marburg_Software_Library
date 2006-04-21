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
  CubeIndex<IBASIS,DIM,CUBEBASIS>::CubeIndex(const unsigned int num,
					     const CUBEBASIS* basis)
    : basis_(basis)
  {
    // to be decreased successively
    unsigned int act_num = num;

    const unsigned int j0 = basis_->j0();
    unsigned int j = j0;

    unsigned int tmp2 = 0;

    bool gen = 0;

    // determine the level of the index
    while (true) {
      unsigned int tmp = 1;
      for (unsigned int i = 0; i < DIM; i++)
	tmp *= (basis_->bases()[i])->Deltasize(j);
      if (tmp > num)
	break;
      tmp2 = tmp;
      j++;
    }
    if (j == j0) {
      j_ = j;
      gen = 1;
    }
    else
      j_ = j-1;

    act_num -= tmp2;

    tmp2 = 0;

    // determine the type of the index 
    unsigned int tmp2_old = 0;
    if (gen)
      e_ = type_type();
    else {
      MultiIndex<int,DIM> type;
      type[DIM-1] = 1;
      bool exit = 0;
      // loop over types
      while (true) {
	unsigned int tmp = 1;
	for (unsigned int i = 0; i < DIM; i++) {
	  if (type[i] == 0)
	    tmp *= (basis_->bases()[i])->Deltasize(j_);
	  else
	    tmp *= (basis_->bases()[i])->Nablasize(j_);
	}
	tmp2_old = tmp2;
	tmp2 += tmp;
	
	// right type found
	if (tmp2 > act_num)
	  break;

	// determine next type
	for (unsigned int i = DIM-1; i >= 0; i--) {
	  if ( type[i] == 1 ) {
	    type[i] = 0;
	    exit = (i == 0);
	  }
	  else {
	    type[i]++;
	    break;
	  }
	}
	if (exit)
	  break;
      }//end while
      e_ = type;
    }// end else

    act_num -= tmp2_old;
    tmp2 = 1;

    // determine the position of the index
    for (unsigned int i = DIM-1; i > 0; i--) {
      if (e_[i] == 0) {
	tmp2 *= (basis_->bases()[i])->Deltasize(j_);
      }
      else {
	tmp2 *= (basis_->bases()[i])->Nablasize(j_);
      }
    }


    act_num += 1;
    for (unsigned int i = 0; i < DIM; i++) {
      unsigned int tmp = 0;
      if (act_num <= tmp2) {
 	if (e_[i] == 0)
	  k_[i] = (basis_->bases()[i])->DeltaLmin();
 	else
 	  k_[i] = (basis_->bases()[i])->Nablamin();
      }
      else {// act_num > tmp
	if (tmp2 == 1)
	  tmp = act_num-1;
	else if ((tmp2 != 1) && ((act_num % tmp2) != 0))
	  tmp = act_num / tmp2;
	else if ( (act_num % tmp2) == 0 ) 
	  tmp = act_num / tmp2 - 1;

 	if (e_[i] == 0)
 	  k_[i] = tmp + (basis_->bases()[i])->DeltaLmin();
 	else
 	  k_[i] = tmp + (basis_->bases()[i])->Nablamin();
	
	if ( (act_num % tmp2) != 0 )
	  act_num = act_num % tmp2;
	else
	  act_num = tmp2;
      }
      if ((i+1) < DIM) {
	if (e_[i+1] == 0) {
	  tmp2 /= (basis_->bases()[i+1])->Deltasize(j_);
	}
	else {
	  tmp2 /= (basis_->bases()[i+1])->Nablasize(j_);
	}
      }
    }
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
  const unsigned int 
  CubeIndex<IBASIS,DIM,CUBEBASIS>::number() const
  {
    unsigned int result = 0;
    bool gen = 1;

    // check if wavelet is a generator
    for (unsigned int i = 0; i < DIM; i++) {
      if (e_[i] == 1) {
	gen = 0;
	break;
      }
    }

    // determine how many wavelets there are on all the levels
    // below the level of this index
    if (! gen) {
      result = 1;
      for (unsigned int i = 0; i < DIM; i++)
	result *= (basis_->bases()[i])->Deltasize(j_);
    }
    
    // now determine how many wavelets there are on the same level
    // that belong to another type less than the type of this wavelet,
    // add the result to res afterwards
    if (! gen) {
      MultiIndex<int,DIM> type;
      type[DIM-1] = 1;
      bool exit = 0;
      // loop over all ''smaller'' types
      while (type < e_) {
	unsigned int tmp = 1;
	for (unsigned int i = 0; i < DIM; i++) {
	  if (type[i] == 0)
	    tmp *= (basis_->bases()[i])->Deltasize(j_);
	  else
	    tmp *= (basis_->bases()[i])->Nablasize(j_);
	}
	result += tmp;
	// determine next type
	for (unsigned int i = DIM-1; i >= 0; i--) {
	  if ( type[i] == 1 ) {
	    type[i] = 0;
	    exit = (i == 0);
	  }
	  else {
	    type[i]++;
	    break;
	  }
	}
	if (exit)
	  break;
      }//end while
    }// end if

    // count the wavelets that belong to the same level and to the same type
    // whose indices are smaller than this index,
    // add the result to res
    for (unsigned int i = 0; i < DIM; i++) {
      unsigned int tmp = 1;
      if (e_[i] == 0) {
	if (k_[i] == (basis_->bases()[i])->DeltaLmin())
	  continue;
      }
      else
	if (k_[i] == (basis_->bases()[i])->Nablamin())
	  continue;
      
      if (e_[i] == 0)
	tmp *= k_[i]-(basis_->bases()[i])->DeltaLmin();
      else
	tmp *= k_[i]-(basis_->bases()[i])->Nablamin();
      
      for (unsigned int l = i+1; l < DIM; l++) {
	if (e_[l] == 0)
	  tmp *= (basis_->bases()[i])->Deltasize(j_);
	else
	  tmp *= (basis_->bases()[i])->Nablasize(j_);
      }
      result += tmp;
    }
    return result;
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
