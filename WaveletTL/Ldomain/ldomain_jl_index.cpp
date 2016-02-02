// implementation for ldomain_jl_index.h

namespace WaveletTL
{
  LDomainJLIndex::LDomainJLIndex()
    : j_(1), k_(-1,-1)
  {
    // e and c are zero by default
  }

  LDomainJLIndex::LDomainJLIndex(const int j,
				 const type_type& e,
				 const component_type& c,
				 const translation_type& k)
    : j_(j), e_(e), c_(c), k_(k)
  {
  }

  LDomainJLIndex::LDomainJLIndex(const LDomainJLIndex& lambda)
    : j_(lambda.j_), e_(lambda.e_), c_(lambda.c_), k_(lambda.k_)
  {
  }

  LDomainJLIndex&
  LDomainJLIndex::operator = (const LDomainJLIndex& lambda)
  {
    j_ = lambda.j();
    e_ = lambda.e();
    c_ = lambda.c();
    k_ = lambda.k();
    return *this;
  }

  bool
  LDomainJLIndex::operator == (const LDomainJLIndex& lambda) const
  {
    return (j_ == lambda.j() &&
	    e_ == lambda.e() &&
	    c_ == lambda.c() &&
	    k_ == lambda.k());
  }
  
  bool
  LDomainJLIndex::operator < (const LDomainJLIndex& lambda) const
  {
    // standard lexicographic order on (j,e,p,k),
    // we assume that e and k are already lexicographically ordered (cf. MultiIndex)
    return (j_ < lambda.j() ||
	    (j_ == lambda.j() &&
	     (e_ < lambda.e() ||
	      (e_ == lambda.e() &&
	       (c_ < lambda.c() ||
		(c_ == lambda.c() && k_ < lambda.k()
		 )
		)
	       )
	      )
	     )
	    );
  }

  LDomainJLIndex&
  LDomainJLIndex::operator ++ ()
  {
    // try to advance k[1]
    if (k_[0] <= -1) { // left boundary (k_[0]==-(1<<j),c_[0]==1) or left half of the L
      if (k_[1] < (c_[1] == 1 ? (1<<j_) : (1<<j_)-1)) {
	k_[1]++;
	return *this;
      }
    } else {
      if (k_[0] == 0) { // central column
	if (k_[1] < (c_[1] == 1 ? (1<<j_) : (1<<j_)-1)) {
	  if (k_[1] >= 0) { // center or upper half (c_[0]==1)
	    k_[1]++;
	    return *this;
	  } else { // k_[1]<0, we may have to skip (0,0)
	    if (k_[1] < -1) {
	      k_[1]++;
	      return *this;
	    } else { // k_[1] == -1
	      if (c_[0] == 1) {
		if (c_[1] == 1) {
		  k_[1]++;
		  return *this;
		} else {
		  k_[1]+=2;
		  return *this;
		}
	      }
	    }
	  }
	}
      } else {
	if (k_[0] >= 1) { // right half of the L or right boundary (k_[0]==(1<<j),c_[0]==1)
	  if (k_[1] < (c_[1] == 1 ? 0 : -1)) {
	    k_[1]++;
	    return *this;
	  }
	}
      }
    }

    // try to advance k[0]
    if (k_[0] < (c_[0] == 1 ? (1<<j_) : (1<<j_)-1)) {
      k_[0]++; // i.e., we are no longer at the left boundary
      k_[1] = (c_[1] == 0 ? 1-(1<<j_) : -(1<<j_));
      return *this;
    }

    // try to advance c[1]
    if (c_[1] == 0) {
      // advance c[1] and set k[0],k[1] minimal
      c_[1] = 1;
      k_[0] = (c_[0] == 0 ? 1-(1<<j_) : -(1<<j_));
      k_[1] = -(1<<j_);
      return *this;
    }

    // try to advance c[0]
    if (c_[0] == 0) {
      // advance c[0] and set c[1],k[0],k[1] minimal
      c_[0] = 1;
      c_[1] = 0;
      k_[0] = -(1<<j_);
      k_[1] = 1-(1<<j_);
      return *this;
    }

    // try to advance e[1]
    if (e_[1] == 0) {
      // advance e[1] and set c[i],k[i] minimal
      e_[1] = 1;
      c_[0] = c_[1] = 0;
      k_[0] = k_[1] = 1-(1<<j_);
      return *this;
    }
    
    // try to advance e[0]
    if (e_[0] == 0) {
      // advance e[0] and set e[1],c[i],k[i] minimal
      e_[0] = 1;
      e_[1] = 0;
      c_[0] = c_[1] = 0;
      k_[0] = k_[1] = 1-(1<<j_);
      return *this;
    }

    // advance j and set e[i],c[i],k[i] minimal
    j_++;
    e_[0] = 0;
    e_[1] = 1;
    c_[0] = c_[1] = 0;
    k_[0] = k_[1] = 1-(1<<j_);
    return *this;
  }
  

  LDomainJLIndex
  first_generator(const int j)
  {
    assert(j >= 1);
    const LDomainJLIndex::type_type e;      // e=(0,0)
    const LDomainJLIndex::component_type c; // c=(0,0)
    const LDomainJLIndex::translation_type k(1-(1<<j), 1-(1<<j));
    return LDomainJLIndex(j, e, c, k);
  }

  LDomainJLIndex
  last_generator(const int j)
  {
    assert(j >= 1);
    const LDomainJLIndex::type_type e; // e=(0,0)
    const LDomainJLIndex::component_type c(1,1);
    const LDomainJLIndex::translation_type k(1<<j,0); // generator sits at upper right corner of the "forefoot"
    return LDomainJLIndex(j, e, c, k);
  }
  
  LDomainJLIndex
  first_wavelet(const int j)
  {
    assert(j >= 1);
    const LDomainJLIndex::type_type e(0,1);
    const LDomainJLIndex::component_type c(0,0);
    const LDomainJLIndex::translation_type k(1-(1<<j),1-(1<<j));
    return LDomainJLIndex(j, e, c, k);
  }

  LDomainJLIndex
  first_wavelet(const int j,
 		const LDomainJLIndex::type_type& e)
  {
    assert(j >= 1);
    const LDomainJLIndex::component_type c(0,0);
    const LDomainJLIndex::translation_type k(1-(1<<j),1-(1<<j));
    return LDomainJLIndex(j, e, c, k);
  }
  

  LDomainJLIndex
  last_wavelet(const int j)
  {
    assert(j >= 1);
    const LDomainJLIndex::type_type e(1,1);
    const LDomainJLIndex::component_type c(1,1);
    const LDomainJLIndex::translation_type k(1<<j,0); // wavelet sits at upper right corner of the "forefoot"
    return LDomainJLIndex(j, e, c, k);
  }

  bool index_is_valid(const int j, const int e0, const int e1,
		      const int c0, const int c1, const int k0, const int k1)
  {
    // return a result as fast as possible
    if (k0 >= 1-(1<<j) && k0 <= -1 && k1 >= 1-(1<<j) && k1 <= (1<<j)-1) return true; // in "shaft"
    if (k0 >= 0 && k0 <= (1<<j)-1 && k1 >= 1-(1<<j) && k1 <= -1) return true;        // in "forefoot"
    if (k0 == -(1<<j) && c0 == 1) {
      if (k1 >= 1-(1<<j) && k1 <= (1<<j)-1) return true; // left edge of "shaft"
      if (c1 == 1) {
	if (k1 == -(1<<j)) return true;         // lower left corner
	if (k1 == (1<<j)) return true;          // upper left corner
      }
    }
    if (k1 == -(1<<j) && c1 == 1) {
      if (k0 >= 1-(1<<j) && k0 <= (1<<j)-1) return true; // bottom
      if (k0 == (1<<j) && c0 == 1) return true;          // lower right corner
    }
    if (k0 == (1<<j) && c0 == 1) {
      if (k1 >= 1-(1<<j) && k1 <= -1) return true;       // right edge of "forefoot"
      if (k1 == 0 && c1 == 1) return true;               // upper right corner of "forefoot"
    }
    if (k1 == 0 && c1 == 1) {
      if (k0 >= 1 && k0 <= (1<<j)-1) return true;        // upper edge of "forefoot"
      if (k0 == 0 && c0 == 1) return true;               // origin
    }
    if (k0 == 0 && c0 == 1) {
      if (k1 >= 1 && k1 <= (1<<j)-1) return true;        // right edge of "shaft"
      if (k1 == (1<<j) && c1 == 1) return true;          // upper right corner of "shaft"
    }
    if (k1 == (1<<j) && k0 >= 1-(1<<j) && k0 <= -1 && c1 == 1) return true; // top

    return false;
  }

}
