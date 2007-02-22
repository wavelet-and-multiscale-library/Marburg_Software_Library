// implementation for jl_index.h

namespace WaveletTL
{
  JLIndex::JLIndex()
    : j_(1), e_(0), c_(0), k_(1)
  {
  }

  JLIndex::JLIndex(const JLIndex& index)
    : j_(index.j_), e_(index.e_), c_(index.c_), k_(index.k_)
  {
  }

  JLIndex::JLIndex(const JLIndex* index)
    : j_(index->j_), e_(index->e_), c_(index->c_), k_(index->k_)
  {
  }



  JLIndex::JLIndex(const int j, const type_type e, const component_type c, const translation_type k)
    : j_(j), e_(e), c_(c), k_(k)
  {
  }
  
  JLIndex&
  JLIndex::operator = (const JLIndex& lambda)
  {
    j_ = lambda.j();
    e_ = lambda.e();
    c_ = lambda.c();
    k_ = lambda.k();

    return *this;
  }

  bool
  JLIndex::operator == (const JLIndex& lambda) const
  {
    return (j_ == lambda.j() &&
	    e_ == lambda.e() &&
	    c_ == lambda.c() &&
	    k_ == lambda.k());
  }
  
  bool
  JLIndex::operator < (const JLIndex& lambda) const
  {
    return (j_ < lambda.j() ||
	    (j_ == lambda.j() && (e_ < lambda.e() ||
				  (e_ == lambda.e() && (c_ < lambda.c() ||
							(c_ == lambda.c() && k_< lambda.k()))))));
  }
  
  JLIndex&
  JLIndex::operator ++ ()
  {
    if (e_==0) {
      if (c_==0) {
	if (k_==(1<<j_)-1) {
	  c_ = 1;
	  k_ = 0;
	}
	else
	  k_++;
      } else {
	// c_==1
	if (k_==1<<j_) {
	  e_ = 1;
	  c_ = 0;
	  k_ = 1;
	}
	else
	  k_++;
      }
    } else {
      // e_==1
      if (c_==0) {
	if (k_==(1<<j_)-1) {
	  c_ = 1;
	  k_ = 0;
	}
	else
	  k_++;
      } else {
	// c_==1
	if (k_==1<<j_) {
	  j_++;
	  c_=0;
	  k_=1;
	}
	else
	  k_++;
      }
    }
    
    return *this;
  }
  
  JLIndex first_generator(const int j)
  {
    assert(j >= 1);
    return JLIndex(j, 0, 0, 1);
  }
  
  JLIndex last_generator(const int j)
  {
    assert(j >= 1);
    return JLIndex(j, 0, 1, 1<<j);
  }
  
  JLIndex first_wavelet(const int j)
  {
    assert(j >= 1);
    return JLIndex(j, 1, 0, 1);
  }
  
  JLIndex last_wavelet(const int j)
  {
    assert(j >= 1);
    return JLIndex(j, 1, 1, 1<<j);
  }

}
