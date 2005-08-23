// implementation for i_index.h

namespace WaveletTL
{
  template <class IBASIS>
  IIndex<IBASIS>::IIndex()
  {
    k_ = IBASIS::DeltaLmin(); // leftmost
    e_ = 0;                   // generator
    j_ = IBASIS::j0();        // on the coarsest level
  }
  
  template <class IBASIS>
  IIndex<IBASIS>::IIndex(const IIndex<IBASIS>& lambda)
  {
    j_ = lambda.j();
    e_ = lambda.e();
    k_ = lambda.k();
  }

  template <class IBASIS>
  IIndex<IBASIS>::IIndex(const int j, const int e, const int k)
  {
    j_ = j;
    e_ = e;
    k_ = k;
  }

  template <class IBASIS>
  IIndex<IBASIS>&
  IIndex<IBASIS>::operator = (const IIndex<IBASIS>& lambda)
  {
    j_ = lambda.j();
    e_ = lambda.e();
    k_ = lambda.k();

    return *this;
  }

  template <class IBASIS>
  bool
  IIndex<IBASIS>::operator == (const IIndex<IBASIS>& lambda) const
  {
    return (j_ == lambda.j() &&
	    e_ == lambda.e() &&
	    k_ == lambda.k());
  }

  template <class IBASIS>
  bool
  IIndex<IBASIS>::operator < (const IIndex<IBASIS>& lambda) const
  {
    return (j_ < lambda.j() ||
	    (j_ == lambda.j() && e_ < lambda.e()) ||
	    (j_ == lambda.j() && e_ == lambda.e() && k_ < lambda.k()));
  }

  template <class IBASIS>
  IIndex<IBASIS>&
  IIndex<IBASIS>::operator ++ ()
  {
    switch (e_) {
    case 0:
      if (k_ == IBASIS::DeltaRmax(j_)) {
	e_ = 1;
	k_ = 0;
      }
      else
	k_++;
      break;
    case 1:
      if (k_ == (1<<j_)-1) {
	j_++;
	k_ = 0;
      }
      else
	k_++;
      break;
    default:
      break;
    }
    
    return *this;
  }

  template <class IBASIS>
  IIndex<IBASIS> first_generator(const int j)
  {
    assert(j >= IBASIS::j0());
    return IIndex<IBASIS>(j, 0, IBASIS::DeltaLmin());
  }
  
  template <class IBASIS>
  IIndex<IBASIS> last_generator(const int j)
  {
    assert(j >= IBASIS::j0());
    return IIndex<IBASIS>(j, 0, IBASIS::DeltaRmax(j));
  }

  template <class IBASIS>
  IIndex<IBASIS> first_wavelet(const int j)
  {
    assert(j >= IBASIS::j0());
    return IIndex<IBASIS>(j, 1, IBASIS::Nablamin());
  }
  
  template <class IBASIS>
  IIndex<IBASIS> last_wavelet(const int j)
  {
    assert(j >= IBASIS::j0());
    return IIndex<IBASIS>(j, 1, IBASIS::Nablamax(j));
  }

}
