// implementation for ring_index.h

#include <interval/periodic.h>
#include <interval/spline_basis.h>

namespace WaveletTL
{
  template <int d, int dt, int s0, int s1>
  RingIndex<d,dt,s0,s1>::RingIndex()
    : j_(RingBasis<d,dt,s0,s1>::j0()),
      k_(PeriodicBasis<CDFBasis<d,dt> >::DeltaLmin(),
	 SplineBasis<d,dt,P_construction,s0,s1,0,0,SplineBasisData_j0<d,dt,P_construction,s0,s1,0,0>::j0>::DeltaLmin())
  {
  }

  template <int d, int dt, int s0, int s1>
  RingIndex<d,dt,s0,s1>::RingIndex(const RingIndex<d,dt,s0,s1>& lambda)
    : j_(lambda.j_), e_(lambda.e_), k_(lambda.k_)
  {
  }

  template <int d, int dt, int s0, int s1>
  RingIndex<d,dt,s0,s1>::RingIndex(const RingIndex<d,dt,s0,s1>* lambda)
    : RingIndex(*lambda)
  {
  }
  
  template <int d, int dt, int s0, int s1>
  RingIndex<d,dt,s0,s1>::RingIndex(const int j,
				   const type_type& e,
				   const translation_type& k)
    : j_(j), e_(e), k_(k)
  {
  }

  template <int d, int dt, int s0, int s1>
  RingIndex<d,dt,s0,s1>&
  RingIndex<d,dt,s0,s1>::operator = (const RingIndex<d,dt,s0,s1>& lambda)
  {
    j_ = lambda.j();
    e_ = lambda.e();
    k_ = lambda.k();
    
    return *this;
  }

  template <int d, int dt, int s0, int s1>
  inline
  bool
  RingIndex<d,dt,s0,s1>::operator == (const RingIndex& lambda) const
  {
    return (j_ == lambda.j() &&
	    e_ == lambda.e() &&
	    k_ == lambda.k());
  }

  template <int d, int dt, int s0, int s1>
  inline
  bool
  RingIndex<d,dt,s0,s1>::operator < (const RingIndex& lambda) const
  {
    return (j_ < lambda.j() ||
	    (j_ == lambda.j() && (e_ < lambda.e() ||
				  (e_ == lambda.e() && k_ < lambda.k()))));
  }

  template <int d, int dt, int s0, int s1>
  RingIndex<d,dt,s0,s1>&
  RingIndex<d,dt,s0,s1>::operator ++ ()
  {
    bool eplusplus = false;

    typedef PeriodicBasis<CDFBasis<d,dt> > Basis0;
    typedef SplineBasis<d,dt,P_construction,s0,s1,0,0,SplineBasisData_j0<d,dt,P_construction,s0,s1,0,0>::j0> Basis1;

    // increment k
    int last_index = (e_[1] == 0 ? Basis1::DeltaRmax(j_) : Basis1::Nablamax(j_));
    if (k_[1] == last_index) {
      k_[1] = (e_[1] == 0 ? Basis1::DeltaLmin() : Basis1::Nablamin());
      last_index = (e_[0] == 0 ? Basis0::DeltaRmax(j_) : Basis0::Nablamax(j_));
      if (k_[0] == last_index) {
	k_[0] = (e_[0] == 0 ? Basis0::DeltaLmin() : Basis0::Nablamin());
	eplusplus = true;
      } else {
	++k_[0];
      }
    } else {
      ++k_[1];
    }
    
    // increment e, if necessary
    bool jplusplus = false;
    if (eplusplus) {
      for (int i = 1; i >= 0; i--) {
 	if (e_[i] == 1) {
 	  e_[i] = 0;
 	  jplusplus = (i == 0);
 	} else {
 	  ++e_[i];
 	  break;
 	}
      }
     
      if (!jplusplus) { // then choose lowest translation index k=k(j,e)
	k_[0] = (e_[0] == 0 ? Basis0::DeltaLmin() : Basis0::Nablamin());
	k_[1] = (e_[1] == 0 ? Basis1::DeltaLmin() : Basis1::Nablamin());
      }
    }

    // increment j, if necessary
    if (jplusplus) {
      ++j_;
      // choose lowest type e=(0,1) and lowest translation index k=k(j,e)
      e_[0] = 0;
      k_[0] = Basis0::DeltaLmin();
      e_[1] = 1;
      k_[1] = Basis1::Nablamin();
    }
    
    return *this;
  }

  template <int d, int dt, int s0, int s1>
  const int
  RingIndex<d,dt,s0,s1>::number() const
  {
    int result = 0;
    
    const int j0 = RingBasis<d,dt,s0,s1>::j0();
    const int ecode = e_[0]*2+e_[1];

    // determine how many wavelets there are with levels j0<=j'<j
    if (j_ > j0 && ecode > 0)
      result += RingBasis<d,dt,s0,s1>::Deltasize(j_);
   
    // determine how many wavelets there are with level j and type 0<=e'<e
    switch(ecode) {
    case 1: // e=(0,1)
      if (j_ == j0)
	result +=
	  RingBasis<d,dt,s0,s1>::Basis0::Deltasize(j_)
	  * RingBasis<d,dt,s0,s1>::Basis1::Deltasize(j_);
      break;
    case 2: // e=(1,0)
      result +=
	RingBasis<d,dt,s0,s1>::Basis0::Deltasize(j_)
	* (j_ == j0
	   ? RingBasis<d,dt,s0,s1>::Basis1::Deltasize(j_+1)
	   : RingBasis<d,dt,s0,s1>::Basis1::Nablasize(j_));
      break;
    case 3: // e=(1,1)
      result += (j_ == j0
		 ? RingBasis<d,dt,s0,s1>::Basis0::Deltasize(j_)
		 * RingBasis<d,dt,s0,s1>::Basis1::Deltasize(j_+1)
		 + RingBasis<d,dt,s0,s1>::Basis0::Nablasize(j_)
		 * RingBasis<d,dt,s0,s1>::Basis1::Deltasize(j_)
		 : RingBasis<d,dt,s0,s1>::Basis0::Deltasize(j_)
		 * RingBasis<d,dt,s0,s1>::Basis1::Nablasize(j_)
		 + RingBasis<d,dt,s0,s1>::Basis0::Nablasize(j_)
		 * RingBasis<d,dt,s0,s1>::Basis1::Deltasize(j_));
      break;
    case 0: // add nothing
      break;
    }

    // determine how many wavelets there are with level j, type e and k'<k
    switch(ecode) {
    case 0: // e=(0,0)
      result +=
	(k_[0]-RingBasis<d,dt,s0,s1>::Basis0::DeltaLmin())
	* RingBasis<d,dt,s0,s1>::Basis1::Deltasize(j_)
	+ k_[1]-RingBasis<d,dt,s0,s1>::Basis1::DeltaLmin();
      break;
    case 1: // e=(0,1)
      result +=
	(k_[0]-RingBasis<d,dt,s0,s1>::Basis0::DeltaLmin())
	* RingBasis<d,dt,s0,s1>::Basis1::Nablasize(j_)
	+ k_[1]-RingBasis<d,dt,s0,s1>::Basis1::Nablamin();
      break;
    case 2:
      result +=
	(k_[0]-RingBasis<d,dt,s0,s1>::Basis0::Nablamin())
	* RingBasis<d,dt,s0,s1>::Basis1::Deltasize(j_)
	+ k_[1]-RingBasis<d,dt,s0,s1>::Basis1::DeltaLmin();	
      break;
    case 3:
      result +=
	(k_[0]-RingBasis<d,dt,s0,s1>::Basis0::Nablamin())
	* RingBasis<d,dt,s0,s1>::Basis1::Nablasize(j_)
	+ k_[1]-RingBasis<d,dt,s0,s1>::Basis1::Nablamin();
      break;
    }

    return result;
  }
 
}
