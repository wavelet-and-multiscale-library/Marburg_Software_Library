// implementation for indexq1D.h

namespace WaveletTL
{
  template <class IFRAME>
  IndexQ1D<IFRAME>::IndexQ1D(const typename IFRAME::Index& ind,
			   const unsigned int patch, const unsigned int dir,
			   const unsigned int der)
    : ind_(ind), patch_(patch), dir_(dir), der_(der)
  {
      num_=ind.number();
  }

  template <class IFRAME>
  bool
  IndexQ1D<IFRAME>::operator < (const IndexQ1D<IFRAME>& lambda) const
  {
    return num_ < lambda.number() ||
      (
       num_ == lambda.number() &&      
       (
	patch_ < lambda.patch() ||
	(
	 patch_ == lambda.patch() &&
	 (
	  dir_ < lambda.direction() ||
	  (
	   dir_ == lambda.direction() && der_ < lambda.derivative()   
	   )
	  )
	 )
	)
       );
  };
  
  template <class IFRAME>
  bool IndexQ1D<IFRAME>::operator == (const IndexQ1D<IFRAME>& lambda) const
  {
    return (num_ == lambda.number()) && (patch_ == lambda.patch()) && (der_ == lambda.derivative());
  };

  template <class IFRAME>
  bool IndexQ1D<IFRAME>::operator != (const IndexQ1D<IFRAME>& lambda) const
  {
    return !(*this == lambda);
  };
  
  template <class IFRAME>
  bool IndexQ1D<IFRAME>::operator <= (const IndexQ1D<IFRAME>& lambda) const
  {
    return (*this < lambda) || (*this == lambda);
 
  };
}

