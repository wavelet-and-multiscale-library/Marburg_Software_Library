// implementation for index1D.h

namespace FrameTL
{
  template <class IBASIS>
  Index1D<IBASIS>::Index1D(const typename IBASIS::Index& ind,
			   const unsigned int p, const unsigned int dir,
			   const unsigned int der)
    : ind_(ind), p_(p), dir_(dir), der_(der)
  {
  }

  template <class IBASIS>
  bool
  Index1D<IBASIS>::operator < (const Index1D<IBASIS>& lambda) const
  {
    return ind_ < lambda.index() ||
      (
       ind_ == lambda.index() &&
       (
	p_ < lambda.p() ||
	(
	 p_ == lambda.p() &&
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
  
  template <class IBASIS>
  bool Index1D<IBASIS>::operator == (const Index1D<IBASIS>& lambda) const
  {
    return (ind_ == lambda.index()) && (p_ == lambda.p()) && (der_ == lambda.derivative());
  };

  template <class IBASIS>
  bool Index1D<IBASIS>::operator != (const Index1D<IBASIS>& lambda) const
  {
    return !(*this == lambda);
  };
  
  template <class IBASIS>
  bool Index1D<IBASIS>::operator <= (const Index1D<IBASIS>& lambda) const
  {
    return (*this < lambda) || (*this == lambda);
 
  };
}

