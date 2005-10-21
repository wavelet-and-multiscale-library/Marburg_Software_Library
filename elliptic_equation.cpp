// implementation for elliptic_equation.h


namespace FrameTL
{

  template <class IBASIS, unsigned int DIM>
  EllipticEquation<IBASIS,DIM>::EllipticEquation(const EllipticBVP<DIM>& ell_bvp,
						 const AggregatedFrame<IBASIS,DIM>* frame)
    : ell_bvp_(ell_bvp), frame_(frame)
  {
  }

  template <class IBASIS, unsigned int DIM>
  double
  EllipticEquation<IBASIS,DIM>::a_same_patches(const typename AggregatedFrame<IBASIS,DIM>::Index& lambda,
					       const typename AggregatedFrame<IBASIS,DIM>::Index& nu,
					       const unsigned int p) const
  {
    return 0.0;
  }

  template <class IBASIS, unsigned int DIM>
  double
  EllipticEquation<IBASIS,DIM>::a_different_patches(const typename AggregatedFrame<IBASIS,DIM>::Index& lambda,
						    const typename AggregatedFrame<IBASIS,DIM>::Index& nu,
						    const unsigned int p) const
  {
    return 0.0;
  }


  template <class IBASIS, unsigned int DIM>
  double
  EllipticEquation<IBASIS,DIM>::a(const typename AggregatedFrame<IBASIS,DIM>::Index& lambda,
				  const typename AggregatedFrame<IBASIS,DIM>::Index& nu,
				  const unsigned int p) const
  {
    //the cases lambda.p() == nu.p() and lambda.p() != nu.p() have
    //to be treated seperately
    return lambda.p() == nu.p() ? a_same_patches(lambda, nu, p) : a_different_patches(lambda, nu, p);
  }

}
