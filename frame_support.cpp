// implementation for frame_support.h

//#include <utils/multiindex.h>

namespace FrameTL
{
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  bool
  intersect_supports(const AffineLinearMapping<DIM_d>& ch1,
		     const AffineLinearMapping<DIM_d> ch2,
		     const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
		     const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& mu)
  {
    ;
  }
  
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  bool
  intersect_supports(const LinearBezierMapping& ch1,
		     const LinearBezierMapping& ch2,
		     const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
		     const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& mu)
  {
    ;
  }

}
