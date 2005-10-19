// implementation for aggregated_frame.h

namespace FrameTL
{
  template<class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  AggregatedFrame<IBASIS,DIM_d,DIM_m>::AggregatedFrame(const Atlas<DIM_d, DIM_m>* atlas,
						       const Array1D<FixedArray1D<int,2*DIM_d> >& bc,
						       const Array1D<FixedArray1D<int,2*DIM_d> >& bcT)
    : atlas_(atlas), delete_atlas(false), bc_(bc)
  {
    lifted_bases.resize((atlas_->charts()).size());
    for (typename Array1D<FixedArray1D<int,2*DIM_d> >::size_type i(0); 
	 i < (atlas_->charts()).size(); ++i)
      {
	lifted_bases[i] = MappedCubeBasis<IBASIS,DIM_d,DIM_m>((atlas_->charts())[i],bc[i],bcT[i]);
      }
  }

//   template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
//   AggregatedFrame<IBASIS,DIM_d,DIM_m>::~AggregatedFrame()
//   {
//     if (delete_atlas)
//       delete atlas_;
//   }

  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  inline
  const MappedCubeBasis<IBASIS, DIM_d, DIM_m>&
  AggregatedFrame<IBASIS,DIM_d,DIM_m>::get_local_basis(const unsigned int i)
  {
    assert((0 <= i) && (i < (atlas_->charts()).size()));
    return lifted_bases[i];
  }
  
  

//   template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
//   inline
//   std::ostream& operator << (std::ostream&,
// 			     const AggregatedFrame<IBASIS, DIM_d, DIM_m>&)
//   {
    
//   }


 }
