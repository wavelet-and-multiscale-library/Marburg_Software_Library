// implementation for aggregated_frame.h

namespace FrameTL
{
  template<class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  AggregatedFrame<IBASIS,DIM_d,DIM_m>::AggregatedFrame(const Atlas<DIM_d, DIM_m>* atlas,
						       const Array1D<FixedArray1D<int,2*DIM_d> >& bc,
						       const Array1D<FixedArray1D<int,2*DIM_d> >& bcT)
    : atlas_(atlas), bc_(bc), bcT_(bcT)
  {
    lifted_bases.resize((atlas_->charts()).size());

#if 0 //FORGET ABOUT IT
   
    //we want to make sure that in case some mapped cube bases
    //fulfill exactly the same boundary conditions, only one instance
    //of such a basis is created
    for (unsigned int  i = 0; i < (atlas_->charts()).size(); ++i)
      {
	MappedCubeBasis<IBASIS,DIM_d,DIM_m>* b = 0;
	for (typename list<IBASIS*>::const_iterator it(bases_infact.begin());
	     it != bases_infact.end(); ++it)
	  {
	    //compare boundary conditions
	    
	  }
      }
#endif

    for (unsigned int  i = 0; i < (atlas_->charts()).size(); ++i)
      lifted_bases[i] = new MappedCubeBasis<IBASIS,DIM_d,DIM_m>((atlas_->charts())[i],bc[i],bcT[i]);
    
  }

  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  AggregatedFrame<IBASIS,DIM_d,DIM_m>::~AggregatedFrame()
  {
    for (unsigned int  i = 0; i < (atlas_->charts()).size(); ++i)
      delete lifted_bases[i];          
  }

//   template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
//   inline
//   const MappedCubeBasis<IBASIS, DIM_d, DIM_m>*
//   AggregatedFrame<IBASIS,DIM_d,DIM_m>::get_local_basis(const unsigned int i)
//   {
//     assert((0 <= i) && (i < (atlas_->charts()).size()));
//     return lifted_bases[i];
//   }  
  

//   template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
//   inline
//   std::ostream& operator << (std::ostream&,
// 			     const AggregatedFrame<IBASIS, DIM_d, DIM_m>&)
//   {
    
//   }


 }
