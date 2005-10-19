// implementation for frame_index.h

#include <math.h>
#include <cmath>

 namespace FrameTL
 {

   template <class IBASIS, unsigned int DIM>
   FrameIndex<IBASIS,DIM>::FrameIndex()
   {
   }

   template <class IBASIS, unsigned int DIM>
   FrameIndex<IBASIS,DIM>::FrameIndex(const FrameIndex& ind)
     : p_(ind.p()), cbI_(ind.get_CubeIndex()), frame_(ind.get_frame())
   {
   }

   template <class IBASIS, unsigned int DIM>
   FrameIndex<IBASIS,DIM>::FrameIndex(const AggregatedFrame<IBASIS,DIM_d,DIM_m>* frame,
				      const CubeIndex<IBASIS,DIM>& c,
				      const int patch,
				      const unsigned int num_patches_)
     : p_(patch), cbI_(c), num_patches(num_patches_),  frame_(frame)
   {
   }

   template <class IBASIS, unsigned int DIM>
   FrameIndex<IBASIS,DIM>::FrameIndex(const AggregatedFrame<IBASIS,DIM_d,DIM_m>* frame,
				      const int j,
				      const type_type& e,
				      const translation_type& k,
				      const unsigned int patch,
				      CubeBasis<IBASIS,DIM_d>* basis,
				      const unsigned int num_patches_)

     : p_(patch), num_patches(num_patches_), frame_(frame)
   {
     cbI_(j, e, k, basis);
   }

   template <class IBASIS, unsigned int DIM>
   inline
   bool
   FrameIndex<IBASIS,DIM>::operator == (const FrameIndex& lambda) const
   {
     return (cbI_ == lambda.get_CubeIndex()) && (p_ == lambda.p());
   }

   template <class IBASIS, unsigned int DIM>
   inline
   FrameIndex<IBASIS,DIM>&
   FrameIndex<IBASIS,DIM>::operator ++ ()
   {
//      if (p_ < num_patches)
//        p_++;
//      //last patch
//      if (p_ == num_patches-1)
//        {
// 	 p_ = 0;
// 	 ++cbI_;
//        }

   }

   template <class IBASIS, unsigned int DIM>
   inline
   bool
   FrameIndex<IBASIS,DIM>::operator < (const FrameIndex& lambda) const
   {
     return (cbI_.j() < lambda.j() ||
      (
       cbI_.j() == lambda.j() &&
       (
	cbI_.e() < lambda.e() ||
	(
	 cbI_.e() == lambda.e() &&
	 (
	  p_ < lambda.p()) ||
	 (
	  p_ == lambda.p() && cbI_.k() < lambda.k()
	  )
	 )
	)
       )
      );
   }
   

   template <class IBASIS, unsigned int DIM>
   inline
   std::ostream&
   operator << (std::ostream& os, const FrameIndex<IBASIS,DIM>& ind)
   {
     os << ind.get_CubeIndex() << endl
	<< "on patch " << ind.p();
     return os;
  }



 }
