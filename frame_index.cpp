// implementation for frame_index.h

#include <math.h>
#include <cmath>
#include <interval/i_index.h>

 namespace FrameTL
 {

   template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
   FrameIndex<IBASIS, DIM_d, DIM_m>::FrameIndex(const AggregatedFrame<IBASIS,DIM_d,DIM_m>* frame)
     : frame_(frame)
   {  
     if (frame_ == 0) {
       j_ = 0; // invalid (e and k are initialized by zero automatically)
       p_ = 0;
       num_ = -1;
     } else {
       j_ = frame_->j0(); // coarsest level;
       // e_ is zero by default: generator
       p_ = 0;
       for (unsigned int i = 0; i < DIM_d; i++)
	 k_[i] = WaveletTL::first_generator<IBASIS>(frame_->bases()[0]->bases()[i], j_).k();
       
       num_ = -1;

     }
   }

   template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
   FrameIndex<IBASIS, DIM_d, DIM_m>::FrameIndex(const FrameIndex& ind)
     : frame_(ind.frame()), j_(ind.j()), e_(ind.e()), p_(ind.p()), k_(ind.k()), num_(ind.number())
   {
   }

   template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
   FrameIndex<IBASIS, DIM_d, DIM_m>::FrameIndex(const AggregatedFrame<IBASIS,DIM_d,DIM_m>* frame,
						const CubeIndex<IBASIS,DIM_d>& c,
						const int p)
     : frame_(frame), j_(c.j()), e_(c.e()), p_(p), k_(c.k()), num_(-1)
   {
   }

   template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
   FrameIndex<IBASIS, DIM_d, DIM_m>::FrameIndex(const int num,
						const AggregatedFrame<IBASIS,DIM_d,DIM_m>* frame)
     : frame_(frame)
   {
     num_ = num; 

     unsigned int n_p = frame_->n_p();

     // to be decreased successively
     int act_num = num_;

     int j = frame_->j0();

     int tmp2 = 0;

     bool gen = 0;

     // determine the level of the index
     while (true) {
       int tmp = 0;
       for (unsigned int p = 0; p < n_p; p++) {
	 int tmp3 = 1;
	 for (unsigned int i = 0; i < DIM_d; i++)
	   tmp3 *= ((frame_->bases()[p])->bases())[i]->Deltasize(j);
	 tmp += tmp3;
       }
       if (tmp > num)
	 break;
       tmp2 = tmp;
       j++;
     }
     if (j == frame_->j0()) {
       j_ = j;
       gen = 1;
     }
     else
       j_ = j-1;

     act_num -= tmp2;

     tmp2 = 0;

     // determine the type of the index 
     int tmp2_old = 0;
     if (gen)
       e_ = type_type();
     else {
       MultiIndex<int,DIM_d> type;
       type[DIM_d-1] = 1;
       bool exit = 0;
       // loop over types
       while (true) {
	 int tmp = 0;
	 for (unsigned int p = 0; p < n_p; p++) {
	   int tmp3 = 1;
	   for (unsigned int i = 0; i < DIM_d; i++) { 
	     if (type[i] == 0)
	       tmp3 *= ((frame_->bases()[p])->bases())[i]->Deltasize(j_);
	     else
	       tmp3 *= ((frame_->bases()[p])->bases())[i]->Nablasize(j_);
	   }
	   tmp += tmp3;
	 }
	 tmp2_old = tmp2;
	 tmp2 += tmp;
	
	 // right type found
	 if (tmp2 > act_num)
	   break;

	 // determine next type
	 for (unsigned int i = DIM_d-1; i >= 0; i--) {
	   if ( type[i] == 1 ) {
	     type[i] = 0;
	     exit = (i == 0);
	   }
	   else {
	     type[i]++;
	     break;
	   }
	 }
	 if (exit)
	   break;
       }//end while
       e_ = type;
     }// end else

     act_num -= tmp2_old;

     tmp2     = 0;
     tmp2_old = 0;

     // determine the patchnumber of this index
     unsigned int p = 0;
     for ( ; p < n_p; p++) {
       int tmp = 1;
       for (unsigned int i = 0; i < DIM_d; i++) { 
	 if (e_[i] == 0)
	   tmp *= ((frame_->bases()[p])->bases())[i]->Deltasize(j_);
	 else
	   tmp *= ((frame_->bases()[p])->bases())[i]->Nablasize(j_);
       }
       
       tmp2_old = tmp2;
       tmp2 += tmp;
       if (tmp2 > act_num)
	 break;
     }

     p_ = p;
     act_num -= tmp2_old;
     
     tmp2 = 1;

     // determine the position of the index
     for (unsigned int i = DIM_d-1; i > 0; i--) {
       if (e_[i] == 0) {
	 tmp2 *= ((frame_->bases()[p_])->bases())[i]->Deltasize(j_);
       }
       else {
	 tmp2 *= ((frame_->bases()[p_])->bases())[i]->Nablasize(j_);
       }
     }
     
     act_num += 1;
     for (unsigned int i = 0; i < DIM_d; i++) {
       int tmp = 0;
       if (act_num <= tmp2) {
	 if (e_[i] == 0)
	   k_[i] = ((frame_->bases()[p_])->bases())[i]->DeltaLmin();
	 else
	   k_[i] = ((frame_->bases()[p_])->bases())[i]->Nablamin();
       }
       else {// act_num > tmp
	 if (tmp2 == 1)
	   tmp = act_num-1;
	 else if ((tmp2 != 1) && ((act_num % tmp2) != 0))
	   tmp = act_num / tmp2;
	 else if ( (act_num % tmp2) == 0 ) 
	   tmp = act_num / tmp2 - 1;

	 if (e_[i] == 0)
	   k_[i] = tmp + ((frame_->bases()[p_])->bases())[i]->DeltaLmin();
	 else
	   k_[i] = tmp + ((frame_->bases()[p_])->bases())[i]->Nablamin();
	
	 if ( (act_num % tmp2) != 0 )
	   act_num = act_num % tmp2;
	 else
	   act_num = tmp2;
       }
       if ((i+1) < DIM_d) {
	 if (e_[i+1] == 0) {
	   tmp2 /= ((frame_->bases()[p_])->bases())[i]->Deltasize(j_);
	 }
	 else {
	   tmp2 /= ((frame_->bases()[p_])->bases())[i]->Nablasize(j_);
	 }
       }
     }
   }

   template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
   FrameIndex<IBASIS, DIM_d, DIM_m>::FrameIndex(const AggregatedFrame<IBASIS,DIM_d,DIM_m>* frame,
				      const int j,
				      const type_type& e,
				      const unsigned int p,
				      const translation_type& k)

     : frame_(frame), j_(j), e_(e), p_(p), k_(k), num_(-1)
   {
   }

   template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
   inline
   bool
   FrameIndex<IBASIS, DIM_d, DIM_m>::operator == (const FrameIndex& lambda) const
   {
     return (j_ == lambda.j() &&
	     e_ == lambda.e() &&
	     p_ == lambda.p() &&
	     k_ == lambda.k());
   }

   template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
   inline
   FrameIndex<IBASIS, DIM_d, DIM_m>&
   FrameIndex<IBASIS, DIM_d, DIM_m>::operator ++ ()
   {

     if (num_ > -1)
       ++num_;
     
     const unsigned int num_patches = (frame_->bases()).size();
     assert(0 <= p_ && p_ < num_patches);

     //bool eplusplus = false;
     bool pplusplus = false;
     for (int i = DIM_d-1; i >= 0; i--) {
       const int last_index = (e_[i] == 0
			       ? ((frame_->bases()[p_])->bases())[i]->DeltaRmax(j_)
			       : ((frame_->bases()[p_])->bases())[i]->Nablamax(j_));
       if (k_[i] == last_index) {
	 pplusplus = (i == 0);
	 if (pplusplus) {
	   k_[i] = (e_[i] == 0
		    ? ((frame_->bases()[(p_ < num_patches-1) ? p_+1 : 0])->bases())[i]->DeltaLmin()
		    : ((frame_->bases()[(p_ < num_patches-1) ? p_+1 : 0])->bases())[i]->Nablamin());
	 }
	 else 
	   k_[i] = (e_[i] == 0
		    ? ((frame_->bases()[p_])->bases())[i]->DeltaLmin()
		    : ((frame_->bases()[p_])->bases())[i]->Nablamin());
	 //eplusplus = (i == 0);
	 //pplusplus = (i == 0);
       } else {
	 ++k_[i];
	 break;
       }
     }

     bool eplusplus = false;
     if (pplusplus) {
       eplusplus = !(p_ < num_patches-1);
       p_ = (p_ < num_patches-1) ? p_+1 : 0;
     }
     else return *this;

     bool jplusplus = false;
     if (eplusplus) {
       for (int i = DIM_d-1; i >= 0; i--) {
	 if (e_[i] == 1) {
	   e_[i] = 0;
	   jplusplus = (i == 0);
	 } else {
	   ++e_[i];
	   break;
	 }
       }
       
       if (!jplusplus) // then choose lowest translation index k=k(j,e)
	 for (unsigned int i = 0; i < DIM_d; i++)
	   k_[i] = (e_[i] == 0
		    ? ((frame_->bases()[p_])->bases())[i]->DeltaLmin()
		    : ((frame_->bases()[p_])->bases())[i]->Nablamin());
     }
     
     if (jplusplus) {
       ++j_;
       // choose lowest type e=(0,...,0,1) and lowest translation index k=k(j,e)
       for (unsigned int i = 0; i < DIM_d-1; i++) {
	 e_[i] = 0;
	 k_[i] = ((frame_->bases()[p_])->bases())[i]->DeltaLmin();
       }
       e_[DIM_d-1] = 1;
       k_[DIM_d-1] = ((frame_->bases()[p_])->bases())[DIM_d-1]->Nablamin();
     }
     
     return *this;
   }

  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  void
  FrameIndex<IBASIS, DIM_d, DIM_m>::set_number()
  {
    type_type gen_type;
    assert( e_ != gen_type || j_ == frame_->j0() );

    unsigned int n_p = frame_->n_p();

    int result = 0;
    bool gen = 1;

    // check if wavelet is a generator
    for (unsigned int i = 0; i < DIM_d; i++) {
      if (e_[i] == 1) {
	gen = 0;
	break;
      }
    }

    // determine how many wavelets there are on all the levels
    // below the level of this index
    if (! gen) {
      result = 0;
      for (unsigned int p = 0; p < n_p; p++) {
	int tmp = 1;
	for (unsigned int i = 0; i < DIM_d; i++)
	  tmp *= ((frame_->bases()[p])->bases())[i]->Deltasize(j_);
	result += tmp;
      }
    }
   
    // now determine how many wavelets there are on the same level
    // that belong to another type less than the type of this wavelet,
    // add the result to res afterwards
    if (! gen) {
      MultiIndex<int,DIM_d> type;
      type[DIM_d-1] = 1;
      bool exit = 0;
      // loop over all ''smaller'' types
      while (type < e_) {
	for (unsigned int p = 0; p < n_p; p++) {
	  int tmp = 1;
	  for (unsigned int i = 0; i < DIM_d; i++) {
	    if (type[i] == 0)
	      tmp *= ((frame_->bases()[p])->bases())[i]->Deltasize(j_);
	    else
	      tmp *= ((frame_->bases()[p])->bases())[i]->Nablasize(j_);
	  }
	  result += tmp;
	}

	// determine next type
	for (unsigned int i = DIM_d-1; i >= 0; i--) {
	  if ( type[i] == 1 ) {
	    type[i] = 0;
	    exit = (i == 0);
	  }
	  else {
	    type[i]++;
	    break;
	  }
	}
	if (exit)
	  break;
      }//end while
    }// end if


    // now we are on the right level and on the right type.
    // next we have to go to right patch
    for (unsigned int p = 0; p < p_ ; p++) {
      int tmp = 1;
      for (unsigned int i = 0; i < DIM_d; i++) {
	if (e_[i] == 0)
	  tmp *= ((frame_->bases()[p])->bases())[i]->Deltasize(j_);
	else
	  tmp *= ((frame_->bases()[p])->bases())[i]->Nablasize(j_);
      }
      result += tmp;
    }

    // count the wavelets that belong to the same level and to the same type
    // on the same patch, whose indices are smaller than this index,
    // add the result to res
    for (unsigned int i = 0; i < DIM_d; i++) {
      int tmp = 1;
      if (e_[i] == 0) {
	if (k_[i] == ((frame_->bases()[p_])->bases())[i]->DeltaLmin())
	  continue;
      }
      else
	if (k_[i] == ((frame_->bases()[p_])->bases())[i]->Nablamin())
	  continue;
      
      if (e_[i] == 0) {
	tmp *= k_[i]-((frame_->bases()[p_])->bases())[i]->DeltaLmin();
      }
      else
	tmp *= k_[i]-((frame_->bases()[p_])->bases())[i]->Nablamin();
      
      for (unsigned int l = i+1; l < DIM_d; l++) {
	if (e_[l] == 0)
	  tmp *= ((frame_->bases()[p_])->bases())[i]->Deltasize(j_);
	else
	  tmp *= ((frame_->bases()[p_])->bases())[i]->Nablasize(j_);
      }
      result += tmp;
    }

    num_ = result;
  }

   template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
   FrameIndex<IBASIS, DIM_d, DIM_m>&
   FrameIndex<IBASIS, DIM_d, DIM_m>::operator = (const FrameIndex& lambda)
   {
     j_ = lambda.j();
     e_ = lambda.e();
     p_ = lambda.p();
     k_ = lambda.k();
     frame_ = lambda.frame();
     num_ = lambda.number();

     return *this;
   }


   template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
   inline
   bool
   FrameIndex<IBASIS, DIM_d, DIM_m>::operator < (const FrameIndex& lambda) const
   {

     return j_ < lambda.j() ||
      (
       j_ == lambda.j() &&
       (
	e_ < lambda.e() ||
	(
	 e_ == lambda.e() &&
	 (
	  p_ < lambda.p() ||
	  (
	   p_ == lambda.p() && k_ < lambda.k()
	  )
	  )
	 )
	)
      );
   }


   template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m, class FRAME>
   FrameIndex<IBASIS,DIM_d,DIM_m>
   first_generator(const FRAME* frame, const int j)
   {
     assert(j >= frame->j0());
     
     typename FrameIndex<IBASIS,DIM_d,DIM_m>::type_type e;//== 0
     typename FrameIndex<IBASIS,DIM_d,DIM_m>::translation_type k;
     for (unsigned int i = 0; i < DIM_d; i++)
       k[i] = WaveletTL::first_generator<IBASIS>(frame->bases()[0]->bases()[i], j).k();
     
     return FrameIndex<IBASIS,DIM_d,DIM_m>(frame, j, e, 0, k);
   }
   
   template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m, class FRAME>
   FrameIndex<IBASIS,DIM_d,DIM_m>
   last_generator(const FRAME* frame, const int j)
   {
     assert(j >= frame->j0());

     typename FrameIndex<IBASIS,DIM_d,DIM_m>::type_type e;//== 0
     typename FrameIndex<IBASIS,DIM_d,DIM_m>::translation_type k;
     for (unsigned int i = 0; i < DIM_d; i++)
       k[i] = WaveletTL::last_generator<IBASIS>(frame->bases()[frame->n_p()-1]->bases()[i], j).k();

     return FrameIndex<IBASIS,DIM_d,DIM_m>(frame, j, e, frame->bases().size()-1, k); 
   }

   template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m, class FRAME>
   FrameIndex<IBASIS,DIM_d,DIM_m>
   first_wavelet(const FRAME* frame, const int j)
   {
     assert(j >= frame->j0());

     typename FrameIndex<IBASIS,DIM_d,DIM_m>::type_type e;//== 0
     typename FrameIndex<IBASIS,DIM_d,DIM_m>::translation_type k; 
     for (unsigned int i = 0; i < DIM_d-1; i++)
       k[i] = WaveletTL::first_generator<IBASIS>(frame->bases()[0]->bases()[i], j).k();

     k[DIM_d-1] = WaveletTL::first_wavelet<IBASIS>(frame->bases()[0]->bases()[DIM_d-1], j).k();
     e[DIM_d-1] = 1;

     return FrameIndex<IBASIS,DIM_d,DIM_m>(frame, j, e, 0, k); 
   }

   template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m, class FRAME>
   FrameIndex<IBASIS,DIM_d,DIM_m>
   last_wavelet(const FRAME* frame, const int j)
   {
     assert(j >= frame->j0());
     
     typename FrameIndex<IBASIS,DIM_d,DIM_m>::type_type e;//== 0
     typename FrameIndex<IBASIS,DIM_d,DIM_m>::translation_type k; 
     for (unsigned int i = 0; i < DIM_d; i++) {
       k[i] = WaveletTL::last_wavelet<IBASIS>(frame->bases()[frame->bases().size()-1]->bases()[i], j).k();
       e[i] = 1;
     }
          
     return FrameIndex<IBASIS,DIM_d,DIM_m>(frame, j, e, frame->bases().size()-1, k); 
   }


   template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m, class FRAME>
   int
   first_generator_num(const FRAME* frame)
   {
     FrameIndex<IBASIS,DIM_d,DIM_m> ind(first_generator<IBASIS,DIM_d,DIM_m>(frame, frame->j0()));
     ind.set_number();
     return ind.number();
   }
   
   template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m, class FRAME>
   int
   last_generator_num(const FRAME* frame)
   {
     FrameIndex<IBASIS,DIM_d,DIM_m> ind(last_generator<IBASIS,DIM_d,DIM_m>(frame, frame->j0()));
     ind.set_number();
     return ind.number();
   }

   template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m, class FRAME>
   int
   first_wavelet_num(const FRAME* frame, const int j)
   {
     FrameIndex<IBASIS,DIM_d,DIM_m> ind(first_wavelet<IBASIS,DIM_d,DIM_m>(frame, j));
     ind.set_number();
     return ind.number();
   }

   template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m, class FRAME>
   int
   last_wavelet_num(const FRAME* frame, const int j)
   {
     FrameIndex<IBASIS,DIM_d,DIM_m> ind(last_wavelet<IBASIS,DIM_d,DIM_m>(frame, j));
     ind.set_number();
     return ind.number();
   }

   
   template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
   inline
   std::ostream&
   operator << (std::ostream& os, const FrameIndex<IBASIS, DIM_d, DIM_m>& lambda)
   {
     os << "("
	<< lambda.j()
	<< ","
	<< lambda.e()
	<< ","
	<< lambda.p()
	<< ","
	<< lambda.k()
	<< ")" << " number = " << lambda.number(); 
    return os;

   }



 }
