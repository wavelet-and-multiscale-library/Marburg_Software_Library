
#include "ldomain_frame_index.h"
#include "ldomain_frame.h"

// implementation for ldomain_index.h

namespace WaveletTL
{
  template <class IFRAME>
  LDomainFrameIndex<IFRAME>::LDomainFrameIndex(const LDomainFrame<IFRAME>* frame)
    : frame_(frame), patch_(0)
  {
    if (frame_ == 0) {
      
      p_[0] = 0; 
      p_[1] = 0;
      j_[0] = 0; 
      j_[1] = 0;// invalid (e and k are initialized by zero automatically)
      num_ = -1;
    } else {
      j_ = frame_->j0(); // coarsest level;
      // e_ is zero by default: generator
      k_[0] = frame->frame1d().DeltaLmin();
      k_[1] = k_[0]+1;
      num_ = 0;
    }
  }

  template <class IFRAME>
  LDomainFrameIndex<IFRAME>::LDomainFrameIndex(const polynomial_type& p,
                                     const level_type& j,
				     const type_type& e,
				     const int patch,
				     const translation_type& k,
                                     const unsigned int number,
				     const LDomainFrame<IFRAME>* frame)
    : frame_(frame), p_(p),j_(j), e_(e), patch_(patch), k_(k), num_(number)
  {
  }

  template <class IFRAME>
  LDomainFrameIndex<IFRAME>::LDomainFrameIndex(const LDomainFrameIndex& lambda)
    : frame_(lambda.frame_), p_(lambda.p_), j_(lambda.j_), e_(lambda.e_), patch_(lambda.patch_), k_(lambda.k_), num_(lambda.num_)
  {
  }

  template <class IFRAME>
  LDomainFrameIndex<IFRAME>::LDomainFrameIndex(const LDomainFrameIndex* lambda)
    : frame_(lambda->frame_), p_(lambda->p_), j_(lambda->j_), e_(lambda->e_), patch_(lambda->patch_), k_(lambda->k_), num_(lambda->num_)
  {
  }


  template <class IFRAME>
  bool
  LDomainFrameIndex<IFRAME>::operator == (const LDomainFrameIndex& lambda) const
  {
    return (p_ == lambda.p() &&
            j_ == lambda.j() &&
	    e_ == lambda.e() &&
	    patch_ == lambda.patch() &&
	    k_ == lambda.k());
  }

  template <class IFRAME>
  LDomainFrameIndex<IFRAME>&
  LDomainFrameIndex<IFRAME>::operator ++ ()
  {
      num_++;
      level_type j0(frame_->j0());
    // decide whether the patch number has to be increased
    bool pplusplus = false;
    for (int i = 1; i >= 0; i--) {
      // determine the highest possible translation index into the i-th direction,
      // this will in general depend on the current patch number
      int last_index = 0;
      if (e_[i] == 0) { // generator in the i-th direction
	switch(patch_) {
	case 0:
	case 1:
	case 2:
	  last_index = frame_->frame1d().DeltaRmax(j_[i])-1;
	  break;
	case 3:
	  last_index = (i == 0
			? frame_->frame1d().DeltaRmax(j_[i])-1
			: 0); // by convention
	  break;
	case 4:
	  last_index = (i == 0
			? 0 // by convention
			: frame_->frame1d().DeltaRmax(j_[i])-1);
	  break;
	}
      } else {
	
          switch(patch_) {//If we want k more quarklets to be mirrored, we have to set last_index-=k for i==0, and last_index=k for i!=0,
                          //also you need to change the first indices in the following routines
            case 0:
            case 1:
            case 2:
              last_index = frame_->frame1d().Nablamax(j_[i]);
              break;
            case 3:
              last_index = (i == 0
                            ? frame_->frame1d().Nablamax(j_[i])
                            : 0); // by convention
              break;
            case 4:
              last_index = (i == 0
                            ? 0 // by convention
                            : frame_->frame1d().Nablamax(j_[i]));
              break;
            }
          
        
	
      }

      if (k_[i] == last_index) {
	// reset k_[i] to the lowest possible translation index
	if (e_[i] == 0) { // generator, 
	  switch(patch_) {
	  case 0:
	  case 1:
	  case 2:
	    k_[i] = frame_->frame1d().DeltaLmin()+1;
	    break;
	  case 3:
	    k_[i] = (i == 0
		     ? frame_->frame1d().DeltaLmin()+1
		     : 0); // by convention
	    break;
	  case 4:
	    k_[i] = (i == 0
		     ? 0 // by convention
		     : frame_->frame1d().DeltaLmin()+1);
	    break;
	  }
	} else { // quarklet, minimal translation index is independent from the patch number
            switch(patch_) {
            case 0:
              k_[i] = (i == 0
                       ? frame_->frame1d().Nablamin()
                       : frame_->frame1d().Nablamin()+1);//here
              break;
            case 1:
              k_[i] = frame_->frame1d().Nablamin(); 
              break;
            case 2:
              k_[i] = (i == 0
                       ? frame_->frame1d().Nablamin()+1//here, and some more. :)
                       : frame_->frame1d().Nablamin());
              break;
            case 3:
              k_[i] = (i == 0
                       ? frame_->frame1d().Nablamin()+1
                       : 0); // by convention
              break;
            case 4:
              k_[i] = (i == 0
                       ? 0 // by convention
                       : frame_->frame1d().Nablamin()+1);
              break;
            }
//	  k_[i] = frame_->frame1d().Nablamin(); // should be 0
	}
	pplusplus = (i == 0);
      } else {
	++k_[i];
	break;
      }
    }
//    cout << k_[0] << k_[1] << endl;
    bool eplusplus = false;
    if (pplusplus) {
      switch (patch_) {
      case 0:
      case 1:
	patch_++;
	break;
      case 2:
	/*if (e_[1] == 1) {
	  if (e_[0] == 0)
	    patch_ = 4; // there are no (0,1) quarklets on the interface 3
	  else
	    eplusplus = true; // there are no (1,1) quarklets on the interfaces
	} else*/ patch_ = 3;
	break;
      case 3:
	/*if (e_[0] == 1)
	  eplusplus = true; // there are no (1,*) quarklets on the interface 4
	else*/
	  patch_ = 4;
	break;
      case 4:
	eplusplus = true; // highest patch number reached
	break;
      }

      if (!eplusplus) { // then choose lowest translation index k=k(j,e,p)
	switch(patch_) { // we know that patch_>0
	case 1:
          k_[0] = (e_[0] == 0
		   ? frame_->frame1d().DeltaLmin()+1
		   : frame_->frame1d().Nablamin());
	  k_[1] = (e_[1] == 0
		   ? frame_->frame1d().DeltaLmin()+1
		   : frame_->frame1d().Nablamin());
//        cout << "2.: " << k_[0] << k_[1] << endl; 
//        cout << "Probe: " << frame_->frame1d().DeltaLmin()+1 << endl; 
	  break;  
	case 2:
	  k_[0] = (e_[0] == 0
		   ? frame_->frame1d().DeltaLmin()+1
		   : frame_->frame1d().Nablamin()+1);
	  k_[1] = (e_[1] == 0
		   ? frame_->frame1d().DeltaLmin()+1
		   : frame_->frame1d().Nablamin());
//        cout << "2.: " << k_[0] << k_[1] << endl; 
//        cout << "Probe: " << frame_->frame1d().DeltaLmin()+1 << endl; 
	  break;
	case 3:
	  k_[0] = (e_[0] == 0
		   ? frame_->frame1d().DeltaLmin()+1
		   : frame_->frame1d().Nablamin());
	  k_[1] = 0; // by convention;
	  break;
	case 4:
	  k_[0] = 0; // by convention
	  k_[1] = (e_[1] == 0
		   ? frame_->frame1d().DeltaLmin()+1
		   : frame_->frame1d().Nablamin());
	}
      }
    } else return *this;

    bool jplusplus = false;
    if (eplusplus) {
      // advance e
      if (e_[0] == 1 && e_[1] == 1)
	jplusplus = true;
      else {
	if (e_[0] == 1)
	  e_[1] = 1;
	else {
	  if (e_[1] == 0)
	    e_[1] = 1;
	  else {
              if(j_[1]==j0[1]){
                e_[0] = 1;
                e_[1] = 0;
              }
              else{
                e_[0] = 1;
                e_[1] = 1;
              }
                
	  }
	}

	// choose lowest patch number ...
	patch_ = 0;

	// ... and lowest translation index k = k(j,e,0)
	k_[0] = (e_[0] == 0
		 ? frame_->frame1d().DeltaLmin()+1
		 : frame_->frame1d().Nablamin());
	k_[1] = (e_[1] == 0
		 ? frame_->frame1d().DeltaLmin()+1
		 : frame_->frame1d().Nablamin()+1);
      }

    } else return *this;

    if (jplusplus) {
      // else: determine next level index
            // "small loop" "e_++" (j_ is not changed)
            // iterate over all combinations of generators/quarklets for all dimensions with j_[i]=j0[i]
            // this looks like binary addition: (in 3 Dim:) gwg is followed by gww (g=0,w=1)
            bool done = true;
            for (int i(1); i >= 0; i--)
            {
                // find first position on level j0
                if (j_[i] == j0[i])
                {
                    if (e_[i] == 1)
                    {
                        e_[i]=0;
                        k_[i]=frame_->frame1d().DeltaLmin()+1;
                    } else
                    {
                        e_[i]=1;
                        k_[i]=frame_->frame1d().Nablamin();
                        done = false;
                        break;
                    }
                }
            }
            // done == true bedeutet, dass alle Komponenten auf level j0() quarklets waren.
            // "big loop" "j_++"
            if (done == true)
            {
//                if (j_[1] != j0[1])
//                        {
//                            // increase left neighbor
//                            j_[0]=j_[0]+1;
//                            e_[0]=1;
//                            k_[0]=frame_->frame1d().Nablamin();
//                            int temp = j_[1]-j0[1];
//                            j_[1]=j0[1]+temp-1;
//                            e_[1]= (temp == 1?0:1);
//                            k_[1]= (temp == 1?frame_->frame1d().DeltaLmin()+1:frame_->frame1d().Nablamin());
//                            
//                        }
//                else{
//                            j_[0]=j0[1]+j_[0]-j0[0]+1;
//                            e_[0]=1;
//                            k_[0]=frame_->frame1d().Nablamin();
////                            j_[0]=j0[0];
////                            e_[0]=0;
////                            k_[0]=frame_->frame1d().DeltaLmin()+1;
//                        
//                             // unnoetig, da i==0 gilt.
//                }
                
                for(int i =1;i>=0;i--){
                    if (i == 1)
                    {
                        if (j_[1] != j0[1])
                        {
                            // increase left neighbor
                            j_[0]=j_[0]+1;
                            e_[0]=1;
                            k_[0]=frame_->frame1d().Nablamin();
                            int temp = j_[1]-j0[1];
                            j_[1]=j0[1]+temp-1;
                            e_[1]= (temp == 1?0:1);
                            k_[1]= (temp == 1?frame_->frame1d().DeltaLmin()+1:frame_->frame1d().Nablamin()+1);
                            break;
                        }
                    }
//                i=0;    
                
                    else // i == 0. "big loop" arrived at the last index. We have to increase the level!
                    {
                        j_[1]=j0[1]+j_[0]-j0[0]+1;
                        e_[1]=1;
                        k_[1]=frame_->frame1d().Nablamin()+1;
                        j_[0]=j0[0];
                        e_[0]=0;
                        k_[0]=frame_->frame1d().DeltaLmin()+1;
                        break; // unnoetig, da i==0 gilt.
                    }
                }
                    
                 
            }
            
        
        
//        ++j_;

      // choose lowest type (we have to advance to a quarklet) ...
//      e_[0] = 0;
//      e_[1] = 1;
      
      // ... lowest patch number ...
      patch_ = 0;
      
      // ... and lowest translation index k = k(j,(0,1),0)
//      k_[0] = frame_->frame1d().DeltaLmin()+1;
//      k_[1] = frame_->frame1d().Nablamin();
    }
    
    return *this;
  }

  template <class IFRAME>
  bool
  LDomainFrameIndex<IFRAME>::operator < (const LDomainFrameIndex& lambda) const
  {
    // standard lexicographic order on (j,e,p,k),
    // we assume that e and k are already lexicographically ordered (cf. MultiIndex)
//      return (num_<lambda.number());
    return (multi_degree(p_) < multi_degree(lambda.p())  ||
                     ((multi_degree(p_) == multi_degree(lambda.p()) && p_ < lambda.p())  ||
                      (p_ == lambda.p() && 
                       ( multi_degree(j_) < multi_degree(lambda.j()) ||
                        ((multi_degree(j_) == multi_degree(lambda.j()) && j_ < lambda.j()) ||
                         (j_ == lambda.j() && 
                          (e_ < lambda.e() ||
                           (e_ == lambda.e() &&
                            (patch_ < lambda.patch() ||
                             (patch_ == lambda.patch() &&  k_.lex(lambda.k())                           
                             )
                            )
                           )
                          ) 
                         )
                        )                                             
                       )                       
                      )                     
                     )
            );
  }

//  template <class IFRAME>
//  const int
//  LDomainFrameIndex<IFRAME>::number() const
//  {
//    const int ecode(e()[0]+2*e()[1]);
//    
//    const int DeltaLmin = frame_->frame1d().DeltaLmin();
//    const int Nablamin  = frame_->frame1d().Nablamin();
//    const int Deltasize = frame_->frame1d().Deltasize(j());
//    const int Nablasize = frame_->frame1d().Nablasize(j());
//
//    switch(ecode) {
//    case 0: // e=(0,0)
//      switch(p()) {
//      case 0:
//	return
//	  (k()[0]-DeltaLmin-1)*(Deltasize-2)
//	  + (k()[1]-DeltaLmin-1);
//	break;
//      case 1:
//	return
//	  (Deltasize-2)*(Deltasize-2)
//	  + (k()[0]-DeltaLmin-1)*(Deltasize-2)
//	  + (k()[1]-DeltaLmin-1);
//	break;
//      case 2:
//	return
//	  2*(Deltasize-2)*(Deltasize-2)
//	  + (k()[0]-DeltaLmin-1)*(Deltasize-2)
//	  + (k()[1]-DeltaLmin-1);
//	break;
//      case 3:
//	return
//	  3*(Deltasize-2)*(Deltasize-2)
//	  + (k()[0]-DeltaLmin-1);
//	break;
//      case 4:
//	return
//	  3*(Deltasize-2)*(Deltasize-2)
//	  + (Deltasize-2)
//	  + (k()[1]-DeltaLmin-1);
//	break;
//      default:
//	return 42; // for the compiler, this will never happen
//      }
//      break;
//    case 1: // e=(1,0)
//      switch(p()) {
//      case 0:
//	return
//	  frame_->Deltasize(j())
//	  + frame_->Nabla01size(j())
//	  + (k()[0]-Nablamin)*(Deltasize-2)
//	  + (k()[1]-DeltaLmin-1);
//	break;
//      case 1:
//	return
//	  frame_->Deltasize(j())
//	  + frame_->Nabla01size(j())
//	  + (Deltasize-2)*Nablasize
//	  + (k()[0]-Nablamin)*(Deltasize-2)
//	  + (k()[1]-DeltaLmin-1);
//	break;
//      case 2:
//	return
//	  frame_->Deltasize(j())
//	  + frame_->Nabla01size(j())
//	  + 2*(Deltasize-2)*Nablasize
//	  + (k()[0]-Nablamin)*(Deltasize-2)
//	  + (k()[1]-DeltaLmin-1);
//	break;
//      case 3:
//	return
//	  frame_->Deltasize(j())
//	  + frame_->Nabla01size(j())
//	  + 3*(Deltasize-2)*Nablasize
//	  + (k()[0]-Nablamin);
//	break;
//      default:
//	return 42; // for the compiler, this will never happen
//      }
//      break;
//    case 2: // e=(0,1)
//      switch(p()) {
//      case 0:
//	return
//	  frame_->Deltasize(j())
//	  + (k()[0]-DeltaLmin-1)*Nablasize
//	  + (k()[1]-Nablamin);
//	break;
//      case 1:
//	return
//	  frame_->Deltasize(j())
//	  + (Deltasize-2)*Nablasize
//	  + (k()[0]-DeltaLmin-1)*Nablasize
//	  + (k()[1]-Nablamin);
//	break;
//      case 2:
//	return
//	  frame_->Deltasize(j())
//	  + 2*(Deltasize-2)*Nablasize
//	  + (k()[0]-DeltaLmin-1)*Nablasize
//	  + (k()[1]-Nablamin);
//	break;
//      case 4:
//	return
//	  frame_->Deltasize(j())
//	  + 3*(Deltasize-2)*Nablasize
//	  + (k()[1]-Nablamin);
//	break;
//      default:
//	return 42; // for the compiler, this will never happen
//      }
//    case 3: // e=(1,1)
//      switch(p()) {
//      case 0:
//	return
//	  frame_->Deltasize(j())
//	  + frame_->Nabla01size(j())
//	  + frame_->Nabla10size(j())
//	  + (k()[0]-Nablamin)*Nablasize
//	  + (k()[1]-Nablamin);
//	break;
//      case 1:
//	return
//	  frame_->Deltasize(j())
//	  + frame_->Nabla01size(j())
//	  + frame_->Nabla10size(j())
//	  + Nablasize*Nablasize
//	  + (k()[0]-Nablamin)*Nablasize
//	  + (k()[1]-Nablamin);
//	break;
//      case 2:
//	return
//	  frame_->Deltasize(j())
//	  + frame_->Nabla01size(j())
//	  + frame_->Nabla10size(j())
//	  + 2*Nablasize*Nablasize
//	  + (k()[0]-Nablamin)*Nablasize
//	  + (k()[1]-Nablamin);
//	break;
//      default:
//	return 42; // for the compiler, this will never happen
//      }
//      break;
//    default:
//      return 42; // for the compiler, this will never happen
//    }
//  }
  
  template <class IFRAME>
  LDomainFrameIndex<IFRAME>
  first_generator(const LDomainFrame<IFRAME>* frame, const typename LDomainFrameIndex<IFRAME>::level_type& j, const typename LDomainFrameIndex<IFRAME>::polynomial_type& p)
  {
    assert(j >= frame->j0());

    typename LDomainFrameIndex<IFRAME>::type_type e;

    // setup lowest translation index for e=(0,0), p=0
    typename LDomainFrameIndex<IFRAME>::translation_type k(frame->frame1d().DeltaLmin()+1,
						      frame->frame1d().DeltaLmin()+1);
    
        return LDomainFrameIndex<IFRAME>(p,j, e, 0, k, p.number()* frame->get_Nablasize(), frame);
    
//        return LDomainFrameIndex<IFRAME>(p,j, e, 0, k, number, frame);
  }

  template <class IFRAME>
  LDomainFrameIndex<IFRAME>
  last_generator(const LDomainFrame<IFRAME>* frame, const typename LDomainFrameIndex<IFRAME>::level_type& j,const typename LDomainFrameIndex<IFRAME>::polynomial_type& p )
  {
    assert(j >= frame->j0() && j[0]==j[1]);

    typename LDomainFrameIndex<IFRAME>::type_type e;

    // setup highest translation index for e=(0,0), p=4
    typename LDomainFrameIndex<IFRAME>::translation_type k(0, frame->frame1d().DeltaRmax(j[1])-1);
    
    return LDomainFrameIndex<IFRAME>(p,j, e, 4, k, p.number()* frame->get_Nablasize()+frame->Deltasize(j[0])-1, frame);
  }

  template <class IFRAME>
  LDomainFrameIndex<IFRAME>
  first_quarklet(const LDomainFrame<IFRAME>* frame, const typename LDomainFrameIndex<IFRAME>::level_type& j, const typename LDomainFrameIndex<IFRAME>::polynomial_type& p)
  {
    assert(j >= frame->j0());

    typename LDomainFrameIndex<IFRAME>::type_type e;
    typename LDomainFrameIndex<IFRAME>::translation_type k;
    typename LDomainFrameIndex<IFRAME>::level_type jdiff;
    /*(frame->frame1d().DeltaLmin()+1,
						      frame->frame1d().Nablamin()+1);*/
    
    bool sofar_only_generators = true;
    
    if (j[0] == frame->j0()[0] )
    {
        e[0] = 0;
        k[0] = frame->frame1d().DeltaLmin()+1;
    } else
    {
        e[0] = 1;
        k[0] = frame->frame1d().Nablamin();
        sofar_only_generators = false;
    }
    
    if ( (sofar_only_generators == true) || (j[1] != frame->j0()[1]) )
    {
        e[1] = 1;
        k[1] = frame->frame1d().Nablamin()+1;
        //sofar_only_generators = false;
    } else
    {
        e[1] = 0;
        k[1] = frame->frame1d().DeltaLmin()+1;
    }
    
    jdiff[0]= j[0]-frame->j0()[0], jdiff[1]= j[1]-frame->j0()[1];
    int level = jdiff.number();
    int number = p.number()* frame->get_Nablasize()+frame->get_first_wavelet_numbers()[level];
    
    
    
    return LDomainFrameIndex<IFRAME>(p, j, e, 0, k, number, frame);
  }

//  template <class IFRAME>
//  LDomainFrameIndex<IFRAME>
//  first_quarklet(const LDomainFrame<IFRAME>* frame,
//		const typename LDomainFrameIndex<IFRAME>::level_type& j,
//		const typename LDomainFrameIndex<IFRAME>::type_type& e,
//                const typename LDomainFrameIndex<IFRAME>::polynomial_type& p)
//  {
//    assert(j >= frame->j0());
//
////    typename LDomainFrameIndex<IFRAME>::type_type e(ewish);
//    
//    // setup lowest translation index appropriately
//    typename LDomainFrameIndex<IFRAME>::translation_type k;
//    const int ecode(e[0]+2*e[1]);
//    if (ecode == 0) {
//      // e = (0,0)
//      k[0] = k[1] = frame->frame1d().DeltaLmin()+1;
//    } else {
//      if (ecode == 1) {
//	// e = (1,0)
//	k[0] = frame->frame1d().Nablamin()+1;
//	k[1] = frame->frame1d().DeltaLmin()+1;
//      } else {
//	if (ecode == 2) {
//	  // e = (0,1)
//	  k[0] = frame->frame1d().DeltaLmin()+1;
//	  k[1] = frame->frame1d().Nablamin()+1;
//	} else {
//	  // e = (1,1)
//	  k[0] = k[1] = frame->frame1d().Nablamin()+1;
//	}
//      }
//    }
//    
//    return LDomainFrameIndex<IFRAME>(p, j, e, 0, k, frame);
//  }
  
  template <class IFRAME>
  LDomainFrameIndex<IFRAME>
  last_quarklet(const LDomainFrame<IFRAME>* frame, const typename LDomainFrameIndex<IFRAME>::level_type& j, const typename LDomainFrameIndex<IFRAME>::polynomial_type& p)
  {
    assert(j >= frame->j0());
    
    typename LDomainFrameIndex<IFRAME>::type_type e(1, 1);
    typename LDomainFrameIndex<IFRAME>::level_type jdiff;
    
    // setup highest translation index for e=(1,1), p=2
    typename LDomainFrameIndex<IFRAME>::translation_type k(0,
						      frame->frame1d().Nablamax(j[1]));
    
    jdiff[0]= j[0]-frame->j0()[0], jdiff[1]= j[1]-frame->j0()[1];
    int level = jdiff.number();
    int number = p.number()* frame->get_Nablasize()+frame->get_last_wavelet_numbers()[level];
    
    return LDomainFrameIndex<IFRAME>(p, j, e, 4, k, number, frame);
  }
}
