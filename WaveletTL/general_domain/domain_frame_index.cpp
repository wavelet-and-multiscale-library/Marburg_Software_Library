
#include "domain_frame_index.h"
#include "domain_frame.h"

// implementation for domain_index.h

namespace WaveletTL
{
  template <class IFRAME, int NPATCHES>
  DomainFrameIndex<IFRAME, NPATCHES>::DomainFrameIndex(const DomainFrame<IFRAME, NPATCHES>* frame)
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
      k_[0] = frame->frame1d()->DeltaLmin();
      k_[1] = k_[0]+1;
      num_ = 0;
    }
  }

  template <class IFRAME, int NPATCHES>
  DomainFrameIndex<IFRAME, NPATCHES>::DomainFrameIndex(const polynomial_type& p,
                                     const level_type& j,
				     const type_type& e,
				     const int patch,
				     const translation_type& k,
                                     const unsigned int number,
				     const DomainFrame<IFRAME, NPATCHES>* frame)
    : frame_(frame), p_(p),j_(j), e_(e), patch_(patch), k_(k), num_(number)
  {
  }

  template <class IFRAME, int NPATCHES>
  DomainFrameIndex<IFRAME, NPATCHES>::DomainFrameIndex(const DomainFrameIndex& lambda)
    : frame_(lambda.frame_), p_(lambda.p_), j_(lambda.j_), e_(lambda.e_), patch_(lambda.patch_), k_(lambda.k_), num_(lambda.num_)
  {
  }

  template <class IFRAME, int NPATCHES>
  DomainFrameIndex<IFRAME, NPATCHES>::DomainFrameIndex(const DomainFrameIndex* lambda)
    : frame_(lambda->frame_), p_(lambda->p_), j_(lambda->j_), e_(lambda->e_), patch_(lambda->patch_), k_(lambda->k_), num_(lambda->num_)
  {
  }


  template <class IFRAME, int NPATCHES>
  bool
  DomainFrameIndex<IFRAME, NPATCHES>::operator == (const DomainFrameIndex& lambda) const
  {
    return (p_ == lambda.p() &&
            j_ == lambda.j() &&
	    e_ == lambda.e() &&
	    patch_ == lambda.patch() &&
	    k_ == lambda.k());
  }

  template <class IFRAME, int NPATCHES>
  DomainFrameIndex<IFRAME, NPATCHES>&
  DomainFrameIndex<IFRAME, NPATCHES>::operator ++ ()
  {
//      cout << "DeltaLmin=" << frame_->frame1d()->).DeltaLmin() <<endl;
//      cout << "DeltaRmax=" << frame_->frame1d()->).DeltaRmax(3) <<endl;
      num_++;
      level_type j0(frame_->j0());
    // decide whether the patch number has to be increased
    bool patchplusplus = false;
    for (int i = 1; i >= 0; i--) { //i: direction y, then x
      // determine the highest possible translation index into the i-th direction,
      // this will in general depend on the current patch number
      int last_index = 0;
      if (e_[i] == 0) { // generator in the i-th direction
        if(patch_<frame_->num_real_patches()){ //current patch is real
            last_index = frame_->frame1d()->DeltaRmax(j_[i])-1;
        }
        else{ //current patch is interface
            int patch=frame_->get_extensions()[patch_-frame_->num_real_patches()][0];
            int target_patch=frame_->get_extensions()[patch_-frame_->num_real_patches()][1];
            if(frame_->get_corners()[target_patch][0]==frame_->get_corners()[patch][0]){ //extension in y-direction
                    if(frame_->get_corners()[target_patch][1]<frame_->get_corners()[patch][1]){ //extension from north to south
                        last_index = (i == 0
			? frame_->frames(patch_,0)->DeltaRmax(j_[i])
			: frame_->frame1d()->DeltaLmin()); 
                    }
                    else{
                         //extension from south to north
                        last_index = (i == 0
			? frame_->frames(patch_,0)->DeltaRmax(j_[i])
			: frame_->frame1d()->DeltaRmax(j_[i])); 
                    }
            }
            else{ //extension in x-direction
                    if(frame_->get_corners()[target_patch][0]<frame_->get_corners()[patch][0]){ //extension from east to west
                        last_index = (i == 0
			? frame_->frame1d()->DeltaLmin() 
			: frame_->frames(patch_,1)->DeltaRmax(j_[i]));
                    }
                    else{
                         //extension from west to east
                        last_index = (i == 0
			? frame_->frame1d()->DeltaRmax(j_[i])
			: frame_->frames(patch_,1)->DeltaRmax(j_[i]));
                    }
            }
        }

      } else { //wavelet
          if(patch_<frame_->num_real_patches()){ //current patch is real
            //last_index = frame_->frame1d()->Nablamax(j_[i]);
            //loop over non-logical patches
            last_index=frame_->frame1d()->Nablamax(j_[i]);  
            for(unsigned int ext1=0;ext1<frame_->get_extensions().size();ext1++){
                if(frame_->get_extensions()[ext1][0]==patch_){ //extension from patch to another one
                    int target_patch=frame_->get_extensions()[ext1][1];
                    if(frame_->get_corners()[target_patch][0]==frame_->get_corners()[patch_][0]){ //extension in y-direction
                        if(frame_->get_corners()[target_patch][1]>frame_->get_corners()[patch_][1] && i==1){ //extension from south to north
                             //quarklet indexing always starts at 0. in case of extension, quarklet 0 lies on the interface
                            last_index=frame_->frame1d()->Nablamax(j_[i])-1;
                        }
                        
                    }
                    if(frame_->get_corners()[target_patch][1]==frame_->get_corners()[patch_][1] ){ //extension in x-direction
                        if(frame_->get_corners()[target_patch][0]>frame_->get_corners()[patch_][0] && i==0){ //extension from west to east
                            last_index=frame_->frame1d()->Nablamax(j_[i])-1; //drop last wavelet index to interface
                        }
                        
                    }
                }
            }
          }
          
          else{ //current patch is interface
            int patch=frame_->get_extensions()[patch_-frame_->num_real_patches()][0];
            int target_patch=frame_->get_extensions()[patch_-frame_->num_real_patches()][1];
            if(frame_->get_corners()[target_patch][0]==frame_->get_corners()[patch][0]){ //extension in y-direction
                    if(frame_->get_corners()[target_patch][1]<frame_->get_corners()[patch][1]){ //extension from north to south
                        last_index = (i == 0
                            ? frame_->frame1d()->Nablamax(j_[i]) 
                            : frame_->frame1d()->Nablamin()); 
                    }
                    else{
                         //extension from south to north
                         last_index = (i == 0
                            ? frame_->frame1d()->Nablamax(j_[i])
                            : frame_->frame1d()->Nablamax(j_[i])); 
                    }
            }
            else{ //extension in x-direction
                    if(frame_->get_corners()[target_patch][0]<frame_->get_corners()[patch][0]){ //extension from east to west
                        last_index = (i == 0
                            ? frame_->frame1d()->Nablamin() 
                            : frame_->frame1d()->Nablamax(j_[i]));
                    }
                    else{
                         //extension from west to east
                        last_index = (i == 0
                            ? frame_->frame1d()->Nablamax(j_[i])
                            : frame_->frame1d()->Nablamax(j_[i]));
                    }
            }
          }
                  
      }

      if (k_[i] == last_index) {
//          cout << "Hier1" << endl;
	// reset k_[i] to the lowest possible translation index
	if (e_[i] == 0) { // generator, 
          if(patch_<frame_->num_real_patches()){ //current patch is real
            k_[i] = frame_->frame1d()->DeltaLmin()+1;
          }
          else{ //interface
            int patch=frame_->get_extensions()[patch_-frame_->num_real_patches()][0];
            int target_patch=frame_->get_extensions()[patch_-frame_->num_real_patches()][1];
            if(frame_->get_corners()[target_patch][0]==frame_->get_corners()[patch][0]){ //extension in y-direction
                    if(frame_->get_corners()[target_patch][1]<frame_->get_corners()[patch][1]){ //extension from north to south
                        k_[i] = (i == 0
                            ? frame_->frame1d()->DeltaLmin()+1
                            : frame_->frame1d()->DeltaLmin()); 
                    }
                    else{
                         //extension from south to north
                         k_[i] = (i == 0
                            ? frame_->frame1d()->DeltaLmin()+1
                            : frame_->frame1d()->DeltaRmax(j_[i]));
                    }
            }
            else{ //extension in x-direction
                    if(frame_->get_corners()[target_patch][0]<frame_->get_corners()[patch][0]){ //extension from east to west
                        k_[i] = (i == 0
                            ? frame_->frame1d()->DeltaLmin() 
                            : frame_->frame1d()->DeltaLmin()+1);
                    }
                    else{
                         //extension from west to east
                        k_[i] = (i == 0
                            ? frame_->frame1d()->DeltaRmax(j_[i]) 
                            : frame_->frame1d()->DeltaLmin()+1);
                    }
            }
          }
	  
	} else { // quarklet
            if(patch_<frame_->num_real_patches()){ //current patch is real
            //last_index = frame_->frame1d()->Nablamax(j_[i]);
            //loop over non-logical patches
            k_[i]=frame_->frame1d()->Nablamin();  
            for(unsigned int ext1=0;ext1<frame_->get_extensions().size();ext1++){
                if(frame_->get_extensions()[ext1][0]==patch_){ //extension from patch to another one
                    int target_patch=frame_->get_extensions()[ext1][1];
                    if(frame_->get_corners()[target_patch][0]==frame_->get_corners()[patch_][0]){ //extension in y-direction
                        if(frame_->get_corners()[target_patch][1]<frame_->get_corners()[patch_][1] && i==1){ //extension from north to south
                             //quarklet indexing always starts at 0. in case of extension, quarklet 0 lies on the interface
                            k_[i]=frame_->frame1d()->Nablamin()+1;
                        }
                        
                    }
                    if(frame_->get_corners()[target_patch][1]==frame_->get_corners()[patch_][1] ){ //extension in x-direction
                        if(frame_->get_corners()[target_patch][0]<frame_->get_corners()[patch_][0] && i==0){ //extension from east to west
                            k_[i]=frame_->frame1d()->Nablamin()+1; //drop 0-th wavelet index to interface
                        }
                        
                    }
                }
            }
            }
            else{ //current patch is interface
            int patch=frame_->get_extensions()[patch_-frame_->num_real_patches()][0];
            int target_patch=frame_->get_extensions()[patch_-frame_->num_real_patches()][1];
            if(frame_->get_corners()[target_patch][0]==frame_->get_corners()[patch][0]){ //extension in y-direction
                    if(frame_->get_corners()[target_patch][1]<frame_->get_corners()[patch][1]){ //extension from north to south
                        k_[i] = (i == 0
                            ? frame_->frame1d()->Nablamin()+1 //really?
                            : frame_->frame1d()->Nablamin()); 
                    }
                    else{
                         //extension from south to north
                         k_[i] = (i == 0
                            ? frame_->frame1d()->Nablamin()+1 //really?
                            : frame_->frame1d()->Nablamax(j_[i]));
                    }
            }
            else{ //extension in x-direction
                    if(frame_->get_corners()[target_patch][0]<frame_->get_corners()[patch][0]){ //extension from east to west
                        k_[i] = (i == 0
                            ? frame_->frame1d()->Nablamin() 
                            : frame_->frame1d()->Nablamin()+1); //really?
                    }
                    else{
                         //extension from west to east
                        k_[i] = (i == 0
                            ? frame_->frame1d()->Nablamax(j_[i])
                            : frame_->frame1d()->Nablamin()+1); //really?
                    }
            }
            }
//	  k_[i] = frame_->frame1d()->Nablamin(); // should be 0
	}
	patchplusplus = (i == 0);
      } else {
	++k_[i];
	break;
      }
    }
//    cout << k_[0] << k_[1] << endl;
    bool eplusplus = false;
    if (patchplusplus) {
//        cout<<"jump to next patch"<<endl;
        if(patch_<frame_->num_real_patches()+frame_->num_logical_patches()-1){
            patch_++;
        }
        else{
            eplusplus=true;
        }
     
        //jump to next patch
      if (!eplusplus) { // then choose lowest translation index k=k(j,e,p)
          if(patch_<frame_->num_real_patches()){ //current patch is real             
            k_[0] = (e_[0] == 0
		   ? frame_->frame1d()->DeltaLmin()+1
		   : frame_->frame1d()->Nablamin());
	    k_[1] = (e_[1] == 0
		   ? frame_->frame1d()->DeltaLmin()+1
		   : frame_->frame1d()->Nablamin()); 
            //loop over non-logical patches
            for(unsigned int ext1=0;ext1<frame_->get_extensions().size();ext1++){
                if(frame_->get_extensions()[ext1][0]==patch_){ //extension from patch to another one
                    int target_patch=frame_->get_extensions()[ext1][1];
                    if(frame_->get_corners()[target_patch][0]==frame_->get_corners()[patch_][0]){ //extension in y-direction
                        if(frame_->get_corners()[target_patch][1]<frame_->get_corners()[patch_][1] ){ //extension from north to south
                             //quarklet indexing always starts at 0. in case of extension, quarklet 0 lies on the interface
                            k_[0] = (e_[0] == 0
                                ? frame_->frame1d()->DeltaLmin()+1
                                : frame_->frame1d()->Nablamin());
                            k_[1] = (e_[1] == 0
                                ? frame_->frame1d()->DeltaLmin()+1
                            : frame_->frame1d()->Nablamin()+1);
                        }
                        
                    }
                    if(frame_->get_corners()[target_patch][1]==frame_->get_corners()[patch_][1] ){ //extension in x-direction
                        if(frame_->get_corners()[target_patch][0]<frame_->get_corners()[patch_][0] ){ //extension from east to west
                            k_[0] = (e_[0] == 0
                                ? frame_->frame1d()->DeltaLmin()+1
                                : frame_->frame1d()->Nablamin()+1);
                            k_[1] = (e_[1] == 0
                                ? frame_->frame1d()->DeltaLmin()+1
                            : frame_->frame1d()->Nablamin());
                        }
                        
                    }
                }
            }
          }
          else{ //current patch is interface
//              cout<<"jump to (logical) patch"<<patch_<<endl;
            int patch=frame_->get_extensions()[patch_-frame_->num_real_patches()][0];
            int target_patch=frame_->get_extensions()[patch_-frame_->num_real_patches()][1];
//                      cout<<"extension from patch "<<patch<<"to patch "<<target_patch<<endl;
            if(frame_->get_corners()[target_patch][0]==frame_->get_corners()[patch][0]){ //extension in y-direction
                    if(frame_->get_corners()[target_patch][1]<frame_->get_corners()[patch][1]){ //extension from north to south
                        k_[0] = (e_[0] == 0
                            ? frame_->frame1d()->DeltaLmin()+1
                            : frame_->frame1d()->Nablamin());
                        k_[1] = (e_[1] == 0
                            ? frame_->frame1d()->DeltaLmin()
                            : frame_->frame1d()->Nablamin());
                    }
                    else{
                         //extension from south to north
                        k_[0] = (e_[0] == 0
                            ? frame_->frame1d()->DeltaLmin()+1
                            : frame_->frame1d()->Nablamin());
                        k_[1] = (e_[1] == 0
                            ? frame_->frame1d()->DeltaRmax(j_[1])
                            : frame_->frame1d()->Nablamax(j_[1]));
                    }
            }
            else{ //extension in x-direction
                    if(frame_->get_corners()[target_patch][0]<frame_->get_corners()[patch][0]){ //extension from east to west
                        k_[0] = (e_[0] == 0
                            ? frame_->frame1d()->DeltaLmin()
                            : frame_->frame1d()->Nablamin()); // by convention
                        k_[1] = (e_[1] == 0
                            ? frame_->frame1d()->DeltaLmin()+1
                            : frame_->frame1d()->Nablamin());
                    }
                    else{
                         //extension from west to east
                        k_[0] = (e_[0] == 0
                            ? frame_->frame1d()->DeltaRmax(j_[0])
                            : frame_->frame1d()->Nablamax(j_[0])); // by convention
                        k_[1] = (e_[1] == 0
                            ? frame_->frame1d()->DeltaLmin()+1
                            : frame_->frame1d()->Nablamin());
                        
                    }
            }
          }
      }
    } else return *this;
//    cout<<"end jump to next patch"<<endl;
    bool jplusplus = false;
    if (eplusplus) {
      // advance e
//        cout<<"bin hier"<<endl;
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
//        cout<<"increase type"<<endl; //to do: extension from south no north, west to east
	// ... and lowest translation index k = k(j,e,0)
        k_[0] = (e_[0] == 0
                ? frame_->frame1d()->DeltaLmin()+1
                : frame_->frame1d()->Nablamin());
        k_[1] = (e_[1] == 0
                ? frame_->frame1d()->DeltaLmin()+1
                : frame_->frame1d()->Nablamin());
	for(unsigned int ext1=0;ext1<frame_->get_extensions().size();ext1++){
                if(frame_->get_extensions()[ext1][0]==patch_){ //extension from patch to another one
                    int target_patch=frame_->get_extensions()[ext1][1];
                    if(frame_->get_corners()[target_patch][0]==frame_->get_corners()[patch_][0]){ //extension in y-direction
                        if(frame_->get_corners()[target_patch][1]<frame_->get_corners()[patch_][1] ){ //extension from north to south
                             //quarklet indexing always starts at 0. in case of extension, quarklet 0 lies on the interface
                            k_[0] = (e_[0] == 0
                                ? frame_->frame1d()->DeltaLmin()+1
                                : frame_->frame1d()->Nablamin());
                            k_[1] = (e_[1] == 0
                                ? frame_->frame1d()->DeltaLmin()+1
                            : frame_->frame1d()->Nablamin()+1);
                        }
                        //new AS
//                        if(frame_->get_corners()[target_patch][1]>frame_->get_corners()[patch_][1] ){ //extension from south to north
//                             //quarklet indexing always starts at 0. in case of extension, quarklet 0 lies on the interface
//                            k_[0] = (e_[0] == 0
//                                ? frame_->frame1d()->DeltaLmin()+1
//                                : frame_->frame1d()->Nablamin());
//                            k_[1] = (e_[1] == 0
//                                ? frame_->frame1d()->DeltaLmin()+1
//                            : frame_->frame1d()->Nablamin()+0);
//                        }
                        
                    }
                    if(frame_->get_corners()[target_patch][1]==frame_->get_corners()[patch_][1] ){ //extension in x-direction
                        if(frame_->get_corners()[target_patch][0]<frame_->get_corners()[patch_][0] ){ //extension from east to west
                            k_[0] = (e_[0] == 0
                                ? frame_->frame1d()->DeltaLmin()+1
                                : frame_->frame1d()->Nablamin()+1);
                            k_[1] = (e_[1] == 0
                                ? frame_->frame1d()->DeltaLmin()+1
                            : frame_->frame1d()->Nablamin());
                        }
                        //new AS
//                        if(frame_->get_corners()[target_patch][0]>frame_->get_corners()[patch_][0] ){ //extension from west to east
//                            k_[0] = (e_[0] == 0
//                                ? frame_->frame1d()->DeltaLmin()+1
//                                : frame_->frame1d()->Nablamin()+0);
//                            k_[1] = (e_[1] == 0
//                                ? frame_->frame1d()->DeltaLmin()+1
//                            : frame_->frame1d()->Nablamin());
//                        }
                        
                    }
                }
        }
      }

    } else return *this;

    if (jplusplus) {
//        cout<<"increase level"<<endl;
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
                        k_[i]=frame_->frame1d()->DeltaLmin()+1;
                    } else
                    {
                        e_[i]=1;
                        k_[i]=frame_->frame1d()->Nablamin();
                        done = false;
                        break;
                    }
                }
            }
            // done == true bedeutet, dass alle Komponenten auf level j0() quarklets waren.
            // "big loop" "j_++"
            if (done == true) //is this loop dependent on patch?
            {
                
                for(int i =1;i>=0;i--){
                    if (i == 1)
                    {
                        if (j_[1] != j0[1])
                        {
                            // increase left neighbor
                            j_[0]=j_[0]+1;
                            e_[0]=1;
                            k_[0]=frame_->frame1d()->Nablamin();
                            int temp = j_[1]-j0[1];
                            j_[1]=j0[1]+temp-1;
                            e_[1]= (temp == 1?0:1);
                            //special treatment needed by extension
                            k_[1]= (temp == 1?frame_->frame1d()->DeltaLmin()+1:frame_->frame1d()->Nablamin()+0);
                            for(unsigned int ext1=0;ext1<frame_->get_extensions().size();ext1++){
                                if(frame_->get_extensions()[ext1][0]==0){ //extension from patch 0 to another one
                                    int target_patch=frame_->get_extensions()[ext1][1];
                                    if(frame_->get_corners()[target_patch][0]==frame_->get_corners()[0][0]){ //extension in y-direction
                                        if(frame_->get_corners()[target_patch][1]<frame_->get_corners()[0][1] ){ //extension from north to south
                                            k_[1]= (temp == 1?frame_->frame1d()->DeltaLmin()+1:frame_->frame1d()->Nablamin()+1);
                                        }
                                    }
                                }
                            }    
                            break;
                        }
                    }
//                i=0;    
                
                    else // i == 0. "big loop" arrived at the last index. We have to increase the level!
                    {
                        j_[1]=j0[1]+j_[0]-j0[0]+1;
                        e_[1]=1;
                        //special treatment needed by extension
                        //k_[1]=frame_->frame1d()->Nablamin()+1;
                        //new AS
                        k_[1]=frame_->frame1d()->Nablamin()+0;
                        for(unsigned int ext1=0;ext1<frame_->get_extensions().size();ext1++){
                                if(frame_->get_extensions()[ext1][0]==0){ //extension from patch 0 to another one
                                    int target_patch=frame_->get_extensions()[ext1][1];
                                    if(frame_->get_corners()[target_patch][0]==frame_->get_corners()[0][0]){ //extension in y-direction
                                        if(frame_->get_corners()[target_patch][1]<frame_->get_corners()[0][1] ){ //extension from north to south
                                            k_[1]=frame_->frame1d()->Nablamin()+1;
                                        }
                                    }
                                }
                        } 
                        j_[0]=j0[0];
                        e_[0]=0;
                        k_[0]=frame_->frame1d()->DeltaLmin()+1;
                        break; // unnoetig, da i==0 gilt.
                    }
                }
                    
                 
            }
            
//            cout<<"k0="<<k_[0]<<endl;
//            cout<<"k1="<<k_[1]<<endl;

               //        ++j_;

      // choose lowest type (we have to advance to a quarklet) ...
//      e_[0] = 0;
//      e_[1] = 1;
      
      // ... lowest patch number ...
      patch_ = 0;
      
      // ... and lowest translation index k = k(j,(0,1),0)
//      k_[0] = frame_->frame1d()->DeltaLmin()+1;
//      k_[1] = frame_->frame1d()->Nablamin();
      
    }
    
    return *this;
  }

  template <class IFRAME, int NPATCHES>
  bool
  DomainFrameIndex<IFRAME, NPATCHES>::operator < (const DomainFrameIndex& lambda) const
  {
    // standard lexicographic order on (j,e,p,k),
    // we assume that e and k are already lexicographically ordered (cf. MultiIndex)
      if(frame_->get_setup_full_collection()){
        return (num_<lambda.number());
      }
      else{
      
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
  }


  
  template <class IFRAME, int NPATCHES>
  DomainFrameIndex<IFRAME, NPATCHES>
  first_generator(const DomainFrame<IFRAME, NPATCHES>* frame, const typename DomainFrameIndex<IFRAME, NPATCHES>::level_type& j, const typename DomainFrameIndex<IFRAME, NPATCHES>::polynomial_type& p)
  {
    assert(j >= frame->j0());

    typename DomainFrameIndex<IFRAME, NPATCHES>::type_type e;

    // setup lowest translation index for e=(0,0), p=0
    typename DomainFrameIndex<IFRAME, NPATCHES>::translation_type k(frame->frame1d()->DeltaLmin()+1,
						      frame->frame1d()->DeltaLmin()+1);
    
        return DomainFrameIndex<IFRAME, NPATCHES>(p,j, e, 0, k, p.number()* frame->get_Nablasize(), frame);
    
//        return DomainFrameIndex<IFRAME>(p,j, e, 0, k, number, frame);
  }

  template <class IFRAME, int NPATCHES>
  DomainFrameIndex<IFRAME, NPATCHES>
  last_generator(const DomainFrame<IFRAME, NPATCHES>* frame, const typename DomainFrameIndex<IFRAME, NPATCHES>::level_type& j,const typename DomainFrameIndex<IFRAME, NPATCHES>::polynomial_type& p )
  {
    assert(j >= frame->j0() && j[0]==j[1]);

    typename DomainFrameIndex<IFRAME, NPATCHES>::type_type e;
    typename DomainFrameIndex<IFRAME, NPATCHES>::translation_type k;

    // setup highest translation index for e=(0,0), p=last patch
    int patch=frame->get_extensions()[frame->num_logical_patches()-1][0];
    int target_patch=frame->get_extensions()[frame->num_logical_patches()-1][1];
    if(frame.get_corners()[target_patch][0]==frame.get_corners()[patch][0]){ //extension in y-direction
                    if(frame->get_corners()[target_patch][1]<frame->get_corners()[patch][1]){ //extension from north to south
                        k[0]=frame->frame1d()->DeltaRmax(j[0])-1, k[1]=frame->frame1d()->DeltaLmin();
                    }
                    else{
                         //extension from south to north
                        k[0]=frame->frame1d()->DeltaRmax(j[0])-1, k[1]=frame->frame1d()->DeltaRmax(j[1]);
                    }
    }
    else{ //extension in x-direction
                    if(frame->get_corners()[target_patch][0]<frame->get_corners()[patch][0]){ //extension from east to west
                        k[0]=frame->frame1d()->DeltaLmin(), k[1]=frame->frame1d()->DeltaRmax(j[1])-1;
                    }
                    else{
                         //extension from west to east
                        k[0]=frame->frame1d()->DeltaRmax(j[0]), k[1]=frame->frame1d()->DeltaRmax(j[1])-1;
                    }
    }
    
    
    return Index(p,j, e, frame->num_real_patches()+frame->num_logical_patches()-1, k, p.number()* frame->get_Nablasize()+frame->Deltasize(j[0])-1, frame);
  }

  template <class IFRAME, int NPATCHES>
  DomainFrameIndex<IFRAME, NPATCHES>
  first_quarklet(const DomainFrame<IFRAME, NPATCHES>* frame, const typename DomainFrameIndex<IFRAME, NPATCHES>::level_type& j, const typename DomainFrameIndex<IFRAME, NPATCHES>::polynomial_type& p)
  {
    assert(j >= frame->j0());

    typename DomainFrameIndex<IFRAME, NPATCHES>::type_type e;
    typename DomainFrameIndex<IFRAME, NPATCHES>::translation_type k;
    typename DomainFrameIndex<IFRAME, NPATCHES>::level_type jdiff;
    /*(frame->frame1d()->DeltaLmin()+1,
						      frame->frame1d()->Nablamin()+1);*/
    
    bool sofar_only_generators = true;
    
    int patch=0; //first quarklet always lives on patch 0
    k[0]=k[1]=frame->frame1d()->Nablamin();  
            for(int ext1=0;ext1<frame->num_logical_patches();ext1++){
                if(frame->get_extensions()[ext1][0]==patch){ //extension from patch to another one
                    int target_patch=frame->get_extensions()[ext1][1];
                    if(frame->get_corners()[target_patch][0]==frame->get_corners()[patch][0]){ //extension in y-direction
                        if(frame->get_corners()[target_patch][1]<frame->get_corners()[patch][1]){ //extension from north to south
                             //quarklet indexing always starts at 0. in case of extension, quarklet 0 lies on the interface
                            if (j[0] == frame->j0()[0])
                            {
                                e[0] = 0;
                                k[0] = frame->frame1d()->DeltaLmin()+1;
                            } else
                            {
                                e[0] = 1;
                                k[0] = frame->frame1d()->Nablamin();
                                sofar_only_generators = false;
                            }
    
                            if ( (sofar_only_generators == true) || (j[1] != frame->j0()[0]) )
                            {
                                e[1] = 1;
                                k[1] = frame->frame1d()->Nablamin()+1;
                                //sofar_only_generators = false;
                            } else
                            {
                                e[1] = 0;
                                k[1] = frame->frame1d()->DeltaLmin()+1;
                            }
                        }
                        
                    }
                    if(frame->get_corners()[target_patch][1]==frame->get_corners()[patch][1] ){ //extension in x-direction
                        if(frame->get_corners()[target_patch][0]<frame->get_corners()[patch][0] ){ //extension from east to west
                             //drop 0-th wavelet index to interface
                            if (j[1] == frame->j0()[1])
                            {
                                e[1] = 0;
                                k[1] = frame->frame1d()->DeltaLmin()+1;
                            } else
                            {
                                e[1] = 1;
                                k[1] = frame->frame1d()->Nablamin();
                                sofar_only_generators = false;
                            }
    
                            if ( (sofar_only_generators == true) || (j[0] != frame->j0()[1]) )
                            {
                                e[0] = 1;
                                k[0] = frame->frame1d()->Nablamin()+1;
                                //sofar_only_generators = false;
                            } else
                            {
                                e[0] = 0;
                                k[0] = frame->frame1d()->DeltaLmin()+1;
                            }
                        }
                        
                    }
                }
            }
    
    jdiff[0]= j[0]-frame->j0()[0], jdiff[1]= j[1]-frame->j0()[1];
    int level = jdiff.number();
    int number = p.number()* frame->get_Nablasize()+frame->get_first_wavelet_numbers()[level];
    
    
    
    return DomainFrameIndex<IFRAME, NPATCHES>(p, j, e, 0, k, number, frame);
  }


  template <class IFRAME, int NPATCHES>
  DomainFrameIndex<IFRAME, NPATCHES>
  last_quarklet(const DomainFrame<IFRAME, NPATCHES>* frame, const typename DomainFrameIndex<IFRAME, NPATCHES>::level_type& j, const typename DomainFrameIndex<IFRAME, NPATCHES>::polynomial_type& p)
  {
    assert(j >= frame->j0());
    
    typename DomainFrameIndex<IFRAME, NPATCHES>::type_type e(1, 1);
    typename DomainFrameIndex<IFRAME, NPATCHES>::level_type jdiff;
    typename DomainFrameIndex<IFRAME, NPATCHES>::translation_type k;
    
    // setup highest translation index for e=(1,1), p=2
//    typename DomainFrameIndex<IFRAME, NPATCHES>::translation_type k(0,
//						      frame->frame1d()->Nablamax(j[1]));
//    
    int patch=frame->get_extensions()[frame->num_logical_patches()-1][0];
            int target_patch=frame->get_extensions()[frame->num_logical_patches()-1][1];
            if(frame->get_corners()[target_patch][0]==frame->get_corners()[patch][0]){ //extension in y-direction
                    if(frame->get_corners()[target_patch][1]<frame->get_corners()[patch][1]){ //extension from north to south
                        
                            k[0]= frame->frame1d()->Nablamax(j[0]) ;
                            k[1]= frame->frame1d()->Nablamin(); 
                    }
                    else{
                         //extension from south to north
                         
                            k[0]= frame->frame1d()->Nablamax(j[0]);
                            k[1]= frame->frame1d()->Nablamax(j[1]); 
                    }
            }
            else{ //extension in x-direction
                    if(frame->get_corners()[target_patch][0]<frame->get_corners()[patch][0]){ //extension from east to west
                        
                            k[0]= frame->frame1d()->Nablamin(); 
                            k[1]= frame->frame1d()->Nablamax(j[1]);
                    }
                    else{
                         //extension from west to east
                        
                            k[0]= frame->frame1d()->Nablamax(j[0]);
                            k[1]= frame->frame1d()->Nablamax(j[1]);
                    }
            }
    
    jdiff[0]= j[0]-frame->j0()[0], jdiff[1]= j[1]-frame->j0()[1];
    int level = jdiff.number();
    int number = p.number()* frame->get_Nablasize()+frame->get_last_wavelet_numbers()[level];
    
    return DomainFrameIndex<IFRAME, NPATCHES>(p, j, e, frame->num_real_patches()+frame->num_logical_patches()-1, k, number, frame);
  }
}
