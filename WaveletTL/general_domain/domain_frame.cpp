// implementation for domain_frame.h

#include <cmath>
#include <time.h>
#include <iostream>

//#include "domain_frame_index.h"

using std::cout;
using std::endl;

namespace WaveletTL
{
  template <class IFRAME, int NPATCHES>
  DomainFrame<IFRAME, NPATCHES>::DomainFrame(const IntervalFrame* frame1d, const IntervalFrame* frame1d_11, 
          const IntervalFrame* frame1d_01, const IntervalFrame* frame1d_10, 
          const Array1D<Point<2, int> >& corners, const Array1D<FixedArray1D<int,2> >& extensions)
    : frame1d_(frame1d), frame1d_11_(frame1d_11), frame1d_01_(frame1d_01), frame1d_10_(frame1d_10), corners_(corners), extensions_(extensions)
  {
      j0_[0] = frame1d_11_->j0();
      j0_[1] = frame1d_11_->j0();

      setup_full_collection_ = false;
      precomputed_supports_ = false;
  }



  template <class IFRAME, int NPATCHES>
  const int
  DomainFrame<IFRAME, NPATCHES>::Deltasize(const int j) const {
    const unsigned int Deltaj = frame1d_->Deltasize(j);
    return num_real_patches()*(Deltaj-2)*(Deltaj-2)+num_logical_patches()*(Deltaj-2);
  }


  
  template <class IFRAME, int NPATCHES>
  typename DomainFrame<IFRAME, NPATCHES>::Index
  DomainFrame<IFRAME, NPATCHES>::first_generator(const level_type& j, const polynomial_type& p) const
  {
    assert(j >= j0_);

    typename Index::type_type e;

    // setup lowest translation index for e=(0,0), p=0
    typename Index::translation_type k(frame1d_->DeltaLmin()+1, frame1d_->DeltaLmin()+1);
    
//   if(number==-1)
        return Index(p,j, e, 0, k, p.number()* Nablasize_, this);
//    else
//        return Index(p,j, e, 0, k, number, this);
  }

  template <class IFRAME, int NPATCHES>
  typename DomainFrame<IFRAME, NPATCHES>::Index
  DomainFrame<IFRAME, NPATCHES>::last_generator(const level_type& j, const polynomial_type& p) const
  {
    assert(j >= j0_);

    typename Index::type_type e;
    typename Index::translation_type k;

    // setup highest translation index for e=(0,0), p=last patch
    int patch=extensions_[num_logical_patches()-1][0];
    int target_patch=extensions_[num_logical_patches()-1][1];
    if(corners_[target_patch][0]==corners_[patch][0]){ //extension in y-direction
                    if(corners_[target_patch][1]<corners_[patch][1]){ //extension from north to south
                        k[0]=frame1d_->DeltaRmax(j[0])-1, k[1]=frame1d_->DeltaLmin();
                    }
                    else{
                         //extension from south to north
                        k[0]=frame1d_->DeltaRmax(j[0])-1, k[1]=frame1d_->DeltaRmax(j[1]);
                    }
    }
    else{ //extension in x-direction
                    if(corners_[target_patch][0]<corners_[patch][0]){ //extension from east to west
                        k[0]=frame1d_->DeltaLmin(), k[1]=frame1d_->DeltaRmax(j[1])-1;
                    }
                    else{
                         //extension from west to east
                        k[0]=frame1d_->DeltaRmax(j[0]), k[1]=frame1d_->DeltaRmax(j[1])-1;
                    }
    }
    
    
    return Index(p,j, e, num_real_patches()+num_logical_patches()-1, k, p.number()* Nablasize_+Deltasize(j[0])-1, this);
//    return Index(p,j, e, num_real_patches()+num_logical_patches()-1, k, 0, this);
  }

  template <class IFRAME, int NPATCHES>
  typename DomainFrame<IFRAME, NPATCHES>::Index
  DomainFrame<IFRAME, NPATCHES>::first_quarklet(const level_type& j, const polynomial_type& p, const int& number) const
  {
    assert(j >= j0());

    typename Index::type_type e;
    typename Index::translation_type k;
    
    bool sofar_only_generators = true;
    
    int patch=0; //first quarklet always lives on patch 0
    k[0]=k[1]=frame1d_->Nablamin();  
            for(int ext1=0;ext1<num_logical_patches();ext1++){
                if(extensions_[ext1][0]==patch){ //extension from patch to another one
                    int target_patch=extensions_[ext1][1];
                    if(corners_[target_patch][0]==corners_[patch][0]){ //extension in y-direction
                        if(corners_[target_patch][1]<corners_[patch][1]){ //extension from north to south
                             //quarklet indexing always starts at 0. in case of extension, quarklet 0 lies on the interface
                            if (j[0] == j0_[0])
                            {
                                e[0] = 0;
                                k[0] = frame1d_->DeltaLmin()+1;
                            } else
                            {
                                e[0] = 1;
                                k[0] = frame1d_->Nablamin();
                                sofar_only_generators = false;
                            }
    
                            if ( (sofar_only_generators == true) || (j[1] != j0_[0]) )
                            {
                                e[1] = 1;
                                k[1] = frame1d_->Nablamin()+1;
                                //sofar_only_generators = false;
                            } else
                            {
                                e[1] = 0;
                                k[1] = frame1d_->DeltaLmin()+1;
                            }
                        }
                        
                    }
                    if(corners_[target_patch][1]==corners_[patch][1] ){ //extension in x-direction
                        if(corners_[target_patch][0]<corners_[patch][0] ){ //extension from east to west
                             //drop 0-th wavelet index to interface
                            if (j[1] == j0_[1])
                            {
                                e[1] = 0;
                                k[1] = frame1d_->DeltaLmin()+1;
                            } else
                            {
                                e[1] = 1;
                                k[1] = frame1d_->Nablamin();
                                sofar_only_generators = false;
                            }
    
                            if ( (sofar_only_generators == true) || (j[0] != j0_[1]) )
                            {
                                e[0] = 1;
                                k[0] = frame1d_->Nablamin()+1;
                                //sofar_only_generators = false;
                            } else
                            {
                                e[0] = 0;
                                k[0] = frame1d_->DeltaLmin()+1;
                            }
                        }
                        
                    }
                }
            }
    
    
    
    if (number==-1)
    {
        level_type jdiff;
        jdiff[0]= j[0]-j0_[0], jdiff[1]= j[1]-j0_[1];
        int level = jdiff.number();
        int altnumber = p.number()* Nablasize_+first_wavelet_numbers[level];
        return Index(p,j, e, 0, k, altnumber, this);
    }
    else
       return Index(p,j, e, 0, k, number, this); 
  }



  template <class IFRAME, int NPATCHES>
  typename DomainFrame<IFRAME, NPATCHES>::Index
  DomainFrame<IFRAME, NPATCHES>::last_quarklet(const level_type& j, const polynomial_type& p, const int& number) const
  {
    assert(j >= j0());
    
    typename Index::type_type e(1, 1);
    typename Index::translation_type k;

    // setup highest translation index for e=(1,1), p=2
//    typename Index::translation_type k(0, frame1d_->Nablamax(j[1]));
            int patch=extensions_[num_logical_patches()-1][0];
            int target_patch=extensions_[num_logical_patches()-1][1];
            if(corners_[target_patch][0]==corners_[patch][0]){ //extension in y-direction
                    if(corners_[target_patch][1]<corners_[patch][1]){ //extension from north to south
                        
                            k[0]= frame1d_->Nablamax(j[0]) ;
                            k[1]= frame1d_->Nablamin(); 
                    }
                    else{
                         //extension from south to north
                         
                            k[0]= frame1d_->Nablamax(j[0]);
                            k[1]= frame1d_->Nablamax(j[1]); 
                    }
            }
            else{ //extension in x-direction
                    if(corners_[target_patch][0]<corners_[patch][0]){ //extension from east to west
                        
                            k[0]= frame1d_->Nablamin(); 
                            k[1]= frame1d_->Nablamax(j[1]);
                    }
                    else{
                         //extension from west to east
                        
                            k[0]= frame1d_->Nablamax(j[0]);
                            k[1]= frame1d_->Nablamax(j[1]);
                    }
            }
//    cout << "in last_quarklet level_type" << endl;
    if (number==-1)
    {
        level_type jdiff;
        jdiff[0]= j[0]-j0_[0], jdiff[1]= j[1]-j0_[1];
        int level = jdiff.number();
        int altnumber = p.number()* Nablasize_+last_wavelet_numbers[level];
        return Index(p,j, e, num_real_patches()+num_logical_patches()-1, k, altnumber, this);
    }
    else
    
       return Index(p,j, e, num_real_patches()+num_logical_patches()-1, k, number, this);
  }
  
  template <class IFRAME, int NPATCHES>
  typename DomainFrame<IFRAME, NPATCHES>::Index
  DomainFrame<IFRAME, NPATCHES>::last_quarklet(const int& levelsum, const polynomial_type& p, const int& number) const
  {
    assert(levelsum >= (int) multi_degree(j0_));
    
    typename Index::type_type e(1, 1);
    typename Index::level_type j(levelsum - j0_[0],  j0_[0]);
    typename Index::translation_type k;

    // setup highest translation index for e=(1,1), p=2
//    typename Index::translation_type k(0, frame1d_->Nablamax(j[1]));
            int patch=extensions_[num_logical_patches()-1][0];
            int target_patch=extensions_[num_logical_patches()-1][1];
            if(corners_[target_patch][0]==corners_[patch][0]){ //extension in y-direction
                    if(corners_[target_patch][1]<corners_[patch][1]){ //extension from north to south
                        
                            k[0]= frame1d_->Nablamax(j[0]) ;
                            k[1]= frame1d_->Nablamin(); 
                    }
                    else{
                         //extension from south to north
                         
                            k[0]= frame1d_->Nablamax(j[0]);
                            k[1]= frame1d_->Nablamax(j[1]); 
                    }
            }
            else{ //extension in x-direction
                    if(corners_[target_patch][0]<corners_[patch][0]){ //extension from east to west
                        
                            k[0]= frame1d_->Nablamin(); 
                            k[1]= frame1d_->Nablamax(j[1]);
                    }
                    else{
                         //extension from west to east
                        
                            k[0]= frame1d_->Nablamax(j[0]);
                            k[1]= frame1d_->Nablamax(j[1]);
                    }
            }
//    cout << "in last_quarklet" << endl;
    if (number==-1)
    {
//        cout << "altnumber" << endl;
//        cout << "j: " << j << ", j.number" << j.number() << endl;
//        cout << "j0_: " << j0_ <<", j0_.number" << j0_.number() << endl;
        level_type jdiff;
        jdiff[0]= j[0]-j0_[0], jdiff[1]= j[1]-j0_[1];
        int level = jdiff.number();
//        cout << "level: " << level << endl;
        int altnumber = p.number()* Nablasize_+last_wavelet_numbers[level];
        return Index(p,j, e, num_real_patches()+num_logical_patches()-1, k, altnumber, this);
    }
    else{
//        cout << "number" << endl;
       return Index(p,j, e, num_real_patches()+num_logical_patches()-1, k, number, this); 
    }
  }
  
  template <class IFRAME, int NPATCHES>
  void
  DomainFrame<IFRAME, NPATCHES>::support(const int& lambda_num, Support& supp) const
  {

     if (precomputed_supports_){
         supp=all_supports_[lambda_num];
     }
     else {
         cout << "DomainFrame<IFRAME>::support(): precompute supports first!" << endl;
         abort();      
     }
  }


  template <class IFRAME, int NPATCHES>
  void
  DomainFrame<IFRAME, NPATCHES>::support(const Index& lambda, Support& supp) const
  {
//      cout<<"computing support"<<endl;

     if (precomputed_supports_){
         supp=all_supports_[lambda.number()];
     }
     else{
//         cout<<"computing new support"<<endl;
	  supp.j[0] = lambda.j()[0]+lambda.e()[0];
          supp.j[1] = lambda.j()[1]+lambda.e()[1];
	  
          if(lambda.patch()<num_real_patches()){ //lambda lives on one patch
              frames(lambda.patch(),0)->support(typename IFRAME::Index(lambda.p()[0],
                                                     lambda.j()[0],
						     lambda.e()[0],
						     lambda.k()[0],
						     frame1d_),
			      supp.xmin[lambda.patch()],
			      supp.xmax[lambda.patch()]);
              
              frames(lambda.patch(),1)->support(typename IFRAME::Index(lambda.p()[1],
                                                     lambda.j()[1],
                                                     lambda.e()[1],
						     lambda.k()[1],
						     frame1d_),
			      supp.ymin[lambda.patch()],
			      supp.ymax[lambda.patch()]);
              for (int i1=0;i1<num_real_patches();i1++){
                if(i1!=lambda.patch()){
                    supp.xmin[i1] = -1;
                }
              }
          }
          else{ //lambda has been extended and lives on two patches
            //int dir;
            int patch=extensions_[lambda.patch()-num_real_patches()][0];
            int target_patch=extensions_[lambda.patch()-num_real_patches()][1];
            if(corners_[target_patch][0]==corners_[patch][0]){ //extension in y-direction
                //dir=0; //direction of the interface (frame))
                //support on interface
                frames(lambda.patch(),0)->support(typename IFRAME::Index(lambda.p()[0],
                                                     lambda.j()[0],
						     lambda.e()[0],
						     lambda.k()[0],
						     frame1d_),
			      supp.xmin[patch],
			      supp.xmax[patch]);	                            
                supp.xmin[target_patch]=supp.xmin[patch];
                supp.xmax[target_patch]=supp.xmax[patch];
                        
                //support on mother patch
                frames(patch,1)->support(typename IFRAME::Index(lambda.p()[1],
                                                     lambda.j()[1],
						     lambda.e()[1],
						     lambda.k()[1],
						     frame1d_),
			      supp.ymin[patch],
			      supp.ymax[patch]);

                //reflect support
                supp.ymax[target_patch]=(1<<supp.j[1])-supp.ymin[patch];
                supp.ymin[target_patch]=(1<<supp.j[1])-supp.ymax[patch];

            }
            else{ //extension in x-direction
                //dir=1;
                //support on interface
                frames(lambda.patch(),1)->support(typename IFRAME::Index(lambda.p()[1],
                                                     lambda.j()[1],
						     lambda.e()[1],
						     lambda.k()[1],
						     frame1d_),
			      supp.ymin[patch],
			      supp.ymax[patch]);	                            
                supp.ymin[target_patch]=supp.ymin[patch];
                supp.ymax[target_patch]=supp.ymax[patch];
                        
                //support on mother patch
                frames(patch,0)->support(typename IFRAME::Index(lambda.p()[0],
                                                     lambda.j()[0],
						     lambda.e()[0],
						     lambda.k()[0],
						     frame1d_),
			      supp.xmin[patch],
			      supp.xmax[patch]);

                //reflect support
                supp.xmax[target_patch]=(1<<supp.j[0])-supp.xmin[patch];
                supp.xmin[target_patch]=(1<<supp.j[0])-supp.xmax[patch];
            }
            
            //no support on all other patches
            for (int i1=0;i1<num_real_patches();i1++){
                if(i1!=patch && i1!=target_patch){
                    supp.xmin[i1] = -1;
                }
              }
            
            
          }
          
	  
     }

	
 
  }
  
  template <class IFRAME, int NPATCHES>
  Array1D<SampledMapping<2> >
  DomainFrame<IFRAME, NPATCHES>::evaluate
  (const typename DomainFrame<IFRAME, NPATCHES>::Index& lambda,
   const int resolution) const
  {
    Array1D<SampledMapping<2> > r(num_real_patches());

    typedef typename DomainFrame<IFRAME, NPATCHES>::Index Index;

    typename Index::type_type zero;
    
      
      FixedArray1D<Array1D<double>,2> values;
      values[0].resize((1<<resolution)+1); // values in x-direction
      values[1].resize((1<<resolution)+1); // values in y-direction

      if(lambda.patch()<num_real_patches()){ //lambda lives on one patch
         values[0] = frames(lambda.patch(),0)->evaluate(typename IFRAME::Index(lambda.p()[0],
                                                              lambda.j()[0],
							      lambda.e()[0],
							      lambda.k()[0],
							      frame1d_),
				       resolution).values();
 	values[1] = frames(lambda.patch(),1)->evaluate(typename IFRAME::Index(lambda.p()[1],
                                                              lambda.j()[1],
							      lambda.e()[1],
							      lambda.k()[1],
							      frame1d_),
				       resolution).values();
 	r[lambda.patch()] = SampledMapping<2>(Point<2>(corners_[lambda.patch()][0],corners_[lambda.patch()][1]), Point<2>(corners_[lambda.patch()][0]+1,corners_[lambda.patch()][1]+1), values);
	
        //zero evaluation on all other patches
 	for (int i = 0; i <= 1<<resolution; i++) {
 	  values[0][i] = values[1][i] = 0;
 	}
        for (int i1=0;i1<num_real_patches();i1++){
            if(i1!=lambda.patch()){
                r[i1] = SampledMapping<2>(Point<2>(corners_[i1][0],corners_[i1][1]), Point<2>(corners_[i1][0]+1,corners_[i1][1]+1), values);
            }
        }
 	

      }
      else{ //lambda has been extended and lives on two patches
          cout<<"evaluate quarklet on interface"<<endl;
          int dir;
            int patch=extensions_[lambda.patch()-num_real_patches()][0];
            int target_patch=extensions_[lambda.patch()-num_real_patches()][1];
            if(corners_[target_patch][0]==corners_[patch][0]){ //extension in y-direction
                dir=0; //direction of the interface (frame))
            }
            else{ //extension in x-direction
                dir=1;
            }     
                //evaluate on interface
                values[dir] = frames(lambda.patch(),dir)->evaluate(typename IFRAME::Index(lambda.p()[dir],
                                                              lambda.j()[dir],
							      lambda.e()[dir],
							      lambda.k()[dir],
							      frame1d_),
				       resolution).values();
                //evaluate on mother patch
                values[1-dir] = frames(patch,1-dir)->evaluate(typename IFRAME::Index(lambda.p()[1-dir],
                                                              lambda.j()[1-dir],
							      lambda.e()[1-dir],
							      lambda.k()[1-dir],
							      frame1d_),
				       resolution).values();

                //evaluation on mother patch
                r[patch] = SampledMapping<2>(Point<2>(corners_[patch][0],corners_[patch][1]), Point<2>(corners_[patch][0]+1,corners_[patch][1]+1), values);
                //reflection and evaluation on target patch
                for (int i = 0; i <= (1<<resolution)/2; i++){           
                    double temp = values[1-dir][i];
                    values[1-dir][i] = values[1-dir][(1<<resolution)-i];
                    values[1-dir][(1<<resolution)-i] = temp;
                }
                r[target_patch] = SampledMapping<2>(Point<2>(corners_[target_patch][0],corners_[target_patch][1]), Point<2>(corners_[target_patch][0]+1,corners_[target_patch][1]+1), values);
                
                //zero evaluation on all other patches
                for (int i = 0; i <= 1<<resolution; i++) {
                    values[0][i] = values[1][i] = 0;
                }
                for (int i1=0;i1<num_real_patches();i1++){
                    if(i1!=patch && i1!=target_patch){
                        r[i1] = SampledMapping<2>(Point<2>(corners_[i1][0],corners_[i1][1]), Point<2>(corners_[i1][0]+1,corners_[i1][1]+1), values);
                    }
                }
              
      }    
    
    return r;
  }
  
  template <class IFRAME, int NPATCHES>
  Array1D<SampledMapping<2> >
  DomainFrame<IFRAME, NPATCHES>::evaluate
  (const InfiniteVector<double, typename DomainFrame<IFRAME, NPATCHES>::Index>& coeffs,
   const int resolution) const
  {
    // first prepare a plot of the zero function
    FixedArray1D<Array1D<double>,2> values;
    values[0].resize((1<<resolution)+1);
    values[1].resize((1<<resolution)+1);
    for (int i = 0; i <= 1<<resolution; i++) {
      values[0][i] = values[1][i] = 0;
    }
    Array1D<SampledMapping<2> > result(num_real_patches());
    for (int i1=0;i1<num_real_patches();i1++){
                result[i1] = SampledMapping<2>(Point<2>(corners_[i1][0],corners_[i1][1]), Point<2>(corners_[i1][0]+1,corners_[i1][1]+1), values);
    }
          
    // add all plots of the single functions
    typedef typename DomainFrame<IFRAME, NPATCHES>::Index Index;
    for (typename InfiniteVector<double,Index>::const_iterator it(coeffs.begin()),
 	   itend(coeffs.end()); it != itend; ++it) {
      Array1D<SampledMapping<2> > temp(evaluate(it.index(), resolution));
      for (int i1=0;i1<num_real_patches();i1++){
           result[i1].add(*it, temp[i1]);     
      }
      
      
    }
    
    return result;
  }
  


  template <class IFRAME, int NPATCHES>
  void
  DomainFrame<IFRAME, NPATCHES>::setup_full_collection()
  {
    if (jmax_ == -1 || jmax_ < (int) multi_degree(j0())) {
      cout << "DomainFrame<IFRAME>::setup_full_collection(): specify a maximal level of resolution first!" << endl;
      abort();
    }   
    
    cout << "setting up collection of quarklet indices..." << endl;

    MultiIndex<int,2> level_it;
    level_it[0] = jmax_ - multi_degree(j0());
    const int numoflevels = level_it.number() + 1;
    level_it[0] = 0;
    first_wavelet_numbers.resize(numoflevels);
    last_wavelet_numbers.resize(numoflevels);
    
    typename Index::polynomial_type p(0,0), pmax(pmax_,0), pmin(0,0);
    set<Index> Lambda;
    int waveletlevel=0;
    Index ind = first_generator(j0_, p);
    for (int k = 0;; k++) {
//        full_collection[k] = ind;
        Lambda.insert(ind);
//        cout << ind << ", " << endl;
        if(ind==first_quarklet(ind.j(),pmin,k)){
                first_wavelet_numbers[waveletlevel]=k;
                
                cout << "first: " << k << endl;
        }
        if(ind==last_quarklet(ind.j(),pmin,k)){
            last_wavelet_numbers[waveletlevel]=k;
            cout << "last: " << k << endl;
            waveletlevel++;
        }
        if(ind==last_quarklet(jmax_, p,k)){
            if(multi_degree(p)==0){
                cout << "Nablasize wird gesetzt auf " << k+1 << endl;
                Nablasize_=k+1;
            }
            if(p==pmax)
                break;
            
            ++p;
            ind=first_generator(j0_, p);
        }
            else
                ++ind;
    }
    cout<<"precomputing supports"<<endl;
    full_collection.resize(Lambda.size());
    all_supports_.resize(Lambda.size());
    int i = 0;
    for (typename set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it, i++){
        full_collection[i] = *it;
        Support supp;
        support(*it,supp);
        all_supports_[i] = supp;
        
//        cout << *it << ", " << (*it).number() << endl;
    }
    precomputed_supports_ = true;
    
    
    
//    cout << "Nablasize in setup: " << Nablasize_ << endl;
    
    cout << "done setting up collection of quarklet indices and precompute supports..." << endl;
    cout << "total degrees of freedom between j0_ = " << j0_ << " and (jmax_= " << jmax_ << ", pmax_= " << pmax_ << ") is " << full_collection.size() << endl;
    setup_full_collection_ = true;
    
    

  }
  
  



}
