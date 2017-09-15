// implementation for ldomain_frame.h

#include <cmath>
#include <time.h>
#include <iostream>

#include "ldomain_frame_index.h"

using std::cout;
using std::endl;

namespace WaveletTL
{
  template <class IFRAME>
  LDomainFrame<IFRAME>::LDomainFrame(const IntervalFrame* frame1d, const IntervalFrame* frame1d_11, 
          const IntervalFrame* frame1d_01, const IntervalFrame* frame1d_10)
    : frame1d_(frame1d), frame1d_11_(frame1d_11), frame1d_01_(frame1d_01), frame1d_10_(frame1d_10)
  {
      j0_[0] = frame1d_11_->j0();
      j0_[1] = frame1d_11_->j0();
  }

//  template <class IFRAME>
//  LDomainFrame<IFRAME>::LDomainFrame(const IntervalFrame& frame1d)
//    : frame1d_(frame1d), frame1d_11_(true, true), frame1d_01_(false, true)
//  {
////    j0_[0] = frame1d.j0();
////    j0_[1] = frame1d.j0();
//    j0_[0] = frame1d_11_.j0();
//    j0_[1] = frame1d_11_.j0();
//  }

  template <class IFRAME>
  const int
  LDomainFrame<IFRAME>::Deltasize(const int j) const {
    const unsigned int Deltaj = frame1d_->Deltasize(j);
    return 3*(Deltaj-2)*(Deltaj-2)+2*(Deltaj-2);
  }

//  template <class IBASIS>
//  const int
//  LDomainBasis<IBASIS>::Nabla01size(const int j) const {
//    return (3*basis1d().Deltasize(j)-5)*(1<<j);
//  }
//  
//  template <class IBASIS>
//  const int
//  LDomainBasis<IBASIS>::Nabla10size(const int j) const {
//    return (3*basis1d().Deltasize(j)-5)*(1<<j);
//  }
//  
//  template <class IBASIS>
//  const int
//  LDomainBasis<IBASIS>::Nabla11size(const int j) const {
//    return 3*(1<<(2*j));
//  }

  
  template <class IFRAME>
  typename LDomainFrame<IFRAME>::Index
  LDomainFrame<IFRAME>::first_generator(const level_type& j, const polynomial_type& p) const
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

  template <class IFRAME>
  typename LDomainFrame<IFRAME>::Index
  LDomainFrame<IFRAME>::last_generator(const level_type& j, const polynomial_type& p) const
  {
    assert(j >= j0_);

    typename Index::type_type e;

    // setup highest translation index for e=(0,0), p=4
    typename Index::translation_type k(frame1d_->DeltaLmin(), frame1d_->DeltaRmax(j[1])-1);
    
    return Index(p,j, e, 4, k, p.number()* Nablasize_+Deltasize(j[0])-1, this);
  }

  template <class IFRAME>
  typename LDomainFrame<IFRAME>::Index
  LDomainFrame<IFRAME>::first_quarklet(const level_type& j, const polynomial_type& p, const int& number) const
  {
    assert(j >= j0());

    typename Index::type_type e;
    typename Index::translation_type k;
    
    bool sofar_only_generators = true;
    
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

//  template <class IFRAME>
//  typename LDomainFrame<IFRAME>::Index
//  LDomainFrame<IFRAME>::first_quarklet(const level_type& j, const type_type& e, const polynomial_type& p) const
//  {
//    assert(j >= j0());
//    
////    typename Index::type_type e(ewish);
//    
//    // setup lowest translation index appropriately 
//    typename Index::translation_type k;
//    const int ecode(e[0]+2*e[1]);
//    if (ecode == 0) {
//      // e = (0,0)
//      k[0] = k[1] = frame1d_->DeltaLmin()+1;
//    } else {
//      if (ecode == 1) {
//	// e = (1,0)
//	k[0] = frame1d_->Nablamin()+1;
//	k[1] = frame1d_->DeltaLmin()+1;
//      } else {
//	if (ecode == 2) {
//	  // e = (0,1)
//	  k[0] = frame1d_->DeltaLmin()+1;
//	  k[1] = frame1d_->Nablamin()+1;
//	} else {
//	  // e = (1,1)
//	  k[0] = k[1] = frame1d_->Nablamin()+1;
//	}
//      }
//    }
//    
//    return Index(p,j, e, 0, k, 0, this);
//  }

  template <class IFRAME>
  typename LDomainFrame<IFRAME>::Index
  LDomainFrame<IFRAME>::last_quarklet(const level_type& j, const polynomial_type& p, const int& number) const
  {
    assert(j >= j0());
    
    typename Index::type_type e(1, 1);

    // setup highest translation index for e=(1,1), p=2
    typename Index::translation_type k(0, frame1d_->Nablamax(j[1]));
//    cout << "in last_quarklet level_type" << endl;
    if (number==-1)
    {
        level_type jdiff;
        jdiff[0]= j[0]-j0_[0], jdiff[1]= j[1]-j0_[1];
        int level = jdiff.number();
        int altnumber = p.number()* Nablasize_+last_wavelet_numbers[level];
        return Index(p,j, e, 4, k, altnumber, this);
    }
    else
    
       return Index(p,j, e, 4, k, number, this);
  }
  
  template <class IFRAME>
  typename LDomainFrame<IFRAME>::Index
  LDomainFrame<IFRAME>::last_quarklet(const int& levelsum, const polynomial_type& p, const int& number) const
  {
    assert(levelsum >= (int) multi_degree(j0_));
    
    typename Index::type_type e(1, 1);
    typename Index::level_type j(levelsum - j0_[0],  j0_[0]);

    // setup highest translation index for e=(1,1), p=2
    typename Index::translation_type k(0, frame1d_->Nablamax(j[1]));
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
        return Index(p,j, e, 4, k, altnumber, this);
    }
    else{
//        cout << "number" << endl;
       return Index(p,j, e, 4, k, number, this); 
    }
  }
  
  template <class IFRAME>
  void
  LDomainFrame<IFRAME>::support(const int& lambda_num, Support& supp) const
  {

     if (precomputed_supports_){
         supp=all_supports_[lambda_num];
     }
     else {
         cout << "LDomainFrame<IFRAME>::support(): precompute supports first!" << endl;
         abort();      
     }
  }


  template <class IFRAME>
  void
  LDomainFrame<IFRAME>::support(const Index& lambda, Support& supp) const
  {
//      int lam_num=lambda.number();
////     check whether the supp(psi_lambda) already exists in the cache
//    typename SupportCache::iterator supp_lb(supp_cache.lower_bound(lam_num));
//    typename SupportCache::iterator supp_it(supp_lb);
//    if (supp_lb == supp_cache.end() ||
//	supp_cache.key_comp()(lam_num, supp_lb->first))
//      {
////        cout << "Neuberechnung" << endl;
////	 compute supp(psi_lambda) and insert it into the cache
//	typedef typename SupportCache::value_type value_type;

     if (precomputed_supports_){
         supp=all_supports_[lambda.number()];
     }
     else{
	  
	  supp.j[0] = lambda.j()[0]+lambda.e()[0];
          supp.j[1] = lambda.j()[1]+lambda.e()[1];
	  
	  switch (lambda.patch()) {
	  case 0:
	    // psi_lambda completely lives on patch 0
	    
	    frame1d_11_->support(typename IFRAME::Index(lambda.p()[0],
                                                     lambda.j()[0],
						     lambda.e()[0],
						     lambda.k()[0],
						     frame1d_),
			      supp.xmin[0],
			      supp.xmax[0]);
	    
	    frame1d_01_->support(typename IFRAME::Index(lambda.p()[1],
                                                     lambda.j()[1],
                                                     lambda.e()[1],
						     lambda.k()[1],
						     frame1d_),
			      supp.ymin[0],
			      supp.ymax[0]);
	    
	    supp.xmin[1] = supp.xmin[2] = -1;
	    
	    break;
	  case 1:
	    // psi_lambda completely lives on patch 1
	    
	    frame1d_11_->support(typename IFRAME::Index(lambda.p()[0],
                                                     lambda.j()[0],
						     lambda.e()[0],
						     lambda.k()[0],
						     frame1d_),
			      supp.xmin[1],
			      supp.xmax[1]);
	    
	    frame1d_11_->support(typename IFRAME::Index(lambda.p()[1],
                                                     lambda.j()[1],
						     lambda.e()[1],
						     lambda.k()[1],
						     frame1d_),
			      supp.ymin[1],
			      supp.ymax[1]);
	    
	    supp.xmin[0] = supp.xmin[2] = -1;
	    
	    break;
	  case 2:
	    // psi_lambda completely lives on patch 2
	    
	    frame1d_01_->support(typename IFRAME::Index(lambda.p()[0],
                                                     lambda.j()[0],
						     lambda.e()[0],
						     lambda.k()[0],
						     frame1d_),
			      supp.xmin[2],
			      supp.xmax[2]);
	    
	    frame1d_11_->support(typename IFRAME::Index(lambda.p()[1],
                                                     lambda.j()[1],
						     lambda.e()[1],
						     lambda.k()[1],
						     frame1d_),
			      supp.ymin[2],
			      supp.ymax[2]);
	    
	    supp.xmin[0] = supp.xmin[1] = -1;
	    
	    break;
	  case 3:
	    // psi_lambda lives on patches 0 and 1
              
	    frame1d_11_->support(typename IFRAME::Index(lambda.p()[0],
                                                     lambda.j()[0],
						     lambda.e()[0],
						     lambda.k()[0],
						     frame1d_),
			      supp.xmin[0],
			      supp.xmax[0]);
	    
	    frame1d_01_->support(typename IFRAME::Index(lambda.p()[1],
                                                     lambda.j()[1],
						     lambda.e()[1],
						     (lambda.e()[1]==0?frame1d_->DeltaLmin()
                                                     :frame1d_->Nablamin()),
						     frame1d_),
			      supp.ymin[0],
			      supp.ymax[0]);
	    
	    frame1d_11_->support(typename IFRAME::Index(lambda.p()[0],
                                                     lambda.j()[0],
						     lambda.e()[0],
						     lambda.k()[0],
						     frame1d_),
			      supp.xmin[1],
			      supp.xmax[1]);
	    
	    //mirrored, therefore we need free boundary conditions in the north
            frame1d_10_->support(typename IFRAME::Index(lambda.p()[1],
                                                     lambda.j()[1],
						     lambda.e()[1],
						     (lambda.e()[1]==0?frame1d_->DeltaRmax(lambda.j()[1])
                                                     //!!Attention: only works, if only the most outter quark is mirrored
                                                     :frame1d_->Nablamax(lambda.j()[1])),
						     frame1d_),
			      supp.ymin[1],
			      supp.ymax[1]);
	    
	    supp.xmin[2] = -1;
	    
	    break;
	  case 4:
	    // psi_lambda lives on patches 1 and 2
	    
	    frame1d_01_->support(typename IFRAME::Index(lambda.p()[0],
                                                     lambda.j()[0],
						     lambda.e()[0],
						     (lambda.e()[0]==0?frame1d_->DeltaRmax(lambda.j()[0])
                                                     :frame1d_->Nablamax(lambda.j()[0])),
						     frame1d_),
			      supp.xmin[1],
			      supp.xmax[1]);
	    
	    frame1d_11_->support(typename IFRAME::Index(lambda.p()[1],
                                                     lambda.j()[1],
						     lambda.e()[1],
						     lambda.k()[1],
						     frame1d_),
			      supp.ymin[1],
			      supp.ymax[1]);
	    
	    //mirrored, therefore we need free boundary conditions in the east
            frame1d_10_->support(typename IFRAME::Index(lambda.p()[0],
                                                     lambda.j()[0],
						     lambda.e()[0],
						     (lambda.e()[0]==0?frame1d_->DeltaLmin()
                                                     :frame1d_->Nablamin()),
						     frame1d_),
			      supp.xmin[2],
			      supp.xmax[2]);
	    
	    frame1d_11_->support(typename IFRAME::Index(lambda.p()[1],
                                                     lambda.j()[1],
						     lambda.e()[1],
						     lambda.k()[1],
						     frame1d_),
			      supp.ymin[2],
			      supp.ymax[2]);
	    
	    supp.xmin[0] = -1;
	    
	    break;
	  }
     }

	
//	supp_it = supp_cache.insert(supp_lb, value_type(lam_num, supp));
//
//
//      }
//    else
//      {
//        
////        typename SupportCache::iterator it(supp_cache.begin());
////        typename SupportCache::iterator itend(supp_cache.end());
////        for(typename SupportCache::iterator it = supp_cache.begin();it!=itend;it++){
////            cout << endl << it->first << endl;
////            cout << "patch 0: [" << it->second.xmin[0]  <<" , " <<it->second.xmax[0] <<"]x["<<it->second.ymin[0] <<" , "<<it->second.ymax[0] <<"]"<<  endl;
////            cout << "patch 1: [" << it->second.xmin[1]  <<" , " <<it->second.xmax[1] <<"]x["<<it->second.ymin[1] <<" , "<<it->second.ymax[1] <<"]"<<  endl;
////            cout << "patch 2: [" << it->second.xmin[2]  <<" , " <<it->second.xmax[2] <<"]x["<<it->second.ymin[2] <<" , "<<it->second.ymax[2] <<"]"<<  endl;
////            
////        }
//	// cache hit, copy the precomputed support
////        cout << "Cache hitted" << supp_it->first << endl;
//        
//  	supp.j[0] = supp_it->second.j[0];
//        supp.j[1] = supp_it->second.j[1];
//  	for (unsigned int i = 0; i < 3; i++) {
//  	  supp.xmin[i] = supp_it->second.xmin[i];
//  	  supp.xmax[i] = supp_it->second.xmax[i];
//  	  supp.ymin[i] = supp_it->second.ymin[i];
//  	  supp.ymax[i] = supp_it->second.ymax[i];
//  	}
//	
//
//      }  
  }
  
  template <class IFRAME>
  Array1D<SampledMapping<2> >
  LDomainFrame<IFRAME>::evaluate
  (const typename LDomainFrame<IFRAME>::Index& lambda,
   const int resolution) const
  {
    Array1D<SampledMapping<2> > r(3);

    typedef typename LDomainFrame<IFRAME>::Index Index;

    typename Index::type_type zero;
    
      
      FixedArray1D<Array1D<double>,2> values;
      values[0].resize((1<<resolution)+1); // values in x-direction
      values[1].resize((1<<resolution)+1); // values in y-direction

      switch (lambda.patch()) {
      case 0:
 	// psi_lambda completely lives on patch 0
 	values[0] = frame1d_11_->evaluate(typename IFRAME::Index(lambda.p()[0],
                                                              lambda.j()[0],
							      lambda.e()[0],
							      lambda.k()[0],
							      frame1d_),
				       resolution).values();
 	values[1] = frame1d_01_->evaluate(typename IFRAME::Index(lambda.p()[1],
                                                              lambda.j()[1],
							      lambda.e()[1],
							      lambda.k()[1],
							      frame1d_),
				       resolution).values();
 	r[0] = SampledMapping<2>(Point<2>(-1, 0), Point<2>(0,1), values);
	
 	for (int i = 0; i <= 1<<resolution; i++) {
 	  values[0][i] = values[1][i] = 0;
 	}
 	r[1] = SampledMapping<2>(Point<2>(-1,-1), Point<2>(0,0), values);
 	r[2] = SampledMapping<2>(Point<2>( 0,-1), Point<2>(1,0), values);
 	break;
      case 1:
 	// psi_lambda completely lives on patch 1
 	values[0] = frame1d_11_->evaluate(typename IFRAME::Index(lambda.p()[0],
                                                              lambda.j()[0],
							      lambda.e()[0],
							      lambda.k()[0],
							      frame1d_),
				       resolution).values();
 	values[1] = frame1d_11_->evaluate(typename IFRAME::Index(lambda.p()[1],
                                                              lambda.j()[1],
							      lambda.e()[1],
							      lambda.k()[1],
							      frame1d_),
				       resolution).values();
 	r[1] = SampledMapping<2>(Point<2>(-1,-1), Point<2>(0,0), values);

 	for (int i = 0; i <= 1<<resolution; i++) {
 	  values[0][i] = values[1][i] = 0;
 	}
 	r[0] = SampledMapping<2>(Point<2>(-1, 0), Point<2>(0,1), values);
 	r[2] = SampledMapping<2>(Point<2>( 0,-1), Point<2>(1,0), values);
	break;
      case 2:
 	// psi_lambda completely lives on patch 2
 	values[0] = frame1d_01_->evaluate(typename IFRAME::Index(lambda.p()[0],
                                                              lambda.j()[0],
							      lambda.e()[0],
							      lambda.k()[0],
							      frame1d_),
				       resolution).values();
 	values[1] = frame1d_11_->evaluate(typename IFRAME::Index(lambda.p()[1],
                                                              lambda.j()[1],
							      lambda.e()[1],
							      lambda.k()[1],
							      frame1d_),
				       resolution).values();
 	r[2] = SampledMapping<2>(Point<2>( 0,-1), Point<2>(1,0), values);
        
        //Test
//        cout << "Test" << endl;
//        frame1d_11()->evaluate(typename IFRAME::Index(lambda.p()[1],
//                                                              lambda.j()[1],
//							      lambda.e()[1],
//							      lambda.k()[1],
//							      frame1d()),
//				       resolution).matlab_output(cout);
//        //TESTEND

 	for (int i = 0; i <= 1<<resolution; i++) {
 	  values[0][i] = values[1][i] = 0;
 	}
 	r[0] = SampledMapping<2>(Point<2>(-1, 0), Point<2>(0,1), values);
 	r[1] = SampledMapping<2>(Point<2>(-1,-1), Point<2>(0,0), values);
 	break;
      case 3:
  	// psi_lambda lives on patches 0 and 1
  	values[0] = frame1d_11_->evaluate(typename IFRAME::Index(lambda.p()[0],
                                                              lambda.j()[0],
							      lambda.e()[0],
							      lambda.k()[0],
							      frame1d_),
				       resolution).values();
 	values[1] = frame1d_01_->evaluate(typename IFRAME::Index(lambda.p()[1],
                                                              lambda.j()[1],
							      lambda.e()[1],
							      lambda.k()[1],
							      frame1d_),
				       resolution).values();
// 	for (int i = 0; i <= 1<<resolution; i++) values[1][i] *= M_SQRT1_2;
 	r[0] = SampledMapping<2>(Point<2>(-1, 0), Point<2>(0,1), values);


 	for (int i = 0; i <= (1<<resolution)/2; i++){           // values[1][i] *= M_SQRT1_2;
            double temp = values[1][i];
            values[1][i] = values[1][(1<<resolution)-i];
            values[1][(1<<resolution)-i] = temp;
        }
 	r[1] = SampledMapping<2>(Point<2>(-1,-1), Point<2>(0,0), values);

 	for (int i = 0; i <= 1<<resolution; i++) {
	  values[0][i] = values[1][i] = 0;
	}
	r[2] = SampledMapping<2>(Point<2>( 0,-1), Point<2>(1,0), values);
	break;
      case 4:
 	// psi_lambda lives on patches 1 and 2

 	values[0] = frame1d_01_->evaluate(typename IFRAME::Index(lambda.p()[0],
                                                              lambda.j()[0],
							      lambda.e()[0],
							      lambda.k()[0],
							      frame1d_),
				       resolution).values();
 	values[1] = frame1d_11_->evaluate(typename IFRAME::Index(lambda.p()[1],
                                                              lambda.j()[1],
							      lambda.e()[1],
							      lambda.k()[1],
							      frame1d_),
				       resolution).values();
        
//	for (int i = 0; i <= 1<<resolution; i++) values[1][i] *= M_SQRT1_2;
	r[2] = SampledMapping<2>(Point<2>( 0,-1), Point<2>(1,0), values);
        
        for (int i = 0; i <= (1<<resolution)/2; i++){           // values[1][i] *= M_SQRT1_2;
            double temp = values[0][i];
            values[0][i] = values[0][(1<<resolution)-i];
            values[0][(1<<resolution)-i] = temp;
        }
        
        r[1] = SampledMapping<2>(Point<2>(-1,-1), Point<2>(0,0), values);
	
	for (int i = 0; i <= 1<<resolution; i++) {
	  values[0][i] = values[1][i] = 0;
	}
	r[0] = SampledMapping<2>(Point<2>(-1, 0), Point<2>(0,1), values);
	break;
      }
     
    
    return r;
  }
  
  template <class IFRAME>
  Array1D<SampledMapping<2> >
  LDomainFrame<IFRAME>::evaluate
  (const InfiniteVector<double, typename LDomainFrame<IFRAME>::Index>& coeffs,
   const int resolution) const
  {
    // first prepare a plot of the zero function
    FixedArray1D<Array1D<double>,2> values;
    values[0].resize((1<<resolution)+1);
    values[1].resize((1<<resolution)+1);
    for (int i = 0; i <= 1<<resolution; i++) {
      values[0][i] = values[1][i] = 0;
    }
    Array1D<SampledMapping<2> > result(3);
    result[0] = SampledMapping<2>(Point<2>(-1, 0), Point<2>(0,1), values);
    result[1] = SampledMapping<2>(Point<2>(-1,-1), Point<2>(0,0), values);
    result[2] = SampledMapping<2>(Point<2>( 0,-1), Point<2>(1,0), values);
          
    // add all plots of the single functions
    typedef typename LDomainFrame<IFRAME>::Index Index;
    for (typename InfiniteVector<double,Index>::const_iterator it(coeffs.begin()),
 	   itend(coeffs.end()); it != itend; ++it) {
      Array1D<SampledMapping<2> > temp(evaluate(it.index(), resolution));
      result[0].add(*it, temp[0]);
      result[1].add(*it, temp[1]);
      result[2].add(*it, temp[2]);
    }
    
    return result;
  }
  


  template <class IFRAME>
  void
  LDomainFrame<IFRAME>::setup_full_collection()
  {
    if (jmax_ == -1 || jmax_ < (int) multi_degree(j0())) {
      cout << "LDomainFrame<IFRAME>::setup_full_collection(): specify a maximal level of resolution first!" << endl;
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
//        cout << ind << ", " << ind.number() << endl;
        if(ind==first_quarklet(ind.j(),pmin,k)){
                first_wavelet_numbers[waveletlevel]=k;
                
//                cout << "first: " << k << endl;
        }
        if(ind==last_quarklet(ind.j(),pmin,k)){
            last_wavelet_numbers[waveletlevel]=k;
//            cout << "last: " << k << endl;
            waveletlevel++;
        }
        if(ind==last_quarklet(jmax_, p,k)){
            if(multi_degree(p)==0){
//                cout << "Nablasize wird gesetzt auf " << k+1 << endl;
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
