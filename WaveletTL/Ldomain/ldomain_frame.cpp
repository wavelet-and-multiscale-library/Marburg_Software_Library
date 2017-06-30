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
  LDomainFrame<IFRAME>::LDomainFrame()
    : frame1d_(false, false)
  {

  }

  template <class IFRAME>
  LDomainFrame<IFRAME>::LDomainFrame(const IntervalFrame& frame1d)
    : frame1d_(frame1d)
  {

    j0_[0] = frame1d.j0();
    j0_[1] = frame1d.j0();
  }


  
  template <class IFRAME>
  typename LDomainFrame<IFRAME>::Index
  LDomainFrame<IFRAME>::first_generator(const level_type& j, const polynomial_type& p, const int& number) const
  {
    assert(j >= j0());

    typename Index::type_type e;

    // setup lowest translation index for e=(0,0), p=0
    typename Index::translation_type k(frame1d().DeltaLmin()+1, frame1d().DeltaLmin()+1);
    
    return Index(p,j, e, 0, k, number, this);
  }

  template <class IFRAME>
  typename LDomainFrame<IFRAME>::Index
  LDomainFrame<IFRAME>::last_generator(const level_type& j, const polynomial_type& p) const
  {
    assert(j >= j0());

    typename Index::type_type e;

    // setup highest translation index for e=(0,0), p=4
    typename Index::translation_type k(0, frame1d().DeltaRmax(j)-1);
    
    return Index(p,j, e, 4, k, 0, this);
  }

  template <class IFRAME>
  typename LDomainFrame<IFRAME>::Index
  LDomainFrame<IFRAME>::first_quarklet(const level_type& j, const polynomial_type& p) const
  {
    assert(j >= j0());

    typename Index::type_type e(0, 1);

    // setup lowest translation index for e=(0,1), p=0
    typename Index::translation_type k(frame1d().DeltaLmin()+1, frame1d().Nablamin()+1);
    
    return Index(p,j, e, 0, k, 0, this);
  }

  template <class IFRAME>
  typename LDomainFrame<IFRAME>::Index
  LDomainFrame<IFRAME>::first_quarklet(const level_type& j, const type_type& e, const polynomial_type& p) const
  {
    assert(j >= j0());
    
//    typename Index::type_type e(ewish);
    
    // setup lowest translation index appropriately 
    typename Index::translation_type k;
    const int ecode(e[0]+2*e[1]);
    if (ecode == 0) {
      // e = (0,0)
      k[0] = k[1] = frame1d().DeltaLmin()+1;
    } else {
      if (ecode == 1) {
	// e = (1,0)
	k[0] = frame1d().Nablamin()+1;
	k[1] = frame1d().DeltaLmin()+1;
      } else {
	if (ecode == 2) {
	  // e = (0,1)
	  k[0] = frame1d().DeltaLmin()+1;
	  k[1] = frame1d().Nablamin()+1;
	} else {
	  // e = (1,1)
	  k[0] = k[1] = frame1d().Nablamin()+1;
	}
      }
    }
    
    return Index(p,j, e, 0, k, 0, this);
  }

  template <class IFRAME>
  typename LDomainFrame<IFRAME>::Index
  LDomainFrame<IFRAME>::last_quarklet(const level_type& j, const polynomial_type& p) const
  {
    assert(j >= j0());
    
    typename Index::type_type e(1, 1);

    // setup highest translation index for e=(1,1), p=2
    typename Index::translation_type k(0, frame1d().Nablamax(j[1]));
    
    return Index(p,j, e, 4, k, 0, this);
  }
  
  template <class IFRAME>
  typename LDomainFrame<IFRAME>::Index
  LDomainFrame<IFRAME>::last_quarklet(const int levelsum, const polynomial_type& p) const
  {
    assert(levelsum >= (int) multi_degree(j0()));
    
    typename Index::type_type e(1, 1);
    typename Index::level_type j(levelsum - j0_[0],  j0_[0]);

    // setup highest translation index for e=(1,1), p=2
    typename Index::translation_type k(0, frame1d().Nablamax(j[1]));
    
    return Index(p,j, e, 4, k, 0, this);
  }
  


  template <class IFRAME>
  void
  LDomainFrame<IFRAME>::support(const Index& lambda, Support& supp) const
  {
    // check whether the supp(psi_lambda) already exists in the cache
    typename SupportCache::iterator supp_lb(supp_cache.lower_bound(lambda));
    typename SupportCache::iterator supp_it(supp_lb);
    if (supp_lb == supp_cache.end() ||
	supp_cache.key_comp()(lambda, supp_lb->first))
      {
	// compute supp(psi_lambda) and insert it into the cache
	typedef typename SupportCache::value_type value_type;


	  
	  supp.j[0] = lambda.j()[0];
          supp.j[1] = lambda.j()[1];
	  
	  switch (lambda.patch()) {
	  case 0:
	    // psi_lambda completely lives on patch 0
	    
	    frame1d().support(typename IFRAME::Index(lambda.p()[0],
                                                     lambda.j()[0],
						     lambda.e()[0],
						     lambda.k()[0],
						     &frame1d()),
			      supp.xmin[0],
			      supp.xmax[0]);
	    
	    frame1d().support(typename IFRAME::Index(lambda.p()[1],
                                                     lambda.j()[1],
                                                     lambda.e()[1],
						     lambda.k()[1],
						     &frame1d()),
			      supp.ymin[0],
			      supp.ymax[0]);
	    
	    supp.xmin[1] = supp.xmin[2] = -1;
	    
	    break;
	  case 1:
	    // psi_lambda completely lives on patch 1
	    
	    frame1d().support(typename IFRAME::Index(lambda.p()[0],
                                                     lambda.j()[0],
						     lambda.e()[0],
						     lambda.k()[0],
						     &frame1d()),
			      supp.xmin[1],
			      supp.xmax[1]);
	    
	    frame1d().support(typename IFRAME::Index(lambda.p()[1],
                                                     lambda.j()[1],
						     lambda.e()[1],
						     lambda.k()[1],
						     &frame1d()),
			      supp.ymin[1],
			      supp.ymax[1]);
	    
	    supp.xmin[0] = supp.xmin[2] = -1;
	    
	    break;
	  case 2:
	    // psi_lambda completely lives on patch 2
	    
	    frame1d().support(typename IFRAME::Index(lambda.p()[0],
                                                     lambda.j()[0],
						     lambda.e()[0],
						     lambda.k()[0],
						     &frame1d()),
			      supp.xmin[2],
			      supp.xmax[2]);
	    
	    frame1d().support(typename IFRAME::Index(lambda.p()[1],
                                                     lambda.j()[1],
						     lambda.e()[1],
						     lambda.k()[1],
						     &frame1d()),
			      supp.ymin[2],
			      supp.ymax[2]);
	    
	    supp.xmin[0] = supp.xmin[1] = -1;
	    
	    break;
	  case 3:
	    // psi_lambda lives on patches 0 and 1
              
	    frame1d().support(typename IFRAME::Index(lambda.p()[0],
                                                     lambda.j()[0],
						     lambda.e()[0],
						     lambda.k()[0],
						     &frame1d()),
			      supp.xmin[0],
			      supp.xmax[0]);
	    
	    frame1d().support(typename IFRAME::Index(lambda.p()[1],
                                                     lambda.j()[1],
						     lambda.e()[1],
						     (lambda.e()[1]==0?frame1d().DeltaLmin()
                                                     :frame1d().Nablamin()),
						     &frame1d()),
			      supp.ymin[0],
			      supp.ymax[0]);
	    
	    frame1d().support(typename IFRAME::Index(lambda.p()[0],
                                                     lambda.j()[0],
						     lambda.e()[0],
						     lambda.k()[0],
						     &frame1d()),
			      supp.xmin[1],
			      supp.xmax[1]);
	    
	    frame1d().support(typename IFRAME::Index(lambda.p()[1],
                                                     lambda.j()[1],
						     lambda.e()[1],
						     (lambda.e()[1]==0?frame1d().DeltaRmax(lambda.j()[1])
                                                     :frame1d().Nablamax(lambda.j()[1])),
						     &frame1d()),
			      supp.ymin[1],
			      supp.ymax[1]);
	    
	    supp.xmin[2] = -1;
	    
	    break;
	  case 4:
	    // psi_lambda lives on patches 1 and 2
	    
	    frame1d().support(typename IFRAME::Index(lambda.p()[0],
                                                     lambda.j()[0],
						     lambda.e()[0],
						     (lambda.e()[0]==0?frame1d().DeltaRmax(lambda.j()[0])
                                                     :frame1d().Nablamax(lambda.j()[0])),
						     &frame1d()),
			      supp.xmin[1],
			      supp.xmax[1]);
	    
	    frame1d().support(typename IFRAME::Index(lambda.p()[1],
                                                     lambda.j()[1],
						     lambda.e()[1],
						     lambda.k()[1],
						     &frame1d()),
			      supp.ymin[1],
			      supp.ymax[1]);
	    
	    frame1d().support(typename IFRAME::Index(lambda.p()[0],
                                                     lambda.j()[0],
						     lambda.e()[0],
						     (lambda.e()[0]==0?frame1d().DeltaLmin()
                                                     :frame1d().Nablamin()),
						     &frame1d()),
			      supp.xmin[2],
			      supp.xmax[2]);
	    
	    frame1d().support(typename IFRAME::Index(lambda.p()[1],
                                                     lambda.j()[1],
						     lambda.e()[1],
						     lambda.k()[1],
						     &frame1d()),
			      supp.ymin[2],
			      supp.ymax[2]);
	    
	    supp.xmin[0] = -1;
	    
	    break;
	  }
        

	
	supp_it = supp_cache.insert(supp_lb, value_type(lambda, supp));


      }
    else
      {
	// cache hit, copy the precomputed support
  	supp.j[0] = supp_it->second.j[0];
        supp.j[1] = supp_it->second.j[1];
  	for (unsigned int i = 0; i < 3; i++) {
  	  supp.xmin[i] = supp_it->second.xmin[i];
  	  supp.xmax[i] = supp_it->second.xmax[i];
  	  supp.ymin[i] = supp_it->second.ymin[i];
  	  supp.ymax[i] = supp_it->second.ymax[i];
  	}
	

      }  
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
 	values[0] = frame1d().evaluate(typename IFRAME::Index(lambda.p()[0],
                                                              lambda.j()[0],
							      lambda.e()[0],
							      lambda.k()[0],
							      &frame1d()),
				       resolution).values();
 	values[1] = frame1d().evaluate(typename IFRAME::Index(lambda.p()[1],
                                                              lambda.j()[1],
							      lambda.e()[1],
							      lambda.k()[1],
							      &frame1d()),
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
 	values[0] = frame1d().evaluate(typename IFRAME::Index(lambda.p()[0],
                                                              lambda.j()[0],
							      lambda.e()[0],
							      lambda.k()[0],
							      &frame1d()),
				       resolution).values();
 	values[1] = frame1d().evaluate(typename IFRAME::Index(lambda.p()[1],
                                                              lambda.j()[1],
							      lambda.e()[1],
							      lambda.k()[1],
							      &frame1d()),
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
 	values[0] = frame1d().evaluate(typename IFRAME::Index(lambda.p()[0],
                                                              lambda.j()[0],
							      lambda.e()[0],
							      lambda.k()[0],
							      &frame1d()),
				       resolution).values();
 	values[1] = frame1d().evaluate(typename IFRAME::Index(lambda.p()[1],
                                                              lambda.j()[1],
							      lambda.e()[1],
							      lambda.k()[1],
							      &frame1d()),
				       resolution).values();
 	r[2] = SampledMapping<2>(Point<2>( 0,-1), Point<2>(1,0), values);

 	for (int i = 0; i <= 1<<resolution; i++) {
 	  values[0][i] = values[1][i] = 0;
 	}
 	r[0] = SampledMapping<2>(Point<2>(-1, 0), Point<2>(0,1), values);
 	r[1] = SampledMapping<2>(Point<2>(-1,-1), Point<2>(0,0), values);
 	break;
      case 3:
  	// psi_lambda lives on patches 0 and 1
  	values[0] = frame1d().evaluate(typename IFRAME::Index(lambda.p()[0],
                                                              lambda.j()[0],
							      lambda.e()[0],
							      lambda.k()[0],
							      &frame1d()),
				       resolution).values();
 	values[1] = frame1d().evaluate(typename IFRAME::Index(lambda.p()[1],
                                                              lambda.j()[1],
							      lambda.e()[1],
							      lambda.k()[1],
							      &frame1d()),
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

 	values[0] = frame1d().evaluate(typename IFRAME::Index(lambda.p()[0],
                                                              lambda.j()[0],
							      lambda.e()[0],
							      lambda.k()[0],
							      &frame1d()),
				       resolution).values();
 	values[1] = frame1d().evaluate(typename IFRAME::Index(lambda.p()[1],
                                                              lambda.j()[1],
							      lambda.e()[1],
							      lambda.k()[1],
							      &frame1d()),
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


    
    typename Index::polynomial_type p(0,0), pmax(pmax_,0);
    set<Index> Lambda;
    Index ind = first_generator(j0_, p);
    for (int k = 0;; k++) {
//        full_collection[k] = ind;
        Lambda.insert(ind);
//        cout << ind << ", " << ind.number() << endl;
 
        if(ind==last_quarklet(jmax_, p)){
            if(p==pmax)
                break;
            ++p;
            ind=first_generator(j0_, p, k+1);
        }
            else
                ++ind;
    }
    
    full_collection.resize(Lambda.size());
    int i = 0;
    for (typename set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it, i++){
        full_collection[i] = *it;
        cout << *it << ", " << (*it).number() << endl;
    }
    
    cout << "done setting up collection of quarklet indices..." << endl;
    cout << "total degrees of freedom between j0_ = " << j0_ << " and (jmax_= " << jmax_ << ", pmax_= " << pmax_ << ") is " << full_collection.size() << endl;
    

  }



}
