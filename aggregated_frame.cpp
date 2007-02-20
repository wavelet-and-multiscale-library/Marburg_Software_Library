// implementation for aggregated_frame.h

namespace FrameTL
{
  template<class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  AggregatedFrame<IBASIS,DIM_d,DIM_m>::AggregatedFrame(const Atlas<DIM_d, DIM_m>* atlas,
						       const Array1D<FixedArray1D<int,2*DIM_d> >& bc,
						       const Array1D<FixedArray1D<int,2*DIM_d> >& bcT,
						       const int jmax)
    : atlas_(atlas), bc_(bc), bcT_(bcT), jmax_(jmax)
  {
    lifted_bases.resize((atlas_->charts()).size());

    for (unsigned int  i = 0; i < (atlas_->charts()).size(); ++i)
      lifted_bases[i] = new MappedCubeBasis<IBASIS,DIM_d,DIM_m>((atlas_->charts())[i],bc[i],bcT[i]);
    
    j0_ = lifted_bases[0]->j0();
    cout << "minimal level = " << j0_ << endl;

    cout << "setting up full collection of wavelet indices level by level..." << endl;
    full_collection_levelwise.resize(jmax - j0_ + 2);
    

    typedef AggregatedFrame<IBASIS,DIM_d,DIM_m> FRAME;

    int count = 0;
    int degrees_of_freedom = 0;

    // allocate memory for all indices on level
    for (int j = j0_-1; j <= jmax; j++) {
      int degees_of_freedom_on_lev = 0;
      Index first;
      Index last;
      if (j == j0_-1) {
	// determine how many functions there are
	for (unsigned int p = 0; p < n_p(); p++) {
	  int tmp = 1;
	  for (unsigned int i = 0; i < DIM_d; i++) {
	    tmp *= bases()[p]->bases()[i]->Deltasize(j0_);
	  }
	  degees_of_freedom_on_lev += tmp;
	}
	
 	first =  FrameTL::first_generator<IBASIS, DIM_d, DIM_m, FRAME>(this, j0_);
 	last = FrameTL::last_generator<IBASIS, DIM_d, DIM_m, FRAME>(this, j0_);
	cout << "degrees of freedom on level = " << degees_of_freedom_on_lev << endl;
      }
      else {
	int degees_of_freedom_on_lev_min_1 = 0;
	// determine how many functions there are
	for (unsigned int p = 0; p < n_p(); p++) {
	  int tmp = 1;
	  int tmp2 = 1;
	  for (unsigned int i = 0; i < DIM_d; i++) {
	    tmp *= bases()[p]->bases()[i]->Deltasize(j+1);
	    tmp2 *= bases()[p]->bases()[i]->Deltasize(j);
	  }
	  degees_of_freedom_on_lev += tmp;
	  degees_of_freedom_on_lev_min_1 += tmp2;
	}
	degees_of_freedom_on_lev -= degees_of_freedom_on_lev_min_1;
 	first = FrameTL::first_wavelet<IBASIS, DIM_d, DIM_m, FRAME>(this, j);
 	last = FrameTL::last_wavelet<IBASIS, DIM_d, DIM_m, FRAME>(this, j);
	cout << "degrees of freedom on level = " << degees_of_freedom_on_lev << endl;
      }
      int k = 0;
      
      full_collection_levelwise[count].resize(degees_of_freedom_on_lev);
      for (Index ind = first; ind <= last; ++ind) {
	full_collection_levelwise[count][k] = ind;
	k++;
      }
      count++;
      degrees_of_freedom += degees_of_freedom_on_lev;
    }
    cout << "done setting up full collection of wavelet indices level by level..." << endl;

    cout << "preprocessing all supports on cubes..." << endl;
    all_supports.resize(degrees_of_freedom);

    cout << "degrees of freedom = " << degrees_of_freedom << endl;

    count = 0;
    for (int j = j0_-1; j <= jmax; j++) {
      for (unsigned int i = 0; i < full_collection_levelwise[j-j0()+1].size(); i++) {
	Index* ind = &full_collection_levelwise[j-j0()+1][i];

	typename WaveletTL::CubeBasis<IBASIS,DIM_d>::Support supp_of_ind;

	typedef typename WaveletTL::CubeBasis<IBASIS,DIM_d>::Index CubeIndex;

	WaveletTL::support<IBASIS,DIM_d>(*bases()[ind->p()], 
					 CubeIndex(ind->j(),
						   ind->e(),
						   ind->k(),
						   bases()[ind->p()]),
					 supp_of_ind);

	all_supports[count].j = supp_of_ind.j;
	for (unsigned int i = 0; i < DIM_d; i++) {
	  all_supports[count].a[i] = supp_of_ind.a[i];
	  all_supports[count].b[i] = supp_of_ind.b[i];
	}
		
	count++;
      }
    }
    cout << "done preprocessing all supports on cubes" << endl;


    cout << "setting up full collection of wavelet indices..." << endl;
    full_collection.resize(degrees_of_freedom);
    int k = 0;
    for (Index ind = FrameTL::first_generator<IBASIS, DIM_d, DIM_m, FRAME>(this, j0_);
	 ind <= FrameTL::last_wavelet<IBASIS, DIM_d, DIM_m, FRAME>(this, jmax); ++ind) {
      full_collection[k] = ind;
      k++;
    }
    cout << "done setting up full collection of wavelet indices..." << endl;
  }

  template<class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  AggregatedFrame<IBASIS,DIM_d,DIM_m>::AggregatedFrame(const Atlas<DIM_d, DIM_m>* atlas,
						       const Array1D<FixedArray1D<int,2*DIM_d> >& bc,
						       const int jmax)
    : atlas_(atlas), bc_(bc), jmax_(jmax)
  {
    
    lifted_bases.resize((atlas_->charts()).size());
    cout << "boundary conditions = " << bc << endl;

    for (unsigned int  i = 0; i < (atlas_->charts()).size(); ++i)
      lifted_bases[i] = new MappedCubeBasis<IBASIS,DIM_d,DIM_m>((atlas_->charts())[i],bc[i]);
    
    j0_ = lifted_bases[0]->j0();
    cout << "minimal level = " << j0_ << endl;

    cout << "setting up collection of wavelet indices level by level..." << endl;
    full_collection_levelwise.resize(jmax - j0_ + 2);
    

    typedef AggregatedFrame<IBASIS,DIM_d,DIM_m> FRAME;

    int count = 0;
    int degrees_of_freedom = 0;

    // allocate memory for all indices on level
    for (int j = j0_-1; j <= jmax; j++) {
      int degrees_of_freedom_on_lev = 0;
      Index first;
      Index last;
      if (j == j0_-1) {
	// determine how many functions there are
	for (unsigned int p = 0; p < n_p(); p++) {
	  int tmp = 1;
	  for (unsigned int i = 0; i < DIM_d; i++) {
	    tmp *= bases()[p]->bases()[i]->Deltasize(j0_);
	  }
	  degrees_of_freedom_on_lev += tmp;
	}
	
 	first =  FrameTL::first_generator<IBASIS, DIM_d, DIM_m, FRAME>(this, j0_);
 	last = FrameTL::last_generator<IBASIS, DIM_d, DIM_m, FRAME>(this, j0_);
	cout << "degrees of freedom on level = " << degrees_of_freedom_on_lev << endl;
      }
      else {
	int degrees_of_freedom_on_lev_min_1 = 0;
	// determine how many functions there are
	for (unsigned int p = 0; p < n_p(); p++) {
	  int tmp = 1;
	  int tmp2 = 1;
	  for (unsigned int i = 0; i < DIM_d; i++) {
	    tmp *= bases()[p]->bases()[i]->Deltasize(j+1);
	    tmp2 *= bases()[p]->bases()[i]->Deltasize(j);
	  }
	  degrees_of_freedom_on_lev += tmp;
	  degrees_of_freedom_on_lev_min_1 += tmp2;
	}
	degrees_of_freedom_on_lev -= degrees_of_freedom_on_lev_min_1;
 	first = FrameTL::first_wavelet<IBASIS, DIM_d, DIM_m, FRAME>(this, j);
 	last = FrameTL::last_wavelet<IBASIS, DIM_d, DIM_m, FRAME>(this, j);
	cout << "degrees of freedom on level = " << degrees_of_freedom_on_lev << endl;
      }
      int k = 0;
      
      full_collection_levelwise[count].resize(degrees_of_freedom_on_lev);
      for (Index ind = first; ind <= last; ++ind) {
	full_collection_levelwise[count][k] = ind;
	//cout << full_collection_levelwise[count][k] << endl;
	k++;
      }
      count++;
      degrees_of_freedom += degrees_of_freedom_on_lev;
    }
    cout << "done setting up collection of wavelet indices level by level..." << endl;

    cout << "preprocessing all supports on cubes..." << endl;
    all_supports.resize(degrees_of_freedom);

    cout << "degrees of freedom = " << degrees_of_freedom << endl;

    count = 0;
    for (int j = j0_-1; j <= jmax; j++) {
      for (unsigned int i = 0; i < full_collection_levelwise[j-j0()+1].size(); i++) {
	Index* ind = &full_collection_levelwise[j-j0()+1][i];

	typename WaveletTL::CubeBasis<IBASIS,DIM_d>::Support supp_of_ind;

	typedef typename WaveletTL::CubeBasis<IBASIS,DIM_d>::Index CubeIndex;

	WaveletTL::support<IBASIS,DIM_d>(*bases()[ind->p()], 
					 CubeIndex(ind->j(),
						   ind->e(),
						   ind->k(),
						   bases()[ind->p()]),
					 supp_of_ind);

	all_supports[count].j = supp_of_ind.j;
	for (unsigned int i = 0; i < DIM_d; i++) {
	  all_supports[count].a[i] = supp_of_ind.a[i];
	  all_supports[count].b[i] = supp_of_ind.b[i];
	}
	count++;
      }
    }
    cout << "done preprocessing all supports on cubes" << endl;

    cout << "setting up collection of wavelet indices..." << endl;
    full_collection.resize(degrees_of_freedom);
    int k = 0;
    for (Index ind = FrameTL::first_generator<IBASIS, DIM_d, DIM_m, FRAME>(this, j0_);
	 ind <= FrameTL::last_wavelet<IBASIS, DIM_d, DIM_m, FRAME>(this, jmax); ++ind) {
      full_collection[k] = ind;
      k++;
    }
    cout << "done setting up collection of wavelet indices..." << endl;
  }


  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  AggregatedFrame<IBASIS,DIM_d,DIM_m>::~AggregatedFrame()
  {
    for (unsigned int  i = 0; i < (atlas_->charts()).size(); ++i)
      delete lifted_bases[i];          
  }  

  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  FrameIndex<IBASIS,DIM_d,DIM_m>
  AggregatedFrame<IBASIS,DIM_d,DIM_m>::first_generator(const int j) const
  {
    assert(j >= j0());
     
    typename FrameIndex<IBASIS,DIM_d,DIM_m>::type_type e;//== 0
    typename FrameIndex<IBASIS,DIM_d,DIM_m>::translation_type k;
    for (unsigned int i = 0; i < DIM_d; i++)
      k[i] = WaveletTL::first_generator<IBASIS>(bases()[0]->bases()[i], j).k();
     
    return FrameIndex<IBASIS,DIM_d,DIM_m>(this, j, e, 0, k);
  }
   
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  FrameIndex<IBASIS,DIM_d,DIM_m>
  AggregatedFrame<IBASIS,DIM_d,DIM_m>::last_generator(const int j) const
  {
    assert(j >= j0());

    typename FrameIndex<IBASIS,DIM_d,DIM_m>::type_type e;//== 0
    typename FrameIndex<IBASIS,DIM_d,DIM_m>::translation_type k;
    for (unsigned int i = 0; i < DIM_d; i++)
      k[i] = WaveletTL::last_generator<IBASIS>(bases()[bases().size()-1]->bases()[i], j).k();

    return FrameIndex<IBASIS,DIM_d,DIM_m>(this, j, e, bases().size()-1, k); 
  }

  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  FrameIndex<IBASIS,DIM_d,DIM_m>
  AggregatedFrame<IBASIS,DIM_d,DIM_m>::first_wavelet(const int j) const
  {
    assert(j >= j0());

    typename FrameIndex<IBASIS,DIM_d,DIM_m>::type_type e;//== 0
    typename FrameIndex<IBASIS,DIM_d,DIM_m>::translation_type k; 
    for (unsigned int i = 0; i < DIM_d-1; i++)
      k[i] = WaveletTL::first_generator<IBASIS>(bases()[0]->bases()[i], j).k();

    k[DIM_d-1] = WaveletTL::first_wavelet<IBASIS>(bases()[0]->bases()[DIM_d-1], j).k();
    e[DIM_d-1] = 1;

    return FrameIndex<IBASIS,DIM_d,DIM_m>(this, j, e, 0, k); 
  }

  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  FrameIndex<IBASIS,DIM_d,DIM_m>
  AggregatedFrame<IBASIS,DIM_d,DIM_m>::last_wavelet(const int j) const
  {
    assert(j >= j0());
     
    typename FrameIndex<IBASIS,DIM_d,DIM_m>::type_type e;//== 0
    typename FrameIndex<IBASIS,DIM_d,DIM_m>::translation_type k; 
    for (unsigned int i = 0; i < DIM_d; i++) {
      k[i] = WaveletTL::last_wavelet<IBASIS>(bases()[bases().size()-1]->bases()[i], j).k();
      e[i] = 1;
    }
          
    return FrameIndex<IBASIS,DIM_d,DIM_m>(this, j, e, bases().size()-1, k); 
  }

}
