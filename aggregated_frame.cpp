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
    
    j0_ = lifted_bases[0]->j0();
    cout << "minimal level = " << j0_ << endl;

//     int degees_of_freedom = 0;
//     for (unsigned int p = 0; p < n_p(); p++) {
//       unsigned int tmp = 1;
//       for (unsigned int i = 0; i < DIM_d; i++) {
// 	tmp *= (bases()[p])->bases()[i]->Deltasize(jmax + 1);
//       }
//       degees_of_freedom += tmp;
//     }

//     cout << "setting up collection of wavelet indices..." << flush;
//     indices_.resize(degees_of_freedom);
    
//     for (int i = 0; i < degees_of_freedom; i++) {
//       Index g(i,this);
//       indices_[i] = g;
//     }
//     cout << " done." << endl;


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

//     int degees_of_freedom = 0;
//     for (unsigned int p = 0; p < n_p(); p++) {
//       unsigned int tmp = 1;
//       for (unsigned int i = 0; i < DIM_d; i++) {
// 	tmp *= (bases()[p])->bases()[i]->Deltasize(jmax + 1);
//       }
//       degees_of_freedom += tmp;
//     }

//     cout << "setting up collection of wavelet indices..." << flush;
//     cout << degees_of_freedom << endl;
//     indices_.resize(degees_of_freedom);
//     for (int i = 0; i < degees_of_freedom; i++) {
//       Index g(i,this);
//       indices_[i] = g;
//     }
//     cout << " done." << endl;
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
      k[i] = WaveletTL::last_generator<IBASIS>(bases()[0]->bases()[i], j).k();

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
