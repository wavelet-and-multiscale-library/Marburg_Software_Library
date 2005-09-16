// implementation for atlas.h



namespace FrameTL
{

//################### Atlas ###################
  
  template <unsigned int DIM_m, unsigned int DIM_d>
  inline
  Atlas<DIM_m,DIM_d>::Atlas ()
  {
  }

  template <unsigned int DIM_m, unsigned int DIM_d>
  inline
  Atlas<DIM_m,DIM_d>::Atlas (const int npatches)
    : num_patches(nPatches)
  {
  }

  template <unsigned int DIM_m, unsigned int DIM_d>
  inline
  Atlas<DIM_m,DIM_d>::Atlas (const Atlas& A)
    : charts(A.get_charts()), adjacency_matrix(A.get_adjacency_matrix()),
      num_patches(A.get_num_patches())
  {
  }

  template <unsigned int DIM_m, unsigned int DIM_d>
  inline
  Atlas<DIM_m,DIM_d>::Atlas (const Array1D< Parametrization<DIM_m,DIM_d>* >& ch)
    : charts(ch)
  {
  }

  template <unsigned int DIM_m, unsigned int DIM_d>
  inline
  const Array1D<Parametrization<DIM_m, DIM_d> *>&
  Atlas<DIM_m,DIM_d>::get_charts() const
  {
    return charts;
  }

  template <unsigned int DIM_m, unsigned int DIM_d>
  inline
  const Parametrization<DIM_m, DIM_d>&
  Atlas<DIM_m,DIM_d>::get_chart(const unsigned int p) const
  {
    assert(0 <= p && p <= num_patches);
    return *(charts[p]);
  }

  template <unsigned int DIM_m, unsigned int DIM_d>
  inline
  void Atlas<DIM_m,DIM_d>::add_chart(Parametrization<DIM_m, DIM_d>& kappa)
  {
    ++num_patches;
    Array1D< Parametrization<DIM_m,DIM_d>* > temporary_param_list(num_patches);
    for (unsigned int i = 0; i < num_patches; i++)
      {
	temporary_param_list[i] = charts[i];
      }
    temporary_param_list[num_patches-1] = &kappa;
    charts = temporary_param_list;
  }

  template <unsigned int DIM_m, unsigned int DIM_d>
  inline
  void Atlas<DIM_m,DIM_d>::set_chart(const int p,
				     const Parametrization<DIM_m, DIM_d>& kappa)
  {
    if ( p > num_patches)
      add_chart(kappa);
    else
      charts[p-1] = kappa;
  }

  template <unsigned int DIM_m, unsigned int DIM_d>
  inline
  const unsigned int
  Atlas<DIM_m,DIM_d>::get_num_patches() const
  {
    return num_patches;
  }
  
  template <unsigned int DIM_m, unsigned int DIM_d>
  inline
  void Atlas<DIM_m,DIM_d>::set_adjacent(const int p1, const int p2)
  {
    adjacency_matrix(p1-1, p2-1) = true;
  }

  template <unsigned int DIM_m, unsigned int DIM_d>
  inline
  const bool
  Atlas<DIM_m,DIM_d>::adjacent(const int p1, const int p2)
  {
    return adjacency_matrix(p1-1,p2-1);
  }

  template <unsigned int DIM_m, unsigned int DIM_d>
  inline
  const Array1D<unsigned int>&
  Atlas<DIM_m,DIM_d>::get_circumjacent (const Point<DIM_d>& p)
  {
    list<unsigned short int> circum_patches;
    //loop over all patches
    for (int i = 0; i < num_patches; i++)
      {
	if( *(charts[i]).point_in_patch(p) )
	  circum_patches.push_back(i+1);
      }
    Array1D<unsigned int> res(circum_patches.size());
    unsigned int k = 0;
    for (list<unsigned int>::const_iterator it = circum_patches.begin();
	 it != circum_patches.end(); ++it)
      {
	res[k] = circum_patches.pop_front();
	k++;
      }
    return res;
  }

  template <unsigned int DIM_m, unsigned int DIM_d>
  inline
  std::ostream&
  operator << (std::ostream& os, const Atlas<DIM_m, DIM_d>& A)
  {
    os << "Atlas consists of " 
       << A.get_num_patches() 
       << " patch(es):" << endl;
    for (unsigned int i = 0; i < A.get_num_patches(); i++)
      {
	os << "patch " << i << ":" << endl
	   << (A.get_chart(i)) << endl;  
      }
    return os;
  }

}

