// implementation for atlas.h

namespace FrameTL
{

//################### Atlas ###################
  
  template <unsigned int DIM_m, unsigned int DIM_d>
  inline
  Atlas<DIM_m,DIM_d>::Atlas () : num_patches(0)
  {
  }

  template <unsigned int DIM_m, unsigned int DIM_d>
  inline
  Atlas<DIM_m,DIM_d>::Atlas (const int nPatches)
    : num_patches(nPatches), charts(nPatches), adjacency_matrix(nPatches)
  {
    for (unsigned int i = 0; i < adjacency_matrix.column_dimension(); i++)
      adjacency_matrix(i,i) = true;
  }

  template <unsigned int DIM_m, unsigned int DIM_d>
  inline
  Atlas<DIM_m,DIM_d>::Atlas (const Atlas<DIM_m,DIM_d>& A)
    : num_patches(A.get_num_patches()), charts(A.get_charts()),
      adjacency_matrix(A.get_adjacency_matrix())
  {
  }

  template <unsigned int DIM_m, unsigned int DIM_d>
  inline
  Atlas<DIM_m,DIM_d>::Atlas (const Array1D< Parametrization<DIM_m,DIM_d>* >& ch)
    : num_patches(ch.size()), charts(ch), adjacency_matrix(ch.size())
  {
    for (unsigned int i = 0; i < adjacency_matrix.column_dimension(); i++)
      adjacency_matrix(i,i) = true;
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
    assert(0 <= p && p < num_patches);
    return *(charts[p]);
  }

  template <unsigned int DIM_m, unsigned int DIM_d>
  inline
  void Atlas<DIM_m,DIM_d>::add_chart(Parametrization<DIM_m, DIM_d>& kappa)
  {
    if (num_patches == 0)
      {
	charts.resize(1);
	charts[0] = &kappa;
	++num_patches;
	return;
      }

    ++num_patches;

    Array1D< Parametrization<DIM_m,DIM_d>* > temporary_param_list(num_patches);
    for (unsigned int i = 0; i < num_patches-1; i++)
      {
	temporary_param_list[i] = charts[i];
      }
    temporary_param_list[num_patches-1] = &kappa;
    charts.resize(num_patches);
    charts = temporary_param_list;
  }

  template <unsigned int DIM_m, unsigned int DIM_d>
  inline
  void Atlas<DIM_m,DIM_d>::set_chart(const unsigned int p,
				     Parametrization<DIM_m, DIM_d>& kappa)
  {
    if ( p > num_patches)
      add_chart(kappa);
    else
      charts[p-1] = &kappa;
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
    assert( (1 <= p1-1 < adjacency_matrix.row_dimension()) &&
	    (1 <= p2-1 < adjacency_matrix.column_dimension()));
    adjacency_matrix(p1-1, p2-1) = true;
  }

  template <unsigned int DIM_m, unsigned int DIM_d>
  inline
  const bool
  Atlas<DIM_m,DIM_d>::adjacent(const int p1, const int p2)
  {
    assert( (1 <= p1-1 < adjacency_matrix.row_dimension()) &&
	    (1 <= p2-1 < adjacency_matrix.column_dimension()));
    return adjacency_matrix(p1-1,p2-1);
  }

  template <unsigned int DIM_m, unsigned int DIM_d>
  inline
  const SymmetricMatrix<bool>&
  Atlas<DIM_m,DIM_d>::get_adjacency_matrix() const
  {
    return adjacency_matrix;
  }

  template <unsigned int DIM_m, unsigned int DIM_d>
  inline
  void Atlas<DIM_m,DIM_d>::set_adjacency_matrix(const SymmetricMatrix<bool>& M)
  {
    adjacency_matrix = M;
  }

  template <unsigned int DIM_m, unsigned int DIM_d>
  inline
  const Array1D<unsigned int>
  Atlas<DIM_m,DIM_d>::get_circumjacent (const Point<DIM_d>& p)
  {
    list<unsigned short int> circum_patches;
    //loop over all patches
    for (unsigned int i = 0; i < num_patches; i++)
      {
	if( (*(charts[i])).point_in_patch(p) )
	  circum_patches.push_back(i+1);
      }
    Array1D<unsigned int> res(circum_patches.size());
    unsigned int k = 0;
    for (list<unsigned short int>::const_iterator it = circum_patches.begin();
	 it != circum_patches.end(); ++it)
      {
	res[k] = circum_patches.front();
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
	os << "patch " << i+1 << ":" << endl
	  << (A.get_chart(i)).toString() << endl;  
      }
    os << "adjacency relation:" << endl
       << A.get_adjacency_matrix();

    return os;
  }

}

