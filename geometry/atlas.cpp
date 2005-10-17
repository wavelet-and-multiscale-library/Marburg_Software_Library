// implementation for atlas.h

namespace MathTL
{
  template <unsigned int DIM_d, unsigned int DIM_m>
  Atlas<DIM_d,DIM_m>::Atlas(const Array1D<Chart<DIM_d,DIM_m>* >& charts,
			    const SymmetricMatrix<bool>& adjacency)
    : charts_(charts), adjacency_matrix(adjacency)
  {
  }

  template <unsigned int DIM_d, unsigned int DIM_m>
  std::ostream&
  operator << (std::ostream& s, const Atlas<DIM_d, DIM_m>& A)
  {
    if (A.charts().size() == 0)
      s << "Atlas: empty atlas!" << endl;
    else {
      s << "Atlas: consists of "
	<< A.charts().size()
	<< " patch(es):" << endl;

      for (unsigned int i = 0; i < A.charts().size(); i++) {
  	s << "patch " << i << ":" << endl
  	  << *(A.charts()[i]);
      }

      s << "adjacency relations:" << endl
 	<< A.get_adjacency_matrix();
    }

    return s;
  }

}
