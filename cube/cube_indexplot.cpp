// implementation for cube_indexplot.h

namespace WaveletTL
{
  template <class IBASIS>
  void plot_indices(const CubeBasis<IBASIS,2>* basis,
		    const InfiniteVector<double, typename CubeBasis<IBASIS,2>::Index>& coeffs,
		    const int jmax,
		    std::ostream& os,
		    const char* colormap = "cool",
		    bool boxed = false,
		    bool colorbar = true,
		    const double a = -6)
  {
  }

}
