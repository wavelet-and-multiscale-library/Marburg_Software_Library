// implementation for i_indexplot.h

namespace WaveletTL
{
  template <class IBASIS>
  void plot_indices(const IBASIS* basis,
		    const InfiniteVector<double, typename IBASIS::Index>& coeffs,
		    const int jmax,
		    std::ostream& os)
  {
    const int j0 = basis->j0();

    // first plot all generators on the coarsest level
    const int n_generators = basis->Deltasize(j0);
    
  }
}
