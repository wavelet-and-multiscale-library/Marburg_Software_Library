// implementation for i_indexplot.h

namespace WaveletTL
{
  template <class IBASIS>
  void plot_indices(const IBASIS* basis,
		    const InfiniteVector<double, typename IBASIS::Index>& coeffs,
		    const int jmax,
		    std::ostream& os,
		    MathTL::MatlabColorMap colormap)
  {
    typedef typename IBASIS::Index Index;
    
    const int j0 = basis->j0();
    const double maxnorm = linfty_norm(coeffs);

    // first plot all generator coefficients on the coarsest level
    const int n_generators = basis->Deltasize(j0);
    const double h0 = 1./n_generators;
    for (int i = 0; i < n_generators; i++) {
      double x = i * h0;
      double y = j0-0.5;
      double red, green, blue;
      MathTL::MatlabColorMap colormap = jet;
      const double c = coeffs.get_coefficient(Index(j0,0,basis->DeltaLmin()+i,basis));
      if (c == 0)
	red = green = blue = 1.0;
      else
	MathTL::get_color(c/maxnorm, colormap, red, green, blue);
      
      os << "rectangle('position',["
	 << x << "," << y << "," << h0 << "," << 1.0 << "],"
  	 << "'LineWidth',0.125,"
	 << "'FaceColor',[" << red << "," << green << "," << blue << "])" << endl;
    }

    // plot the wavelet coefficients
    for (int j = j0; j <= jmax; j++) {
      const int n_wavelets = basis->Nablasize(j);
      const double hj = 1./n_wavelets;
      for (int i = 0; i < n_wavelets; i++) {
	double x = i * hj;
	double y = j+0.5;
	double red, green, blue;
	const double c = coeffs.get_coefficient(Index(j,1,basis->Nablamin()+i,basis));
	if (c == 0)
	  red = green = blue = 1.0;
	else
	  MathTL::get_color(c/maxnorm, colormap, red, green, blue);
	
	os << "rectangle('position',["
	   << x << "," << y << "," << hj << "," << 1.0 << "],"
 	   << "'LineWidth',0.125,"
	   << "'FaceColor',[" << red << "," << green << "," << blue << "])" << endl;
      }
    }
  }
}
