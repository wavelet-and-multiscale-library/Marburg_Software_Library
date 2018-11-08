// implementation for i_indexplot.h

namespace WaveletTL
{
  template <class IBASIS>
  void plot_indices(const IBASIS* basis,
		    const InfiniteVector<double, typename IBASIS::Index>& coeffs,
		    const int jmax,
		    std::ostream& os,
		    const char* colormap,
		    bool boxed,
		    bool colorbar,
		    const double a)
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
      const double c = coeffs.get_coefficient(Index(j0,0,basis->DeltaLmin()+i,basis));
      if (c == 0) {
// 	// draw an empty rectangle
// 	os << "rectangle('position',["
// 	   << x << "," << y << "," << h0 << "," << 1.0 << "],"
// 	   << "'LineWidth',0.125,"
// 	   << "'FaceColor',[1.0,1.0,1.0])" << endl;
      } else {
	// draw a patch (to have Matlab manage the colormap)
	os << "patch("
	   << "[" << x << ";" << x+h0 << ";" << x+h0 << ";" << x << "],"
	   << "[" << y << ";" << y << ";" << y+1 << ";" << y+1 << "],"
	   << std::max(log10(fabs(c)/maxnorm),a);
	if (!boxed)
	  os << ",'EdgeColor','none'";
	os 
// 	   << ",'LineWidth',0.125"
// 	   << ",'LineStyle','none'"
// 	   << ",'EraseMode','background'"
	  << ")" << endl;
      }
    }

//     // plot some empty boxes above level 7
//     for (int j = 7; j <= jmax; j++) {
//       os << "rectangle('position',["
// 	 << 0.0 << "," << j+0.5 << "," << 1.0 << "," << 1.0 << "],"
// 	 << "'LineWidth',0.125,"
// 	 << "'FaceColor',[1.0,1.0,1.0])" << endl;
//     }
    
    // plot the wavelet coefficients
    for (int j = j0; j <= jmax; j++) {
      const int n_wavelets = basis->Nablasize(j);
      const double hj = 1./n_wavelets;
      for (int i = 0; i < n_wavelets; i++) {
 	double x = i * hj;
 	double y = j+0.5;
  	const double c = coeffs.get_coefficient(Index(j,1,basis->Nablamin()+i,basis));
 	if (c == 0) {
	  if (j <= 6) {
// 	    // draw an empty rectangle
// 	    os << "rectangle('position',["
// 	       << x << "," << y << "," << hj << "," << 1.0 << "],"
// 	       << "'LineWidth',0.125,"
// 	       << "'FaceColor',[1.0,1.0,1.0])" << endl;
	  }
	} else {
	  // draw a patch
	  os << "patch("
	     << "[" << x << ";" << x+hj << ";" << x+hj << ";" << x << "],"
	     << "[" << y << ";" << y << ";" << y+1 << ";" << y+1 << "],"
	     << std::max(log10(fabs(c)/maxnorm),a);
	  if (!boxed)
	    os << ",'EdgeColor','none'";
	  os
//   	     << ",'LineStyle','none'"
//   	     << ",'LineWidth',0.125"
	    << ")" << endl;
	}
      }
    }

    // draw boxes and axes on top of figure
    os << "box on" << endl;
    os << "set(gca,'Layer','top')" << endl;
    
//     // set x axis limits
//     os << "set(gca,'XLim',[" << -0.01 << " " << 1.01 << "])" << endl;

    // set y axis limits
    os << "set(gca,'YLim',[" << j0-0.5 << " " << jmax+1.5 << "])" << endl;

    // set axis labels
    os << "xlabel 'k'" << endl
       << "ylabel 'level j'" << endl;
    
    // turn off x ticks
    os << "set(gca,'XTick',[])" << endl;

    // set colormap
    os << "colormap " << colormap << endl;
    
    // set color axis limits
    os << "set(gca,'CLim',[" << a << " 0])" << endl;

    // set the y-tick marks
    os << "set(gca,'YTick',[";
    for (int j = j0; j <= jmax+1; j++)
      os << j-1 << " ";
    os << "])" << endl;

    // plot a colorbar (or not)
    if (colorbar) {
      os << "colorbar" << endl;
    }
  }

  template <class IBASIS>
  void plot_indices2(const IBASIS* basis,
		     const InfiniteVector<double, typename IBASIS::Index>& coeffs,
		     const int jmax,
		     std::ostream& os,
		     const char* colormap,
		     bool boxed,
		     bool colorbar,
		     const double a)
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
      const double c = coeffs.get_coefficient(Index(j0,0,basis->DeltaLmin()+i));
      if (c == 0) {
// 	// draw an empty rectangle
// 	os << "rectangle('position',["
// 	   << x << "," << y << "," << h0 << "," << 1.0 << "],"
// 	   << "'LineWidth',0.125,"
// 	   << "'FaceColor',[1.0,1.0,1.0])" << endl;
      } else {
	// draw a patch (to have Matlab manage the colormap)
	os << "patch("
	   << "[" << x << ";" << x+h0 << ";" << x+h0 << ";" << x << "],"
	   << "[" << y << ";" << y << ";" << y+1 << ";" << y+1 << "],"
	   << std::max(log10(fabs(c)/maxnorm),a);
	if (!boxed)
	  os << ",'EdgeColor','none'";
	os 
// 	   << ",'LineWidth',0.125"
// 	   << ",'LineStyle','none'"
// 	   << ",'EraseMode','background'"
	  << ")" << endl;
      }
    }

//     // plot some empty boxes above level 7
//     for (int j = 7; j <= jmax; j++) {
//       os << "rectangle('position',["
// 	 << 0.0 << "," << j+0.5 << "," << 1.0 << "," << 1.0 << "],"
// 	 << "'LineWidth',0.125,"
// 	 << "'FaceColor',[1.0,1.0,1.0])" << endl;
//     }
    
    // plot the wavelet coefficients
    for (int j = j0; j <= jmax; j++) {
      const int n_wavelets = basis->Nablasize(j);
      const double hj = 1./n_wavelets;
      for (int i = 0; i < n_wavelets; i++) {
 	double x = i * hj;
 	double y = j+0.5;
 	const double c = coeffs.get_coefficient(Index(j,1,basis->Nablamin()+i));
 	if (c == 0) {
	  if (j <= 6) {
// 	    // draw an empty rectangle
// 	    os << "rectangle('position',["
// 	       << x << "," << y << "," << hj << "," << 1.0 << "],"
// 	       << "'LineWidth',0.125,"
// 	       << "'FaceColor',[1.0,1.0,1.0])" << endl;
	  }
	} else {
	  // draw a patch
	  os << "patch("
	     << "[" << x << ";" << x+hj << ";" << x+hj << ";" << x << "],"
	     << "[" << y << ";" << y << ";" << y+1 << ";" << y+1 << "],"
	     << std::max(log10(fabs(c)/maxnorm),a);
	  if (!boxed)
	    os << ",'EdgeColor','none'";
	  os
//   	     << ",'LineStyle','none'"
//   	     << ",'LineWidth',0.125"
	    << ")" << endl;
	}
      }
    }

    // draw boxes and axes on top of figure
    os << "box on" << endl;
    os << "set(gca,'Layer','top')" << endl;
    
//     // set x axis limits
//     os << "set(gca,'XLim',[" << -0.01 << " " << 1.01 << "])" << endl;

    // set y axis limits
    os << "set(gca,'YLim',[" << j0-0.5 << " " << jmax+1.5 << "])" << endl;

    // set axis labels
    os << "xlabel 'k'" << endl
       << "ylabel 'level j'" << endl;
    
    // turn off x ticks
    os << "set(gca,'XTick',[])" << endl;

    // set colormap
    os << "colormap " << colormap << endl;
    
    // set color axis limits
    os << "set(gca,'CLim',[" << a << " 0])" << endl;

    // set the y-tick marks
    os << "set(gca,'YTick',[";
    for (int j = j0; j <= jmax+1; j++)
      os << j-1 << " ";
    os << "])" << endl;

    // plot a colorbar (or not)
    if (colorbar) {
      os << "colorbar" << endl;
    }
  }
}
