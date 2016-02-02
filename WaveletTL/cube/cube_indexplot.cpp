// implementation for cube_indexplot.h

namespace WaveletTL
{
  template <class CUBEBASIS>
  void plot_indices(const CUBEBASIS* basis,
		    const InfiniteVector<double, typename CUBEBASIS::Index>& coeffs,
		    const int jmax,
		    std::ostream& os,
		    const char* colormap = "cool",
		    bool boxed = false,
		    bool colorbar = true,
		    const double aa = -6)
  {
    typedef typename CUBEBASIS::Index Index;
    typedef typename CUBEBASIS::IntervalBasis Basis1D;

    const Basis1D* xbasis = basis->bases()[0];
    const Basis1D* ybasis = basis->bases()[1];
    
    const int j0 = basis->j0();
    const double maxnorm = linfty_norm(coeffs);

    // determine smallest coefficient
    typename InfiniteVector<double, typename CUBEBASIS::Index>::const_iterator it = coeffs.begin();

    double a = 0.;

#if 0
    a = log10(fabs(*it));
    for (; it != coeffs.end();++it) {
      if (*it < a )
	a = log10(fabs(*it));
    }
#endif

    a = aa;

    //const double threshold = pow(10,aa);
    const double threshold = 1e-15;

    cout << "linfinity norm of coefficients = " << maxnorm << endl;
    cout << "number of coefficients" << coeffs.size() << endl;
    int count = 0;

    const int num_plots_col = colorbar ? jmax-j0+1 : jmax-j0+2;

    // first plot all generator coefficients on the coarsest level

    const int n_generators_x = xbasis->Deltasize(j0);
    const int n_generators_y = ybasis->Deltasize(j0);

    const double h0_x = 1./n_generators_x;
    const double h0_y = 1./n_generators_y;
    
    os << "colormap(" << colormap << ");" << endl;
    os << "X=" << colormap << ";" << endl;
    os << "subplot(4," << jmax-j0+1 << ",1);" << endl;
    for (int j = 0; j < n_generators_y; j++) {
      double y = j * h0_y;
      for (int i = 0; i < n_generators_x; i++) {
 	double x = i * h0_x;
	
	const double c = coeffs.get_coefficient(Index(j0, MultiIndex<int,2>(0,0),
						      MultiIndex<int,2>(xbasis->DeltaLmin()+i,ybasis->DeltaLmin()+j),
						      basis));

	if (fabs(c) < threshold) {
// 	  // draw an empty rectangle
// 	  os << "rectangle('position',["
// 	     << x << "," << y << "," << h0_x << "," << h0_y << "],"
// 	     << "'LineWidth',0.125,"
// 	     << "'FaceColor',[1.0,1.0,1.0]";
// 	  if (!boxed)
// 	    os << ",'EdgeColor','none'";
// 	  os << ");" << endl;
	} 
	else {
	  count++;	  
	  double col = (std::max(log10(fabs(c)/maxnorm),a)+(-a))/(-a);
	  os << "Y=ceil(" << col << " * length(colormap));" << endl;
	  os << "if Y==0" << endl << "Y=Y+1;" << endl << "end;" << endl;
	  os << "patch("
	     << "[" << x << ";" << x+h0_x << ";" << x+h0_x << ";" << x << "],"
	     << "[" << y << ";" << y << ";" << y+h0_y << ";" << y+h0_y << "],"
	     << "X(Y,:)";

	  if (!boxed)
	    os << ",'EdgeColor','none'";
	  os 
 	    // 	   << ",'LineWidth',0.125"
 	    // 	   << ",'LineStyle','none'"
 	    // 	   << ",'EraseMode','background'"
 	    << ")" << endl;
 	}
      }
      // draw boxes and axes on top of figure
      os << "box on" << endl;
      os << "set(gca,'Layer','top')" << endl;

      // turn off x ticks
      os << "set(gca,'XTick',[])" << endl;
      // turn off y ticks
      os << "set(gca,'YTick',[])" << endl;

      os << "ylabel('gen-gen:               ','rotation',0)" << endl;
      os << "title 'level " << j0 << "'" << endl;
      os << "axis([0,1,0,1]);" << endl;

    }
//     // plot a colorbar (or not)
//     if (colorbar) {
//       //os << "subplot(4," << jmax-j0+1 << "," << "[1 4]);" << endl;
//       //os << "colorbar('location','NorthOutside')" << endl;
//       os << "colorbar" << endl;
//     }


    // plot the generator-wavelet coefficients
    for (int j = j0; j <= jmax; j++) {
      os << "subplot(4," << jmax-j0+1 << "," << num_plots_col+(j-j0)+1  << ");" << endl;
      const int n_generators_x = xbasis->Deltasize(j);
      const int n_wavelets_y = ybasis->Nablasize(j);

      const double hj_x = 1./n_generators_x;
      const double hj_y = 1./n_wavelets_y;


      for (int k = 0; k < n_wavelets_y; k++) {
	double y = k * hj_y;
	for (int i = 0; i < n_generators_x; i++) {
 	double x = i * hj_x;

	const double c = coeffs.get_coefficient(Index(j, MultiIndex<int,2>(0,1),
						      MultiIndex<int,2>(xbasis->DeltaLmin()+i,ybasis->Nablamin()+k),
						      basis));
	if (fabs(c) < threshold) {
// 	  if (j <= 6) {
// 	    // draw an empty rectangle
// 	    os << "rectangle('position',["
// 	       << x << "," << y << "," << hj_x << "," << hj_y << "],"
// 	       << "'LineWidth',0.125,"
// 	       << "'FaceColor',[1.0,1.0,1.0]";
// 	    if (!boxed)
// 	      os << ",'EdgeColor','none'";
// 	    os << ")" << endl;
// 	  }
	} else {
	  count++;
	  double col = (std::max(log10(fabs(c)/maxnorm),a)+(-a))/(-a);
	  os << "Y=ceil(" << col << " * length(colormap));" << endl;
	  os << "if Y==0" << endl << "Y=Y+1;" << endl << "end;" << endl;
	  os << "patch("
	     << "[" << x << ";" << x+hj_x << ";" << x+hj_x << ";" << x << "],"
	     << "[" << y << ";" << y << ";" << y+hj_y << ";" << y+hj_y << "],"
	     << "X(Y,:)";
	  if (!boxed)
	    os << ",'EdgeColor','none'";
	  os
// 	    << ",'LineStyle','none'"
// 	    << ",'LineWidth',0.125"
	    << ")" << endl;
	}
	}
      }
      // draw boxes and axes on top of figure
      os << "box on" << endl;
      os << "set(gca,'Layer','top')" << endl;

      // turn off x ticks
      os << "set(gca,'XTick',[])" << endl;
      // turn off y ticks
      os << "set(gca,'YTick',[])" << endl;

      if (j == j0)
	os << "ylabel('gen-wav:               ','rotation',0)" << endl;

      os << "title 'level " << j << "'" << endl;
      os << "axis([0,1,0,1]);" << endl;

    }  

    // plot the wavelet-generator coefficients
    for (int j = j0; j <= jmax; j++) {
      os << "subplot(4," << jmax-j0+1 << "," << 2*num_plots_col+(j-j0)+1  << ");" << endl;
      const int n_wavelets_x = xbasis->Nablasize(j);
      const int n_generators_y = ybasis->Deltasize(j);

      const double hj_x = 1./n_wavelets_x;
      const double hj_y = 1./n_generators_y;


      for (int k = 0; k < n_generators_y; k++) {
	double y = k * hj_y;
	for (int i = 0; i < n_wavelets_x; i++) {
 	double x = i * hj_x;

	const double c = coeffs.get_coefficient(Index(j, MultiIndex<int,2>(1,0),
						      MultiIndex<int,2>(xbasis->Nablamin()+i,ybasis->DeltaLmin()+k),
						      basis));

	if (fabs(c) < threshold) {
// 	  if (j <= 6) {
// 	    // draw an empty rectangle
// 	    os << "rectangle('position',["
// 	       << x << "," << y << "," << hj_x << "," << hj_y << "],"
// 	       << "'LineWidth',0.125,"
// 	       << "'FaceColor',[1.0,1.0,1.0]";
// 	    if (!boxed)
// 	      os << ",'EdgeColor','none'";
// 	    os << ")" << endl;
// 	  }
	} else {
	  count++;
	  double col = (std::max(log10(fabs(c)/maxnorm),a)+(-a))/(-a);
	  os << "Y=ceil(" << col << " * length(colormap));" << endl;
	  os << "if Y==0" << endl << "Y=Y+1;" << endl << "end;" << endl; 
	  os << "patch("
	     << "[" << x << ";" << x+hj_x << ";" << x+hj_x << ";" << x << "],"
	     << "[" << y << ";" << y << ";" << y+hj_y << ";" << y+hj_y << "],"
	     << "X(Y,:)";
	  if (!boxed)
	    os << ",'EdgeColor','none'";
	  os
// 	    << ",'LineStyle','none'"
// 	    << ",'LineWidth',0.125"
	    << ")" << endl;
	}
	}
      }
      // draw boxes and axes on top of figure
      os << "box on" << endl;
      os << "set(gca,'Layer','top')" << endl;
      
      // turn off x ticks
      os << "set(gca,'XTick',[])" << endl;
      // turn off y ticks
      os << "set(gca,'YTick',[])" << endl;
      
      if (j == j0)
	os << "ylabel('wav-gen:               ','rotation',0)" << endl;

      os << "title 'level " << j << "'" << endl;
      os << "axis([0,1,0,1]);" << endl;
      
    }

    // plot the wavelet-wavelet coefficients
    for (int j = j0; j <= jmax; j++) {
      os << "subplot(4," << jmax-j0+1 << "," << 3*num_plots_col+(j-j0)+1  << ");" << endl;
      const int n_wavelets_x = xbasis->Nablasize(j);
      const int n_wavelets_y = ybasis->Nablasize(j);
      const double hj_y = 1./n_wavelets_y;
      const double hj_x = 1./n_wavelets_x;

      for (int k = 0; k < n_wavelets_y; k++) {
	double y = k * hj_y;
	for (int i = 0; i < n_wavelets_x; i++) {
 	double x = i * hj_x;

	const double c = coeffs.get_coefficient(Index(j, MultiIndex<int,2>(1,1),
						      MultiIndex<int,2>(xbasis->Nablamin()+i,ybasis->Nablamin()+k),
						      basis));
	if (fabs(c) < threshold) {
// 	  if (j <= 6) {
// 	    // draw an empty rectangle
// 	    os << "rectangle('position',["
// 	       << x << "," << y << "," << hj_x << "," << hj_y << "],"
// 	       << "'LineWidth',0.125,"
// 	       << "'FaceColor',[1.0,1.0,1.0]";
// 	    if (!boxed)
// 	      os << ",'EdgeColor','none'";
// 	    os << ")" << endl;
// 	  }
	} else {
	  count++;
	  double col = (std::max(log10(fabs(c)/maxnorm),a)+(-a))/(-a);
	  os << "Y=ceil(" << col << " * length(colormap));" << endl;
	  os << "if Y==0" << endl << "Y=Y+1;" << endl << "end;" << endl;
	  os << "patch("
	     << "[" << x << ";" << x+hj_x << ";" << x+hj_x << ";" << x << "],"
	     << "[" << y << ";" << y << ";" << y+hj_y << ";" << y+hj_y << "],"
	     << "X(Y,:)";
	  if (!boxed)
	    os << ",'EdgeColor','none'";
	  os
// 	    << ",'LineStyle','none'"
// 	    << ",'LineWidth',0.125"
	    << ")" << endl;
	}
	}
      }
      // draw boxes and axes on top of figure
      os << "box on" << endl;
      os << "set(gca,'Layer','top')" << endl;
      
      // turn off x ticks
      os << "set(gca,'XTick',[])" << endl;
      // turn off y ticks
      os << "set(gca,'YTick',[])" << endl;

      if (j == j0)
	os << "ylabel('wav-wav:               ','rotation',0)" << endl;
      os << "title 'level " << j << "'" << endl;    
      os << "axis([0,1,0,1]);" << endl;
    }
    cout << "considered coefficients = " << count << endl;
  }
}
