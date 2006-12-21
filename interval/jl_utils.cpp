// implementation for jl_tools.h

#include <numerics/bezier.h>

namespace WaveletTL
{
  double evaluate(const unsigned int derivative,
		  const int j, const int e, const int c, const int k,
 		  const double x)
  {
    assert(derivative <= 1); // currently we only support derivatives up to the first order

    double r = 0;

    if (e == 0) {
      // generator
      const double factor = (c==0 ? sqrt(35./26.) : sqrt(105./2.)); // 1/||phi_c||_2
      r = factor * (derivative == 0
		    ? MathTL::EvaluateHermiteSpline_td  (c, j, k, x)
		    : MathTL::EvaluateHermiteSpline_td_x(c, j, k, x));
    } else {
      // wavelet
      const double phi0factor = sqrt(35./26.);
      const double phi1factor = sqrt(105./2.);
      if (c == 0) {
 	// type psi_0
 	// psi_0(x) = -2*phi_0(2*x+1)+4*phi_0(2*x)-2*phi_0(2*x-1)-21*phi_1(2*x+1)+21*phi_1(2*x-1)
	const double factor = sqrt(35./352.);
	
 	int m = 2*k-1; // m-2k=-1
 	// phi_0(2x+1)
 	r += (-2.0)*M_SQRT1_2 * factor/phi0factor
 	  * evaluate(derivative, j+1, 0, 0, m, x);
	// phi_1(2x+1)
	r += (-21.0)*M_SQRT1_2 * factor/phi1factor
	  * evaluate(derivative, j+1, 0, 1, m, x);
	
 	// m=2k <-> m-2k=0
 	m++;
 	// phi_0(2x)
	r += 4.0*M_SQRT1_2 * factor/phi0factor
	  * evaluate(derivative, j+1, 0, 0, m, x);
	
 	// m=2k+1 <-> m-2k=1
 	m++;
 	// phi_0(2x-1)
 	r += (-2.0)*M_SQRT1_2 * factor/phi0factor
 	  * evaluate(derivative, j+1, 0, 0, m, x);
	// phi_1(2x-1)
	r += 21.0*M_SQRT1_2 * factor/phi1factor
	  * evaluate(derivative, j+1, 0, 1, m, x);
      } else { // c == 1
 	// type psi_1
 	// psi_1(x) = phi_0(2*x+1)-phi_0(2*x-1)+ 9*phi_1(2*x+1)+12*phi_1(2*x)+ 9*phi_1(2*x-1)
	const double factor = sqrt(35./48.);
	
 	int m = 2*k-1; // m-2k=-1
 	// phi_0(2x+1)
 	r += M_SQRT1_2 * factor/phi0factor
 	  * evaluate(derivative, j+1, 0, 0, m, x);
	// phi_1(2x+1)
	r += 9.0*M_SQRT1_2 * factor/phi1factor
	  * evaluate(derivative, j+1, 0, 1, m, x);
	
 	// m=2k <-> m-2k=0
 	m++;
 	// phi_1(2x)
 	r += 12.0*M_SQRT1_2 * factor/phi1factor
 	  * evaluate(derivative, j+1, 0, 1, m, x);
	
 	// m=2k+1 <-> m-2k=1
 	m++;
 	// phi_0(2x-1)
 	r += (-1.0)*M_SQRT1_2 * factor/phi0factor
 	  * evaluate(derivative, j+1, 0, 0, m, x);
	// phi_1(2x-1)
	r += 9.0*M_SQRT1_2 * factor/phi1factor
	  * evaluate(derivative, j+1, 0, 1, m, x);
      }
    }

    return r;
  }
  

  void evaluate(const int j, const int e, const int c, const int k,
		const Array1D<double>& points, Array1D<double>& funcvalues, Array1D<double>& dervalues)
  {
    const unsigned int npoints(points.size());
    funcvalues.resize(npoints);
    dervalues.resize(npoints);
    if (e == 0) {
      // generator
      const double factor = (c==0 ? sqrt(35./26.) : sqrt(105./2.)); // 1/||phi_c||_2
      for (unsigned int m(0); m < npoints; m++) {
	funcvalues[m] = factor * MathTL::EvaluateHermiteSpline_td  (c, j, k, points[m]);
	dervalues[m]  = factor * MathTL::EvaluateHermiteSpline_td_x(c, j, k, points[m]);
      }
    } else {
      // wavelet
      for (unsigned int i(0); i < npoints; i++) {
	funcvalues[i] = 0;
	dervalues[i] = 0;
      }
      
      const double phi0factor = sqrt(35./26.);
      const double phi1factor = sqrt(105./2.);
      
      Array1D<double> help1, help2;
      if (c == 0) {
	// type psi_0
	// psi_0(x) = -2*phi_0(2*x+1)+4*phi_0(2*x)-2*phi_0(2*x-1)-21*phi_1(2*x+1)+21*phi_1(2*x-1)
	const double factor = sqrt(35./352.);
	
	int m = 2*k-1; // m-2k=-1
	// phi_0(2x+1)
	evaluate(j+1, 0, 0, m, points, help1, help2);
	for (unsigned int i = 0; i < npoints; i++) {
	  funcvalues[i] += (-2.0)*M_SQRT1_2 * factor/phi0factor * help1[i];
	  dervalues[i]  += (-2.0)*M_SQRT1_2 * factor/phi0factor * help2[i];
	}
	// phi_1(2x+1)
	evaluate(j+1, 0, 1, m, points, help1, help2);
	for (unsigned int i = 0; i < npoints; i++) {
	  funcvalues[i] += (-21.0)*M_SQRT1_2 * factor/phi1factor * help1[i];
	  dervalues[i]  += (-21.0)*M_SQRT1_2 * factor/phi1factor * help2[i];
	}
	
	// m=2k <-> m-2k=0
	m++;
	// phi_0(2x)
	evaluate(j+1, 0, 0, m, points, help1, help2);
	for (unsigned int i = 0; i < npoints; i++) {
	  funcvalues[i] += 4.0*M_SQRT1_2 * factor/phi0factor * help1[i];
	  dervalues[i]  += 4.0*M_SQRT1_2 * factor/phi0factor * help2[i];
	}
	
	// m=2k+1 <-> m-2k=1
	m++;
	// phi_0(2x-1)
	evaluate(j+1, 0, 0, m, points, help1, help2);
	for (unsigned int i = 0; i < npoints; i++) {
	  funcvalues[i] += (-2.0)*M_SQRT1_2 * factor/phi0factor * help1[i];
	  dervalues[i]  += (-2.0)*M_SQRT1_2 * factor/phi0factor * help2[i];
	}
	// phi_1(2x-1)
	evaluate(j+1, 0, 1, m, points, help1, help2);
	for (unsigned int i = 0; i < npoints; i++) {
	  funcvalues[i] += 21.0*M_SQRT1_2 * factor/phi1factor * help1[i];
	  dervalues[i]  += 21.0*M_SQRT1_2 * factor/phi1factor * help2[i];
	}
      } else { // c == 1
	// type psi_1
	// psi_1(x) = phi_0(2*x+1)-phi_0(2*x-1)+ 9*phi_1(2*x+1)+12*phi_1(2*x)+ 9*phi_1(2*x-1)
	const double factor = sqrt(35./48.);
	
	int m = 2*k-1; // m-2k=-1
	// phi_0(2x+1)
	evaluate(j+1, 0, 0, m, points, help1, help2);
	for (unsigned int i = 0; i < npoints; i++) {
	  funcvalues[i] += M_SQRT1_2 * factor/phi0factor * help1[i];
	  dervalues[i]  += M_SQRT1_2 * factor/phi0factor * help2[i];
	}
	// phi_1(2x+1)
	evaluate(j+1, 0, 1, m, points, help1, help2);
	for (unsigned int i = 0; i < npoints; i++) {
	  funcvalues[i] += 9.0*M_SQRT1_2 * factor/phi1factor * help1[i];
	  dervalues[i]  += 9.0*M_SQRT1_2 * factor/phi1factor * help2[i];
	}
	
	// m=2k <-> m-2k=0
	m++;
	// phi_1(2x)
	evaluate(j+1, 0, 1, m, points, help1, help2);
	for (unsigned int i = 0; i < npoints; i++) {
	  funcvalues[i] += 12.0*M_SQRT1_2 * factor/phi1factor * help1[i];
	  dervalues[i]  += 12.0*M_SQRT1_2 * factor/phi1factor * help2[i];
	}
	
	// m=2k+1 <-> m-2k=1
	m++;
	// phi_0(2x-1)
	evaluate(j+1, 0, 0, m, points, help1, help2);
	for (unsigned int i = 0; i < npoints; i++) {
	  funcvalues[i] += (-1.0)*M_SQRT1_2 * factor/phi0factor * help1[i];
	  dervalues[i]  += (-1.0)*M_SQRT1_2 * factor/phi0factor * help2[i];
	}
	// phi_1(2x-1)
	evaluate(j+1, 0, 1, m, points, help1, help2);
	for (unsigned int i = 0; i < npoints; i++) {
	  funcvalues[i] += 9.0*M_SQRT1_2 * factor/phi1factor * help1[i];
	  dervalues[i]  += 9.0*M_SQRT1_2 * factor/phi1factor * help2[i];
	}
      }
    }
  }  

}
