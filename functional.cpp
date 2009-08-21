// implementation for functional.h

#include <numerics/gauss_data.h>

using namespace MathTL;


namespace FrameTL
{
  template<class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  double
  Functional<IBASIS,DIM_d,DIM_m>::evaluate(const typename AggregatedFrame<IBASIS, DIM_d, DIM_m>::Index& lambda) const
  {
    // calculate: \int_\Omega f(Kappa(x)) \psi^\Box (x) \sqrt(\sqrt(det ((D Kappa)^T(x) (D Kappa)(x)) ))
    // recall: \sqrt(\sqrt(...)) = Gram_factor

    // number of intervals in summed quadrature rule
#ifdef RHS_QUADRATURE_GRANULARITY
    const int N = RHS_QUADRATURE_GRANULARITY;
#else
    const int N = 1;
#endif
    
    
    double r = 0;
    
    const unsigned int p = lambda.p();
 
    typedef WaveletTL::CubeBasis<IBASIS,DIM_d> CUBEBASIS;
 
    // first compute supp(psi_lambda)
 
    typename CUBEBASIS::Support supp;
    
    typename CubeBasis<IBASIS,DIM_d>::Index lambda_c(lambda.j(),
						   lambda.e(),
						   lambda.k(),frame_->bases()[p]);

    
    WaveletTL::support<IBASIS,DIM_d>(*(frame_->bases()[p]), lambda_c, supp); 
    const Chart<DIM_d>* chart = frame_->atlas()->charts()[p];
    // setup Gauss points and weights for a composite quadrature formula:
    const int N_Gauss = 6;
    //const double h = ldexp(1.0, -supp.j); // granularity for the quadrature
    const double h = 1.0 / (1 << supp.j); // granularity for the quadrature
    FixedArray1D<Array1D<double>,DIM_d> gauss_points, gauss_weights, v_values;
    for (unsigned int i = 0; i < DIM_d; i++) {
      gauss_points[i].resize(N * N_Gauss * (supp.b[i]-supp.a[i]));
      gauss_weights[i].resize(N * N_Gauss * (supp.b[i]-supp.a[i]));
      for (int patch = supp.a[i]; patch < supp.b[i]; patch++)
	for(int m=0;m<N;m++)
	  for (int n = 0; n < N_Gauss; n++) {
	    gauss_points[i][N*(patch-supp.a[i])*N_Gauss+m*N_Gauss+n]
	      = h*(N*2*patch+1+GaussPoints[N_Gauss-1][n]+2.0*m)/(2.*N);
	    gauss_weights[i][N*(patch-supp.a[i])*N_Gauss+m*N_Gauss+n]
	      = h*GaussWeights[N_Gauss-1][n]/N;
	  }
    }

    // compute the point values of the integrand (where we use that it is a tensor product)
    for (unsigned int i = 0; i < DIM_d; i++) {
       
      WaveletTL::evaluate(*(frame_->bases()[p]->bases()[i]), 0,
			  typename IBASIS::Index(lambda.j(),
						 lambda.e()[i],
						 lambda.k()[i],
						 frame_->bases()[p]->bases()[i]),
			  gauss_points[i], v_values[i]);
      
    }
    // iterate over all points and sum up the integral shares
    int index[DIM_d]; // current multiindex for the point values
    for (unsigned int i = 0; i < DIM_d; i++)
      index[i] = 0;
    
    Point<DIM_d> x;
    Point<DIM_d> x_patch;
    
    bool exit = false;
    while (!exit) {
      for (unsigned int i = 0; i < DIM_d; i++)
	x[i] = gauss_points[i][index[i]];

      chart->map_point(x,x_patch);

      double share = g_->value(x_patch) * chart->Gram_factor(x);
      for (unsigned int i = 0; i < DIM_d; i++)
	share *= gauss_weights[i][index[i]] * v_values[i][index[i]];
      r += share;

      // "++index"
      for (unsigned int i = 0; i < DIM_d; i++) {
	if (index[i] == N * N_Gauss * (supp.b[i]-supp.a[i])-1) {
	  index[i] = 0;
	  exit = (i == DIM_d-1);
	} else {
	  index[i]++;
	  break;
	}
      }
    }
    return r;
  }
}
