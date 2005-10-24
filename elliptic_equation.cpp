// implementation for elliptic_equation.h

#include <numerics/gauss_data.h>
#include <cube/cube_basis.h>

namespace FrameTL
{

  template <class IBASIS, unsigned int DIM>
  EllipticEquation<IBASIS,DIM>::EllipticEquation(const EllipticBVP<DIM>& ell_bvp,
						 const AggregatedFrame<IBASIS,DIM>* frame)
    : ell_bvp_(ell_bvp), frame_(frame)
  {
  }

  template <class IBASIS, unsigned int DIM>
  double
  EllipticEquation<IBASIS,DIM>::a_same_patches(const typename AggregatedFrame<IBASIS,DIM>::Index& lambda,
					       const typename AggregatedFrame<IBASIS,DIM>::Index& nu,
					       const unsigned int p) const
  {
    return 0.0;
  }

  template <class IBASIS, unsigned int DIM>
  double
  EllipticEquation<IBASIS,DIM>::a_different_patches(const typename AggregatedFrame<IBASIS,DIM>::Index& lambda,
						    const typename AggregatedFrame<IBASIS,DIM>::Index& nu,
						    const unsigned int p) const
  {
    return 0.0;
  }


  template <class IBASIS, unsigned int DIM>
  double
  EllipticEquation<IBASIS,DIM>::a(const typename AggregatedFrame<IBASIS,DIM>::Index& lambda,
				  const typename AggregatedFrame<IBASIS,DIM>::Index& nu,
				  const unsigned int p) const
  {
    //the cases lambda.p() == nu.p() and lambda.p() != nu.p() have
    //to be treated seperately
    return lambda.p() == nu.p() ? a_same_patches(lambda, nu, p) : a_different_patches(lambda, nu, p);
  }

  template <class IBASIS, unsigned int DIM>
  double
  EllipticEquation<IBASIS,DIM>::f(const typename AggregatedFrame<IBASIS,DIM>::Index& lambda) const
  {
//     double r = 0.0;
    
//     const unsigned int p = lambda.p();
    
// //     FixedArray1D<int,DIM> j;
// //     for (unsigned int i = 0; i < DIM; i++) {
// //       j[i] = lambda.j() + lambda.e()[i];
// //     }
    
//     typedef typename IBASIS::Index Index1D;
      
//     FixedArray1D<int,DIM> k1;
//     FixedArray1D<int,DIM> k2;
    
//     //get all supports of 1D functions
//     //of the corresponding cube wavelet
//     for (unsigned int i = 0; i < DIM; i++) {
//       IBASIS* basis_p = frame_->bases()[p]->bases()[i];
//       support(*basis_p,
// 	      Index1D(lambda.j()+lambda.e()[i],lambda.e()[i],lambda.k()[i],basis_p),
// 	      k1[i], k2[i]);
//     }

//     // Set up Gauss points and weights for a composite quadrature formula:
//     const unsigned int N_Gauss = 5;

//     FixedArray1D<Array1D<double>*,DIM> gauss_points;
//     for (unsigned int i = 0; i < DIM; i++) {
//       //allocate memory for gaussian knots
//       gauss_points[i] = new Array1D<double>(N_Gauss*(k2[i]-k1[i]));
//       const double h = ldexp(1.0, -(lambda.j()+lambda.e()[i]));
//       for (int patch = k1[i]; patch < k2[i]; patch++) // refers to 2^{-j}[patch,patch+1]
// 	for (unsigned int n = 0; n < N_Gauss; n++)
// 	  (*gauss_points[i])[(patch-(k1[i]))*N_Gauss+n] = h*(2*patch+1+GaussPoints[N_Gauss-1][n])/2;
//     }



//     FixedArray1D<Array1D<double>*,DIM> wav_values;
//     //compute point values of the
//     //cube wavelet ebcoded by index lambda
//     //evaluate(frame_->bases()[p], wav_values );

//     Point<DIM> x;
//     //calculate approximate integral value
//     for (unsigned int i = 0; i < DIM; i++) {
//       const double h = ldexp(1.0, -(lambda.j()+lambda.e()[i]));
//       // - add all integral shares
//       for (int patch = k1[i], id = 0; patch < k2[i]; patch++)
// 	for (unsigned int n = 0; n < N_Gauss; n++, id++) {
// 	  const double t = (*gauss_points[i])[id];
// 	  const double gauss_weight = GaussWeights[N_Gauss-1][n] * h;
// 	}
//     }



//     //deleting temporary
//     for (unsigned int i = 0; i < DIM; i++) {
//       delete gauss_points[i];
//     }

//#######################################################################

    double r = 0;

    const unsigned int p = lambda.p();

    //typedef WaveletTL::MappedCubeBasis<IBASIS,DIM> MAPPEDCUBEBASIS;
    typedef WaveletTL::CubeBasis<IBASIS,DIM> CUBEBASIS;
    // first compute supp(psi_lambda)
    //typename MAPPEDCUBEBASIS::Support supp;
    typename CUBEBASIS::Support supp;
    
    typename CubeBasis<IBASIS,DIM>::Index lambda_c(lambda.j(),
						   lambda.e(),
						   lambda.k(),frame_->bases()[p]);
    

    //    WaveletTL::support<IBASIS,DIM,MAPPEDCUBEBASIS>(*(frame_->bases()[p]), lambda_c, supp);
    WaveletTL::support<IBASIS,DIM,CUBEBASIS>(*(frame_->bases()[p]), lambda_c, supp);
    cout << supp.j << endl;
    for (int i = 0; i < DIM; i++)
      cout << (supp.a)[i] << " " << (supp.b)[i] << endl;



    // setup Gauss points and weights for a composite quadrature formula:
    const int N_Gauss = 2;
    const double h = ldexp(1.0, -supp.j); // granularity for the quadrature
    FixedArray1D<Array1D<double>,DIM> gauss_points, gauss_weights, v_values;
    for (unsigned int i = 0; i < DIM; i++) {
      gauss_points[i].resize(N_Gauss*(supp.b[i]-supp.a[i]));
      gauss_weights[i].resize(N_Gauss*(supp.b[i]-supp.a[i]));
      for (int patch = supp.a[i]; patch < supp.b[i]; patch++)
	for (int n = 0; n < N_Gauss; n++) {
	  gauss_points[i][(patch-supp.a[i])*N_Gauss+n]
	    = h*(2*patch+1+GaussPoints[N_Gauss-1][n])/2.;
	  gauss_weights[i][(patch-supp.a[i])*N_Gauss+n]
	    = GaussWeights[N_Gauss-1][n];
	}
    }
//     cout << "f(), Gauss points and weights set up:" << endl;
//     for (unsigned int i = 0; i < DIM; i++)
//       cout << "i=" << i << " with points "
// 	   << gauss_points[i] << " and weights " << gauss_weights[i] << endl;

    // compute the point values of the integrand (where we use that it is a tensor product)
    for (unsigned int i = 0; i < DIM; i++) {
      WaveletTL::evaluate(*(frame_->bases()[p]->bases()[i]), 0,
	       typename IBASIS::Index(lambda.j(),
				      lambda.e()[i],
				      lambda.k()[i],
				      frame_->bases()[p]->bases()[i]),
	       gauss_points[i], v_values[i]);
      cout << v_values[i] << endl;
      cout << "#########################" << endl;
    }

    // iterate over all points and sum up the integral shares
    int index[DIM]; // current multiindex for the point values
    for (unsigned int i = 0; i < DIM; i++)
      index[i] = 0;
    
    Point<DIM> x;
    Point<DIM> x_patch;
    
    while (true) {
      for (unsigned int i = 0; i < DIM; i++)
	x[i] = gauss_points[i][index[i]];
      //cout << x << endl;
      frame_->atlas()->charts()[p]->map_point(x,x_patch);
      //cout << x_patch << endl;
      double share = ell_bvp_.f(x_patch) * frame_->atlas()->charts()[p]->Gram_factor(x);
      //cout << "Gram = " << frame_->atlas()->charts()[p]->Gram_factor(x) << endl;
      for (unsigned int i = 0; i < DIM; i++)
	share *= gauss_weights[i][index[i]] * v_values[i][index[i]];
      r += share;

      // "++index"
      bool exit = false;
      for (unsigned int i = 0; i < DIM; i++) {
	if (index[i] == N_Gauss*(supp.b[i]-supp.a[i])-1) {
	  index[i] = 0;
	  exit = (i == DIM-1);
	} else {
	  index[i]++;
	  break;
	}
      }
      if (exit) break;
    }
    
    return r;
  }
  
}
