// implementation for ds_expansion.h

namespace WaveletTL
{
  template <int d, int dT>
  double integrate(const Function<1>* f,
		   const DSBasis<d,dT>& basis,
		   const typename DSBasis<d,dT>::Index& lambda)
  {
    double r = 0;
    
    const int j = lambda.j()+lambda.e();
    int k1, k2;
    support(basis, lambda, k1, k2);
    
    // Set up Gauss points and weights for a composite quadrature formula:
    const unsigned int N_Gauss = 5;
    const double h = ldexp(1.0, -j);
    Array1D<double> gauss_points (N_Gauss*(k2-k1));
    for (int patch = k1; patch < k2; patch++) // refers to 2^{-j}[patch,patch+1]
      for (unsigned int n = 0; n < N_Gauss; n++)
	gauss_points[(patch-k1)*N_Gauss+n] = h*(2*patch+1+GaussPoints[N_Gauss-1][n])/2;
    
    //   cout << "gauss points: " << gauss_points << endl;
    
    // - add all integral shares
    for (unsigned int n = 0; n < N_Gauss; n++)
      {
	const double gauss_weight = GaussWeights[N_Gauss-1][n] * h;
	for (int patch = k1; patch < k2; patch++)
	  {
	    const double t = gauss_points[(patch-k1)*N_Gauss+n];
	    
	    const double ft = f->value(Point<1>(t));
	    if (ft != 0)
	      r += ft
		* evaluate(basis, 0, lambda, t)
		* gauss_weight;
	  }
      }
    
    return r;
  }
  
  template <int d, int dT>
  void
  expand(const Function<1>* f,
	 const DSBasis<d,dT>& basis,
	 const bool primal,
	 const int jmax,
	 InfiniteVector<double, typename DSBasis<d,dT>::Index>& coeffs)
  {
    typedef typename DSBasis<d,dT>::Index Index;
    const int j0 = basis.j0();

    for (Index lambda = first_generator(&basis, j0);;++lambda)
      {
	coeffs.set_coefficient(lambda, integrate(f, basis, lambda));
	if (lambda == last_wavelet(&basis, jmax))
	  break;
      }

    if (!primal) {
    }
  }
}
