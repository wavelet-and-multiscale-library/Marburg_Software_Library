#include <iostream>

#include <utils/function.h>
#include <algebra/vector.h>
#include <numerics/gauss_data.h>

#include <Rd/cdf_basis.h>
#include <Rd/cdf_utils.h>

using namespace std;
using namespace WaveletTL;

class TestPolynomial
  : public Function<1>
{
  public:
  inline double value(const Point<1>& p,
		      const unsigned int component = 0) const
  {
    return -p[0]*p[0]+1;
  }
  
  void vector_value(const Point<1> &p,
		    Vector<double>& values) const
  {
    values.resize(1, false);
    values[0] = value(p);
  }
};

// integrate a (smooth) function f against a primal CDF generator or wavelet
template <int d, int dT>
double integrate(const Function<1>& f, const RIndex& lambda)
{
  CDFBasis<d,dT> basis;
  
  double r = 0;

  const int j = lambda.j() + lambda.e();

  const int k1 = psi_supp_left<d,dT>() + lambda.k();
  const int k2 = psi_supp_right<d,dT>() + lambda.k();

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
	  
	  const double ft = f.value(Point<1>(t));
 	  if (ft != 0)
 	    r += ft
 	      * basis.evaluate(0, lambda, t)
 	      * gauss_weight;
	}
    }
  
  return r;
}

int main()
{
  cout << "Test expansion of a polynomial in a dual CDF basis..." << endl;

  TestPolynomial p;
  cout << "- sample values of a test polynomial:" << endl;
  SampledMapping<1> s(Grid<1>(-1.0, 1.0, 10), p);
  s.matlab_output(cout);

//   cout << "- integration of p against phi<2,2> yields:" << endl
//        << integrate<2,2>(p, RIndex(0,0,0)) << endl;
//   cout << "- integration of p against psi<2,2> yields:" << endl
//        << integrate<2,2>(p, RIndex(0,1,0)) << endl;

  const int d = 3;
  const int dT = 3; // be sure to use a continuous dual here!
  CDFBasis<d,dT> basis;
  typedef CDFBasis<d,dT>::Index Index;

  InfiniteVector<double,Index> coeffs;
  const int j = 1;
  for (int k = -12*(1<<j); k <= 12*(1<<j); k++)
    coeffs[Index(j,0,k)] = integrate<d,dT>(p,Index(j,0,k));

   cout << "- integrals of p against some primal generators:" << endl
        << coeffs << endl;

  cout << "- evaluating this linear combination of dual generators yields the pointwise error on [-8,8]" << endl;
  SampledMapping<1> s2(basis.evaluate(0, coeffs, false, -8, 8, 3));
  Vector<double> error(s2.points().size());
  for (unsigned int i = 0; i < error.size(); i++)
    error[i] = fabs(s2.values()[i]-p.value(Point<1>(s2.points()[i])));
  cout << error << endl;
  cout << "(max. error: " << linfty_norm(error) << ")" << endl;
  
  return 0;
}
