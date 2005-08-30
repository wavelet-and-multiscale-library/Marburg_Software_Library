#include <iostream>

#include <utils/function.h>
#include <algebra/vector.h>
#include <numerics/gauss_data.h>
#include <geometry/sampled_mapping.h>

#include <interval/ds_basis.h>

using namespace std;
using namespace WaveletTL;

class TestPolynomial
  : public Function<1>
{
  public:
  inline double value(const Point<1>& p,
		      const unsigned int component = 0) const
  {
    return 1;
//     return -p[0]*p[0]+1;
  }
  
  void vector_value(const Point<1> &p,
		    Vector<double>& values) const
  {
    values.resize(1, false);
    values[0] = value(p);
  }
};

// integrate a (smooth) function f against a primal DS generator or wavelet
template <int d, int dT>
double integrate(const Function<1>& f, const DSBasis<d,dT>& basis, const typename DSBasis<d,dT>::Index& lambda)
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
	  
	  const double ft = f.value(Point<1>(t));
	  if (ft != 0)
	    r += ft
	      * evaluate(basis, 0, lambda, t)
	      * gauss_weight;
	}
    }
  
  return r;
}

int main()
{
  cout << "Test expansion of a polynomial in a DS basis..." << endl;

  TestPolynomial p;
  cout << "- sample values of a test polynomial:" << endl;
  SampledMapping<1> s(Grid<1>(-1.0, 1.0, 10), p);
  s.matlab_output(cout);

  const int d = 2;
  const int dT = 2;

  typedef DSBasis<d,dT> Basis;
  typedef Basis::Index Index;

  Basis basis;
  InfiniteVector<double,Index> coeffs;

  const int j0 = basis.j0();
//   const int jmax = j0;
  
  for (Index lambda = first_generator(&basis, j0);;++lambda)
    {
      cout << "calculating integral against psi_lambda with lambda=" << lambda << endl;
      coeffs.set_coefficient(lambda, integrate(p,basis,lambda));
      if (lambda == last_generator(&basis, j0))
//   	if (lambda == last_wavelet(&basis, jmax))
	break;
    }
  
  cout << "- integrals of p against all primal generators on level j0:" << endl
       << coeffs << endl;

  cout << "- evaluation of this linear combination of dual generators yields the pointwise error on [-1,1]:" << endl;
  SampledMapping<1> s2(evaluate(basis, coeffs, false, 5));
  Vector<double> error(s2.points().size());
  for (unsigned int i = 0; i < error.size(); i++)
    error[i] = fabs(s2.values()[i]-p.value(Point<1>(s2.points()[i])));
  cout << error << endl;
  cout << "(max. error: " << linfty_norm(error) << ")" << endl;
  
  return 0;
}
