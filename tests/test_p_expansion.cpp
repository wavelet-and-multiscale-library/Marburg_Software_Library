#include <iostream>

#include <utils/function.h>
#include <algebra/vector.h>
#include <numerics/gauss_data.h>
#include <geometry/sampled_mapping.h>

#include <interval/p_basis.h>
#include <interval/p_expansion.h>

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

class Hat
  : public Function<1>
{
  public:
  inline double value(const Point<1>& p,
		      const unsigned int component = 0) const
  {
    return std::max(0.0,0.5-abs(p[0]-0.5));
  }
  
  void vector_value(const Point<1> &p,
		    Vector<double>& values) const
  {
    values.resize(1, false);
    values[0] = value(p);
  }
};

int main()
{
  cout << "Test expansion of a polynomial in a DS basis..." << endl;

  TestPolynomial p;
  cout << "- sample values of a test polynomial:" << endl;
  SampledMapping<1> s(Grid<1>(-1.0, 1.0, 10), p);
  s.matlab_output(cout);

  const int d = 2;
  const int dT = 6; // be sure to use a continuous dual here!

  typedef PBasis<d,dT> Basis;
  typedef Basis::Index Index;

   Basis basis(0, 0); // no b.c.'s
//   Basis basis(1, 0); // complementary b.c. at x=0
//   Basis basis(0, 1); // complementary b.c. at x=1
// Basis basis(1, 1); // complementary b.c.'s

  InfiniteVector<double,Index> coeffs;

  const int j0 = basis.j0();
  const int jmax = 8;

  expand(&p, basis, true, j0, coeffs);
  cout << "- integrals of p against all primal generators on level j0:" << endl
       << coeffs << endl;

  cout << "- evaluation of this linear combination of dual generators yields the pointwise error on [-1,1]:" << endl;
  SampledMapping<1> s2(evaluate(basis, coeffs, false, 6));
  Vector<double> error(s2.points().size());
  for (unsigned int i = 0; i < error.size(); i++)
    error[i] = fabs(s2.values()[i]-p.value(Point<1>(s2.points()[i])));
  cout << error << endl;
  cout << "(max. error: " << linfty_norm(error) << ")" << endl;
  
  cout << endl;

  Hat hat;
  cout << "- sample values of the hat function:" << endl;
  SampledMapping<1> shat(Grid<1>(0.0, 1.0, 10), hat);
  shat.matlab_output(cout);
  
  InfiniteVector<double,Index> dual_coeffs;
  expand(&hat, basis, false, jmax, dual_coeffs);
//   cout << "- (approx.) expansion coefficients of the hat function in the primal basis:" << endl
//        << dual_coeffs;
  
  cout << "- pointwise error:" << endl;
  SampledMapping<1> s3(evaluate(basis, dual_coeffs, true, jmax+1));
  error.resize(s3.points().size());
  for (unsigned int i = 0; i < error.size(); i++)
    error[i] = fabs(s3.values()[i]-hat.value(Point<1>(s3.points()[i])));
  cout << error << endl;
  cout << "(max. error: " << linfty_norm(error) << ")" << endl;
  
  return 0;
}
