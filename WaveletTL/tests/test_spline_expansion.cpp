#include <iostream>
#include <time.h>

#include <utils/function.h>
#include <algebra/vector.h>
#include <numerics/gauss_data.h>
#include <geometry/sampled_mapping.h>

#include <interval/spline_basis.h>
#include <interval/spline_expansion.h>

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

class Gaussian
  : public Function<1> {
public:
  inline double value(const Point<1>& p, const unsigned int component = 0) const {
    const double x = p[0];
    const double t = get_time();
    return exp(-300*(x-0.6+0.2*t)*(x-0.6+0.2*t));
  }
  
  void vector_value(const Point<1> &p, Vector<double>& values) const {
    values.resize(1, false);
    values[0] = value(p);
  }
};

class Quadratic
  : public Function<1> {
public:
  inline double value(const Point<1>& p, const unsigned int component = 0) const {
    const double x = p[0];
    return x*(1-x);
  }
  
  void vector_value(const Point<1> &p, Vector<double>& values) const {
    values.resize(1, false);
    values[0] = value(p);
  }
};

int main()
{
  cout << "Test expansion of a polynomial in a SplineBasis..." << endl;

  TestPolynomial p;
  cout << "- sample values of a test polynomial:" << endl;
  SampledMapping<1> s(Grid<1>(-1.0, 1.0, 10), p);
  s.matlab_output(cout);

#if 1
//   const int d = 2, dt = 2, s0 = 0, s1 = 0, sT0 = 0, sT1 = 0; // PBasis, no b.c.'s
//   const int d = 2, dt = 2, s0 = 1, s1 = 1, sT0 = 0, sT1 = 0; // PBasis, homogeneous b.c.'s
//   const int d = 2, dt = 2, s0 = 1, s1 = 0, sT0 = 0, sT1 = 0; // PBasis, homogeneous b.c. at x=0
//   const int d = 2, dt = 2, s0 = 0, s1 = 1, sT0 = 0, sT1 = 0; // PBasis, homogeneous b.c. at x=1
  const int d = 3, dt = 3, s0 = 1, s1 = 1, sT0 = 0, sT1 = 0; // PBasis, homogeneous b.c.'s
  typedef SplineBasis<d,dt,P_construction,s0,s1,sT0,sT1,SplineBasisData_j0<d,dt,P_construction,s0,s1,sT0,sT1>::j0> Basis; 
#else
//   typedef SplineBasis<3,3,DS_construction_bio5,0,0,0,0> Basis; // DSBasis, no b.c.'s
  typedef SplineBasis<3,3,DS_construction_bio5e,0,0,0,0> Basis; // DSBasis, no b.c.'s
#endif

  Basis basis;
 
  typedef Basis::Index Index;

  InfiniteVector<double,Index> coeffs;
  
  const int j0 = basis.j0();
  const int jmax = 10;
  
  expand(&p, basis, true, j0, coeffs);
  cout << "- integrals of p against all primal generators on level j0:" << endl
       << coeffs << endl;
 
  Hat hat;
  cout << "- sample values of the hat function:" << endl;
  SampledMapping<1> shat(Grid<1>(0.0, 1.0, 10), hat);
  shat.matlab_output(cout);
  
  InfiniteVector<double,Index> dual_coeffs;
  expand(&hat, basis, false, jmax, dual_coeffs);
  cout << "- (approx.) expansion coefficients of the hat fct. in the primal basis:" << endl
       << dual_coeffs;
  
  cout << "- pointwise error:" << endl;
  SampledMapping<1> s3(basis.evaluate(dual_coeffs, jmax+1));
  Vector<double> error(s3.points().size());
  for (unsigned int i = 0; i < error.size(); i++)
    error[i] = fabs(s3.values()[i]-hat.value(Point<1>(s3.points()[i])));
  //   cout << error << endl;
  cout << "(max. error: " << linfty_norm(error) << ")" << endl;

  cout << "- expand a quadratic polynomial..." << endl;
  Quadratic q;
  InfiniteVector<double,Index> q_coeffs;
  basis.expand(&q, false, jmax, q_coeffs);
  cout << "- pointwise error:" << endl;
  SampledMapping<1> qs(basis.evaluate(q_coeffs, jmax+1));
  error.resize(qs.points().size());
  for (unsigned int i = 0; i < error.size(); i++)
    error[i] = fabs(qs.values()[i]-q.value(Point<1>(qs.points()[i])));
  cout << "(max. error: " << linfty_norm(error) << ")" << endl;
  
  cout << "- expand a Gaussian..." << endl;
  clock_t tstart, tend;
  Gaussian g;
  
  g.set_time(0);
  InfiniteVector<double,Index> g_coeffs;
  tstart = clock();
  expand(&g, basis, false, jmax, g_coeffs);
  tend = clock();
  cout << "  ... done, time needed: "
       << (double)(tend-tstart)/CLOCKS_PER_SEC
       << "s" << endl;
  cout << "- pointwise error:" << endl;
  SampledMapping<1> gs(basis.evaluate(g_coeffs, jmax+1));
  error.resize(gs.points().size());
  for (unsigned int i = 0; i < error.size(); i++)
    error[i] = fabs(gs.values()[i]-g.value(Point<1>(gs.points()[i])));
//   cout << error << endl;
  cout << "(max. error: " << linfty_norm(error) << ")" << endl;
  

  return 0;
}
