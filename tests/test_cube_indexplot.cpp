#include <iostream>
#include <fstream>

#include <geometry/sampled_mapping.h>
#include <utils/fixed_array1d.h>

#include <interval/spline_basis.h>
#include <cube/cube_basis.h>
#include <cube/cube_support.h>
#include <cube/cube_expansion.h>

using namespace std;
using namespace WaveletTL;

class Plateau
  : public Function<2>
{
  public:
  inline double value(const Point<2>& p,
		      const unsigned int component = 0) const
  {
    return (p[0] < 0.25 || p[0] > 0.75 ? 0 : 1)
      * (p[1] < 0.25 || p[1] > 0.75 ? 0 : 1);
  }
  
  void vector_value(const Point<2> &p,
		    Vector<double>& values) const
  {
    values.resize(1, false);
    values[0] = value(p);
  }
};


int main()
{
  cout << "* expand some interesting functions on the square into their wavelet series..." << endl;

  Function<2>* f = new Plateau();

#if 0
  Grid<2> grid(Point<2>(0.0, 0.0), Point<2>(1.0, 1.0), 1<<5, 1<<5);
  SampledMapping<2> mapping(grid, *f);
  ofstream file("zweidfunction.m");
  mapping.matlab_output(file);
  file.close();
#endif

  typedef SplineBasis<2,2,P_construction> Basis1D;
  Basis1D basis1d("",1,1,0,0); // PBasis, complementary b.c.'s
  
  FixedArray1D<Basis1D*,2> bases1d;
  bases1d[0] = bases1d[1] = &basis1d;

  typedef CubeBasis<Basis1D> Basis;
  typedef Basis::Index Index;
  
  Basis basis(bases1d);

  const int jmax = 6;

  InfiniteVector<double,Index> coeffs;
  expand(f, basis, false, jmax, coeffs);
  cout << "- (approx.) expansion coefficients in the primal basis:" << endl
       << coeffs;

  if (f) delete f;

  return 0;
}
