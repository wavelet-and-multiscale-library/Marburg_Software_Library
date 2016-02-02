#include <iostream>
#include <fstream>

#include <geometry/sampled_mapping.h>
#include <utils/fixed_array1d.h>

#include <interval/p_basis.h>
#include <cube/cube_basis.h>
#include <cube/cube_support.h>
#include <cube/cube_expansion.h>
#include <cube/cube_indexplot.h>

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

class Jump_and_Singularity
  : public Function<2>
{
  public:
  inline double value(const Point<2>& p,
		      const unsigned int component = 0) const
  {
    Point<2> x;
    x[0] = 0.5;
    x[1] = 0.5;
    
    return atan2((p-x)[0],(p-x)[1]);
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

  Function<2>* f  = new Plateau();
  Function<2>* f2 = new Jump_and_Singularity();

#if 1
  Grid<2> grid(Point<2>(0.0, 0.0), Point<2>(1.0, 1.0), 1<<5, 1<<5);
  SampledMapping<2> mapping(grid, *f);
  ofstream file("plateau_function.m");
  mapping.matlab_output(file);
  file.close();

  Grid<2> grid2(Point<2>(0.0, 0.0), Point<2>(1.0, 1.0), 1<<5, 1<<5);
  SampledMapping<2> mapping2(grid, *f2);
  ofstream file2("Jump_and_Singulartity_function.m");
  mapping2.matlab_output(file2);
  file2.close();


#endif

  typedef PBasis<2,2> Basis1D;
  Basis1D basis1d(1,1); // PBasis, complementary b.c.'s
  
  FixedArray1D<Basis1D*,2> bases1d;
  bases1d[0] = bases1d[1] = &basis1d;

  typedef CubeBasis<Basis1D> Basis;
  typedef Basis::Index Index;
  
  Basis basis(bases1d);

  const int jmax = 4;

  InfiniteVector<double,Index> coeffs;
  expand(f, basis, false, jmax, coeffs);
  cout << "- (approx.) expansion coefficients in the primal basis:" << endl
       << coeffs;

  InfiniteVector<double,Index> coeffs2;
  expand(f2, basis, false, jmax, coeffs2);
  cout << "- (approx.) expansion coefficients in the primal basis:" << endl
       << coeffs2;


  std::ofstream plotstream;
  plotstream.open("coefficient_plot.m");

  plot_indices(&basis,coeffs, jmax, plotstream, "jet", true, true);
  plotstream.close();
  if (f) delete f;


  std::ofstream plotstream2;
  plotstream2.open("coefficient_plot2.m");

  plot_indices(&basis,coeffs2, jmax, plotstream2, "jet", false, true);
  plotstream2.close();
  if (f2) delete f2;


  return 0;
}
