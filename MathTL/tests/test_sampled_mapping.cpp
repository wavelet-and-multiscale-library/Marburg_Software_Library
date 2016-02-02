#include <iostream>
#include <functional>
#include <utils/function.h>
#include <geometry/grid.h>
#include <geometry/sampled_mapping.h>

using std::cout;
using std::endl;

using namespace MathTL;

class TestFunction : public Function<1>
{
public:
  TestFunction() : Function<1>(1) {}
  virtual ~TestFunction(){}
  double value(const Point<1>& p,
	       const unsigned int component = 0) const
  {
    return p[0]*p[0];
  }
  
  void vector_value(const Point<1> &p,
		    Vector<double>& values) const
  {
    values.resize(1, false);
    values[0] = value(p);
  }
};

class TestFunction2 : public Function<2>
{
public:
  TestFunction2() : Function<2>(1) {}
  virtual ~TestFunction2(){}
  double value(const Point<2>& p,
	       const unsigned int component = 0) const
  {
    return p[0]*p[1];
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
  cout << "Testing SampledMapping<1> class ..." << endl;
  SampledMapping<1> f;
  cout << "- Matlab output of empty 1D mapping:" << endl;
  f.matlab_output(cout);

  cout << "- Matlab output of a sampled function on a 1D grid:" << endl;
  SampledMapping<1> g(Grid<1>(0.0, 1.0, 5), TestFunction());
  g.matlab_output(cout);

  cout << "- Gnuplot output for the same parameters:" << endl;
  g.gnuplot_output(cout);

  cout << "- testing SampledMapping<1>::add():" << endl;
  g.add(SampledMapping<1>(Grid<1>(0.0, 1.0, 5), TestFunction()));
  g.matlab_output(cout);

  cout << "- Matlab output of a sampled function on a 2D grid:" << endl;
  Grid<2> grid(Point<2>(0.0, 0.0), Point<2>(1.0, 1.0), 4, 4);
  SampledMapping<2> h(grid, TestFunction2());
  h.matlab_output(cout);

  cout << "- testing SampledMapping<2>::add():" << endl;
  h.add(-2.0, SampledMapping<2>(grid, TestFunction2()));
  h.matlab_output(cout);

  return 0;
}
