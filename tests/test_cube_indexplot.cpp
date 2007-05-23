#include <cube/cube_basis.h>
#include <cube/cube_support.h>

#include <geometry/sampled_mapping.h>

using namespace std;
using namespace WaveletTL;

class Plateau
  : public Function<2>
{
  public:
  inline double value(const Point<2>& p,
		      const unsigned int component = 0) const
  {
    return std::max(0.0,0.5-abs(p[0]-0.5))
      * std::max(0.0,0.5-abs(p[1]-0.5));
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

#if 1
  Grid<2> grid(Point<2>(0.0, 0.0), Point<2>(1.0, 1.0), 4, 4);
  SampledMapping<2> mapping(grid, *f);
  mapping.matlab_output(cout);
#endif

  if (f) delete f;

  return 0;
}
