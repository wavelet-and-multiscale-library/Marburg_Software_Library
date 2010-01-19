#include <iostream>
#include <utils/function.h>
#include <interval/p_basis.h>
#include <interval/p_expansion.h>
#include <interval/i_indexplot.h>

using namespace std;
using namespace MathTL;
using namespace WaveletTL;

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

class CutSqrt
  : public Function<1>
{
  public:
  inline double value(const Point<1>& p,
		      const unsigned int component = 0) const
  {
    return (p[0] < 0.1 || p[0] > 0.75
	    ? 0
	    : sqrt(p[0]-0.1));
  }
  
  void vector_value(const Point<1> &p,
		    Vector<double>& values) const
  {
    values.resize(1, false);
    values[0] = value(p);
  }
};

class CutPwr
  : public Function<1>
{
  public:
  inline double value(const Point<1>& p,
		      const unsigned int component = 0) const
  {
    return (p[0] < 0.1 || p[0] > 0.75
	    ? 0
	    : pow(p[0]-0.1,0.75));
  }
  
  void vector_value(const Point<1> &p,
		    Vector<double>& values) const
  {
    values.resize(1, false);
    values[0] = value(p);
  }
};

class CutLin
  : public Function<1>
{
  public:
  inline double value(const Point<1>& p,
		      const unsigned int component = 0) const
  {
    return (p[0] < 0.1 || p[0] > 0.75
	    ? 0
	    : 2*(p[0]-0.1));
  }
  
  void vector_value(const Point<1> &p,
		    Vector<double>& values) const
  {
    values.resize(1, false);
    values[0] = value(p);
  }
};

// kink at 0<a<1
class Kink : public Function<1> {
 public:
  Kink(const double a = 0.5) : a_(a) {}
  
  inline double value(const Point<1>& p, const unsigned int component = 0) const {
    if (0. <= p[0] && p[0] < a_)
      return 1/(2*a_*a_)*p[0]*p[0];
    
    if (a_ <= p[0] && p[0] <= 1.0)
      return 0.5*(1-(p[0]-a_)/(1-a_))*(1-(p[0]-a_)/(1-a_));
    
    return 0.;
  }
  
  void vector_value(const Point<1> &p, Vector<double>& values) const {
    values.resize(1, false);
    values[0] = value(p);
  }
  
 protected:
  double a_;
};

int main() {
  cout << "* expand some interesting functions into their wavelet series..." << endl;

//   Function<1>* f = new Hat();
  Function<1>* f = new CutLin();
//   Function<1>* f = new CutSqrt();
//   Function<1>* f = new CutPwr();
//   Function<1>* f = new Kink(0.5);

  const int d = 2;
  const int dt = 2;

  typedef PBasis<d,dt> Basis;
  typedef Basis::Index Index;
  Basis basis(0,0); // PBasis, no b.c.'s

//   const int j0 = basis.j0();
  const int jmax = 9;

  InfiniteVector<double,Index> coeffs;
  expand(f, basis, false, jmax, coeffs);
  cout << "- (approx.) expansion coefficients in the primal basis:" << endl
       << coeffs;

  cout << "* plotting the coefficient set..." << endl;
  std::ofstream plotstream;
  plotstream.open("coefficient_plot.m");
//   plot_indices(&basis, coeffs, 10, plotstream, "cool", true, true, -5);
//   plot_indices(&basis, coeffs, 9, plotstream, "invgray", true, true, -8);
  plot_indices(&basis, coeffs, 9, plotstream, "invgray", false, true, -8);
  plotstream.close();
  cout << "  ...done!" << endl;

  delete f;

  return 0;
}
