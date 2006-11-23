#include <iostream>
#include <fstream>
#include <sstream>
#include <algebra/vector.h>
#include <utils/function.h>
#include <interval/spline_basis.h>
#include <galerkin/full_gramian.h>
#include <Rd/cdf_utils.h>
#include <numerics/iteratsolv.h>
#include <numerics/quadrature.h>
#include <numerics/schoenberg_splines.h>
#include <interval/interval_bspline.h>

using namespace std;
using namespace MathTL;
using namespace WaveletTL;

class Function1 : public Function<1> {
public:
  inline double value(const Point<1>& p, const unsigned int component = 0) const {
    return 0.5*p[0]*(1-p[0]);
  }
  
  void vector_value(const Point<1> &p, Vector<double>& values) const {
    values.resize(1, false);
    values[0] = value(p);
  }
};

int main()
{
  cout << "Testing FullGramian ..." << endl;

  const unsigned int d = 2;
  const unsigned int dT = 2;

  SplineBasis<d,dT> basis("P","",1,1,0,0); // PBasis, complementary b.c.'s
  FullGramian<d,dT> G(basis);

  cout << "* Gramian matrix on coarsest level j0=" << basis.j0() << ":" << endl
       << G;
  
  G.set_level(basis.j0()+1);
  cout << "* Gramian matrix on next level j0+1=" << basis.j0()+1 << ":" << endl
       << G;

  const unsigned int testcase=1;
  Function<1>* f = 0;

  switch(testcase) {
  case 1:
    f = new Function1();
    break;
  default:
    break;
  }

  cout << "* compute approximate expansions of the test function for several levels..." << endl;
  const int jmin = basis.j0();
  const int jmax = jmin;
//   const int jmax = 18;
  Vector<double> js(jmax-jmin+1);
  Vector<double> Linfty_errors(jmax-jmin+1), L2_errors(jmax-jmin+1);

  for (int j = jmin; j <= jmax; j++) {
    cout << "  j=" << j << ":" << endl;
    js[j-jmin] = j;
    
    G.set_level(j);

//     // setup rhs in the phi_{j,k} basis,
//     Vector<double> rhs_phijk(delta.row_dimension());
//     if (solution == 1) {
//       if (d == 2) {
//  	// exact right-hand side is known
//  	rhs_phijk = sqrt(ldexp(1.0, -j));
//       } else {
// 	// perform quadrature with a composite rule on [0,1]
// 	cout << "  solution 1, quadrature for rhs..." << endl;
// 	SimpsonRule simpson;
// 	CompositeRule<1> composite(simpson, 12);
// 	SchoenbergIntervalBSpline_td<d> sbs(j,0);
// 	for (int k = basis.DeltaLmin(); k <= basis.DeltaRmax(j); k++) {
// 	  sbs.set_k(k);
// 	  rhs_phijk[k-basis.DeltaLmin()]
//  	    = composite.integrate(sbs, // against f=1
//  				  Point<1>(std::max(0.0, (k+ell1<d>())*ldexp(1.0, -j))),
//  				  Point<1>(std::min(1.0, (k+ell2<d>())*ldexp(1.0, -j))));
//  	}
//       }


  }
  
  if (f) delete f;

  return 0;
}
