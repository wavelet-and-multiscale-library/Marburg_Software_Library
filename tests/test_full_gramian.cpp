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

// polynomial
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

// function with kink at x=0.5
class Function2 : public Function<1> {
public:
  inline double value(const Point<1>& p, const unsigned int component = 0) const {
    if (0. <= p[0] && p[0] < 0.5)
      return -sin(3.*M_PI*p[0]) + 2.*p[0]*p[0];
    
    if (0.5 <= p[0] && p[0] <= 1.0)
      return -sin(3.*M_PI*p[0]) + 2.*(1-p[0])*(1-p[0]);
    
    return 0.;
  }
  
  void vector_value(const Point<1> &p, Vector<double>& values) const {
    values.resize(1, false);
    values[0] = value(p);
  }
};

// smooth part of function 2
class Function3 : public Function<1> {
public:
  inline double value(const Point<1>& p, const unsigned int component = 0) const {
    return -sin(3.*M_PI*p[0]);
  }
  
  void vector_value(const Point<1> &p, Vector<double>& values) const {
    values.resize(1, false);
    values[0] = value(p);
  }
};


int main()
{
  cout << "Testing FullGramian ..." << endl;

  const unsigned int d = 3;
  const unsigned int dT = 3;

  SplineBasis<d,dT> basis("P","",1,1,0,0); // PBasis, complementary b.c.'s
  FullGramian<d,dT> G(basis);

  cout << "* Gramian matrix on coarsest level j0=" << basis.j0() << ":" << endl
       << G;
  
  G.set_level(basis.j0()+1);
  cout << "* Gramian matrix on next level j0+1=" << basis.j0()+1 << ":" << endl
       << G;

  const unsigned int testcase=3;
  Function<1>* u = 0;

  switch(testcase) {
  case 1:
    u = new Function1();
    break;
  case 2:
    u = new Function2();
    break;
  case 3:
    u = new Function3();
    break;
  default:
    break;
  }

  cout << "* compute approximate expansions of the test function for several levels..." << endl;
  const int jmin = basis.j0();
//   const int jmax = jmin;
  const int jmax = 20;
  Vector<double> js(jmax-jmin+1);
  Vector<double> Linfty_errors(jmax-jmin+1), L2_errors(jmax-jmin+1);

  for (int j = jmin; j <= jmax; j++) {
    cout << "  j=" << j << ":" << endl;
    js[j-jmin] = j;
    
    G.set_level(j);

    // compute integrals w.r.t. the primal generators on level j
    Vector<double> coeffs_phijk(G.row_dimension());
    SimpsonRule simpson;
    CompositeRule<1> composite(simpson, 12);
    SchoenbergIntervalBSpline_td<d> sbs(j,0);
    for (int k = basis.DeltaLmin(); k <= basis.DeltaRmax(j); k++) {
      sbs.set_k(k);
      ProductFunction<1> integrand(u, &sbs);
      coeffs_phijk[k-basis.DeltaLmin()]
 	= composite.integrate(integrand,
 			      Point<1>(std::max(0.0, (k+ell1<d>())*ldexp(1.0, -j))),
 			      Point<1>(std::min(1.0, (k+ell2<d>())*ldexp(1.0, -j))));
    }
//     cout << "  inner products against phi_{j,k} basis: " << coeffs_phijk << endl;

    // transform rhs into that of psi_{j,k} basis: apply T_{j-1}^T
    Vector<double> rhs(G.row_dimension());
    if (j == basis.j0())
      rhs = coeffs_phijk;
    else
      basis.apply_Tj_transposed(j-1, coeffs_phijk, rhs);
//     cout << "  inner products against psi_{j,k} basis: " << rhs << endl;

    // solve Gramian system
    Vector<double> uj(G.row_dimension()), residual(G.row_dimension()); uj = 0;
    unsigned int iterations;
    CG(G, rhs, uj, 1e-15, 250, iterations);
//     cout << "  solution coefficients: " << uj << endl;
    cout << "  Galerkin system solved with residual (infinity) norm ";
    G.apply(uj, residual);
    residual -= rhs;
    cout << linfty_norm(residual) << endl;

    // compute coefficients of uj in the phi_{j,k} basis: apply T_{j-1}
    Vector<double> uj_phijk(G.row_dimension());
    if (j == basis.j0())
      uj_phijk = uj;
    else
      basis.apply_Tj(j-1, uj, uj_phijk);
//     cout << "  solution coefficients in phi_{j,k} basis: " << uj_phijk << endl;

    // evaluate linear combination of Schoenberg B-splines on a grid
    const unsigned int N = 100;
    const double h = 1./N;
    Vector<double> uj_values(N+1);
    for (unsigned int i = 0; i <= N; i++) {
      const double x = i*h;
      for (unsigned int k = 0; k < G.row_dimension(); k++) {
	sbs.set_k(basis.DeltaLmin()+k);
	uj_values[i] +=
	  uj_phijk[k] * sbs.value(Point<1>(x));
      }
    }
//     cout << "  point values of Galerkin solution: " << uj_values << endl;

    // evaluate exact solution
    Vector<double> uexact_values(N+1);
    for (unsigned int i = 0; i <= N; i++) {
      const double x = i*h;
      uexact_values[i] = u->value(Point<1>(x));
    }
//     cout << "  point values of exact solution: " << uexact_values << endl;

    // compute some errors
    const double Linfty_error = linfty_norm(uj_values-uexact_values);
    cout << "  L_infinity error on a subgrid: " << Linfty_error << endl;
    Linfty_errors[j-jmin] = Linfty_error;

    const double L2_error = sqrt(l2_norm_sqr(uj_values-uexact_values)*h);
    cout << "  L_2 error on a subgrid: " << L2_error << endl;
    L2_errors[j-jmin] = L2_error;
  }

#if 1
  // write Galerkin errors to a file
  ostringstream filename;
  filename << "full_gramian_errors_" << d << "_" << testcase << ".m";
  ofstream galerkin_stream(filename.str().c_str());
  galerkin_stream << "js=" << js << ";" << endl
 		  << "Linfty_errors=" << Linfty_errors << ";" << endl
 		  << "L2_errors=" << L2_errors << ";" << endl;
  galerkin_stream.close();
#endif
  
  if (u) delete u;

  return 0;
}
