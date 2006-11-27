#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <algebra/vector.h>
#include <utils/function.h>
#include <interval/spline_basis.h>
#include <galerkin/full_laplacian.h>
#include <galerkin/full_gramian.h>
#include <Rd/cdf_utils.h>
#include <numerics/iteratsolv.h>
#include <numerics/quadrature.h>
#include <numerics/schoenberg_splines.h>
#include <interval/interval_bspline.h>

using namespace std;
using namespace MathTL;
using namespace WaveletTL;

// exact solution of -u''=1, u(0)=u(1)=0
class Solution1 : public Function<1> {
public:
  inline double value(const Point<1>& p, const unsigned int component = 0) const {
    return 0.5*p[0]*(1-p[0]);
  }
  
  void vector_value(const Point<1> &p, Vector<double>& values) const {
    values.resize(1, false);
    values[0] = value(p);
  }
};

// exact solution of -u''=f, u(0)=u(1)=0 with kink at x=0.5
class Solution2 : public Function<1> {
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

// smooth part of corresponding right-hand side
class RHS2_part : public Function<1> {
  inline double value(const Point<1>& p, const unsigned int component = 0) const {
    return -sin(3.*M_PI*p[0])*9.*M_PI*M_PI - 4.;
  }
  
  void vector_value(const Point<1> &p, Vector<double>& values) const {
    values.resize(1, false);
    values[0] = value(p);
  }
};

// smooth part of solution 2
class Solution3 : public Function<1> {
public:
  inline double value(const Point<1>& p, const unsigned int component = 0) const {
    return -sin(3.*M_PI*p[0]);
  }
  
  void vector_value(const Point<1> &p, Vector<double>& values) const {
    values.resize(1, false);
    values[0] = value(p);
  }
};

// smooth part of corresponding right-hand side
class RHS3 : public Function<1> {
  inline double value(const Point<1>& p, const unsigned int component = 0) const {
    return -sin(3.*M_PI*p[0])*9.*M_PI*M_PI;
  }
  
  void vector_value(const Point<1> &p, Vector<double>& values) const {
    values.resize(1, false);
    values[0] = value(p);
  }
};

// kink at 0<a<1
class Solution4 : public Function<1> {
 public:
  Solution4(const double a = 0.5) : a_(a) {}
  
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

// smooth part of corresponding right-hand side
class RHS4_part : public Function<1> {
public:
  RHS4_part(const double a = 0.5) : a_(a) {}
  
  inline double value(const Point<1>& p, const unsigned int component = 0) const {
    if (0. <= p[0] && p[0] < a_)
      return -1/(a_*a_);
    
    if (a_ <= p[0] && p[0] <= 1.0)
      return -1/((1-a_)*(1-a_));
    
    return 0.;
  }
  
  void vector_value(const Point<1> &p, Vector<double>& values) const {
    values.resize(1, false);
    values[0] = value(p);
  }
  
protected:
  double a_;
};

int main()
{
  cout << "Testing FullLaplacian ..." << endl;

  const unsigned int d = 3;
  const unsigned int dT = 3;

  SplineBasis<d,dT> basis("P","",1,1,0,0); // PBasis, complementary b.c.'s
  FullLaplacian<d,dT> delta(basis);
  FullGramian<d,dT> G(basis);

  cout << "* stiffness matrix on coarsest level j0=" << basis.j0() << ":" << endl
       << delta;
  
  delta.set_level(basis.j0()+1);
  cout << "* stiffness matrix on next level j0+1=" << basis.j0()+1 << ":" << endl
       << delta;
  
  const unsigned int solution = 2;
  double kink = 0; // for Solution4;

  Function<1> *uexact = 0;
  switch(solution) {
  case 1:
    uexact = new Solution1();
    break;
  case 2:
    uexact = new Solution2();
    break;
  case 3:
    uexact = new Solution3();
    break;
  case 4:
    kink = 0.5;
    uexact = new Solution4(kink);
    break;
  case 5:
    kink = 5./7.;
    uexact = new Solution4(kink);
    break;
  default:
    break;
  }

  // setup (approximate) coefficients of u in the primal basis on a sufficiently high level
  const int jref = 16;
  G.set_level(jref);
  
  // 1. compute integrals w.r.t. the primal generators on level jref
  Vector<double> ucoeffs_phijk(G.row_dimension(), false);
  SimpsonRule simpson;
  CompositeRule<1> composite(simpson, 12);
  SchoenbergIntervalBSpline_td<d> sbs(jref,0);
  for (int k = basis.DeltaLmin(); k <= basis.DeltaRmax(jref); k++) {
    sbs.set_k(k);
    ProductFunction<1> integrand(uexact, &sbs);
    ucoeffs_phijk[k-basis.DeltaLmin()]
      = composite.integrate(integrand,
			    Point<1>(std::max(0.0, (k+ell1<d>())*ldexp(1.0, -jref))),
			    Point<1>(std::min(1.0, (k+ell2<d>())*ldexp(1.0, -jref))));
  }
  // 2. transform rhs into that of psi_{j,k} basis: apply T_{j-1}^T
  Vector<double> urhs(G.row_dimension(), false);
  if (jref == basis.j0())
    urhs = ucoeffs_phijk;
  else
    basis.apply_Tj_transposed(jref-1, ucoeffs_phijk, urhs);
  
  // 3. solve Gramian system, yields approximate L_2 wavelet coefficients of u on level jref
  Vector<double> uj_psijk(G.row_dimension());
  unsigned int iterations;
  CG(G, urhs, uj_psijk, 1e-15, 250, iterations);


  cout << "* compute wavelet-Galerkin approximations for several levels..." << endl;
  const int jmin = basis.j0();
//   const int jmax = jmin+2;
  const int jmax = 16;
  Vector<double> js(jmax-jmin+1);
  Vector<double> Linfty_errors(jmax-jmin+1), L2_errors(jmax-jmin+1),
    discr_L2_errors(jmax-jmin+1), discr_H1_errors(jmax-jmin+1);

  for (int j = jmin; j <= jmax; j++) {
    cout << "  j=" << j << ":" << endl;
    js[j-jmin] = j;
    
    delta.set_level(j);

    // setup rhs in the phi_{j,k} basis,
    Vector<double> rhs_phijk(delta.row_dimension(), false);
    if (solution == 1) {
      if (d == 2) {
 	// exact right-hand side is known
 	rhs_phijk = sqrt(ldexp(1.0, -j));
      } else {
	// perform quadrature with a composite rule on [0,1]
	cout << "  solution 1, quadrature for rhs..." << endl;
	SimpsonRule simpson;
	CompositeRule<1> composite(simpson, 12);
	SchoenbergIntervalBSpline_td<d> sbs(j,0);
	for (int k = basis.DeltaLmin(); k <= basis.DeltaRmax(j); k++) {
	  sbs.set_k(k);
	  rhs_phijk[k-basis.DeltaLmin()]
 	    = composite.integrate(sbs, // against f=1
 				  Point<1>(std::max(0.0, (k+ell1<d>())*ldexp(1.0, -j))),
 				  Point<1>(std::min(1.0, (k+ell2<d>())*ldexp(1.0, -j))));
 	}
      }
    } else {
      if (solution == 2) {
	// perform quadrature with a composite rule on [0,1]
	cout << "  solution 2, quadrature for rhs..." << endl;
	RHS2_part f_part;
	SimpsonRule simpson;
	CompositeRule<1> composite(simpson, 72);
	SchoenbergIntervalBSpline_td<d> sbs(j,0);
	for (int k = basis.DeltaLmin(); k <= basis.DeltaRmax(j); k++) {
	  sbs.set_k(k);
	  ProductFunction<1> integrand(&f_part, &sbs);
	  rhs_phijk[k-basis.DeltaLmin()]
 	    = composite.integrate(integrand,
 				  Point<1>(std::max(0.0, (k+ell1<d>())*ldexp(1.0, -j))),
 				  Point<1>(std::min(1.0, (k+ell2<d>())*ldexp(1.0, -j))))
 	    + 4*sbs.value(Point<1>(0.5));
	}
      } else {
	if (solution == 3) {
	  // perform quadrature with a composite rule on [0,1]
	  cout << "  solution 3, quadrature for rhs..." << endl;
	  RHS3 f;
	  SimpsonRule simpson;
	  CompositeRule<1> composite(simpson, 72);
	  SchoenbergIntervalBSpline_td<d> sbs(j,0);
	  for (int k = basis.DeltaLmin(); k <= basis.DeltaRmax(j); k++) {
	    sbs.set_k(k);
	    ProductFunction<1> integrand(&f, &sbs);
	    rhs_phijk[k-basis.DeltaLmin()]
	      = composite.integrate(integrand,
				    Point<1>(std::max(0.0, (k+ell1<d>())*ldexp(1.0, -j))),
				    Point<1>(std::min(1.0, (k+ell2<d>())*ldexp(1.0, -j))));
	  }
	} else {
	  if (solution == 4 || solution == 5) {
	    // perform quadrature with a composite rule on [0,1]
	    cout << "  solution " << solution << ", quadrature for rhs..." << endl;
	    RHS4_part f_part(kink);
	    SimpsonRule simpson;
	    CompositeRule<1> composite(simpson, 72);
	    SchoenbergIntervalBSpline_td<d> sbs(j,0);
	    for (int k = basis.DeltaLmin(); k <= basis.DeltaRmax(j); k++) {
	      sbs.set_k(k);
	      ProductFunction<1> integrand(&f_part, &sbs);
	      // f is piecewise smooth with (potential) jump at x=a
	      rhs_phijk[k-basis.DeltaLmin()] = (1/kink + 1/(1-kink))*sbs.value(Point<1>(kink));
	      if (std::max(0.0, (k+ell1<d>())*ldexp(1.0, -j)) < kink) // generator intersects left half of the interval
		rhs_phijk[k-basis.DeltaLmin()]
		  += composite.integrate(integrand,
					 Point<1>(std::max(0.0, (k+ell1<d>())*ldexp(1.0, -j))),
					 Point<1>(std::min(kink, (k+ell2<d>())*ldexp(1.0, -j))));
	      
	      if (std::min(1.0, (k+ell2<d>())*ldexp(1.0, -j)) > kink) // generator intersects right half of the interval
		rhs_phijk[k-basis.DeltaLmin()]
		  += composite.integrate(integrand,
					 Point<1>(std::max(kink, (k+ell1<d>())*ldexp(1.0, -j))),
					 Point<1>(std::min(1.0, (k+ell2<d>())*ldexp(1.0, -j))));
	    }
	  }
	}
      }
    }
//     cout << "  rhs in phi_{j,k} basis: " << rhs_phijk << endl;

    // transform rhs into that of psi_{j,k} basis:
    // 1. apply T_{j-1}^T
    Vector<double> rhs(delta.row_dimension(), false);
    if (j == basis.j0())
      rhs = rhs_phijk;
    else
      basis.apply_Tj_transposed(j-1, rhs_phijk, rhs);
    // 2. apply D^{-1}
    for (int k(0); k < basis.Deltasize(basis.j0()); k++)
      rhs[k] /= (1<<basis.j0());
    for (int level = basis.j0(); level < j; level++) {
      for (int k(basis.Deltasize(level)); k < basis.Deltasize(level+1); k++)
	rhs[k] /= (1<<level);
    }
//     cout << "  rhs in psi_{j,k} basis: " << rhs << endl;

    // solve Galerkin system
    Vector<double> ulambda(delta.row_dimension()), residual(delta.row_dimension(), false);
    unsigned int iterations;
    CG(delta, rhs, ulambda, 1e-15, 500, iterations);

//     cout << "  solution coefficients: " << ulambda;
    cout << "  Galerkin system solved with residual (infinity) norm ";
    delta.apply(ulambda, residual);
    residual -= rhs;
    cout << linfty_norm(residual) << endl;

    // compute coefficients of ulambda in the phi_{j,k} basis
    // 1. apply D^{-1}
    for (int k(0); k < basis.Deltasize(basis.j0()); k++)
      ulambda[k] /= (1<<basis.j0());
    for (int level = basis.j0(); level < j; level++) {
      for (int k(basis.Deltasize(level)); k < basis.Deltasize(level+1); k++)
	ulambda[k] /= (1<<level);
    }

    // save the L_2 wavelet coefficients
    Vector<double> ulambda_prolong(basis.Deltasize(jref));
    std::copy(ulambda.begin(), ulambda.end(), ulambda_prolong.begin());

    // 2. apply T_{j-1}
    Vector<double> ulambda_phijk(delta.row_dimension(), false);
    if (j == basis.j0())
      ulambda_phijk = ulambda;
    else
      basis.apply_Tj(j-1, ulambda, ulambda_phijk);
//     cout << "  solution coefficients in phi_{j,k} basis: " << ulambda_phijk << endl;

    // evaluate linear combination of Schoenberg B-splines on a grid
    const unsigned int N = 100;
    const double h = 1./N;
    Vector<double> ulambda_values(N+1);
    for (unsigned int i = 0; i <= N; i++) {
      const double x = i*h;
      SchoenbergIntervalBSpline_td<d> sbs(j,0);
      for (unsigned int k = 0; k < delta.row_dimension(); k++) {
	sbs.set_k(basis.DeltaLmin()+k);
	ulambda_values[i] +=
	  ulambda_phijk[k] * sbs.value(Point<1>(x));
      }
    }
//     cout << "  point values of Galerkin solution: " << ulambda_values << endl;

    // evaluate exact solution
    Vector<double> uexact_values(N+1);
    for (unsigned int i = 0; i <= N; i++) {
      const double x = i*h;
      uexact_values[i] = uexact->value(Point<1>(x));
    }
//     cout << "  point values of exact solution: " << uexact_values << endl;

    // compute some errors
    const double Linfty_error = linfty_norm(ulambda_values-uexact_values);
    cout << "  L_infinity error on a subgrid: " << Linfty_error << endl;
    Linfty_errors[j-jmin] = Linfty_error;

    const double L2_error = sqrt(l2_norm_sqr(ulambda_values-uexact_values)*h);
    cout << "  L_2 error on a subgrid: " << L2_error << endl;
    L2_errors[j-jmin] = L2_error;

    Vector<double> coeff_error = ulambda_prolong-uj_psijk;
    const double discr_L2_error = l2_norm(coeff_error);
    cout << "  discr. L_2 error: " << discr_L2_error << endl;
    discr_L2_errors[j-jmin] = discr_L2_error;

    for (int k(0); k < basis.Deltasize(basis.j0()); k++)
      coeff_error[k] *= (1<<basis.j0());
    for (int level = basis.j0(); level < jref; level++) {
      for (int k(basis.Deltasize(level)); k < basis.Deltasize(level+1); k++)
 	coeff_error[k] *= (1<<level);
    }
    const double discr_H1_error = l2_norm(coeff_error);
    cout << "  discr. H^1 error: " << discr_H1_error << endl;
    discr_H1_errors[j-jmin] = discr_H1_error;
  }

#if 0
  // write Galerkin errors to a file
  ostringstream filename;
  filename << "full_galerkin_errors_" << d << "_" << solution << ".m";
  ofstream galerkin_stream(filename.str().c_str());
  galerkin_stream << "js=" << js << ";" << endl
 		  << "Linfty_errors=" << Linfty_errors << ";" << endl
 		  << "L2_errors=" << L2_errors << ";" << endl
		  << "discr_L2_errors=" << discr_L2_errors << ";" << endl
		  << "discr_H1_errors=" << discr_H1_errors << ";" << endl;
  galerkin_stream.close();
#endif

  delete uexact;

  return 0;
}
