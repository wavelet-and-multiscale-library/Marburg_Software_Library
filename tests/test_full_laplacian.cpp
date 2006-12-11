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
#include <interval/p_expansion.h>
#include <interval/spline_expansion.h>

#include "poisson_1d_solutions.h"

using namespace std;
using namespace MathTL;
using namespace WaveletTL;

int main()
{
  cout << "Testing FullLaplacian ..." << endl;

  const unsigned int d = 2;
  const unsigned int dT = 2;

  SplineBasis<d,dT> basis("P","",1,1,0,0); // PBasis, complementary b.c.'s
//   FullLaplacian<d,dT> delta(basis, dyadic);
  FullLaplacian<d,dT> delta(basis, energy);

  cout << "* stiffness matrix on coarsest level j0=" << basis.j0() << ":" << endl
       << delta;

  Vector<double> diagonal(delta.row_dimension());
  for (unsigned int k = 0; k < diagonal.size(); k++)
    diagonal[k] = delta.diagonal(k);
  cout << "* main diagonal of unpreconditioned stiffness matrix:" << endl
       << "  " << diagonal << endl;
  
  delta.set_level(basis.j0()+1);
  cout << "* stiffness matrix on next level j0+1=" << basis.j0()+1 << ":" << endl
       << delta;

  diagonal.resize(delta.row_dimension());
  for (unsigned int k = 0; k < diagonal.size(); k++)
    diagonal[k] = delta.diagonal(k);
  cout << "* main diagonal of unpreconditioned stiffness matrix:" << endl
       << "  " << diagonal << endl;
  
#if 1
  const unsigned int solution = 1;
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
  const int jref = 15;
  Vector<double> uj_psijk;
  expand(uexact, basis, false, jref-1, uj_psijk);

  cout << "* compute wavelet-Galerkin approximations for several levels..." << endl;
  const int jmin = basis.j0();
//   const int jmax = jmin+5;
  const int jmax = 12;
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
    for (unsigned int k(0); k < rhs.size(); k++)
      rhs[k] /= delta.D(k);
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
    for (unsigned int k(0); k < ulambda.size(); k++)
      ulambda[k] /= delta.D(k);

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

    delta.set_level(jref);
    for (unsigned int k(0); k < coeff_error.size(); k++)
      coeff_error[k] *= delta.D(k);
    const double discr_H1_error = l2_norm(coeff_error);
    cout << "  discr. H^1 error: " << discr_H1_error << endl;
    discr_H1_errors[j-jmin] = discr_H1_error;
  }

  delete uexact;

#endif

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

  return 0;
}
