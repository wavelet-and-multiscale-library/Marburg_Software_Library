// implementation for spline_expansion.h

#include <numerics/quadrature.h>
#include <galerkin/full_gramian.h>

namespace WaveletTL
{
  template <int d, int dT>
  void
  expand(const Function<1>* f,
	 const SplineBasis<d,dT,P_construction>& basis,
	 const bool primal,
	 const int jmax,
	 Vector<double>& coeffs)
  {
    assert(jmax >= basis.j0());
    coeffs.resize(basis.Deltasize(jmax+1));

    // 1. compute integrals w.r.t. the primal generators on level jmax
    Vector<double> coeffs_phijk(coeffs.size(), false);
    SimpsonRule simpson;
    CompositeRule<1> composite(simpson, 12); // should be sufficient for many cases
    SchoenbergIntervalBSpline_td<d> sbs(jmax+1,0);
    for (int k = basis.DeltaLmin(); k <= basis.DeltaRmax(jmax+1); k++) {
      sbs.set_k(k);
      ProductFunction<1> integrand(f, &sbs);
      coeffs_phijk[k-basis.DeltaLmin()]
	= composite.integrate(integrand,
			      Point<1>(std::max(0.0, (k+ell1<d>())*ldexp(1.0, -jmax-1))),
			      Point<1>(std::min(1.0, (k+ell2<d>())*ldexp(1.0, -jmax-1))));
    }
    // 2. transform rhs into that of psi_{j,k} basis: apply T_{j-1}^T
    Vector<double> rhs(coeffs.size(), false);
//     if (jmax == basis.j0())
//       rhs = coeffs_phijk;
//     else
    basis.apply_Tj_transposed(jmax, coeffs_phijk, rhs);
    
    if (!primal) {
      FullGramian<d,dT> G(basis);
      G.set_level(jmax+1);
      unsigned int iterations;
      CG(G, rhs, coeffs, 1e-15, 250, iterations);
    } else {
      coeffs.swap(rhs);
    }
  }

  template <int d, int dT>
  void
  expand(const Function<1>* f,
	 const SplineBasis<d,dT,P_construction>& basis,
	 const bool primal,
	 const int jmax,
	 InfiniteVector<double, typename SplineBasis<d,dT,P_construction>::Index>& coeffs)
  {
    Vector<double> coeffs_vector;
    expand(f, basis, primal, jmax, coeffs_vector);
    typedef typename Vector<double>::size_type size_type;
    typedef typename SplineBasis<d,dT,P_construction>::Index Index;
    size_type i(0);
    for (Index lambda(basis.first_generator(basis.j0()));
	 i < coeffs_vector.size(); ++lambda, i++)
      {
	const double coeff = coeffs_vector[i];
	if (fabs(coeff)>1e-15)
	  coeffs.set_coefficient(lambda, coeff);
      } 
  }

}
