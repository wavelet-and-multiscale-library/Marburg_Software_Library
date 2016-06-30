// implementation for biharmonic_equation.h

#include <cmath>
#include <list>
#include <utils/array1d.h>
#include <algebra/vector.h>
#include <algebra/sparse_matrix.h>
#include <numerics/eigenvalues.h>
#include <numerics/gauss_data.h>
#include <numerics/quadrature.h>
#include <geometry/point.h>

using namespace MathTL;

namespace WaveletTL
{
  template <class IBasis>
  BiharmonicEquation1D<IBasis>::BiharmonicEquation1D
  (const WaveletBasis& basis, const Function<1>* g)
    : basis_(basis), g_(g)
  {
  }
  
  template <class IBasis>
  inline
  double
  BiharmonicEquation1D<IBasis>::a(const typename WaveletBasis::Index& lambda,
                             const typename WaveletBasis::Index& nu) const
  {
    return a(lambda, nu, WaveletBasis::primal_polynomial_degree()*WaveletBasis::primal_polynomial_degree());
  }
  
  template <class IBasis>
  double
  BiharmonicEquation1D<IBasis>::a(const typename WaveletBasis::Index& lambda,
                             const typename WaveletBasis::Index& nu,
                             const unsigned int p) const
  {
    double r = 0;

    // first compute the support intersection of \psi_\lambda and \psi_\nu:
    typedef typename WaveletBasis::Support Support;
    Support supp;
    if (intersect_supports(basis_, lambda, nu, supp))
      {
        // Set up Gauss points and weights for a composite quadrature formula:
        // (TODO: maybe use an instance of MathTL::QuadratureRule instead of computing
        // the Gauss points and weights)
        const unsigned int N_Gauss = (p+1)/2;
        const double h = ldexp(1.0, -supp.j);
        Array1D<double> gauss_points (N_Gauss*(supp.k2-supp.k1)), der1values, der2values;
        for (int patch = supp.k1, id = 0; patch < supp.k2; patch++) // refers to 2^{-j}[patch,patch+1]
          for (unsigned int n = 0; n < N_Gauss; n++, id++)
            gauss_points[id] = h*(2*patch+1+GaussPoints[N_Gauss-1][n])/2.;
        
        // - compute point values of the integrands
        evaluate(basis_, 2, lambda, gauss_points, der1values);
        evaluate(basis_, 2, nu, gauss_points, der2values);

        // - add all integral shares
        for (int patch = supp.k1, id = 0; patch < supp.k2; patch++)
          for (unsigned int n = 0; n < N_Gauss; n++, id++) {
            const double gauss_weight = GaussWeights[N_Gauss-1][n] * h;
            r += der1values[id] * der2values[id] * gauss_weight;
          }
      }
    
    return r;
  }

  template <class IBasis>
  double
  BiharmonicEquation1D<IBasis>::f(const typename WaveletBasis::Index& lambda) const
  {
    // f(v) = \int_0^1 g(t)v(t) dt

    double r = 0;

    const int j = lambda.j()+lambda.e();
    int k1, k2;
    support(basis_, lambda, k1, k2);

    // Set up Gauss points and weights for a composite quadrature formula:
    const unsigned int N_Gauss = 7;
    const double h = ldexp(1.0, -j);
    Array1D<double> gauss_points (N_Gauss*(k2-k1)), vvalues;
    for (int patch = k1; patch < k2; patch++) // refers to 2^{-j}[patch,patch+1]
      for (unsigned int n = 0; n < N_Gauss; n++)
        gauss_points[(patch-k1)*N_Gauss+n] = h*(2*patch+1+GaussPoints[N_Gauss-1][n])/2;

    // - compute point values of the integrand
    evaluate(basis_, 0, lambda, gauss_points, vvalues);
    
    // - add all integral shares
    for (int patch = k1, id = 0; patch < k2; patch++)
      for (unsigned int n = 0; n < N_Gauss; n++, id++) {
        const double t = gauss_points[id];
        const double gauss_weight = GaussWeights[N_Gauss-1][n] * h;
            
        const double gt = g_->value(Point<1,double>(t));
        if (gt != 0)
          r += gt
            * vvalues[id]
            * gauss_weight;
      }
    
    return r;
  }
  
  template <class IBasis>
  inline
  double
  BiharmonicEquation1D<IBasis>::D(const typename WaveletBasis::Index& lambda) const
  {
    return sqrt(a(lambda, lambda));
  }
}
