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
    precompute_rhs();
  }


  template <class WBASIS>
  void
  BiharmonicEquation1D<WBASIS>::precompute_rhs()
  {
    typedef typename WaveletBasis::Index Index;
    cout << "precompute rhs.." << endl;
    // precompute the right-hand side on a fine level
    InfiniteVector<double,Index> fhelp;

    for (int i=0; i<basis_.degrees_of_freedom();i++) {
//        cout << "bin hier: " << i << endl;
//        cout << D(*(basis_.get_wavelet(i))) << endl;
//        cout << *(basis_.get_wavelet(i)) << endl;
        const double coeff = f(*(basis_.get_wavelet(i)))/D(*(basis_.get_wavelet(i)));
//        cout << f(*(basis_.get_wavelet(i))) << endl;
//        cout << coeff << endl;
        fhelp.set_coefficient(*(basis_.get_wavelet(i)), coeff);
//        cout << *(basis_.get_wavelet(i)) << endl;
    }

    fnorm_sqr = l2_norm_sqr(fhelp);

    // sort the coefficients into fcoeffs
    fcoeffs.resize(fhelp.size());
    unsigned int id(0), id2(0);
    for (typename InfiniteVector<double,Index>::const_iterator it(fhelp.begin()), itend(fhelp.end());
         it != itend; ++it, ++id)
      fcoeffs[id] = std::pair<Index,double>(it.index(), *it);
    sort(fcoeffs.begin(), fcoeffs.end(), typename InfiniteVector<double,Index>::decreasing_order());

    cout << "end precompute rhs.." << endl;
//    cout << fhelp << endl;
//    cout << fcoeffs << endl;
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


  template <class WBASIS>
  double
  BiharmonicEquation1D<WBASIS>::norm_A() const
  {
    if (normA == 0.0)
    {
      typedef typename WaveletBasis::Index Index;
      std::set<Index> Lambda;
      const int j0 = basis().j0();
      const int jmax = j0+3;

      for (Index lambda = first_generator(&basis(), j0);; ++lambda) {
        Lambda.insert(lambda);
        if (lambda == last_wavelet(&basis(), jmax)) break;
      }

      SparseMatrix<double> A_Lambda;
      setup_stiffness_matrix(*this, Lambda, A_Lambda);

#if 1
      double help;
      unsigned int iterations;
      LanczosIteration(A_Lambda, 1e-6, help, normA, 200, iterations);
      normAinv = 1./help;
#else
      Vector<double> xk(Lambda.size(), false);
      xk = 1;
      unsigned int iterations;
      normA = PowerIteration(A_Lambda, xk, 1e-6, 100, iterations);
#endif
    }

    return normA;
  }


  template <class WBASIS>
  double
  BiharmonicEquation1D<WBASIS>::norm_Ainv() const
  {
    if (normAinv == 0.0)
    {
      typedef typename WaveletBasis::Index Index;
      std::set<Index> Lambda;
      const int j0 = basis().j0();
      const int jmax = j0+3;

      for (Index lambda = first_generator(&basis(), j0);; ++lambda) {
        Lambda.insert(lambda);
        if (lambda == last_wavelet(&basis(), jmax)) break;
      }

      SparseMatrix<double> A_Lambda;
      setup_stiffness_matrix(*this, Lambda, A_Lambda);

#if 1
      double help;
      unsigned int iterations;
      LanczosIteration(A_Lambda, 1e-6, help, normA, 200, iterations);
      normAinv = 1./help;
#else
      Vector<double> xk(Lambda.size(), false);
      xk = 1;
      unsigned int iterations;
      normAinv = InversePowerIteration(A_Lambda, xk, 1e-6, 200, iterations);
#endif
    }

    return normAinv;
  }


  template <class WBASIS>
  inline
  void
  BiharmonicEquation1D<WBASIS>::RHS(const double eta,
                                    InfiniteVector<double, typename WBASIS::Index>& coeffs) const
  {
    coeffs.clear();
    double coarsenorm(0);
    double bound(fnorm_sqr - eta*eta);
    typedef typename WBASIS::Index Index;
    typename Array1D<std::pair<Index, double> >::const_iterator it(fcoeffs.begin());
    do {
      coarsenorm += it->second * it->second;
      coeffs.set_coefficient(it->first, it->second);
      ++it;
    } while (it != fcoeffs.end() && coarsenorm < bound);
  }

}
