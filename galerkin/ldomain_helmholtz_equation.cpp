// implementation for ldomain_helmholtz_equation.h

#include <cmath>
#include <time.h>
#include <utils/fixed_array1d.h>
#include <numerics/gauss_data.h>

namespace WaveletTL
{
  template <int d, int dT>
  LDomainHelmholtzEquation<d,dT>::LDomainHelmholtzEquation
  (const WaveletBasis& basis,
   const char* G_file,
   const char* A_file,
   const int jmax,
   const double alpha,
   const InfiniteVector<double,Index>& y)
    : basis_(basis),
      alpha_(alpha),
      G_(basis, InfiniteVector<double,Index>()),
      GC_(&G_, G_file, jmax),
      A_(basis, InfiniteVector<double,Index>()),
      AC_(&A_, A_file, jmax),
      normA(0.0), normAinv(0.0)
  {
  }

  template <int d, int dT>
  inline
  double
  LDomainHelmholtzEquation<d,dT>::a
  (const Index& lambda,
   const Index& nu) const
  {
    return alpha_ * GC_.a(lambda, nu) + A_.a(lambda, nu);
  }
  
  template <int d, int dT>
  inline
  double
  LDomainHelmholtzEquation<d,dT>::D(const Index& lambda) const
  {
    return sqrt(a(lambda, lambda));
  }

//   template <class IBASIS>
//   double
//   LDomainHelmholtzEquation<IBASIS>::norm_A() const
//   {
//     if (normA == 0.0) {
//       typedef typename WaveletBasis::Index Index;
//       std::set<Index> Lambda;
//       const int j0 = basis().j0();
//       const int jmax = j0+1; // the bigger the more precise...
//       for (Index lambda = basis().first_generator(j0);; ++lambda) {
// 	Lambda.insert(lambda);
// 	if (lambda == basis().last_wavelet(jmax)) break;
//       }
//       SparseMatrix<double> A_Lambda;
//       setup_stiffness_matrix(*this, Lambda, A_Lambda);
      
// #if 1
//       double help;
//       unsigned int iterations;
//       LanczosIteration(A_Lambda, 1e-6, help, normA, 200, iterations);
//       normAinv = 1./help;
// #else
//       Vector<double> xk(Lambda.size(), false);
//       xk = 1;
//       unsigned int iterations;
//       normA = PowerIteration(A_Lambda, xk, 1e-6, 100, iterations);
// #endif
//     }

//     return normA;
//   }
   
//   template <class IBASIS>
//   double
//   LDomainHelmholtzEquation<IBASIS>::norm_Ainv() const
//   {
//     if (normAinv == 0.0) {
//       typedef typename WaveletBasis::Index Index;
//       std::set<Index> Lambda;
//       const int j0 = basis().j0();
//       const int jmax = j0+1;  // the bigger the more precise...
//       for (Index lambda = basis().first_generator(j0);; ++lambda) {
// 	Lambda.insert(lambda);
// 	if (lambda == basis().last_wavelet(jmax)) break;
//       }
//       SparseMatrix<double> A_Lambda;
//       setup_stiffness_matrix(*this, Lambda, A_Lambda);
      
// #if 1
//       double help;
//       unsigned int iterations;
//       LanczosIteration(A_Lambda, 1e-6, help, normA, 200, iterations);
//       normAinv = 1./help;
// #else
//       Vector<double> xk(Lambda.size(), false);
//       xk = 1;
//       unsigned int iterations;
//       normAinv = InversePowerIteration(A_Lambda, xk, 1e-6, 200, iterations);
// #endif
//     }

//     return normAinv;
//   }
  


//   template <class IBASIS>
//   double
//   LDomainHelmholtzEquation<IBASIS>::s_star() const
//   {
//     // notation from [St04a]
//     const double t = operator_order();
//     const int n = 2;
//     const int dT = WaveletBasis::primal_vanishing_moments();
//     const double gamma = WaveletBasis::primal_regularity();
    
//     return std::min((t+dT)/(double)n, (gamma-t)/(n-1.)); // [St04a, Th. 2.3]
//   }


}
