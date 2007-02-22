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
      y_(y),
      G_(basis, InfiniteVector<double,Index>()),
//       GC_(&G_, G_file, jmax),
      GC_(&G_, G_file, jmax, 1.0, 1.0), // dirty
      A_(basis, InfiniteVector<double,Index>()),
//       AC_(&A_, A_file, jmax),
      AC_(&A_, A_file, jmax, 1.0, 1.0), // dirty
//       normA(0.0), normAinv(0.0)
      normA(1.0), normAinv(1.0) // dirty
  {
    y_precond = y_;
    y_precond.scale(this, -1);
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

  template <int d, int dT>
  void
  LDomainHelmholtzEquation<d,dT>::set_alpha(const double alpha) const
  {
    assert(alpha >= 0);
    alpha_ = alpha;
    y_precond = y_;
    y_precond.scale(this, -1);
  }
  
  template <int d, int dT>
  void
  LDomainHelmholtzEquation<d,dT>::set_rhs(const InfiniteVector<double,Index>& y) const
  {
    y_ = y;
    y_precond = y_;
    y_precond.scale(this, -1);
  }


  template <int d, int dT>
  double
  LDomainHelmholtzEquation<d,dT>::norm_A() const
  {
    if (normA == 0.0) {
      typedef typename WaveletBasis::Index Index;
      std::set<Index> Lambda;
      const int j0 = basis().j0();
      const int jmax = j0+1; // the bigger the more precise...
      for (Index lambda = basis().first_generator(j0);; ++lambda) {
	Lambda.insert(lambda);
	if (lambda == basis().last_wavelet(jmax)) break;
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
   
  template <int d, int dT>
  double
  LDomainHelmholtzEquation<d,dT>::norm_Ainv() const
  {
    if (normAinv == 0.0) {
      typedef typename WaveletBasis::Index Index;
      std::set<Index> Lambda;
      const int j0 = basis().j0();
      const int jmax = j0+1;  // the bigger the more precise...
      for (Index lambda = basis().first_generator(j0);; ++lambda) {
	Lambda.insert(lambda);
	if (lambda == basis().last_wavelet(jmax)) break;
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
  
  template <int d, int dT>
  double
  LDomainHelmholtzEquation<d,dT>::s_star() const
  {
    // notation from [St04a]
    const double t = operator_order();
    const int n = 2;
    const int mt = WaveletBasis::primal_vanishing_moments();
    const double gamma = WaveletBasis::primal_regularity();
    
    return std::min((t+mt)/(double)n, (gamma-t)/(n-1.)); // [St04a, Th. 2.3]
  }

  template <int d, int dT>
  void
  LDomainHelmholtzEquation<d,dT>::add_level
  (const Index& lambda,
   //InfiniteVector<double, Index>& w,
   Vector<double>& w,
   const int j,
   const double factor,
   const int J,
   const CompressionStrategy strategy) const
  {
    // We have to compute a (level-)part of the lambda-th column of 
    //   D_alpha^{-1}<(alpha*I-A)Psi,Psi>^T D_alpha^{-1}
    //   = alpha*D_alpha^{-1}<Psi,Psi>^T D_alpha^{-1}
    //     - D_alpha^{-1}D ( D^{-1}<APsi,Psi>^T D^{-1} ) DD_alpha^{-1}

    // Gramian part, help1 = column block of alpha*<Psi,Psi>^T D_alpha^{-1}
    
    InfiniteVector<double,Index> help1, help2;
    Vector<double> help1_full(w.size()), help2_full(w.size());

    GC_.add_level(lambda, help1_full, j, factor * alpha_/D(lambda), J, strategy);
    
    // elliptic part, help2 = column block of D(D^{-1}<-APsi,Psi>^T D^{-1})DD_alpha^{-1}
    AC_.add_level(lambda, help2_full, j, factor*AC_.D(lambda)/D(lambda), J, strategy);
    
    // hack: copy full vectors to sparse ones
    for (unsigned int i = 0; i < help1_full.size(); i++) {
      if (help1_full[i] != 0.) {
	Index ind(basis_.get_wavelet(i));
	help1.set_coefficient(ind, help1_full[i]);
      }
      if (help2_full[i] != 0.) {
	Index ind(basis_.get_wavelet(i));
	help2.set_coefficient(ind, help2_full[i]);
      }
    }

    help2.scale(&AC_, 1); // help2 *= D
    help1.add(help2);
    help1.scale(this, -1); // help1 *= D_alpha^{-1}

    // hack: copy help1 to full vector
    
    Vector<double> help1_full_new(w.size());
    for (typename InfiniteVector<double,Index>::const_iterator it(help1.begin());
 	   it != help1.end(); ++it) {
      help1_full_new[it.index().number()] = *it;
    }
    
    //w.add(help1);
    
    w += help1_full_new;
  }
  
  

}
