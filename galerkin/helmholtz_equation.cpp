// implementation for helmholtz_equation.h

#include <numerics/quadrature.h>
#include <numerics/schoenberg_splines.h>
#include <interval/interval_bspline.h>
#include <numerics/eigenvalues.h>

using namespace MathTL;

namespace WaveletTL
{
  template <int d, int dT>
  HelmholtzEquation1D<d,dT>::HelmholtzEquation1D
  (const WaveletBasis& basis,
   const double alpha,
   const InfiniteVector<double,Index>& y)
    : basis_(basis), alpha_(alpha), y_(y),
      H_(basis_, alpha, no_precond),
      G_(basis, InfiniteVector<double,Index>()),
      GC_(&G_),
      A_(basis, InfiniteVector<double,Index>()),
      AC_(&A_),
      normA(0.0), normAinv(0.0)
  {
    y_precond = y_;
    y_precond.scale(this, -1);
  }

  template <int d, int dT>
  void
  HelmholtzEquation1D<d,dT>::set_alpha(const double alpha) const
  {
    assert(alpha >= 0);
    alpha_ = alpha;
    H_.set_alpha(alpha);
    y_precond = y_;
    y_precond.scale(this, -1);
  }
  
  template <int d, int dT>
  void
  HelmholtzEquation1D<d,dT>::set_rhs(const InfiniteVector<double,Index>& y) const
  {
    y_ = y;
    y_precond = y_;
    y_precond.scale(this, -1);
  }

  template <int d, int dT>
  inline
  double
  HelmholtzEquation1D<d,dT>::D(const typename WaveletBasis::Index& lambda) const
  {
#if 1
    // determine number of index lambda
    size_type number = 0;
    if (lambda.e() == 0) {
      number = lambda.k()-basis_.DeltaLmin();
    } else {
      number = basis_.Deltasize(lambda.j())+lambda.k()-basis_.Nablamin();
    }
    
    H_.set_level(lambda.j()+lambda.e());
    return sqrt(H_.diagonal(number));
#else
    return sqrt(a(lambda, lambda));
#endif
  }
  
  template <int d, int dT>
  inline
  double
  HelmholtzEquation1D<d,dT>::a(const typename WaveletBasis::Index& lambda,
			       const typename WaveletBasis::Index& nu) const
  {
    return alpha_ * GC_.a(lambda, nu) + A_.a(lambda, nu);
  }
  
  template <int d, int dT>
  double
  HelmholtzEquation1D<d,dT>::norm_A() const
  {
    if (normA == 0.0) {
      FullHelmholtz<d,dT> A(basis_, alpha_, energy);
      A.set_level(basis().j0()+4);
      double help;
      unsigned int iterations;
      LanczosIteration(A, 1e-6, help, normA, 200, iterations);
      normAinv = 1./help;
    }

    return normA;
  }
   
  template <int d, int dT>
  double
  HelmholtzEquation1D<d,dT>::norm_Ainv() const
  {
    if (normAinv == 0.0) {
      FullHelmholtz<d,dT> A(basis_, alpha_, energy);
      A.set_level(basis().j0()+4);
      double help;
      unsigned int iterations;
      LanczosIteration(A, 1e-6, help, normA, 200, iterations);
      normAinv = 1./help;
    }

    return normAinv;
  }

  template <int d, int dT>
  void
  HelmholtzEquation1D<d,dT>::add_level (const Index& lambda,
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
   
    //GC_.add_level(lambda, help1, j, factor * alpha_/D(lambda), J, strategy);
    GC_.add_level(lambda, help1_full, j, factor * alpha_/D(lambda), J, strategy);
   
    // elliptic part, help2 = column block of D(D^{-1}<-APsi,Psi>^T D^{-1})DD_alpha^{-1}
    //AC_.add_level(lambda, help2, j, factor*AC_.D(lambda)/D(lambda), J, strategy);
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

    Vector<double> help1_full_new(w.size());
    for (typename InfiniteVector<double,Index>::const_iterator it(help1.begin());
 	   it != help1.end(); ++it) {
      help1_full_new[it.index().number()] = *it;
    }

    //    w.add(help1);
    w += help1_full_new;
  }
  
  
}
