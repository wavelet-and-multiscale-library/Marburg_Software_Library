// implementation for r_q_mask.h

#include <algebra/vector.h>
#include <algebra/matrix.h>
#include <numerics/matrix_decomp.h>

using MathTL::Matrix;
using MathTL::Vector;
using MathTL::QUDecomposition;

namespace WaveletTL
{
  template <unsigned int L, int BEGIN>
  inline
  InfiniteVector<double, int>
  RQRefinementMask<L,BEGIN>::evaluate(const int resolution) const
  {
    return evaluate(0, resolution);
  }
  
  template <unsigned int L, int BEGIN>
  InfiniteVector<double, int>
  RQRefinementMask<L,BEGIN>::evaluate(const int mu,  const int resolution) const
  {
    assert(resolution >= 0);
    
    InfiniteVector<double, int> r;
    
    // First we calculate the values on Z

    // abbreviate begin and end of the support of phi
    int suppleft  = begin();
    int suppright = end();
    
    // exclude the special case of \chi_{[0,1)}
    if (suppleft == suppright-1) {
      r.set_coefficient(0, 1);
    } else {
      // In the following, we set up the eigenvalue problem for the values
      //   V_\alpha := D^\mu\phi(\alpha), \alpha\in\mathbb Z
      // The eigenvector is determined uniquely by the following equations,
      // see [DM] for details (so we don't need an iterative scheme):
      //
      // (3.22) eigenvalue condition
      //   2^{-\mu}V_\alpha = \sum_\beta a_{2\alpha-\beta}V_\beta, \alpha\in\mathbb Z
      //
      // (3.23) orthogonality condition
      //   \sum_\alpha (-\alpha)^\nu V_\alpha = \mu!\delta_{\mu,\nu}, 0\le\nu\le\mu
	
      Matrix<double> A(suppright-suppleft+mu, suppright-suppleft-1);
      Vector<double> b(suppright-suppleft+mu);

      for (int m = 0; m < suppright-suppleft-1; m++) {
	// (3.22)
	for (int n = 0; n < suppright-suppleft-1; n++) {
	  A(m, n) = a(2*(suppleft+1+m)-(suppleft+1+n));
	}
	A(m, m) -= ldexp(1.0, -mu);
	// b[m] = 0;
      }
      
      for (int nu = 0; nu <= mu; nu++) {
	// (3.23)
	for (int n = 0; n < suppright-suppleft-1; n++)
	  A(suppright-suppleft-1+nu, n) = minus1power(nu) * intpower(suppleft+1+n, nu);
	b[suppright-suppleft-1+nu] = (nu == mu ? factorial(mu) : 0);
      }

      // the system matrix is rectangular, but Ax=b is solvable via a QU decomposition
      QUDecomposition<double> qu(A);
      assert(qu.hasFullRank());
      Vector<double> x;
      qu.solve(b, x);
      x.compress(1e-15);

      // reinterpret the entries of x
      for (int n = 0; n < suppright-suppleft-1; n++)
 	r.set_coefficient(suppleft+1+n, x[n]);
    }

    // For the remaining non-integer points we use the refinement relation of phi
    if (resolution > 0) {
      for (int newres(1); newres <= resolution; newres++) {
	// copy the coarse values \phi(2^{-j}m) = \phi(2^{-(j+1)}2m), newres=j+1
	InfiniteVector<double, int> coarse;
	coarse.swap(r);
	for (typename InfiniteVector<double, int>::const_iterator it(coarse.begin());
	     it != coarse.end(); ++it) {
	  r.set_coefficient(2*it.index(), *it);
	}

	// For m=2i+1 odd, we have to compute and insert the new values
	//   \phi(2^{-(j+1)}m) = \sum_k a_k\phi(2^{-j}(m-2^jk)), newres=j+1
	// The requirement suppleft < 2^{-(j+1)m} < suppright is equivalent to
	//   2^j*suppleft <= i <= 2^j*suppright-1
	for (int i = (1<<(newres-1))*suppleft; i <= (1<<(newres-1))*suppright-1; i++) {
	  int m = 2*i+1;
	  for (int k = begin(); k <= end(); k++)
	    r[m] += a(k) * coarse.get_coefficient(m - (1<<(newres-1))*k);
	}
      }   
    }
    
    return r;
  }

  template <unsigned int L, int BEGIN>
  SampledMapping<1>
  RQRefinementMask<L,BEGIN>::evaluate(const int mu,
				      const int p,
				     const int j,
				     const int k,
				     const int a,
				     const int b,
				     const int resolution) const
  {
    // take "relative resolution" into account
    InfiniteVector<double, int> help(evaluate(mu, resolution-j)), v;
    int i = 1;
    //left end of support of \phi_{p,j,k}
    int lsupp((-(int) L+1)/2);
    const double factor(twotothejhalf(j) * ldexp(1.0, j*mu));
    for (typename InfiniteVector<double, int>::const_iterator it(help.begin());
	 it != help.end(); ++it) {
      double factor2 =  (lsupp + i * pow(2, -(resolution-j)))/(L/2);
      i++;
      v.set_coefficient(it.index()+(1<<(resolution-j))*k, *it * factor * pow(factor2,p));//NOTLÖSUNG !!! Funktioniert nicht für Ableitungen
    }
   
    // clipping is done in the SampledMapping constructor
    return SampledMapping<1>(a, b, v, resolution);
  }

  template <unsigned int L, int BEGIN>
  const double
  RQRefinementMask<L,BEGIN>::moment(const unsigned int k) const
  {
    double r(1.0);
	
    if (k > 0)
      {
	r = 0.0;

	for (unsigned int m = 0; m < k; m++) {
	  double alsum = 0;
	  for (int l = begin(); l <= end(); l++) {
	    alsum += a(l) * intpower(l, k-m);
	  }
	  r += binomial(k, m) * moment(m) * alsum;
	}

	r /= (1<<(k+1))-2;
      }
    
    return r;
  }
  
}
