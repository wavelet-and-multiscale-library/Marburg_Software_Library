// implementation for refinable.h

#include <algorithm>
#include <cmath>
#include <iostream>
#include <utils/array1d.h>
#include <utils/tiny_tools.h>
#include <algebra/matrix.h>
#include <algebra/triangular_matrix.h>
#include <numerics/matrix_decomp.h>
#include <numerics/eigenvalues.h>
#include <numerics/differences.h>

namespace WaveletTL
{
  template <class MASK>
  template <unsigned int DERIVATIVE>
  SampledMapping<1> RefinableFunction<MASK>::evaluate(const int j, const int k,
						      const int a, const int b,
						      const int resolution) const
  {
    Grid<1> grid(a, b, (b-a)*(1<<resolution));
    Array1D<double> values((b-a)*(1<<resolution)+1);
    for (unsigned int i(0); i < values.size(); i++) values[i] = 0;

    if (DERIVATIVE == 0u)
      {
	// First we sample the generator on a fine dyadic grid,
	// the wavelet values can then be derived very easily
	
	const int kbegin(LaurentPolynomial<double>::begin().index()),
	  kend(LaurentPolynomial<double>::rbegin().index());
	
	if (kend - kbegin == 1)
	  {
	    // The corresponding generator is the Haar generator,
	    // which we know explicitly.
	    // (the following for loops can be sped up a bit, since we
	    // know the grid points exactly)
	    for (unsigned int i(0); i < grid.size(); i++)
	      {
		if (k*ldexp(1.0,-j) <= grid.points()[i]
		    && grid.points()[i] < (k+1)*ldexp(1.0,-j))
		  values[i] = sqrt(ldexp(1.0,j));
	      }
	  }
	else
	  {
	    // First we compute the values of the original refinable function
	    // (without the j and k parameter). In the end, we will need the
	    // following dyadic resolution:
	    int reshelp = std::max(0, resolution-j);
	    Array1D<double> svalues((kend-kbegin)*(1<<reshelp)+1);
	    for (unsigned int i(0); i < svalues.size(); i++) svalues[i] = 0;
	    
	    // Compute the exact values of phi on the integer numbers.
	    // We assume that phi is continuous, so that
	    //   phi(kbegin) = phi(kend) = 0
	    // For the remaining point values we set up an eigenvalue problem.
	    Matrix<double> A(kend-kbegin-1, kend-kbegin-1);
	    for (int m(kbegin+1); m < kend; m++)
	      for (int n(kbegin+1); n < kend; n++)
		A(m-kbegin-1, n-kbegin-1) = LaurentPolynomial<double>::get_coefficient(2*m-n);
	    
	    // We assume that the largest eigenvalue of A is 1,
	    // so a power iteration should converge to the values of
	    // phi at the integer points.
	    Vector<double> eigenvalues(A.row_dimension(), true);
	    eigenvalues = 1;
	    unsigned int iterations;
	    PowerIteration(A, eigenvalues, 1e-6, 100, iterations); // neglects return value
	    
	    // normalize eigenvector to sum 1 (Strang-Fix condition 0)
	    double sum(0);
	    for (unsigned int m(0); m < eigenvalues.size(); m++)
	      sum += eigenvalues[m];
	    eigenvalues.scale(1.0 / sum);
	    
	    // insert integer point values 
	    for (unsigned int m(0); m < eigenvalues.size(); m++)
	      svalues[(1<<reshelp)*(m+1)] = eigenvalues[m];
	    
	    // for the remaining points we use the refinement equation of phi
	    for (int j1(1); j1 <= reshelp; j1++)
	      {
		// we have to compute the values of phi at the odd multiples of 2^j1
		for (int l((1<<j1)*kbegin+1); l < (1<<j1)*kend; l += 2)
		  {
		    for (int m(kbegin); m <= kend; m++)
		      {
			// the following code can be sped up a bit
			// (by skipping the trivial entries)
			if ((1<<(reshelp-(j1-1)))*(l-(1<<(j1-1))*m)-(1<<reshelp)*kbegin >= 0
			    && (1<<(reshelp-(j1-1)))*(l-(1<<(j1-1))*m)-(1<<reshelp)*kbegin < (int) svalues.size())
			  svalues[(1<<(reshelp-j1))*l-(1<<reshelp)*kbegin]
			    += LaurentPolynomial<double>::get_coefficient(m)
			    * svalues[(1<<(reshelp-(j1-1)))*(l-(1<<(j1-1))*m)-(1<<reshelp)*kbegin];
		      }
		  }
	      }
	    
	    // Now we can calculate the values of \phi_{j,k} at the target resolution
	    // from those of \phi
	    if (j > resolution)
	      {
		// This is a special case, since phi has been sampled at a resolution
		// which we wouldn't have needed: resolution = 0
		for (unsigned int m(0); m < values.size(); m++)
		  {
		    int help = (1<<j)*a+m*(1<<(j-resolution))-k-kbegin;
		    if (help >= 0 && help < (int) svalues.size())
		      values[m] = sqrt(ldexp(1.0,j)) * svalues[help];
		  }
	      }
	    else
	      {
		for (unsigned int m(0); m < values.size(); m++)
		  {
		    int help = (1<<resolution)*a+m-(1<<reshelp)*k-(1<<reshelp)*kbegin;
		    if (help >= 0 && help < (int) svalues.size())
		      values[m] = sqrt(ldexp(1.0,j)) * svalues[help];
		  }
	      }
	  }
      }
    else
      {
	/*
	  To evaluate the n-th derivative of phi, we will use the following
	  important fact:
	  
	  1) If a^*(z) = a(z)*2/(1+z), then there exists an a^*-refinable function
	     \phi^* such that
        
             \phi'(k) = \phi^*(k) - \phi^*(k-1).

	  2) If a^*(z) = a(z)*(1+z)/2, then there exists an a^*-refinable function 
             \phi^* such that
	 
             (\phi^*)'(k) = \phi(k+1) - \phi(k).

	  References:
	  [R92] P.G. Lemari\'e-Rieusset: Analyses multi-r\'esolutions non orthogonales (...),
	  Revista Mat. Iberoamericana, 8 (1992), 221-236
	*/
	LaurentPolynomial<double> mask(*this);
	LaurentPolynomial<double> factor;
	factor.set_coefficient(0, 0.5);
	factor.set_coefficient(1, 0.5); // factor(z) = (1+z)/2
	this->divide(factor.power(DERIVATIVE), mask);

	const int kbegin(LaurentPolynomial<double>::begin().index()),
	  kend(LaurentPolynomial<double>::rbegin().index());
 	const int kstarend(kend - (int) DERIVATIVE);

  	// we compute the values of the corresponding refinable function phi^*
  	// (without the j and k parameter). In the end, we will need the
  	// following dyadic resolution:
  	int reshelp = std::max(0, resolution-j);
 	Array1D<double> svalues((kend-kbegin)*(1<<reshelp)+1);
 	for (unsigned int i(0); i < svalues.size(); i++) svalues[i] = 0;
	
  	// Compute the exact values of phi^* on the integer numbers.
  	// We assume that either phi^* is the Haar function or it is
	// at least continuous, so that
  	//   phi^*(kbegin) = phi^*(kend) = 0
 	InfiniteVector<double, int> phistar;
 	if (kstarend-kbegin-1>0)
 	  {
 	    // For the remaining point values we set up an eigenvalue problem.
 	    Matrix<double> A(kstarend-kbegin-1, kstarend-kbegin-1);
 	    for (int m1(kbegin+1); m1 < kstarend; m1++)
 	      for (int m2(kbegin+1); m2 < kstarend; m2++)
 		A(m1-kbegin-1, m2-kbegin-1) = mask.get_coefficient(2*m1-m2);
	    
 	    // We assume that the largest eigenvalue of A is 1,
 	    // so a power iteration should converge to the values of
 	    // phi at the integer points.
 	    Vector<double> eigenvalues(A.row_dimension(), true);
 	    eigenvalues = 1;
 	    unsigned int iterations;
 	    PowerIteration(A, eigenvalues, 1e-6, 100, iterations); // neglect return value
	    
 	    // normalize eigenvector to sum 1 (Strang-Fix condition 0)
 	    double sum(0);
 	    for (unsigned int m(0); m < eigenvalues.size(); m++)
 	      sum += eigenvalues[m];
 	    eigenvalues.scale(1.0 / sum);
	    
 	    // extract integer point values
 	    for (unsigned int m(0); m < eigenvalues.size(); m++)
 	      phistar[kbegin+m+1] = eigenvalues[m];
 	  }
 	else // derivative is Haar function
 	  {
 	    phistar[kbegin] = 1.0;
 	  }

  	InfiniteVector<double, int> phi(backward_difference<DERIVATIVE>(phistar));

	for (InfiniteVector<double, int>::const_iterator it(phi.begin());
	     it != phi.end(); ++it)
	  {
	    svalues[(1<<reshelp)*(it.index()-kbegin)] = *it;
	  }

  	// for the remaining points we use the refinement equation of phi
  	for (int j1(1); j1 <= reshelp; j1++)
  	  {
 	    // we have to compute the values of phi at the odd multiples of 2^j1
 	    for (int l((1<<j1)*kbegin+1); l < (1<<j1)*kend; l += 2)
 	      {
 		for (int m(kbegin); m <= kend; m++)
 		  {
 		    // the following code can be sped up a bit
 		    // (by skipping the trivial entries)
 		    if ((1<<(reshelp-(j1-1)))*(l-(1<<(j1-1))*m)-(1<<reshelp)*kbegin >= 0
 			&& (1<<(reshelp-(j1-1)))*(l-(1<<(j1-1))*m)-(1<<reshelp)*kbegin < (int) svalues.size())
 		      svalues[(1<<(reshelp-j1))*l-(1<<reshelp)*kbegin]
 			+= (1<<DERIVATIVE) * LaurentPolynomial<double>::get_coefficient(m)
 			* svalues[(1<<(reshelp-(j1-1)))*(l-(1<<(j1-1))*m)-(1<<reshelp)*kbegin];
 		  }
 	      }
 	  }

 	// Now we can calculate the values of \phi_{j,k} at the target resolution
 	// from those of \phi
 	if (j > resolution)
 	  {
 	    // This is a special case, since phi has been sampled at a resolution
 	    // which we wouldn't have needed: resolution = 0
 	    for (unsigned int m(0); m < values.size(); m++)
 	      {
 		int help = (1<<j)*a+m*(1<<(j-resolution))-k-kbegin;
 		if (help >= 0 && help < (int) svalues.size())
 		  values[m] = ldexp(1.0, j*DERIVATIVE) * sqrt(ldexp(1.0,j)) * svalues[help];
 	      }
 	  }
 	else
 	  {
 	    for (unsigned int m(0); m < values.size(); m++)
 	      {
 		int help = (1<<resolution)*a+m-(1<<reshelp)*k-(1<<reshelp)*kbegin;
 		if (help >= 0 && help < (int) svalues.size())
 		  values[m] = ldexp(1.0, j*DERIVATIVE) * sqrt(ldexp(1.0,j)) * svalues[help];
 	      }
 	  }
      }

    return SampledMapping<1>(grid, values);
  }

  template <class MASK>
  InfiniteVector<double, int> RefinableFunction<MASK>::evaluate(const unsigned int mu) const
  {
    InfiniteVector<double, int> r;

    // compute the support cube
    int suppleft(LaurentPolynomial<double>::begin().index());
    int suppright(suppleft);
    
    for (typename MASK::const_iterator it(LaurentPolynomial<double>::begin());
 	 it != LaurentPolynomial<double>::end(); ++it)
      {
	suppleft = std::min(suppleft, it.index());
	suppright = std::max(suppright, it.index());
      }

    int alpha = suppleft+1;
    int beta = suppright-1;

//     // for convenience, collect all integer points from the interior of the support cube
//     MultiIndex<int, DIMENSION> alpha, beta;
//     for (unsigned int i(0); i < DIMENSION; i++)
//       {
// 	alpha[i] = suppleft+1;
// 	beta[i] = suppright-1;
//       }
//     std::set<MultiIndex<int, DIMENSION> > indices
//       (cuboid_indices<int, DIMENSION>(alpha, beta));

    // In the following, we set up the eigenvalue problem for the values
    //   V_\alpha := D^\mu\phi(\alpha), \alpha\in\mathbb Z
    // The eigenvector is determined uniquely by the following equations,
    // see [DM] for details (so we don't need an iterative scheme):
    //
    // (3.22) eigenvalue condition
    //   2^{-\mu}V_\alpha = \sum_\beta a_{2\alpha-\beta}V_\beta, \alpha\in\mathbb Z
    //
    // (3.23) orthogonality condition
    //   \sum_\alpha (-\alpha)^\nu V_\alpha = \mu!\delta_{\mu,\nu}, \nu\le\mu

    unsigned int facmu(faculty(mu));

//     // we also collect all \nu\in\mathbb N^d, such that |\nu|\le\mu|
//     std::set<MultiIndex<unsigned int, DIMENSION> > nus;
//     for (unsigned int degnu(0); degnu <= degmu; degnu++)
//       {
// 	std::set<MultiIndex<unsigned int, DIMENSION> > sofar(nus);
// 	std::set<MultiIndex<unsigned int, DIMENSION> > plus(degree_indices<DIMENSION>(degnu));
// 	std::set_union(sofar.begin(), sofar.end(),
// 		       plus.begin(), plus.end(),
// 		       inserter(nus, nus.begin()));
//       }

    Matrix<double> A((beta-alpha+1) + (mu+1), beta-alpha+1);
    Vector<double> b((beta-alpha+1) + (mu+1));

    for (unsigned int m(0); m < A.column_dimension(); m++)
      {
	// (3.22)
	for (unsigned int n(0); n < A.column_dimension(); n++)
	  {
 	    A(m, n) = LaurentPolynomial<double>::get_coefficient(2*(alpha+m)-(alpha+n));
	  }
	A(m, m) -= ldexp(1.0, -mu);
	// b[m] = 0;
      }

    for (unsigned int m(A.column_dimension()), nu(0); m < A.row_dimension(); m++, nu++)
      {
 	// (3.23)
 	for (unsigned int n(0); n < A.column_dimension(); n++)
 	  {
 	    A(m, n) = minus1power(nu) * intpower(alpha+n, nu);
	  }
	b[m] = (nu == mu ? facmu : 0);
      }

    // the system matrix is rectangular, but Ax=b is solvable via a QR decomposition
    QRDecomposition<double> qr(A);
    assert(qr.hasFullRank());
    Vector<double> x;
    qr.solve(b, x);
    x.compress(1e-15);

    // reinterpret the entries of x
    for (unsigned int n(0); n < A.column_dimension(); n++)
      r.set_coefficient(alpha+n, x[n]);

    return r;
  }


  template <class MASK>
  double RefinableFunction<MASK>::cnk(const unsigned int n, const unsigned int k) const
  {
    double r(0);

    for (LaurentPolynomial<double>::const_iterator it(LaurentPolynomial<double>::begin());
	 it != LaurentPolynomial<double>::end(); ++it)
      r += pow(it.index(), (double)k) * *it;
    
    return r * binomial(n, k) * ldexp(1.0, -((int)n+1));
  }

  template <class MASK>
  double RefinableFunction<MASK>::moment(const unsigned int n) const
  {
    double r(1.0);
    
    if (n > 0)
      {
	r = 0.0;
	
	for (unsigned int k(0); k < n; k++)
	  r += cnk(n, n-k) * moment(k);
	
	r /= 1.0 - ldexp(1.0, -(int)n); // != 0 !
      }
    
    return r;
  }
}
