// implementation for refinable.h

#include <algorithm>
#include <cmath>
#include <iostream>
#include <utils/array1d.h>
#include <algebra/matrix.h>
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

}
