// implementation for multi_refinable.h

#include <algorithm>
#include <cmath>
#include <iostream>
#include <utils/fixed_array1d.h>
#include <utils/tiny_tools.h>
#include <algebra/matrix.h>
#include <numerics/eigenvalues.h>
#include <numerics/differences.h>


using std::cout;
using std::endl;

namespace WaveletTL
{
  template <class MASK, unsigned int DIMENSION>
  InfiniteVector<double, MultiIndex<int, DIMENSION> >
  MultivariateRefinableFunction<MASK, DIMENSION>::evaluate() const
  {
    InfiniteVector<double, MultiIndex<int, DIMENSION> > r;

    // compute a support cube
    int suppleft(begin().index()[0]);
    int suppright(suppleft);
    
    for (typename MASK::const_iterator it(begin()); it != end(); ++it)
      {
	for (unsigned int i(0); i < DIMENSION; i++)
	  {
	    suppleft = std::min(suppleft, it.index()[i]);
	    suppright = std::min(suppright, it.index()[i]);
	  }
      }
    cout << "support cube: ["
	 << suppleft << ","
	 << suppright << "]^"
	 << DIMENSION << endl;

    // set up eigenvalue problem for the support cube
    Matrix<double> A(intpower(suppright-suppleft-1, DIMENSION));
    // HIERGEHTSWEITER

// 	    Matrix<double> A(kend-kbegin-1, kend-kbegin-1);
// 	    for (int m(kbegin+1); m < kend; m++)
// 	      for (int n(kbegin+1); n < kend; n++)
// 		A(m-kbegin-1, n-kbegin-1) = LaurentPolynomial<double>::get_coefficient(2*m-n);
    cout << "matrix for the eigenvalue problem:" << endl << A << endl;

    return r;
  }
}
