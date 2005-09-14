// some helper routines
#include <math.h>
#include <cmath>

#include <algebra/matrix.h>
#include <geometry/point.h>

using MathTL::Matrix;

namespace FrameTL
{

  template <unsigned int DIM, class VALUE>
  inline
  const Point<DIM> apply (const Matrix<VALUE>& M, const Point<DIM>& p)
  {
    
    assert(DIM == M.column_dimension());
    Point<DIM> res;
    for (typename Matrix<VALUE>::size_type i(0); i < M.row_dimension(); i++) 
      {
	for (typename Matrix<VALUE>::size_type j(0); j < M.column_dimension(); j++)
	  {
	    res(i) += M(i,j) * p(j);
	  } 
      }
    return res;
  }
}
