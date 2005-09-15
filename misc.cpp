// some helper routines

#ifndef _MISC_CPP
#define _MISC_CPP

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

  template <unsigned int DIM, class VALUE>
  inline
  void apply (const Matrix<VALUE>& M, const Point<DIM>& p, Point<DIM>& res)
  {    
    assert(DIM == M.column_dimension());
    for (typename Matrix<VALUE>::size_type i(0); i < M.row_dimension(); i++) 
      {
	for (typename Matrix<VALUE>::size_type j(0); j < M.column_dimension(); j++)
	  {
	    res(i) += M(i,j) * p(j);
	  } 
      }
  }

  template <unsigned int DIM>
  inline
  Point<DIM> add (const Point<DIM>& p, const Vector<double>& b)
  {    
    Point<DIM> res;
    for (typename Point<DIM>::size_type i(0); i < p.size(); i++) 
      {
	res(i) = p(i) + b(i);
      }
    return res;
  }
  
  template <unsigned int DIM>
  inline
  Point<DIM> sub (const Point<DIM>& p, const Vector<double>& b)
  {    
    Point<DIM> res;
    for (typename Point<DIM>::size_type i(0); i < p.size(); i++) 
      {
	res(i) = p(i) - b(i);
      }
    return res;
  }
}
#endif
