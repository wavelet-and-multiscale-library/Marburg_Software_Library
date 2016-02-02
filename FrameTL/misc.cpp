// some helper routines

#ifndef _MISC_CPP
#define _MISC_CPP

#include <math.h>
#include <cmath>

#include <algebra/matrix.h>
#include <geometry/point.h>
#include <algebra/tensor.h>

using MathTL::Matrix;
using MathTL::Tensor;

namespace FrameTL
{
  /*!
    apply matrix M to point p (p viewed as a vector)
   */
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

  /*!
    apply matrix M to point p (p viewed as a vector)
   */
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

  /*!
    add Point p, viewed as vector, to Vector b
   */
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
  

  /*!
    difference between Point p, viewed as vector, and Vector b
   */
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
  
//   /*!
//     Checks whether p lies left of right of or on the line specified by the
//     Points p1 and p2. This line is oriented by the vector starting in p1 and
//     ending in p2.
//     returning 0 means RIGHT OF LINE
//     returning 1 means LEFT OF LINE
//     returning 2 means ON LINE
//    */
//   template <unsigned int DIM>
//   inline
//   const unsigned short int pos_wrt_line (const Point<DIM> p,
// 					 const Point<DIM>& p1, const Point<DIM>&  p2)
//   {
//     double d = (p(1)-p1(1)) * (p2(0)-p1(0)) - (p(0)-p1(0)) * (p2(1)-p1(1));
    
//     if( d > 0.0 )
//       return 1;
//     else if (d < 0.0)
//       return 0;
//     else
//       return 2;
//   }

  /*!
    computes l_2 inner product
    will be placed in class Point in the future or replaced
    by appropriate function in class tensor
   */
  template <unsigned int DIM>
  inline
  const double inner_prod (const Point<DIM>& p1, const Point<DIM>& p2)
  {
    double res = .0;
    for (unsigned int i = 0; i < DIM; i++) 
      {
	res += p1(i)*p2(i);
      }
    return res;
  }

  /*!
    computes l_2 inner product
    will be placed in class Point in the future or replaced
    by appropriate function in class tensor
   */
  template <unsigned int DIM>
  inline
  const double  norm2 (Point<DIM>& p)
  {
    double norm = .0;
    for (unsigned int i = 0; i < DIM; i++) 
      {
	norm += p(i)*p(i);
      }
    return sqrt(norm);
  }


  template <unsigned int DIM>
  inline
  void swap (const Point<DIM>& p1, const Point<DIM>&  p2)
  {
    const Point<DIM>& tmp;
    tmp = p1;
    p1 = p2;
    p2 = tmp;
  }

}
#endif
