// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_MATRIX_NORMS_H
#define _MATHTL_MATRIX_NORMS_H

#include <algorithm>

namespace MathTL
{
  /*
    computation of different norms for generic MATRIX classes

    You can use template specialization to speed up the norm computation
    for specific MATRIX classes.
  */

  /*!
    maximum row sum \|M\|_\infty
    (the corresponding operator norm for the l_\infty norm on the vector space)
   */
  template <class MATRIX>
  const double row_sum_norm(const MATRIX& M)
  {
    double r(0.0);

    for (typename MATRIX::size_type i(0); i < M.row_dimension(); i++)
      {
	double rowsum(0.0);
	for (typename MATRIX::size_type j(0); j < M.column_dimension(); j++)
	  rowsum += fabs(M.get_entry(i,j));
	r = std::max(r, rowsum);
      }

    return r;
  }

  /*!
    maximum column sum \|M\|_1
    (the corresponding operator norm for the l_1 norm on the vector space)
   */
  template <class MATRIX>
  const double column_sum_norm(const MATRIX& M)
  {
    double r(0.0);

    for (typename MATRIX::size_type j(0); j < M.column_dimension(); j++)
      {
	double colsum(0.0);
	for (typename MATRIX::size_type i(0); i < M.column_dimension(); i++)
	  colsum += fabs(M.get_entry(i,j));
	r = std::max(r, colsum);
      }

    return r;
  }

  /*!
    Frobenius norm \|M\|_F=(\sum_{i,j}|m_{i,j}|^2)^{1/2}
    (this is compatible to the l_2 norm on the vector space, but not the
    corresponding operator norm)
   */
  template <class MATRIX>
  const double frobenius_norm(const MATRIX& M)
  {
    double r(0.0);
    
    for (typename MATRIX::size_type j(0); j < M.column_dimension(); j++)
      for (typename MATRIX::size_type i(0); i < M.column_dimension(); i++)
	r += M.get_entry(i,j)*M.get_entry(i,j);
    
    return sqrt(r);
  }

  /*!
    maximum absolute entry value of a matrix
   */
  template <class MATRIX>
  const double maximum_norm(const MATRIX& M)
  {
    double r(0.0);

    for (typename MATRIX::size_type i(0); i < M.row_dimension(); i++)
      for (typename MATRIX::size_type j(0); j < M.column_dimension(); j++)
	r = std::max(r, M.get_entry(i,j));
    
    return r;
  }
  
}

#endif
