// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2004                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_MATRIX_IO_H
#define _MATHTL_MATRIX_IO_H

#include <iostream>

namespace MathTL
{
  /*
    generic input/output routines for MATRIX classes with a standard
    signature
    We use the row_dimension() and column_dimension() routines
    for determining the dimensions of the matrix.
  */
  
  /*!
    generic Matlab-style matrix stream output of the form
    [x_{1,1} x_{1,2} ... x_{1,n}; x_{2,1} ... x_{m,n}]
  */
  template <class MATRIX>
  void print_matrix(const MATRIX& M, std::ostream& os)
  {
    os << "[";
    for (unsigned int row(0); row < M.row_dimension(); row++)
      {
	for (unsigned int column(0); column < M.column_dimension(); column++)
	  {
	    os << M(row, column);
	    if (column < M.column_dimension()-1)
	      os << " ";
	  }
	if (row < M.row_dimension()-1)
	  os << "; ";
      }
    os << "]";
  }
}

#endif
