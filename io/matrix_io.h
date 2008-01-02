// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2008                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_MATRIX_IO_H
#define _MATHTL_MATRIX_IO_H

#include <iostream>
#include <iomanip>

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
      [x_{1,1} x_{1,2} ... x_{1,n}; x_{2,1} ... x_{m,n}];
  */
  template <class MATRIX>
  void print_matrix(const MATRIX& M, std::ostream& os)
  {
    os << "[";
    unsigned int precision=15, tabwidth=10;
    unsigned int old_precision = os.precision(precision);
    for (unsigned int row(0); row < M.row_dimension(); row++)
      {
	for (unsigned int column(0); column < M.column_dimension(); column++)
	  {
	    os << std::setw(tabwidth) << std::setprecision(precision)
	       << M.get_entry(row, column);
	    if (column < M.column_dimension()-1)
	      os << " ";
	  }
	if (row < M.row_dimension()-1)
	  os << "; ";
      }
    os << "];" << std::endl;
    os.precision(old_precision);
  }
}

#endif
