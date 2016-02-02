// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_EXTRAPOLATION_H
#define _MATHTL_EXTRAPOLATION_H

#include <algebra/triangular_matrix.h>

namespace MathTL
{
  /*!
    Romberg sequence, n_k=2^{k-1}
   */
  class RombergSequence
  {
  public:
    unsigned int n(const unsigned int k) const { return 1<<(k-1); }
  };

  /*!
    Bulirsch sequence, n_k is a power of two, alternating with 1.5*2^k
   */
  class BulirschSequence
  {
  public:
    unsigned int n(const unsigned int k) const
    {
      return (k == 1
	      ? 1
	      : (k%2 == 0
		 ? 1<<(k/2)
		 : 3*(1<<((k/2)-1))));
    }
  };

  /*!
    harmonic sequence, n_k=k
  */
  class HarmonicSequence
  {
  public:
    unsigned int n(const unsigned int k) const { return k; }
  };
  
  /*!
    For a given generic numerical algorithm T with basic "stepsize" H
    and a given sequence of positive integers
      n_1 < n_2 < ... < n_k < ...
    setup the corresponding extrapolation table.
    The order p has to be specified (asymptotic expansion in h^p).

    The algorithm T is handed over as an instance of the template parameter
    ALGORITHM and should provide a method

      void approximate(const unsigned int n, RESULT& r) const

    to compute the approximation at t+H with n "intermediate steps".
    For concrete problems, the interpretation of H and n can look
    quite differently. The Romberg quadrature, e.g., uses H as the
    interval length b-a, n is the number of subintervals for the trapezoidal sum.
    When applying extrapolation to initial value problems, H is the (time) step
    size and n the number of substeps.

    The template parameter RESULT should be a number or vector class,
    since the entries of the extrapolation table are of this type.
   */
  template <class ALGORITHM, class RESULT = double, class SEQUENCE = RombergSequence>
  class ExtrapolationTable
  {
  public:
    /*!
      constructor with a priori fixed extrapolation table size
    */
    ExtrapolationTable(const ALGORITHM& T,
		       const unsigned int size,
		       const int p = 1);

    /*!
      read access to the extrapolation table
    */
    const LowerTriangularMatrix<RESULT>& table() const { return table_; }
    
  protected:
    //! the extrapolation table
    LowerTriangularMatrix<RESULT> table_;
  };

}

// include implementation of inline functions
#include <numerics/extrapolation.cpp>

#endif
