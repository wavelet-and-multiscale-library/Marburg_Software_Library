// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_COMPRESSION_H
#define _WAVELETTL_COMPRESSION_H

#include <map>
#include <list>
#include <algebra/vector.h>

namespace MathTL 
{
    template <class C>
    class Vector;
}

namespace WaveletTL
{
    //! enum for the different compression strategies, see below
    enum CompressionStrategy {
        CDD1,
        St04a,
        tensor_simple
    };
    
    /*!
      For a given compressible matrix A stemming from the discretization of
      a differential or integral operator in wavelet coordinates, realize the column
      compression stragegies from [CDD1], [St04a] or [DSS09].

      Isotropic wavelets: (eg. cubeequation)
      The lambda-th column of A will be compressed by the J-th truncation rule
      from [CDD1, Prop. 3.4] or [St04a, Th. 2.3 and Th. 3.3]
      It is also possible to specify a fixed maximal level jmax.

      (commented code:)
      Anisotropic wavelets: (eg. tbasis_equation)
      The lambda-th column of A will be compressed by the rule from [DSS09],
      modified for the biorthogonal setting of tbasis/qtbasis. We use that in this case the
      matrix is compressible with s^*=\infty. The meaning of alphak(int x) is assumed
      to be constant.
      You can specify jmax, but observe the meaning ||level_of_lambda||_1 <= jmax
      in this setting.
  
      For integral equations, both a compression in scale and space has to be performed, while
      for differential equations, only a compression in scale is needed,
      see [CDD1, Prop. 3.4].
  
      The compressible matrix A should provide the unpreconditioned bilinear form a
      and the diagonal preconditioner D.
      The problem class PROBLEM should indicate via a function local(), whether the
      matrix A is induced by a local operator or not.
  
      (Remark: this routine could be extended for other compression strategies or
      alternative methods to compute single columns of the stiffness matrix.)
  
      References:
      [CDD1]  Cohen/Dahmen/DeVore:
              Adaptive Wavelet Methods for Elliptic Operator Equations - Convergence Rates
      [St04a] Stevenson:
              On the compressibility of operators in wavelet coordinates
      [DSS09] Dijkema/Schwab/Stevenson:
              An Adaptive Wavelet Method for Solving High-Dimensional Elliptic PDES
    */
    
  
    using MathTL::Vector;

    template <class PROBLEM>
    void add_compressed_column(const PROBLEM& P,
                               const double factor,
                               const typename PROBLEM::Index& lambda,
			       const int J,
			       Vector<double>& w,
			       const int jmax = 999,
			       const CompressionStrategy strategy = St04a,
                               const bool preconditioning = true); // only relevant for the anisotropic case
    template <class PROBLEM>
    void add_compressed_column(const PROBLEM& P,
			       const int p,
			       const double factor,
			       const typename PROBLEM::Index& lambda,
			       const int J,
			       Vector<double>& w,
			       const int jmax = 999,
			       const CompressionStrategy strategy = St04a,
                               const bool preconditioning = true); // only relevant for the anisotropic case);


}

#include <adaptive/compression.cpp>


#endif
