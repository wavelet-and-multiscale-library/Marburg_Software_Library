// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2007                                            |
// | Thorsten Raasch, Manuel Werner, Andreas Schneider                  |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_BIHARMONIC_EQUATION_H
#define _WAVELETTL_BIHARMONIC_EQUATION_H

#include <utils/function.h>
#include <utils/array1d.h>
#include <algebra/infinite_vector.h>
#include <interval/spline_basis.h>
#include <galerkin/galerkin_utils.h>
#include <galerkin/infinite_preconditioner.h>
#include <adaptive/compression.h>

using namespace MathTL;

namespace WaveletTL
{
  /*!
    The biharmonic equation -u^(4)=f, u(0)=u(1)=u'(0)=u'(1)=0 in one space dimension.
  */
  template <class IBasis>
  class BiharmonicEquation1D
    :  public FullyDiagonalEnergyNormPreconditioner<typename IBasis::Index>
  {
  public:
    /*!
      type of the wavelet basis
    */
    typedef IBasis WaveletBasis;

    /*!
      constructor from a given wavelet basis and a right-hand side g
    */
    BiharmonicEquation1D(const WaveletBasis& basis, const Function<1>* g);

    /*!
      read access to the basis
    */
    const WaveletBasis& basis() const { return basis_; }
    
    /*!
      wavelet index class
    */
    typedef typename WaveletBasis::Index Index;

    /*!
      index type of vectors and matrices
    */
    typedef typename Vector<double>::size_type size_type;

    /*!
      space dimension of the problem
    */
    static const int space_dimension = 1;

    /*!
      differential operators are local
    */
    static bool local_operator() { return true; }

    /*!
      (half) order t of the operator
      (inherited from FullyDiagonalDyadicPreconditioner)
    */
    double operator_order() const { return 2.; }
    
    /*!
      evaluate the diagonal preconditioner D
    */
    double D(const Index& lambda) const;

    /*!
      evaluate the (unpreconditioned) bilinear form a
      (inherited from FullyDiagonalEnergyNormPreconditioner)
     */
    double a(const Index& lambda,
             const Index& nu) const;

    /*!
      evaluate the (unpreconditioned) bilinear form a;
      you can specify the order p of the quadrature rule, i.e.,
      (piecewise) polynomials of maximal degree p will be integrated exactly.
      Internally, we use an m-point composite Gauss quadrature rule adapted
      to the singular supports of the spline wavelets involved,
      so that m = (p+1)/2;
    */
    double a(const Index& lambda,
             const Index& nu,
             const unsigned int p) const;

    /*!
      evaluate the (unpreconditioned) right-hand side f
    */
    double f(const typename WaveletBasis::Index& lambda) const;

    /*!
      estimate compressibility exponent s^*
    */
    double s_star() const {
      return 1.0 + WaveletBasis::primal_vanishing_moments(); // [St04a], Th. 2.3 for n=1
    }

    /*!
      w += factor * (stiffness matrix entries in column lambda on level j)
    */
    void add_level (const Index& lambda,
                    InfiniteVector<double, Index>& w, const int j,
                    const double factor,
                    const int J,
                    const CompressionStrategy strategy = St04a) const;

  protected:
    const WaveletBasis& basis_;
    const Function<1>* g_;
  };

}

#include <galerkin/biharmonic_equation.cpp>

#endif // _WAVELETTL_BIHARMONIC_EQUATION_H
