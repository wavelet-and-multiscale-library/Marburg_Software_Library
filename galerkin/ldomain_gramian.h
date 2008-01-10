// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2008                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_LDOMAIN_GRAMIAN_H
#define _WAVELETTL_LDOMAIN_GRAMIAN_H

#include <set>
#include <utils/fixed_array1d.h>
#include <utils/array1d.h>
#include <numerics/bvp.h>

#include <interval/spline_basis.h>
#include <galerkin/galerkin_utils.h>
#include <galerkin/infinite_preconditioner.h>
#include <galerkin/ldomain_equation.h>

using MathTL::FixedArray1D;
using MathTL::EllipticBVP;

namespace WaveletTL
{
  template <int d, int dT, SplineBasisFlavor flavor, int s0, int s1, int sT0, int sT1> class SplineBasis;
 
  template <class IBASIS> class LDomainBasis;

  /*!
    This class models the Gramian matrix of a composite wavelet basis
    on the L--shaped domain in R^2.
  */
  template <class IBASIS>
  class LDomainGramian
    : public LDomainEquation<IBASIS>
  {
  public:
    /*!
      make template argument accessible
    */
    typedef LDomainBasis<IBASIS> WaveletBasis;

    /*!
      wavelet index class
    */
    typedef typename WaveletBasis::Index Index;

    /*!
      constructor from an identity bvp
    */
    LDomainGramian(IdentityBVP<2>* bvp)
      : LDomainEquation<IBASIS>(bvp, false)
    {
    }

    /*!
      (half) order t of the operator
    */
    double operator_order() const { return 0.; }
    
    /*!
      evaluate the diagonal preconditioner D (we don't have any)
    */
    double D(const Index& lambda) const { return 1; }
  };
  
//   // template specialization for the case IBASIS==SplineBasis<d,dT,DS_construction>
//   template <int d, int dT> class LDomainGramian<SplineBasis<d,dT,DS_construction> >;
// 
//   template <int d, int dT>
//   class LDomainGramian<SplineBasis<d,dT,DS_construction> >
// //     : public FullyDiagonalEnergyNormPreconditioner<typename LDomainBasis<SplineBasis<d,dT,DS_construction> >::Index>
// // note: energy norm preconditioning (and the corresponding branch in D()) is commented out
// //       since this is incompatible with the current implementation of LDomainHelmholtzEquation!!!
//   {
//   public:
//     /*!
//       the 1D wavelet basis class
//     */
//     typedef SplineBasis<d,dT,DS_construction> IntervalBasis;
//     
//     /*!
//       the wavelet basis class
//     */
//     typedef LDomainBasis<IntervalBasis> WaveletBasis;
//     
//     /*!
//       read access to the basis
//     */
//     const WaveletBasis& basis() const { return basis_; }
// 
//     /*!
//       wavelet index class
//     */
//     typedef typename WaveletBasis::Index Index;
// 
//     /*!
//       constructor from a given wavelet basis and a given right-hand side y
//     */
//     LDomainGramian(const WaveletBasis& basis,
//  		   const InfiniteVector<double,Index>& y);
//     
//     /*!
//       set right-hand side y
//     */
//     void set_rhs(const InfiniteVector<double,Index>& y) const {
//       y_ = y;
//     }
// 
//     /*!
//       space dimension of the problem
//     */
//     static const int space_dimension = 2;
// 
//     /*!
//       differential operators are local
//     */
//     static bool local_operator() { return true; }
// 
//     /*!
//       (half) order t of the operator
//     */
//     double operator_order() const { return 0.; }
//     
//     /*!
//       estimate compressibility exponent s^*
//     */
//     double s_star() const {
//       return std::min(WaveletBasis::primal_regularity(),
// 		      WaveletBasis::primal_vanishing_moments()/2.); // [St04a], Th. 2.3 for n=2,t=0,gamma
//     }
//     
//     /*!
//       evaluate the diagonal preconditioner D
//     */
//     double D(const Index& lambda) const {
// //       return sqrt(a(lambda,lambda));
//       return 1.0;
//     }
// 
//     /*!
//       evaluate the (unpreconditioned) bilinear form a
//     */
//     double a(const Index& lambda,
//  	     const Index& nu) const;
// 
//     /*!
//       evaluate the (unpreconditioned) right-hand side f
//     */
//     double f(const Index& lambda) const {
//       return y_.get_coefficient(lambda);
//     }
// 
//   protected:
//     const WaveletBasis& basis_;
//     
//     // rhs, mutable to have 'const' method
//     mutable InfiniteVector<double,Index> y_;
//     
//     // estimates for ||A|| and ||A^{-1}||
//     mutable double normA, normAinv;
//   };
}

#include <galerkin/ldomain_gramian.cpp>

#endif
