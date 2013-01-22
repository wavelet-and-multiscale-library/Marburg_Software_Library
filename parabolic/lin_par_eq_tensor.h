// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner, Ulrich Friedrich                   |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_LIN_PAR_EQ_TENSOR_H
#define _WAVELETTL_LIN_PAR_EQ_TENSOR_H

#include <algebra/infinite_vector.h>
#include <numerics/ivp.h>
#include <numerics/w_method.h>
#include <utils/function.h>

//#include <galerkin/gramian.h>
#include <galerkin/cached_tproblem.h>

#include <galerkin/infinite_preconditioner.h>
//#include <adaptive/apply.h>
#include <adaptive/apply_tensor.h>
#include <adaptive/cdd1.h>

using MathTL::InfiniteVector;
using MathTL::AbstractIVP;
using MathTL::Function;
using MathTL::WMethodPreprocessRHSHelper;

namespace WaveletTL
{
  /*!
   * ROW stage equation helper class for linear parabolic equations (see below).
   * The helper class models the preconditioned stage equation
   *
   *       D^{-1}(alpha*I-T)D^{-1}*Dx = D^{-1}y
   *
   * (after multiplication with the Gramian, see implementation)
   *
   * The class has the minimal signature to be used within the APPLY routine
   * or within adaptive solvers like CDD1.
  */
  template <class CACHEDTPROBLEM>
  class LinParEqTenROWStageEquationHelper
    : public FullyDiagonalEnergyNormPreconditioner<typename CACHEDTPROBLEM::Index>
  {
  public:
    /*!
      constructor from alpha, T, the Gramian and y
    */
    LinParEqTenROWStageEquationHelper
    (const double alpha,
     const CACHEDTPROBLEM* T,
     const CACHEDTPROBLEM* G,
     const InfiniteVector<double, typename CACHEDTPROBLEM::Index>& y);

    /*!
      make wavelet basis type accessible
    */
    typedef typename CACHEDTPROBLEM::WaveletBasis WaveletBasis;

    /*!
      wavelet index class
    */
    typedef typename CACHEDTPROBLEM::Index Index;

    /*!
     * Type of wavelet level
     */
    typedef typename Index::level_type index_lt;

    /*!
      read access to the basis
    */
    const WaveletBasis& basis() const { return T->basis(); }

    /*!
      space dimension of the problem
    */
    static const int space_dimension = CACHEDTPROBLEM::space_dimension;

    /*!
      locality of the operator
    */
    static bool local_operator()
    {
        return true; //CACHEDTPROBLEM::local_operator();
    }

    /*!
      (half) order t of the operator
    */
    double operator_order() const
    {
        return T->operator_order();
    }

    /*!
      evaluate the diagonal preconditioner D defined by T
    */
    double D(const Index& lambda) const
    {
        return sqrt(a(lambda, lambda));
    }

    /*!
      evaluate the (unpreconditioned) bilinear form a
    */
    double a(const Index& lambda,
	     const Index& nu) const
    {
// TODO: for local operators cache misshits will happen for G and T at the same time. Introduce a_save that assumes value is already stored in the cache
      return alpha * G->a(lambda, nu) + T->a(lambda, nu);
    }

    /*!
      estimate the spectral norm ||alpha*D^{-2}-A|| from above
    */
    double norm_A() const
    {
      return T->norm_A(); // dirty!
    }

    /*!
      estimate the spectral norm ||(alpha*D^{-2}-A)^{-1}|| from above
    */
    double norm_Ainv() const
    {
      return T->norm_Ainv(); // dirty!
    }

    /*!
      estimate compressibility exponent s^*
    */
//    double s_star() const {
//      return T->s_star();
//    }

    /*!
     * alpha_k is used differently for tensorbasis
    */
    double alphak(const unsigned int k) const {
// TODO check this line
      return T->alphak(10); // argument does not influence the value
    }

    /*!
      evaluate the (unpreconditioned) right-hand side f
    */
    double f(const Index& lambda) const {
      return y.get_coefficient(lambda);
    }

    /*!
      approximate the wavelet coefficient set of the preconditioned right-hand side F
      within a prescribed \ell_2 error tolerance
    */
    void RHS(const double eta,
	     InfiniteVector<double, Index>& coeffs) const
    {
        // dirty
        coeffs = y_scaled;
    }

    /*!
      compute (or estimate) ||F||_2
    */
    double F_norm() const
    {
        return l2_norm(y_scaled);
    }

    /*
     * Called by add_compresed_column
     * w += factor * (stiffness matrix entries in column lambda with ||nu-lambda|| <= range && ||nu|| <= maxlevel)
     * Hack? : works with non dynamic vector w of size = degrees_of_freedom
     *
     * "strategy" argument ist needed for compatibility with compression.h .
     * There is only 1 strategy for TBasis at the moment
     */
    void add_ball(const Index& lambda,
                 //InfiniteVector<double, Index>& w,
                 Vector<double>& w,
                 const int radius,
                 const double factor,
                 const int maxlevel,
                 const CompressionStrategy strategy = tensor_simple,
                 const bool precond = true) const;
   

    /*! Set the rhs vector y. compute y_scaled for RHS method*/
    void set_y(const InfiniteVector<double, Index> ynew);
    /*! Set the parameter alpha*/
    void set_alpha (const double alpha);

  protected:
    double alpha;
    const CACHEDTPROBLEM* T;
    const CACHEDTPROBLEM* G;
    InfiniteVector<double, typename CACHEDTPROBLEM::Index> y;
    InfiniteVector<double, typename CACHEDTPROBLEM::Index> y_scaled;
  };

  /*!
    This class models a linear parabolic equation of the form

      u'(t) = Au(t) + f(t) =: F(t,u(t)),  0 < t <= T
       u(0) = u_0

    where -A: \ell_{2,D} -> \ell_{2,D^{-1}} is a positive isomorphism
    and f: (0,T] -> \ell_{2,D^{-1}}, D being a diagonal preconditioner (see below).
    An equation of this form arises when equivalently reformulating
    a problem of the analogous form

      v'(t) = Bv(t) + g(t) =: G(t,v(t)), 0 < t <= T
       v(0) = v_0

    with isomorphism -B: H -> H' and right-hand side g: (0,T] -> H'
    by means of a biorthogonal wavelet basis Psi=\{\psi_\lambda\},
    setting F(t,v):=<G(t,v^\top\Psi),\tilde\Psi> for all v in \ell_{2,D}.

    The (unpreconditioned) stiffness matrix

      -A = (a(\psi_\lambda',\psi_\lambda))_{\lambda,\lambda'}

    is modeled in the template parameter class CACHEDTPROBLEM.

  */
  template <class CACHEDTPROBLEM>
  class LinearParabolicEquationTensor
    : public AbstractIVP<InfiniteVector<double, typename CACHEDTPROBLEM::Index> >,
      public WMethodPreprocessRHSHelper<InfiniteVector<double, typename CACHEDTPROBLEM::Index> >
  {
  public:
    /*!
      the index class, for convenience
    */
    typedef typename CACHEDTPROBLEM::Index Index;

    /*!
      space dimension of the problem
    */
    static const int space_dimension = CACHEDTPROBLEM::space_dimension;

    /*!
      constructor from a helper object for the stiffness matrix,
      a given initial value u0 in ell_2
      and a driving term f (coeffs. w.r.t. the dual basis) which is constant in time
    */
    LinearParabolicEquationTensor(const CACHEDTPROBLEM* elliptic,
                                  CACHEDTPROBLEM* identity,
                                  const InfiniteVector<double,Index>& u0,
                                  const InfiniteVector<double,Index>& f,
                                  const int jmax = 10);

    /*!
      constructor from a helper object for the stiffness matrix,
      a given initial value u0 in ell_2
      and a time-dependent driving term f

      (gcc 2.95 requires implementation of this constructor already in the *.h file!?)
    */
    LinearParabolicEquationTensor(const CACHEDTPROBLEM* ellipt,
                                  CACHEDTPROBLEM* identi,
  			          const MathTL::InfiniteVector<double,Index>& initial,
                                  MathTL::Function<CACHEDTPROBLEM::space_dimension,double>* f = 0,
                                  const int jmax = 10)
      : elliptic(ellipt), identity(identi),
 	constant_f_(), f_(f), jmax_(jmax)
    {
      AbstractIVP<InfiniteVector<double,typename CACHEDTPROBLEM::Index> >::u0 = initial;
    }

    /*!
      evaluate the right-hand side F(t,v)=Av+f(t) up to a prescribed tolerance
    */
    void evaluate_f(const double t,
		    const InfiniteVector<double,Index>& v,
		    const double tolerance,
		    InfiniteVector<double,Index>& result) const;

    /*!
      evaluate F_t(t,v)=f'(t) up to a prescribed tolerance
    */
    void evaluate_ft(const double t,
		     const InfiniteVector<double,Index>& v,
		     const double tolerance,
		     InfiniteVector<double,Index>& result) const;

    /*!
      solve (alpha*I-A)x=y up to a prescribed tolerance
    */
    void solve_ROW_stage_equation(const double t,
				  const InfiniteVector<double,Index>& v,
				  const double alpha,
				  const InfiniteVector<double,Index>& y,
				  const double tolerance,
				  InfiniteVector<double,Index>& result) const;

    /*!
      preprocess the right-hand side share from the previous steps
      (cf. w_method.h for details)
     Transforms a primal coeficient vector, ie (<u,\tilde\psi_\lambda>) into a dual one, ie (<u,\psi_\lambda>)
    */
    void preprocess_rhs_share(InfiniteVector<double,Index>& wbeforeandafter,
			      const double tolerance) const
    {
      InfiniteVector<double,Index> help(wbeforeandafter);
      //help.scale(identity,1); // should be 1 anyways ... although it isn't ...
      APPLY(*identity, help, tolerance, wbeforeandafter, jmax_, tensor_simple, false);
      //wbeforeandafter.scale(identity,1);
    }

  protected:
    //! pointer to the elliptic subproblem helper
    const CACHEDTPROBLEM* elliptic;

    //! pointer to the cached gramian
    CACHEDTPROBLEM* identity;

    //! Gramian
    //IntervalGramian<typename CACHEDTPROBLEM::WaveletBasis> G;

    //! cached Gramian
    //CachedProblem<IntervalGramian<typename CACHEDTPROBLEM::WaveletBasis> > GC;

    InfiniteVector<double,Index> constant_f_;
    Function<space_dimension>* f_;

    //! maximal level
    const int jmax_;
  };

}

#include <parabolic/lin_par_eq_tensor.cpp>

#endif
