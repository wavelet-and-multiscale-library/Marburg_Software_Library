// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of FrameTL - the Frame Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Manuel Werner                                                      |
// +--------------------------------------------------------------------+

#ifndef _FRAMETL_SIMPLE_ELLIPTIC_EQUATION_H
#define _FRAMETL_SIMPLE_ELLIPTIC_EQUATION_H

#include <aggregated_frame.h>
#include <numerics/bvp.h>
#include <adaptive/compression.h>
#include <interval/i_index.h>
#include <galerkin/infinite_preconditioner.h>
#include <frame_support.h>
#include <index1D.h>

using FrameTL::AggregatedFrame;
using MathTL::EllipticBVP;
using WaveletTL::CompressionStrategy;
using WaveletTL::IntervalIndex;
using WaveletTL::FullyDiagonalEnergyNormPreconditioner;

namespace FrameTL
{

  /*!
    This class models the (preconditioned) infinite-dimensional matrix problem
    
    \f$Au = D^{-1}LD^{-1}u = D^{-1}F\f$

    when reformulating a symmetric, second-order elliptic
    boundary value problem in divergence form over some domain
    Omega in \f$R^d\f$ with boundary \f$\Gamma=\partial \Omega\f$,
    with homogeneous Dirichlet boundary conditions

    \f$-\mbox{div}(a(x)\nabla u(x)) + q(x)u(x) = f(x)\f$ in \f$\Omega\f$<br>
                             \f$u(x) = 0\f$ on \f$\Gamma\f$<br>
                             
    The corresponding bilinear form in

    \f$L = (a(\psi_\nu,\psi_\lambda))_{\lambda,\nu}\f$

    is

    \f$a(u,v) = \int_\Omega \langle a(x) \nabla u(x), \nabla v(x)\rangle  dx +
              \int_\Omega q(x) u(x) v(x) dx\f$
     
    and the right-hand side is
     
    \f$f(v) = \int_\Omega f(x) v(x)  dx\f$.

    The evaluation of \f$a(.,.)\f$ and \f$f\f$ is possible for arguments \f$\psi_\lambda\f$
    which stem from an aggregated wavelet frame \f$\Psi=\{\psi_\lambda\}\f$ of the corresponding
    function space over \f$\Omega\f$.

    WE ASSUME THAT THE COEFFICIENTS OF THE ELLIPTIC PDE ARE SEPERABLE AND SMOOTH AND THAT THE PATCHES
    OF THE UNDERLYING DOMAIN DECOMPOSITION ARE RECATANGULAR AND ALIGNED WITH THE CARTESIAN
    COORDINATES.
    FOR THIS SPECIAL CASE, a(.,.) CAN BE EXACTLY COMPUTED AT UNIT COST AND TENSOR PRODUCT
    STRUCTURE CAN BE EXPLOITED.

    @tparam IBASIS The type of interval basis underlying the construction of the aggregated frame.
    @tparam DIM The dimension of the underlying domain.
   */
  template <class IBASIS, unsigned int DIM>
  class SimpleEllipticEquation
  //: public FullyDiagonalDyadicPreconditioner<typename AggregatedFrame<IBASIS,DIM>::Index>
      : public FullyDiagonalEnergyNormPreconditioner<typename AggregatedFrame<IBASIS,DIM>::Index>
  {
  public:

    /*!
      Constructor. The diagonal of the stiffness matrix and the coefficients of the right-hand side
      are precomputed between minimal and maximal level.
      
      @param ell_bvp The elliptic boundary value problem that is modeled.
      @param frame Pointer to the aggregated frame that is used for discretization.
      @param jmax The maximal level of resolution that is considered.
     */
    SimpleEllipticEquation(const EllipticBVP<DIM>* ell_bvp,
			   const AggregatedFrame<IBASIS,DIM>* frame,
			   const int jmax);

    /*!
      The frame type.
     */
    typedef AggregatedFrame<IBASIS,DIM> Frame;

    /*!
      Dummy typedef to be compatible with WaveletTL
      routines.
     */
    typedef AggregatedFrame<IBASIS,DIM> WaveletBasis;
    
    /*!
      The index type.
    */
    typedef typename Frame::Index Index;

    /*!
      Read access to the frame.
    */
    const AggregatedFrame<IBASIS,DIM>& frame() const { return *frame_; }

    /*!
      Read access to the boundary value problem.
    */
    const EllipticBVP<DIM>&  get_bvp() const { return *ell_bvp_; }

    /*!
      Read access to the frame. The routine is called basis() to be
      compatible with the routines in WaveletTL's compression.h.
    */
    const AggregatedFrame<IBASIS,DIM>& basis() const { return *frame_; }  

    /*!
      Space dimension of the problem.
    */
    static const int space_dimension = DIM;

    /*!
      Differential operators are local.
    */
    static bool local_operator() { return true; }

    /*!
      Order of the operator.
    */
    double operator_order() const { return 1; }
    
    /*!
      Evaluate the diagonal preconditioner D.
    */
    double D(const typename AggregatedFrame<IBASIS,DIM>::Index& lambda) const;

    /*!
      Rescale a coefficient vector by an integer power of D, \f$c \mapsto D^{n}c\f$.
    */
    void rescale(InfiniteVector<double, typename AggregatedFrame<IBASIS,DIM>::Index>& coeffs,
		 const int n) const;

    /*!
      Evaluate the (unpreconditioned) bilinear form a.
    */
    double a(const typename AggregatedFrame<IBASIS,DIM>::Index& lambda,
	     const typename AggregatedFrame<IBASIS,DIM>::Index& nu) const;

    /*!
      Estimate the spectral norm \f$\|A\|\f$.
    */
    double norm_A() const;

    /*!
      Returns spectral norm \f$\|A^{-1}\|\f$.
      An estimate for \f$\|A^{-1}\|\f$ has to be
      externally computed and to be set
      during initialization of the program.
    */
    double norm_Ainv() const { return normAinv; };

    /*!
      Sets estimate for \f$\|A\|\f$.
    */
    void set_norm_A(const double _normA) { normA = _normA; }

    /*!
      Sets estimate for \f$\|A^{-1}\|\f$.
    */
    void set_Ainv(const double nAinv) { normAinv = nAinv; };

    /*!
      Estimate compressibility exponent \f$s^\ast\f$.
    */
    double s_star() const;

    /*!
      Estimate the compression constants alpha_k in
      \f$\|A-A_k\| \leq \alpha_k  2^{-sk}\f$
    */
    double alphak(const unsigned int k) const {
      //return pow(2,(-k))*norm_A(); // suboptimal
      return 2.*norm_A(); // suboptimal
    }

    /*!
      Evaluate the (unpreconditioned) right-hand side f.
    */
    double f(const typename AggregatedFrame<IBASIS,DIM>::Index& lambda) const;

    /*!
      Approximate the wavelet coefficient set of the preconditioned right-hand side
      within a prescribed \f$\ell_2\f$ error tolerance.
    */
    void RHS(const double eta, InfiniteVector<double, 
	     typename AggregatedFrame<IBASIS,DIM>::Index>& coeffs) const;

    /*!
      Approximate the wavelet coefficient set of the preconditioned right-hand side restricted
      to patch p within a prescribed \f$\ell_2\f$ error tolerance.
    */
    void RHS(const double eta, const int p,
	     InfiniteVector<double, 
	     typename AggregatedFrame<IBASIS,DIM>::Index>& coeffs) const;
    
    /*!
      Compute (or estimate) the \f$\ell_2\f$ norm of the right-hand side.
    */
    double F_norm() const { return sqrt(fnorm_sqr); }
    
    /*!
      Compute the \f$\ell_2\f$ norm of the right-hand side coefficients on corresponding to a fixed patch.
    */
    double F_norm_local(const int patch) const { return sqrt(fnorms_sqr_patch[patch]); }

    /*!
      Set the boundary value problem.
    */
    void set_bvp(const EllipticBVP<DIM>*);


    /*!
      Multiplies the stiffness matrix entries of column lambda on level j of the compressed martrix A_J
      by factor and adds the result to w. 
    */
    void add_level (const Index& lambda,
		    InfiniteVector<double, Index>& w, const int j,
		    const double factor,
		    const int J,
		    const CompressionStrategy strategy) const;

   protected:
    
    /*!
      The elliptic boundary value problem.
     */
    const EllipticBVP<DIM>* ell_bvp_;

    /*!
      The underlying aggregated frame.
    */
    const AggregatedFrame<IBASIS,DIM>* frame_;

    // #####################################################################################
    // Caching of appearing 1D integrals when making use of the tensor product structure of the wavelets
    // during the evaluation of the bilinear form.
    typedef std::map<Index1D<IBASIS>,double > Column1D;
    typedef std::map<Index1D<IBASIS>,Column1D> One_D_IntegralCache;
    
    mutable One_D_IntegralCache one_d_integrals;
    // #####################################################################################

  private:

    /*!
      Calculate integral of the product of the interval wavelets or generators
      given by lambda and mu. The class Index1D<IBASIS> also includes
      generators on levels above the coarsest one. This routine is used
      during the evaluation of the bilinear form a(.,.) when use of the tensor
      product structure of the multi-dimensional integrand is made.

      @param lambda Index of first wavelet or generator.
      @param mu Index of second wavelet or generator.
      @param N_Gauss Number of Gauss quadrature knots to be used. This should be cosen equal
      to the spline order in the constant coefficient case to be sure to integrate exactly.
      @param dir The spatial direction under considerattion.
      @param supp_lambda Support of the reference wavelet on the cube having the function given by lambda
      as component in the direction dir.
      @param supp_mu Support of the reference wavelet on the cube having the function given by mu
      as component in the direction dir.
     */
    double integrate(const Index1D<IBASIS>& lambda,
		     const Index1D<IBASIS>& mu,
		     const int N_Gauss,
		     const int dir,
		     const typename CubeBasis<IBASIS,DIM>::Support* supp_lambda,
		     const typename CubeBasis<IBASIS,DIM>::Support* supp_mu) const;

    //! Precompute the right-hand side between minimal and maximal level.
    void compute_rhs();

    //! Precompute the diagonal of the stiffness matrix between minimal and maximal level.
    void compute_diagonal();

    //! The maximal level of resolution.
    const int jmax_;

    //! Right-hand side coefficients up to a fine level, sorted by modulus.
    Array1D<std::pair<typename AggregatedFrame<IBASIS,DIM>::Index,double> > fcoeffs;

    //! Patchwise right-hand side coefficients on a fine level, sorted by modulus.
    Array1D<Array1D<std::pair<typename AggregatedFrame<IBASIS,DIM>::Index,double> > > fcoeffs_patch;

    //! Square root of coefficients on diagonal of stiffness matrix.
    Array1D<double> stiff_diagonal;

    //! (Squared) \f$\ell_2\f$ norm of the precomputed right-hand side.
    double fnorm_sqr;

    /*! 
      (Squared) \f$\ell_2\f$ norm of the respective precomputed right-hand side coefficients
      on each patch.
     */
    Array1D<double> fnorms_sqr_patch;

    // reminder: The keyword mutable can only be applied to non-static
    // and non-const data members of a class. If a data member is declared mutable,
    // then it is legal to assign a value to this data member from
    // a const member function.
    mutable double normA, normAinv;
  };
}

#include <simple_elliptic_equation.cpp>

#endif
