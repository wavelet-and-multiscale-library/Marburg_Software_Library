// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of FrameTL - the Frame Template Library          |
// |                                                                    |
// | Copyright (c) 2016-2019                                            |
// | Christoph Hartmann                                                 |
// +--------------------------------------------------------------------+

#ifndef _FRAMETL_PPOISSON_EQUATION_H
#define _FRAMETL_PPOISSON_EQUATION_H

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

//! global variables for time measurement: (can be deleted after finished optimization)
extern double time_consumption_of_compute_diagonal;
extern double time_consumption_of_compute_rhs;
extern double time_consumption_of_a_same_patches;
extern double time_consumption_of_a_different_patches;


namespace FrameTL
{

  /*!
    This class models the (preconditioned) infinite-dimensional matrix problem (i.e., in particular the p-Poisson equation)

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

    @tparam IBASIS The type of interval basis underlying the construction of the aggregated frame.
    @tparam DIM The dimension of the underlying domain.
  */
  template <class IBASIS, unsigned int DIM>
  class pPoissonEquation
  //: public FullyDiagonalDyadicPreconditioner<typename AggregatedFrame<IBASIS,DIM>::Index>
      : public FullyDiagonalEnergyNormPreconditioner<typename AggregatedFrame<IBASIS,DIM>::Index>
  {
  public:

    /*!
      Constructor. The diagonal of the stiffness matrix and the coefficients of the right-hand side
      are precomputed between minimal and maximal level.

      @param ell_bvp The elliptic boundary value problem that is modeled.
      @param frame Pointer to the aggragated frame that is used for discretization.
      @param max_lev_wav_basis The maximal level of resolution that is considered.
     */
    pPoissonEquation(const PoissonBVP_Coeff<DIM>* ell_bvp,
		             const AggregatedFrame<IBASIS,DIM>* frame,
		             const double parameter_p,
                     const int max_lev_wav_rhs = 4,
                     const int doe_quad_a = 0,
                     const int min_res_quad_a = 0,
                     const int doe_quad_f = 7,
                     const int min_res_quad_f = 0);

    /*!
      The frame type.
     */
    typedef AggregatedFrame<IBASIS,DIM> Frame;

    /*!
      Intverval basis type.
    */
    typedef IBASIS IntervalBasis;

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
    const PoissonBVP_Coeff<DIM>&  get_bvp() const { return *ell_bvp_; }

    /*!
      Read access to the frame. The routine is called basis() to be
      compatible with the routines in WaveletTL's compression.h.
    */
    const AggregatedFrame<IBASIS,DIM>& basis() const { return *frame_; }

    /*!
      maximal level of the underlying wavelet basis
    */
    int jmax() const { return max_level_wavelet_basis; }

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
      the parameter p of the p-Poisson equation
    */
    double param_p;

    /*!
      the maximal level of the wavelet basis
    */
    int max_level_wavelet_basis;

    /*!
      the maximal level of the wavelets considered by the member 'compute_rhs()'
    */
    int max_level_wavelets_rhs;


    /*!
      parameter for the evaluation of the bilinear form a(,), i.e. for the member 'a'.
      degree of exactness of the Gauss quadrature rule used in 'a'
    */
    int degree_of_exactness_quadrature_a;

    /*!
      parameter for the evaluation of the bilinear form a(,), i.e. for the member 'a'.
      maximal side length (=2^-min_res_quadrature_a) of the cubes
      used for the composite Gauss quadrature, i.e., all cubes will have
      side lengths less than or equal to 2^-min_res_quadrature_a
    */
    int min_res_quadrature_a;

    /*!
      parameter for the evaluation of the functional f, i.e. for the member 'f'.
      degree of exactness of the Gauss quadrature rule used in 'f'
    */
    int degree_of_exactness_quadrature_f;

    /*!
      parameter for the evaluation of the functional f, i.e. for the member 'f'.
      maximal side length (=2^-min_res_quadrature_f) of the cubes
      used for the composite Gauss quadrature, i.e., all cubes will have
      side lengths less than or equal to 2^-min_res_quadrature_f
    */
    int min_res_quadrature_f;


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
    void set_normA(const double _normA) const { normA = _normA; }

    /*!
      Sets estimate for \f$\|A^{-1}\|\f$.
    */
    void set_normAinv(const double nAinv) const { normAinv = nAinv; };

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
    void RHS(const double eta,
             InfiniteVector<double, typename AggregatedFrame<IBASIS,DIM>::Index>& coeffs) const;

    /*!
      Approximate the wavelet coefficient set of the preconditioned right-hand side restricted
      to patch p within a prescribed \f$\ell_2\f$ error tolerance.
    */
    void RHS(const double eta,
             const int p,
             InfiniteVector<double, typename AggregatedFrame<IBASIS,DIM>::Index>& coeffs) const;

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
    void set_bvp(const PoissonBVP_Coeff<DIM>*);


    /*!
      Multiplies the stiffness matrix entries of column lambda on level j of the compressed martrix A_J
      by factor and adds the result to w.
    */
    void add_level (const Index& lambda,
		    InfiniteVector<double, Index>& w, const int j,
		    const double factor,
		    const int J,
		    const CompressionStrategy strategy) const;


    /*!
      (cached!) point evaluation of the coefficient function bvp_->a
    */
    double cached_coeff_a(const int global_patch, const int current_patch, const int gauss_point_in_patch, const int j, const int N_Gauss) const;


    /*!
      clear cache 'coeff_cache'
    */
    void clear_coeff_cache()
    {
      coeff_cache.clear();
    }




//! Right-hand side coefficients up to a fine level, sorted by modulus.
    Array1D<std::pair<typename AggregatedFrame<IBASIS,DIM>::Index,double> > fcoeffs;


   protected:

    /*!
      The elliptic boundary value problem.
     */
    const PoissonBVP_Coeff<DIM>* ell_bvp_;
    //const EllipticBVP<DIM>* ell_bvp_;

    /*!
      The underlying aggregated frame.
    */
    const AggregatedFrame<IBASIS,DIM>* frame_;

   private:

    /*!
      Helper routine for a (...,...). Entries in diagonal and non-diagonal
      blocks of the stiffness matrix have to be treated differently.
      This routine is responsible for entries in the diagonal blocks.

      @param n_Gauss_knots The number of Gauss knots used in the Gauss quadrature rule.
     */
    double a_same_patches(const typename AggregatedFrame<IBASIS,DIM>::Index& lambda,
			  const typename AggregatedFrame<IBASIS,DIM>::Index& nu,
			  const unsigned int n_Gauss_knots = 3) const;

    /*!
      Helper routine for a (...,...). Entries in diagonal and non-diagonal
      blocks of the stiffness matrix have to be treated differently.
      This routine is responsible for entries in the non-diagonal blocks.
      We perform a composite Gaussian qudrature rule of fixed order and rank.

      @param n_Gauss_knots The number of Gauss knots used in the Gauss quadrature rule.
     */
    double a_different_patches(const typename AggregatedFrame<IBASIS,DIM>::Index& lambda,
			       const typename AggregatedFrame<IBASIS,DIM>::Index& nu,
			       const unsigned int n_Gauss_knots = 3, const unsigned int rank = 1) const;



//! +++++++++++ tunded version of a_same_patches which caches the point evaluations
//! +++++++++++ of the coefficient function bvp_->a
//! +++++++++++ ONLY for L-shaped domain !!!

    double a_same_patches2(const typename AggregatedFrame<IBASIS,DIM>::Index& lambda,
			  const typename AggregatedFrame<IBASIS,DIM>::Index& nu,
			  const unsigned int doe_quadrature = 5, const unsigned int min_resolution = 0) const;

//! +++++++++++ tunded version of a_different_patches which caches the point evaluations
//! +++++++++++ of the coefficient function bvp_->a
//! +++++++++++ ONLY for L-shaped domain !!! 'rank' must be a power of 2 !!!

    double a_different_patches2(const typename AggregatedFrame<IBASIS,DIM>::Index& lambda,
			       const typename AggregatedFrame<IBASIS,DIM>::Index& nu,
			       const unsigned int doe_quadrature = 5, const unsigned int min_resolution = 0,
			       const unsigned int rank = 1) const;


    //! Precompute the right-hand side between minimal and maximal level.
    void compute_rhs();

    //! Precompute the diagonal of the stiffness matrix between minimal and maximal level.
    void compute_diagonal();



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


     /*!
      cache for coefficient function bvp_->a
    */

    // type of one Func_values.
    // the key codes the gauss point inside the patch
    // the data are the point evaluations of the coefficient function bvp_->a at the gauss points
    typedef std::map<int, double> Func_values;

    // type of one Patch.
    // the key codes the number of Gauss points in one dimension
    // the data are the Func_values
    typedef std::map<int, Func_values> Patch;

    // type of one Partition.
    // the key codes the position, that data are the Patches
    typedef std::map<int, Patch> Partition;

    // the key codes the resolution, the data are the Partitions,
    // i.e., the corresponding partition has 2^(Resolution*DIM) Patches
    typedef std::map<int, Partition> Resolution;

    // the key codes the global patch (of the frame), the data are the Resolutions
    typedef std::map<int, Resolution> Global_Patch;

    // entries cache for bvp_->a
    mutable Global_Patch coeff_cache;



  };
}

#include <pPoisson_equation.cpp>

#endif
