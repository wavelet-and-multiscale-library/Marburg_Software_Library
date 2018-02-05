// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_CUBE_EQUATION_PPOISSON_H
#define _WAVELETTL_CUBE_EQUATION_PPOISSON_H

#include <set>
#include <utils/fixed_array1d.h>
#include <utils/array1d.h>
#include <numerics/bvp.h>

#include <galerkin/galerkin_utils.h>
#include <galerkin/infinite_preconditioner.h>

using MathTL::FixedArray1D;
using MathTL::EllipticBVP;
extern double time_consumption_of_a;

namespace WaveletTL
{
  template <class IBASIS, unsigned int DIM> class CubeBasis;

  /*!
    This class models the (preconditioned) infinite-dimensional matrix problem (tailored for the p-Poisson operator!)

      Au = D^{-1}LD^{-1}u = D^{-1}F

    when reformulating an elliptic boundary value problem on the cube [0,1]^d

      -div(a(x)grad u(x)) + q(x)u(x) = f(x)

    with first (Dirichlet) or second (Neumann) order b.c.'s (modeled in
    the class EllipticBVP) as an equivalent operator equation
    within \ell_2 by means of a wavelet basis.

    The corresponding bilinear form in

      L = (a(\psi_\nu,\psi_\lambda))_{\lambda,\nu}

    is

      a(u,v) = \int_Omega [a(x)grad u(x)grad v(x)+q(x)u(x)v(x)] dx

    and the right-hand side is

      F(v) = \int_0^1 f(x)v(x) dx

    The evaluation of a(.,.) and f is possible for arguments \psi_\lambda
    which stem from a wavelet basis \Psi=\{\psi_\lambda\} of the corresponding
    function space over Omega.
    To achieve independence from the concrete choice of \Psi, the wavelet basis
    class is given as a template parameter CUBEBASIS. It should have a constructor of
    the form

       CUBEBASIS::CUBEBASIS(const FixedArray1D<bool,2*DIM>& bc);

    where bc indicates the enforcement of homogeneous Dirichlet boundary conditions (true).
    A natural concrete value for CUBEBASIS is the CubeBasis<DSBasis<d,dT> >.
  */
  template <class IBASIS, unsigned int DIM, class CUBEBASIS = CubeBasis<IBASIS,DIM> >
  class CubeEquationpPoisson
     : public FullyDiagonalDyadicPreconditioner<typename CUBEBASIS::Index>
  //    : public FullyDiagonalEnergyNormPreconditioner<typename CUBEBASIS::Index>
  {
  public:

    /*!
      constructor from a boundary value problem and specified b.c.'s
    */
    CubeEquationpPoisson(const PoissonBVP_Coeff<DIM>* bvp,
                         const FixedArray1D<bool,2*DIM>& bc,
                         const double parameter_p,
                         const int max_lev_wav_basis = 4,
                         const int max_lev_wav_rhs = 4,
                         const int doe_quad_a = 0,
                         const int min_res_quad_a = 0,
                         const int doe_quad_f = 9,
                         const int min_res_quad_f = 0);

    /*!
      constructor from a boundary value problem and specified b.c.'s
    */
    CubeEquationpPoisson(const PoissonBVP_Coeff<DIM>* bvp,
                         const FixedArray1D<int,2*DIM>& bc,
                         const double parameter_p,
                         const int max_lev_wav_basis = 4,
                         const int max_lev_wav_rhs = 4,
                         const int doe_quad_a = 0,
                         const int min_res_quad_a = 0,
                         const int doe_quad_f = 9,
                         const int min_res_quad_f = 0);

    /*!
      copy constructor
    */
    CubeEquationpPoisson(const CubeEquationpPoisson&);

    /*!
      make template argument accessible
    */
    typedef CUBEBASIS WaveletBasis;

    /*!
      wavelet index class
    */
    typedef typename WaveletBasis::Index Index;

    /*!
      read access to the basis
    */
    const CUBEBASIS& basis() const { return basis_; }

    /*!
      space dimension of the problem
    */
    static const int space_dimension = DIM;

    /*!
      differential operators are local
    */
    static bool local_operator() { return true; }

    /*!
      (half) order t of the operator
      (inherited from FullyDiagonalEnergyNormPreconditioner)
    */
    double operator_order() const { return 1.; }

    /*!
      the parameter p of the p-Poisson equation
    */
    double param_p;

    /*!
      the maximal level of the wavelet basis (CUBEBASIS)
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
      evaluate the diagonal preconditioner D
    */
    double D(const typename WaveletBasis::Index& lambda) const;

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
      Internally, we use an m-point composite tensor product Gauss rule adapted
      to the singular supports of the spline wavelets involved,
      so that m = (p+1)/2;
      you can specify the maximal side length (=2^-min_resolution) of the cubes
      used for the composite Gauss quadrature, i.e., all cubes will have
      side lengths less than or equal to 2^-min_resolution.
      This version caches the point evaluations of the coefficient function bvp_->a
    */
    double a(const typename WaveletBasis::Index& lambda,
	     const typename WaveletBasis::Index& nu,
	     const unsigned int p, const unsigned int min_resolution = 0) const;

    /*!
      The uncached version of a (Christoph: renamed it to 'a2')
    */
    double a2(const typename WaveletBasis::Index& lambda,
	     const typename WaveletBasis::Index& nu,
	     const unsigned int p, const unsigned int min_resolution = 0) const;

    /*!
      estimate the spectral norm ||A||
    */
    double norm_A() const;

    /*!
      estimate the spectral norm ||A^{-1}||
    */
    double norm_Ainv() const;

    /*!
      estimate compressibility exponent s^*
      (we assume that the coefficients a(x),q(x) are smooth)
    */
    double s_star() const;

    /*!
      estimate the compression constants alpha_k in
        ||A-A_k|| <= alpha_k * 2^{-s*k}
    */
    double alphak(const unsigned int k) const {
      return 2*norm_A(); // suboptimal
    }

    /*!
      evaluate the (unpreconditioned) right-hand side f, i.e.: f(lambda) = \int_{\Omega} f(x) * Psi_{\lambda}(x) dx
    */
    double f(const typename WaveletBasis::Index& lambda) const;



    /*!
      approximate the wavelet coefficient set of the preconditioned right-hand side F
      within a prescribed \ell_2 error tolerance
    */
    void RHS(const double eta,
	     InfiniteVector<double,typename WaveletBasis::Index>& coeffs) const;

    /*!
      compute (or estimate) ||F||_2
    */
    double F_norm() const { return sqrt(fnorm_sqr); }

    /*!
      set the boundary value problem
    */
    void set_bvp(const PoissonBVP_Coeff<DIM>*);

    void set_normA(const double norm_A_new);
    void set_normAinv(const double norm_Ainv_new);

    /*!
      (cached!) point evaluation of the coefficient function bvp_->a
    */
    double cached_coeff_a(const int current_patch, const int gauss_point_in_patch, const int j, const int N_Gauss) const;

    /*!
      clear cache 'coeff_cache'
    */
    void clear_coeff_cache()
    {
      coeff_cache.clear();
    }


  //protected:
    //const EllipticBVP<DIM>* bvp_;
    const PoissonBVP_Coeff<DIM>* bvp_;
    CUBEBASIS basis_;

    // right-hand side coefficients on a fine level, sorted by modulus
    Array1D<std::pair<typename WaveletBasis::Index,double> > fcoeffs;

    // precompute the right-hand side
    void compute_rhs();

    // (squared) \ell_2 norm of the precomputed right-hand side
    double fnorm_sqr;

    // estimates for ||A|| and ||A^{-1}||
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

    // entries cache for bvp_->a
    mutable Resolution coeff_cache;


  };
}

#include <galerkin/cube_equation_pPoisson.cpp>

#endif
