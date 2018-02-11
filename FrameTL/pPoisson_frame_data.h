// -*- c++ -*-
// +--------------------------------------------------------------------+
// | This file is part of FrameTL - the Frame Template Library      |
// |                                                                    |
// | Copyright (c) 2016-2019                                            |
// | Christoph Hartmann                                                      |
// +--------------------------------------------------------------------+
/*!
    Class for the computation, output und storage of the following data after each Kacanov iteration:

    double epsilon_MS;              // error bound for CDD1: 'CDD1_SOLVE' terminates if ||u - u_epsilon||_{l_2} < epsilon_CDD1
    double epsilon_stabilization;   // stabilization parameter in Kachanov iteration
    int max_lev_wav_CDD1;           // the maximal level of the wavelets considered by 'CDD1_SOLVE'
    int max_lev_wav_rhs;            // the maximal level of the wavelets considered by 'compute_rhs()' for precomputation of the rhs
    int max_lev_wav_u_exact;        // the maximal level of the wavelets considered for the expansion of the exact solution u_exact
    int min_res_quad_a;             // parameter for the evaluation of the bilinear form a(,). maximal side length (=2^-min_resolution_quadrature) of the cubes
                                    // used for the composite Gauss quadrature
    int doe_quad_a;                 // degree of exactness of the Gauss quadrature rule used in 'a' (Note: only odd doe possible. For even value, doe will be value-1 !!)
    int min_res_quad_f;             // parameter for the evaluation of the functional f(). maximal side length (=2^-min_resolution_quadrature_f) of the cubes
                                    // used for the composite Gauss quadrature. Set value > 0 if rhs is not very smooth.
    int doe_quad_f;                 // degree of exactness of the Gauss quadrature rule used in 'cube_equation_pPoisson::f' (Note: only odd doe possible. For even value, doe will be value-1 !!)
    double gamma_CDD1;
    double c1_patch0;
    double c2_patch0;
    double c1_patch1;
    double c2_patch1;
    double norm_A;
    double kappa;                   // (estimate of the) condition number of the stiffness matrix A
    double F;
    double delta_start;             // (upper bound of) the initial l2 approximation error ||u_exact - u_0||_l2
    int DOF;                        // Degrees of freedom (= number of nonzero wavelet coefficients)
    double energy;                  // the energy J(u_epsilon)
    double error_energy;            // the difference J(u_exact) - J(u_epsilon)
    double constrained_energy;      // the constrained energy J_epsilon(u_epsilon)
    double energy1;         //! test
    double energy2;         //! test
    double error_energy1;   //! test
    double error_energy2;   //! test
    double error_Lp;                // the L_p error ||u_exact - u_epsilon||_Lp
    double error_Lp_deriv[2];       // the L_p errors of the first order partial derivatives, i.e., ||d/dx u_exact - d/dx u_epsilon||_Lp and ||d/dy u_exact - d/dy u_epsilon||_Lp
    double error_W1p;               // the W^1_p error ||u_exact - u_epsilon||_{W^1_p} = ||u_exact - u_epsilon||_Lp + ||d/dx u_exact - d/dx u_epsilon||_Lp + ||d/dy u_exact - d/dy u_epsilon||_Lp
    double error_L_infty;           // the L_infty error ||u_exact - u_epsilon||_L_infty = ess sup |u_exact(x) - u_epsilon(x)|
    double error_L_infty_deriv[2];  // the L_infty errors of the first order partial derivatives, i.e., ||d/dx u_exact - d/dx u_epsilon||_L_infty and ||d/dy u_exact - d/dy u_epsilon||_L_infty
    double error_W1_infty;          // the W^1_infty error ||u_exact - u_epsilon||_{W^1_infty} = ||u_exact - u_epsilon||_L_infty + ||d/dx u_exact - d/dx u_epsilon||_L_infty + ||d/dy u_exact - d/dy u_epsilon||_L_infty
    double rel_error_L_infty;           // relative L_infty error, i.e., ess sup { |u_exact(x) - u_epsilon(x)| / |u_exact(x)| }
    Array1D<SampledMapping<2,double> > u_approx_sm;       // the approximate solution u_epsilon as a sampled mapping with grid resolution 'grid_resolution_sm'
    Array1D<SampledMapping<2,double> > abs_error_sm;      // the absolute error u_epsilon - u_exact as a sampled mapping with grid resolution 'grid_resolution_sm'
    SampledMapping<2,double> rel_error_sm;      // the relative error (u_epsilon - u_exact) / u_exact as a sampled mapping with grid resolution 'grid_resolution_sm'

  */
//
//! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//!  Note: For the storage of the computed data (approximation errors, plots, parameter values, etc.), the subfolder structure './results/plots' must exist in the folder of the executable
//! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//

#ifndef _WAVELETTL_P_POISSON_DATA_H
#define _WAVELETTL_P_POISSON_DATA_H


#ifndef _FRAMETL_AGGREGATED_FRAME_H
#include <aggregated_frame.h>
#endif

#ifndef _FRAMETL_FRAME_EVALUATE_H
#include <frame_evaluate.h>
#endif
#include <cube/cube_indexplot.h>

//! global variables for time measurement: (can be deleted after finished optimization)
extern double time_consumption_of_eval_u_eps;
extern double time_consumption_of_eval_deriv_u_eps;


namespace FrameTL
{

template <class IBASIS, class PROBLEM>
class pPoissonDataFrame
{
public:

  /*! constructor from a problem*/
  pPoissonDataFrame(const AggregatedFrame<IBASIS,2,2>* frame_param, const PROBLEM* prob,
                    const Function<2, double>* exact_solution, const Function<2, double>* derivative_exact_solution,
                    const int grid_resolution);


  //! the type of data that is saved at the end of each Kachanov iteration
  struct Data_structure
  {
    double epsilon_MS;              // error bound for CDD1: 'CDD1_SOLVE' terminates if ||u - u_epsilon||_{l_2} < epsilon_CDD1
    double epsilon_stabilization;   // stabilization parameter in Kachanov iteration
    int max_lev_wav_CDD1;           // the maximal level of the wavelets considered by 'CDD1_SOLVE'
    int max_lev_wav_rhs;            // the maximal level of the wavelets considered by 'compute_rhs()' for precomputation of the rhs
    int max_lev_wav_u_exact;        // the maximal level of the wavelets considered for the expansion of the exact solution u_exact
    int min_res_quad_a;             // parameter for the evaluation of the bilinear form a(,). maximal side length (=2^-min_resolution_quadrature) of the cubes
                                    // used for the composite Gauss quadrature
    int doe_quad_a;                 // degree of exactness of the Gauss quadrature rule used in 'a' (Note: only odd doe possible. For even value, doe will be value-1 !!)
    int min_res_quad_f;             // parameter for the evaluation of the functional f(). maximal side length (=2^-min_resolution_quadrature_f) of the cubes
                                    // used for the composite Gauss quadrature. Set value > 0 if rhs is not very smooth.
    int doe_quad_f;                 // degree of exactness of the Gauss quadrature rule used in 'cube_equation_pPoisson::f' (Note: only odd doe possible. For even value, doe will be value-1 !!)
    double gamma_CDD1;
    double c1_patch0;
    double c2_patch0;
    double c1_patch1;
    double c2_patch1;
    double norm_A;
    double kappa;                   // (estimate of the) condition number of the stiffness matrix A
    double F;
    double delta_start;             // (upper bound of) the initial l2 approximation error ||u_exact - u_0||_l2
    int DOF;                        // Degrees of freedom (= number of nonzero wavelet coefficients)
    double energy;                  // the energy J(u_epsilon)
    double error_energy;            // the difference J(u_exact) - J(u_epsilon)
    double constrained_energy;      // the constrained energy J_epsilon(u_epsilon)
    double energy1;         //! test
    double energy2;         //! test
    double error_energy1;   //! test
    double error_energy2;   //! test
    double error_Lp;                // the L_p error ||u_exact - u_epsilon||_Lp
    double error_Lp_deriv[2];       // the L_p errors of the first order partial derivatives, i.e., ||d/dx u_exact - d/dx u_epsilon||_Lp and ||d/dy u_exact - d/dy u_epsilon||_Lp
    double error_W1p;               // the W^1_p error ||u_exact - u_epsilon||_{W^1_p} = ||u_exact - u_epsilon||_Lp + ||d/dx u_exact - d/dx u_epsilon||_Lp + ||d/dy u_exact - d/dy u_epsilon||_Lp
    double error_L_infty;           // the L_infty error ||u_exact - u_epsilon||_L_infty = ess sup |u_exact(x) - u_epsilon(x)|
    double error_L_infty_deriv[2];  // the L_infty errors of the first order partial derivatives, i.e., ||d/dx u_exact - d/dx u_epsilon||_L_infty and ||d/dy u_exact - d/dy u_epsilon||_L_infty
    double error_W1_infty;          // the W^1_infty error ||u_exact - u_epsilon||_{W^1_infty} = ||u_exact - u_epsilon||_L_infty + ||d/dx u_exact - d/dx u_epsilon||_L_infty + ||d/dy u_exact - d/dy u_epsilon||_L_infty
    double rel_error_L_infty;           // relative L_infty error, i.e., ess sup { |u_exact(x) - u_epsilon(x)| / |u_exact(x)| }
    Array1D<SampledMapping<2,double> > u_approx_sm;       // the approximate solution u_epsilon as a sampled mapping with grid resolution 'grid_resolution_sm'
    Array1D<SampledMapping<2,double> > abs_error_sm;      // the absolute error u_epsilon - u_exact as a sampled mapping with grid resolution 'grid_resolution_sm'
    SampledMapping<2,double> rel_error_sm;      // the relative error (u_epsilon - u_exact) / u_exact as a sampled mapping with grid resolution 'grid_resolution_sm'
  };

  int grid_resolution_sm;      // the grid resolution for sampled mappings, i.e., grid width = 2^-grid_resolution_sm


  std::vector<Data_structure> Data;        //! the actual data


  /*! 'precompute_exact_solution':
      computes:  - the exact solution as a sampled mapping  (with resolution 'grid_resolution_sm')
                 - a Matlab plot of the exact solution,
                 - the energy of the exact solution (for integration, a composite Gauss quadrature formula with degree of exactness = 'doe_quad' and
                   resolution = 'res_quad' is used, i.e., the sidelength of the patches is 2^-res_quad)
  */
  void precompute_exact_solution(const int doe_quad, const int res_quad);


  //! 'compute_data':
  //! Computes the energy J(u_epsilon) (and energy-error) of the (global) approximation u_epsilon (given by the wavelet frame coefficients 'u_epsilon' = last element of 'approximations')
  //! to the p-Poisson equation for p='param_p', as well as the approximation error in various norms (L_p, L_inf, W^1_p ...).
  //!
  //! For integration, a composite gauss quadrature rule with degree of exactness 'doe_quad' is used.
  //! For the resolution of the composite gauss quadrature there are two strategies at choice:
  //!   1) If 'force_fixed_quad_res' == 0, then the resolution of the composite gauss quadrature is set to
  //!      the maximum of: parameter 'min_res_quad' and 'highest wavelet level of u_epsilon' + 1.
  //!      (u_epsilon is piecewise smooth on dyadic squares with sidelength 2^-j, j='highest wavelet level of u_epsilon' + 1)
  //!   2) If 'force_fixed_quad_res' > 0, then the resolution of the composite gauss quadrature is set to this value.
  //!
  //! If 'compute_constrained_energy' == true, also the constrained energy J_epsilon(u_epsilon) is computed and stored in 'constrained energy'.
  //! The values are stored in 'Data'.
  void compute_data(const Array1D<InfiniteVector<double, typename PROBLEM::Index> >& approximations,
                    double param_p, int doe_quad, int min_res_quad, bool compute_constrained_energy,
                    double epsilon_lower, double epsilon_upper, int force_fixed_quad_res);


  //! 'print_data':
  //!  prints the saved data for iteration i, i.e., all data contained in 'Data[i]' (standard output)
  void print_data(const int i, std::ostream& ostr);


  //! 'save_error_plots':
  //! Saves plots (Matlab) of energy, W^1_p and W^1_infty error. x-axis contains the iteration, y-axis the corresponding error */
  void save_error_plots(std::ostream&  ofs_plot_error_energy,
                        std::ostream&  ofs_plot_error_W1p,
                        std::ostream&  ofs_plot_error_W1_infty);


  //! returns the number of saved Kachanov iterations
  int iterations()
  {
    return Data.size();
  }

double dummy_integrate_u(const typename PROBLEM::Index& lambda) const;

//private:
  const AggregatedFrame<IBASIS,2>* frame;
  const PROBLEM* P;     // pointer to the underlying problem
  EvaluateFrame<IBASIS, 2, 2> _frame_eval;

  const Function<2, double>* u_exact;           // pointer to the exact solution
  const Function<2, double>* u_exact_deriv;     // pointer to the derivative of the exact solution

  int current_iteration;

  Array1D<SampledMapping<2> > u_exact_sm;               // the exact solution as a sampled mapping (with grid_resolution 'grid_resolution_sm')

  double u_exact_energy;    // energy of u, i.e., \int_{\Omega} 1/p \lvert \nabla u \rvert^p dx - \int_{\Omega} f u dx

  double u_exact_energy1;   // 1st term of energy, i.e., \int_{\Omega} 1/p \lvert \nabla u \rvert^p dx
  double u_exact_energy2;   // 2nd term of energy, i.e., \int_{\Omega} f u dx

};

}

#include <pPoisson_frame_data.cpp>

#endif
