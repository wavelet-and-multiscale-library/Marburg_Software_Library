// -*- c++ -*-
// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2016-2019                                            |
// | Christoph Hartmann                                                      |
// +--------------------------------------------------------------------+
/*!
    Class for the computation, output und storage of the following data after each Kacanov iteration:

    double epsilon_CDD1;            // error bound for CDD1: 'CDD1_SOLVE' terminates if ||u - u_epsilon||_{l_2} < epsilon_CDD1
    double epsilon_stabilization;   // stabilization parameter
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
    double c1;
    double c2;
    double kappa;                   // (estimate of the) condition number of the stiffness matrix A
    double F;
    double delta_start;             // (upper bound of) the initial l2 approximation error ||u_exact - u_0||_l2
    int DOF;                        // Degrees of freedom (= number of nonzero wavelet coefficients)
    double error_l2;                // \ell_2 error of the (approximate) wavelet coefficients u_epsilon, i.e., ||u_exact - u_epsilon||_l2
    double error_l_infty;           // \ell_infty error of the (approximate) wavelet coefficients u_epsilon, i.e., ||u_exact - u_epsilon||_l_infty
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
    double rel_error_l2;                // relative \ell_2 error of the (approximate) wavelet coefficients u_epsilon, i.e., ||u_exact - u_epsilon||_l2 / ||u_exact||_l2
    double rel_error_L_infty;           // relative L_infty error, i.e., ess sup { |u_exact(x) - u_epsilon(x)| / |u_exact(x)| }
    SampledMapping<2,double> u_approx_sm;       // the approximate solution u_epsilon as a sampled mapping with grid resolution 'grid_resolution_sm'
    SampledMapping<2,double> abs_error_sm;      // the absolute error u_epsilon - u_exact as a sampled mapping with grid resolution 'grid_resolution_sm'
    SampledMapping<2,double> rel_error_sm;      // the relative error (u_epsilon - u_exact) / u_exact as a sampled mapping with grid resolution 'grid_resolution_sm'
    double abs_error_sm_max_nrm;                // L_infty error of abs_error_sm DUE TO PERFORMANCE REASONS: Faster than computation of 'error_L_infty' if 'grid_resolution_sm' << "max. wavelet level of u_epsilon"
    double rel_error_sm_max_nrm;                // L_infty error of rel_error_sm DUE TO PERFORMANCE REASONS: Faster than computation of 'rel_error_L_infty' if 'grid_resolution_sm' << "max. wavelet level of u_epsilon"

  */
//
//! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//!  Note: For the storage of the computed data (approximation errors, plots, parameter values, etc.), the subfolder structure './results/plots' must exist in the folder of the executable
//! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
#ifndef _WAVELETTL_P_POISSON_DATA_H
#define _WAVELETTL_P_POISSON_DATA_H

//! global variables for time measurement: (can be deleted after finished optimization)
extern double time_consumption_of_eval_u_eps;
extern double time_consumption_of_eval_deriv_u_eps;


namespace WaveletTL
{

template <class PROBLEM>
class pPoissonData
{
public:

  /*! constructor from a problem*/
  pPoissonData(const PROBLEM* prob, const Function<2, double>* exact_solution, const Function<2, double>* derivative_exact_solution);


  // the type of data that is saved at the end of one Kachanov iteration
  struct Data_structure
  {
    double epsilon_CDD1;            // error bound for CDD1: 'CDD1_SOLVE' terminates if ||u - u_epsilon||_{l_2} < epsilon_CDD1
    double epsilon_stabilization;   // stabilization parameter
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
    double c1;
    double c2;
    double kappa;                   // (estimate of the) condition number of the stiffness matrix A
    double F;
    double delta_start;             // (upper bound of) the initial l2 approximation error ||u_exact - u_0||_l2
    int DOF;                        // Degrees of freedom (= number of nonzero wavelet coefficients)
    double error_l2;                // \ell_2 error of the (approximate) wavelet coefficients u_epsilon, i.e., ||u_exact - u_epsilon||_l2
    double error_l_infty;           // \ell_infty error of the (approximate) wavelet coefficients u_epsilon, i.e., ||u_exact - u_epsilon||_l_infty
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
    double rel_error_l2;                // relative \ell_2 error of the (approximate) wavelet coefficients u_epsilon, i.e., ||u_exact - u_epsilon||_l2 / ||u_exact||_l2
    double rel_error_L_infty;           // relative L_infty error, i.e., ess sup { |u_exact(x) - u_epsilon(x)| / |u_exact(x)| }
    SampledMapping<2,double> u_approx_sm;       // the approximate solution u_epsilon as a sampled mapping with grid resolution 'grid_resolution_sm'
    SampledMapping<2,double> abs_error_sm;      // the absolute error u_epsilon - u_exact as a sampled mapping with grid resolution 'grid_resolution_sm'
    SampledMapping<2,double> rel_error_sm;      // the relative error (u_epsilon - u_exact) / u_exact as a sampled mapping with grid resolution 'grid_resolution_sm'
    double abs_error_sm_max_nrm;                // L_infty error of abs_error_sm DUE TO PERFORMANCE REASONS: Faster than computation of 'error_L_infty' if 'grid_resolution_sm' << "max. wavelet level of u_epsilon"
    double rel_error_sm_max_nrm;                // L_infty error of rel_error_sm DUE TO PERFORMANCE REASONS: Faster than computation of 'rel_error_L_infty' if 'grid_resolution_sm' << "max. wavelet level of u_epsilon"
  };

  int grid_resolution_sm;      // the grid resolution for sampled mappings, i.e., grid width = 2^-grid_resolution_sm

  // the actual data
  std::vector<Data_structure> Data;

  // computes all terms contained in 'Data_structure'. The values are stored in 'Data'.
  void compute_data(const InfiniteVector<double, typename PROBLEM::Index>& u_epsilon,
                    double param_p, int doe_quad, int min_res_quad, bool compute_constrained_energy,
                    double epsilon_lower, double epsilon_upper);

  /*! computes:  - the wavelet expansion of the exact solution (up to wavelet level 'max_lev_wav_u_exact'),
                 - the exact solution as a sampled mapping and as a matrix,
                 - a Matlab plot of the exact solution,
                 - the energy of the exact solution (for integration, a composite Gauss quadrature formula with degree of exactness = 'doe_quad' and
                   resolution = 'min_res_quad' is used, i.e., the sidelength of the patches is 2^-min_res_quad)
  */
  void precompute_exact_solution(const int max_lev_wav_u_exact, const int doe_quad, const int min_res_quad);

  // prints the saved data for iteration i, i.e., all data contained in 'Data[i]' (standard output)
  void print_data(const int i, std::ostream& ostr);

  // save some error plots */
  void save_error_plots(std::ostream&  ofs_plot_error_energy,
                        std::ostream&  ofs_plot_error_l2,
                        std::ostream&  ofs_plot_error_W1p,
                        std::ostream&  ofs_plot_error_W1_infty);

  // returns the number of saved Kachanov iterations
  int iterations()
  {
    return Data.size();
  }

double dummy_integrate_u(const typename PROBLEM::Index& lambda) const;

//private:

  const PROBLEM* P;     // pointer to the underlying problem

  const Function<2, double>* u_exact;           // pointer to the exact solution
  const Function<2, double>* u_exact_deriv;     // pointer to the derivative of the exact solution

  int current_iteration;

  SampledMapping<2, double> u_exact_sm;                  // the exact solution as a sampled mapping (with grid_resolution 'grid_resolution_sm')
  Matrix<double> u_exact_matrix;                         // the exact solution as a matrix (essentially a copy of 'u_exact_sm'), used to compute the relative pointwise error
  InfiniteVector<double, typename PROBLEM::Index> u_exact_wav_coeff;              // the exact solution as a wavelet expansion

#if 1
//! test
InfiniteVector<double, typename PROBLEM::Index> dummy_u_exact_wav_coeff;
#endif

  double u_exact_energy;

  //! test
  double u_exact_energy1;   // 1st term
  double u_exact_energy2;   // 2nd term

  // energy of the exact wavelet expansion of u
  double energy1_exact_wav;  //! test
  double energy2_exact_wav;  //! test
  double energy_exact_wav;  //! test
  //! test

};

}

#include <galerkin/p_poisson_data_cube.cpp>

#endif
