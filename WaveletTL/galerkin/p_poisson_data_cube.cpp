// -*- c++ -*-
// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2016-2019                                            |
// | Christoph Hartmann                                                      |
// +--------------------------------------------------------------------+
//
// implementation for p_poisson_data_cube.h
namespace WaveletTL
{

/*! constructor */
template <class PROBLEM>
pPoissonData<PROBLEM>::pPoissonData(const PROBLEM* prob, const Function<2, double>* exact_solution, const Function<2, double>* derivative_exact_solution)
  : grid_resolution_sm(7), P(prob), u_exact(exact_solution), u_exact_deriv(derivative_exact_solution), current_iteration(-1), u_exact_energy(0.0)
{
}

/*! computes:   - the wavelet expansion of the exact solution (up to wavelet level 'max_lev_wav_u_exact'),
                - the exact solution as a sampled mapping and as a matrix,
                - a Matlab plot of the exact solution,
                - the energy of the exact solution
*/
template <class PROBLEM>
void
pPoissonData<PROBLEM>::precompute_exact_solution(const int max_lev_wav_u_exact, const int doe_quad, const int min_res_quad)
{
  //! exact solution as wavelet expansion

  cout << "'precompute_exact_solution': -> expand u_exact in primal wavelet basis up to level " << max_lev_wav_u_exact << " .." << endl;
  P->basis().expand(u_exact, true, max_lev_wav_u_exact, u_exact_wav_coeff);       // compute wavelet coefficients of function u_exact. 'expand' uses
                                                                                  // composite Gauss quadrature with N_Gauss=5 points, i.e., degree of exactness is 9.
  cout << "'precompute_exact_solution': .. finished!" << endl;
 // P->basis().expand(u_exact, false, max_lev_wav_u_exact, dummy_u_exact_wav_coeff);


#if 0
Point<2> y;
y[0] = 0.501;
y[1] = 0.501;

double value_u_exact_in_y = u_exact->value(y);
double value_wav_expansion_in_y = evaluate_tuned_dim2_4(P->basis(), u_exact_wav_coeff, y);
double dummy_value_wav_expansion_in_y = evaluate_tuned_dim2_4(P->basis(), dummy_u_exact_wav_coeff, y);


cout << std::setprecision(10);
cout << endl << "u_exact_wav_coeff:" << endl << u_exact_wav_coeff << endl << endl;

cout << endl << "dummy_u_exact_wav_coeff:" << endl << dummy_u_exact_wav_coeff << endl << endl;

cout << "length of u_exact_wav_coeff: " << u_exact_wav_coeff.size() << endl << endl;
cout << "length of dummy_u_exact_wav_coeff: " << dummy_u_exact_wav_coeff.size() << endl << endl;

cout << "u_exact(x) = " << value_u_exact_in_y << endl;
cout << "wavlet_expansion(x) = " << value_wav_expansion_in_y << endl;
cout << "abs_error(x) = " << abs(value_u_exact_in_y - value_wav_expansion_in_y) << endl << endl;

cout << "dummy_wavlet_expansion(x) = " << dummy_value_wav_expansion_in_y << endl;
cout << "abs_error(x) = " << abs(value_u_exact_in_y - dummy_value_wav_expansion_in_y) << endl << endl;

cout << std::setprecision(6);

exit(1);
#endif

#if 0
//! test

typedef typename PROBLEM::Index Index;

const int j0 = P->basis().j0();
double coeff;

for (Index lambda(P->basis().first_generator(j0));; ++lambda)
{

  coeff = dummy_integrate_u(lambda);

  if (fabs(coeff)>1e-15)
  {
	dummy_u_exact_wav_coeff.set_coefficient(lambda, coeff);
  }

  if (lambda == P->basis().last_wavelet(max_lev_wav_u_exact))
  {
    break;
  }
}

//! test
#endif


  //! exact solution as sampled mapping
  Point<2, double> a(0.,0.), b(1., 1.);
  int number_gridpoints( 1 << grid_resolution_sm);
  Grid<2> my_grid(a,b, number_gridpoints);
  SampledMapping<2> dummy_sm(my_grid, *u_exact);
  u_exact_sm = dummy_sm;

  //! exact solution as Matlab plot
  std::ofstream u_exact_stream("./results/plots/u_exact.m");
  u_exact_sm.matlab_output(u_exact_stream);
  u_exact_stream.close();

  //! exact solution as matrix
  Matrix<double> dummy_matrix(u_exact_sm.values());
  u_exact_matrix = dummy_matrix;


  //! compute energy of exact solution
  Vector<double> deriv_values_u_exact;                      // to store one point evaluation of nabla u_exact
  deriv_values_u_exact.resize(2);
  double u_exact_in_x;
  double f_in_x;

  int j = min_res_quad;
  const int dim = 2;                    // space dimension
  double weights;                       // the current Gauss quadrature weight
  double param_p = P->get_problem()->param_p;

  // setup Gauss points and weights for a composite quadrature formula:
  const int N_Gauss = (doe_quad+1)/2;
  const double h = ldexp(1.0, -j);        // granularity for the quadrature (h=2^-j)
  const int total_patches = (1 << j);
  FixedArray1D<Array1D<double>,dim> gauss_points, gauss_weights;

  for (unsigned int i = 0; i < dim; i++)
  {
    gauss_points[i].resize(N_Gauss*total_patches);
	gauss_weights[i].resize(N_Gauss*total_patches);
	for (int patch = 0; patch < total_patches; patch++)
	{
	  for (int n = 0; n < N_Gauss; n++)
	  {
	    gauss_points[i][patch*N_Gauss+n] = h*(2*patch+1+GaussPoints[N_Gauss-1][n])/2.;      //! Attention: GaussPoints refer to intervall [-1,1] !!
	    gauss_weights[i][patch*N_Gauss+n] = h*GaussWeights[N_Gauss-1][n];                   //! Attention: GaussWeights refer to intervall [0,1] !!
	  }
    }
  }

#if 1
 // evaluate u_exact_wav_coeff on entire grid of gauss points on [0,1]^2. 'evaluate' makes use of the tensor product structure of the cube basis.
  Matrix<double> evaluation_u_exact_wav_coeff_grid =  evaluate(P->basis(), u_exact_wav_coeff, gauss_points);

  // evaluate first derivatives u_epsilon on entire grid of gauss points on [0,1]^2. 'evaluate' makes use of the tensor product structure of the cube basis.
  Matrix<double> deriv_x_u_exact_wav_coeff_grid;
  Matrix<double> deriv_y_u_exact_wav_coeff_grid;
  evaluate(P->basis(), u_exact_wav_coeff, gauss_points, deriv_x_u_exact_wav_coeff_grid, deriv_y_u_exact_wav_coeff_grid);
  energy1_exact_wav = 0.0;
  energy2_exact_wav = 0.0;
  FixedArray1D<double,2> deriv_values_u_exact_wav_coeff;
#endif

   // iterate over all points and sum up the integral shares
  int index[dim];                               // current multiindex for the Gauss points/weights
  for (unsigned int i = 0; i < dim; i++)
  {
	index[i] = 0;
  }

  Point<dim> x;

  double integral_f_u = 0.0;

  while (true)
  {
    // the current Gauss point x
    for (unsigned int i = 0; i < dim; i++)
	{
	  x[i] = gauss_points[i][index[i]];
    }
	// product of current Gauss weights
	weights = 1.0;
	for (unsigned int i = 0; i < dim; i++)
	{
	  weights *= gauss_weights[i][index[i]];
    }

	// compute the shares

    // evaluate u_exact (given as Function) in x
    u_exact_in_x =  u_exact->value(x);

    // evaluate f (given as Function) in x
    f_in_x = (P->get_problem())->bvp_->f(x);

    integral_f_u += (weights * u_exact_in_x * f_in_x);

    // evaluate first derivative of u_exact (given as Function) in x
    u_exact_deriv->vector_value(x, deriv_values_u_exact);

#if 1
    // evaluate function value and derivative of the exact wavelet expansion u_exact_wav_coeff
    double value_u_exact_wav_coeff_in_x = evaluation_u_exact_wav_coeff_grid(index[0], index[1]);
    deriv_values_u_exact_wav_coeff[0] = deriv_x_u_exact_wav_coeff_grid(index[0], index[1]);     // fastest version. the above call of 'evaluate' makes use of tensor product structure of the cube basis.
    deriv_values_u_exact_wav_coeff[1] = deriv_y_u_exact_wav_coeff_grid(index[0], index[1]);
    double dummy_share_ewc = 0.0;
	for (unsigned int i = 0; i < dim; i++)
	{
	  dummy_share_ewc += deriv_values_u_exact_wav_coeff[i]*deriv_values_u_exact_wav_coeff[i];
    }

    // compute share for energy
    double share_ewc = ( pow(dummy_share_ewc, param_p / 2.0) / param_p );
 	energy1_exact_wav += (weights * share_ewc);

    energy2_exact_wav += (weights * value_u_exact_wav_coeff_in_x * f_in_x);
#endif

    // compute |nabla u_exact|^2
	double share_u_exact = 0.0;
	for (unsigned int i = 0; i < dim; i++)
	{
	  share_u_exact += deriv_values_u_exact[i]*deriv_values_u_exact[i];
    }

 	// compute share for energy of u_exact
 	double share_energy_u_exact = ( pow(share_u_exact, param_p / 2.0) / param_p );
 	u_exact_energy += (weights * share_energy_u_exact);

 	// "++index"
	bool exit = false;
    for (unsigned int i = 0; i < dim; i++)
    {
      if (index[i] == N_Gauss*total_patches -1)
      {
	    index[i] = 0;
	    exit = (i == dim-1);
      }
      else
      {
	    index[i]++;
	    break;
      }
    }
    if (exit) break;
  }

  //! test
  u_exact_energy1 = u_exact_energy;
  //! test

  InfiniteVector<double, typename PROBLEM::Index> dummy;

  // copy P.fcoeffs (an Array1D) into an InfiniteVector (P.fcoeffs contains preconditioned wavelet coefficients of f!)
  for (typename Array1D<std::pair<typename PROBLEM::Index,double> >::const_iterator it((P->get_problem())->fcoeffs.begin()), itend((P->get_problem())->fcoeffs.end()); it != itend; ++it)
  {
    dummy.set_coefficient(it->first, it->second);
  }


   // compute <f,u_exact> = sum_{lambda \in u_epsilon} u_epsilon_lambda * f_lambda
  double r = 0.0;
#if 0
  for (typename InfiniteVector<double, typename PROBLEM::Index>::const_iterator it(u_exact_wav_coeff.begin()), itend(u_exact_wav_coeff.end()); it != itend; ++it)
  {
    typename PROBLEM::Index lambda;
    lambda = it.index();
    r += ( *it * dummy[lambda] * P->D(lambda) );
  }
#endif

#if 1

  r = integral_f_u;

#endif

  u_exact_energy -= r;

  //! test
  u_exact_energy2 = r;
  //! test

#if 1
energy_exact_wav = energy1_exact_wav - energy2_exact_wav;
#endif

}


//! computes the energy J(u_epsilon) of the approximation u_epsilon (given by the wavelet coefficients 'u_epsilon') to the p-Poisson equation for p='param_p'.
//! For integration, a composite gauss quadrature rule with degree of exactness 'doe_quad' is used.
//! If 'compute_constrained_energy' == true, also the constrained energy J_epsilon(u_epsilon) is computed and stored in 'constrained energy'.
template <class PROBLEM>
void
pPoissonData<PROBLEM>::compute_data(const InfiniteVector<double, typename PROBLEM::Index>& u_epsilon,
                                double param_p, int doe_quad, int min_res_quad, bool compute_constrained_energy,
                                double epsilon_lower, double epsilon_upper)
{
  current_iteration++;
  Data_structure current_data;         // storage for computed values

  double energy = 0.0;                 // the energy J(u_epsilon)
  double constrained_energy = 0.0;     // the constrained energy J_epsilon(u_epsilon)

  double error_Lp = 0.0;               // the error ||u_exact - u_epsilon||_Lp
  double error_Lp_deriv[2];            // the errors ||d/dx u_exact - d/dx u_epsilon||_Lp and ||d/dy u_exact - d/dy u_epsilon||_Lp
  error_Lp_deriv[0] = 0.0;
  error_Lp_deriv[1] = 0.0;

  double error_L_infty = 0.0;           // the error ||u_exact - u_epsilon||_L_infty
  double error_L_infty_deriv[2];        // the errors ||d/dx u_exact - d/dx u_epsilon||_L_infty and ||d/dy u_exact - d/dy u_epsilon||_L_infty
  error_L_infty_deriv[0] = 0.0;
  error_L_infty_deriv[1] = 0.0;

  double rel_error_L_infty = 0.0;       // the relative error ||(u_exact - u_epsilon) / u_exact||_L_infty

  int j;                                // resolution of the composite Gauss quadrature
  const int dim = 2;                    // space dimension

  double weights;                       // the current Gauss quadrature weight

  //! compute and save the values of c1, c2, kappa, F, delta_start and DOF
  current_data.c1 = 1.0/P->norm_Ainv();
  current_data.c2 = P->norm_A();
  current_data.kappa = current_data.c2 / current_data.c1;
  current_data.F = P->F_norm();
  current_data.delta_start = current_data.F / current_data.c1;
  current_data.DOF = u_epsilon.size();


  //! compute and save the l_2 error of the wavelet coefficients u_epsilon
  InfiniteVector<double, typename PROBLEM::Index> error_wav_coeff;
  error_wav_coeff.clear();
  error_wav_coeff += u_exact_wav_coeff;
  error_wav_coeff -= u_epsilon;
  current_data.error_l2 = l2_norm(error_wav_coeff);

  //! compute and save the l_infty error of the wavelet coefficients u_epsilon
  current_data.error_l_infty = linfty_norm(error_wav_coeff);


  //! compute energy, constrained energy, L_p errors and L_infty errors
  FixedArray1D<double,2> deriv_values_u_epsilon;            // to store one point evaluation of nabla u_epsilon
  Vector<double> deriv_values_u_exact;                      // to store one point evaluation of nabla u_exact
  deriv_values_u_exact.resize(2);


  // set minimal resolution of composite gauss quadrature to 'highest wavelet level of u_epsilon' + 1   (u_epsilon is piecewise smoothness on dyadic squares with sidelength 2^-j)
  typename std::set<typename PROBLEM::Index> supp;
  u_epsilon.support(supp);
  typename std::set<typename PROBLEM::Index>::const_iterator it = supp.end();


  if (supp.empty())
  {
    //j = frame->j0() + 1;
    j = 3 + 1;
  }

  else
  {
    --it;
    typename PROBLEM::Index dummy_index = *it;
    j = dummy_index.j() + 1;
  }


  if (j < min_res_quad)
  {
    j = min_res_quad;
  }

#if 0
//! test
j = 5;
//! test
#endif

  cout << "'compute_energy': min_res_quad = " << j << endl;

  // setup Gauss points and weights for a composite quadrature formula:
  const int N_Gauss = (doe_quad+1)/2;
  const double h = ldexp(1.0, -j);        // granularity for the quadrature (h=2^-j)
  const int total_patches = (1 << j);
  FixedArray1D<Array1D<double>,dim> gauss_points, gauss_weights;

  for (unsigned int i = 0; i < dim; i++)
  {
    gauss_points[i].resize(N_Gauss*total_patches);
	gauss_weights[i].resize(N_Gauss*total_patches);
	for (int patch = 0; patch < total_patches; patch++)
	{
	  for (int n = 0; n < N_Gauss; n++)
	  {
	    gauss_points[i][patch*N_Gauss+n] = h*(2*patch+1+GaussPoints[N_Gauss-1][n])/2.;      //! Attention: GaussPoints refer to intervall [-1,1] !!
	    gauss_weights[i][patch*N_Gauss+n] = h*GaussWeights[N_Gauss-1][n];                   //! Attention: GaussWeights refer to intervall [0,1] !!
	  }
    }
  }


clock_t begin_eval_u_eps = clock();
  // evaluate u_epsilon on entire grid of gauss points on [0,1]^2. 'evaluate' makes use of the tensor product structure of the cube basis.
  Matrix<double> evaluation_u_epsilon_grid =  evaluate(P->basis(), u_epsilon, gauss_points);

clock_t end_eval_u_eps = clock();
time_consumption_of_eval_u_eps += double(end_eval_u_eps - begin_eval_u_eps) / CLOCKS_PER_SEC;


clock_t begin_eval_deriv_u_eps = clock();
  // evaluate first derivatives u_epsilon on entire grid of gauss points on [0,1]^2. 'evaluate' makes use of the tensor product structure of the cube basis.
  Matrix<double> deriv_x_u_epsilon_grid;
  Matrix<double> deriv_y_u_epsilon_grid;
  evaluate(P->basis(), u_epsilon, gauss_points, deriv_x_u_epsilon_grid, deriv_y_u_epsilon_grid);

clock_t end_eval_deriv_u_eps = clock();
time_consumption_of_eval_deriv_u_eps += double(end_eval_deriv_u_eps - begin_eval_deriv_u_eps) / CLOCKS_PER_SEC;


  // iterate over all points and sum up the integral shares
  int index[dim];                               // current multiindex for the Gauss points/weights
  for (unsigned int i = 0; i < dim; i++)
  {
	index[i] = 0;
  }

  Point<dim> x;

  double share = 0.0;
  double integral_f_u_epsilon = 0.0;
//  double value_f_in_x;

#if 1
//! test
//double Lp_error_wav_expansion = 0.0;
//double L_infty_error_wav_expansion = 0.0;
//double value_wav_expansion_in_x;
//! test
#endif

  while (true)
  {
    // the current Gauss point x
    for (unsigned int i = 0; i < dim; i++)
	{
	  x[i] = gauss_points[i][index[i]];
    }
	// product of current Gauss weights
	weights = 1.0;
	for (unsigned int i = 0; i < dim; i++)
	{
	  weights *= gauss_weights[i][index[i]];
    }

	// compute the shares

    // evaluate u_epsilon in x
    //double value_u_epsilon_in_x = evaluate(P->basis(), u_epsilon, x);
    //double value_u_epsilon_in_x = evaluate_tuned_dim2_4(P->basis(), u_epsilon, x);       //! tuned Version (only for dimension=2, linear (2,2) Primbs basis with homog. bc's)
    double value_u_epsilon_in_x = evaluation_u_epsilon_grid(index[0], index[1]);        // fastest version. the above call of 'evaluate' makes use of tensor product structure of the cube basis.

    // evaluate u_exact (given as Function) in x
    double value_u_exact_in_x = u_exact->value(x);

    // evaluate f (given as Function) in x
    double value_f_in_x = (P->get_problem())->bvp_->f(x);

#if 0
//! test
value_wav_expansion_in_x = evaluate_tuned_dim2_4(P->basis(), u_exact_wav_coeff, x);
// Lp error
Lp_error_wav_expansion += (weights * pow( abs(value_wav_expansion_in_x - value_u_exact_in_x), param_p));
// L_infty error
if (abs(value_wav_expansion_in_x - value_u_exact_in_x) > L_infty_error_wav_expansion)
{
  L_infty_error_wav_expansion = abs(value_wav_expansion_in_x - value_u_exact_in_x);
}
//!test
#endif

    // compute share for <f,u>
    integral_f_u_epsilon += (weights * value_f_in_x * value_u_epsilon_in_x);

    // compute share for L_p-error
    double share_error_Lp =  pow( abs(value_u_epsilon_in_x - value_u_exact_in_x), param_p);
    error_Lp += ( weights * share_error_Lp );

    // L_infty error
    if (abs(value_u_epsilon_in_x - value_u_exact_in_x) > error_L_infty)
    {
      error_L_infty = abs(value_u_epsilon_in_x - value_u_exact_in_x);
    }

    // relative L_infty error
    if (abs(value_u_exact_in_x) > 1e-14)
    {
      if (abs((value_u_epsilon_in_x - value_u_exact_in_x) / value_u_exact_in_x) > rel_error_L_infty)
      {
        rel_error_L_infty = abs((value_u_epsilon_in_x - value_u_exact_in_x) / value_u_exact_in_x);
      }
    }


    // evaluate first derivative of u_epsilon in x
    //evaluate_tuned_dim2_4(P->basis(), u_epsilon, x, deriv_values_u_epsilon);  //! tuned Version (only for dimension=2, linear (2,2) Primbs basis with homog. bc's)
    deriv_values_u_epsilon[0] = deriv_x_u_epsilon_grid(index[0], index[1]);     // fastest version. the above call of 'evaluate' makes use of tensor product structure of the cube basis.
    deriv_values_u_epsilon[1] = deriv_y_u_epsilon_grid(index[0], index[1]);


    // evaluate first derivative of u_exact (given as Function) in x
    u_exact_deriv->vector_value(x, deriv_values_u_exact);


    // compute share for L_p-error of derivative
    for (unsigned int i = 0; i < dim; i++)
    {
       error_Lp_deriv[i] += ( weights * pow( abs(deriv_values_u_exact[i] - deriv_values_u_epsilon[i]), param_p) );
    }

    // L_infty-error of derivative
    for (unsigned int i = 0; i < dim; i++)
    {
      if (error_L_infty_deriv[i] < abs(deriv_values_u_exact[i] - deriv_values_u_epsilon[i]) )
      {
        error_L_infty_deriv[i] = abs(deriv_values_u_exact[i] - deriv_values_u_epsilon[i]);

    #if 0
        if (error_L_infty_deriv[i] > 1.3)
        {
          cout << "'p_poisson_data': error_L_infty_deriv[" << i << "] = " << error_L_infty_deriv[i] << endl;
          cout << "x = " << x[0] << ", " << x[1] << endl;
          cout << "deriv_values_u_exact[i] = " << deriv_values_u_exact[i] << endl;
          cout << "deriv_values_u_epsilon[i] = " << deriv_values_u_epsilon[i] << endl;
          exit(1);
        }
    #endif
      }
    }


    // compute |nabla u_epsilon|^2
	share = 0.0;
	for (unsigned int i = 0; i < dim; i++)
	{
	  share += deriv_values_u_epsilon[i]*deriv_values_u_epsilon[i];
    }

    // compute share for energy
    double share_energy = ( pow(share, param_p / 2.0) / param_p );
 	energy += (weights * share_energy);


 	// compute share for constrained energy of u_epsilon
 	double share_constr_energy = sqrt(share);
 	if (share_constr_energy < epsilon_lower)
 	{
 	  share_constr_energy = epsilon_lower;
 	}
 	else if (share_constr_energy > epsilon_upper)
 	{
 	  share_constr_energy = epsilon_upper;
 	}
 	constrained_energy += ( weights * ( ( 0.5 * pow(share_constr_energy, param_p - 2.0) * share ) + ( ((1.0/param_p) - 0.5) * pow(share_constr_energy, param_p) ) ) );


	// "++index"
	bool exit = false;
    for (unsigned int i = 0; i < dim; i++)
    {
      if (index[i] == N_Gauss*total_patches -1)
      {
	    index[i] = 0;
	    exit = (i == dim-1);
      }
      else
      {
	    index[i]++;
	    break;
      }
    }
    if (exit) break;
  }


#if 0
double  cp_Lp_error_wav_expansion =  Lp_error_wav_expansion;
Lp_error_wav_expansion = pow(cp_Lp_error_wav_expansion, 1.0/param_p);

x[0] = 0.501;
x[1] = 0.501;

double value_u_exact_in_x = u_exact->value(x);


value_wav_expansion_in_x = evaluate_tuned_dim2_4(P->basis(), u_exact_wav_coeff, x);

double orig_value_wav_expansion_in_x = evaluate(P->basis(), u_exact_wav_coeff, x);

double dummy_value_wav_expansion_in_x = evaluate_tuned_dim2_4(P->basis(), dummy_u_exact_wav_coeff, x);


cout << std::setprecision(10);
cout << u_exact_wav_coeff << endl << endl;
cout << "length of u_exact_wav_coeff: " << u_exact_wav_coeff.size() << endl << endl;

cout << dummy_u_exact_wav_coeff << endl << endl;
cout << "length of dummy_u_exact_wav_coeff: " << dummy_u_exact_wav_coeff.size() << endl << endl;

cout << "u_exact(x) = " << value_u_exact_in_x << endl;
cout << "evaluate_tuned_dim2_4: wavlet_expansion(x) = " << value_wav_expansion_in_x << endl;
cout << "abs_error(x) = " << abs(value_u_exact_in_x - value_wav_expansion_in_x) << endl << endl;

cout << "evaluate: wavlet_expansion(x) = " << orig_value_wav_expansion_in_x << endl << endl;

cout << "dummy_wavlet_expansion(x) = " << dummy_value_wav_expansion_in_x << endl << endl;

cout << "Lp_error_wav_expansion = " << Lp_error_wav_expansion << endl;
cout << "L_infty_error_wav_expansion = " << L_infty_error_wav_expansion << endl;
cout << std::setprecision(6);

exit(1);
#endif

  //! test
  current_data.energy1 = energy;
  //! test


  // compute <f,u_epsilon> = int_{[0,1]^2} f(x)*u_epsilon(x) dx = int_{[0,1]^2} f(x) * [ sum_{lambda \in u_epsilon} u_epsilon_lambda * Psi_lambda(x) ] dx
  //                       = sum_{lambda \in u_epsilon} u_epsilon_lambda * [ int_{[0,1]^2} f(x) * Psi_lambda(x) dx ]
  //                       = sum_{lambda \in u_epsilon} u_epsilon_lambda * f_lambda
  double r2 = 0.0;
  InfiniteVector<double, typename PROBLEM::Index> dummy;

  // copy P.fcoeffs (an Array1D) into an InfiniteVector (P.fcoeffs contains preconditioned wavelet coefficients of f!)
  for (typename Array1D<std::pair<typename PROBLEM::Index,double> >::const_iterator it((P->get_problem())->fcoeffs.begin()), itend((P->get_problem())->fcoeffs.end()); it != itend; ++it)
  {
    dummy.set_coefficient(it->first, it->second);
  }

  for (typename InfiniteVector<double, typename PROBLEM::Index>::const_iterator it(u_epsilon.begin()), itend(u_epsilon.end()); it != itend; ++it)
  {
    typename PROBLEM::Index lambda;
    lambda = it.index();
    r2 += ( *it * dummy[lambda] * P->D(lambda) );
  }

#if 1
  r2 = integral_f_u_epsilon;

#endif

  energy -= r2;
  constrained_energy -= r2;

  double error_Lp_copy = error_Lp;
  error_Lp = pow(error_Lp_copy, 1.0/param_p);

  double error_Lp_deriv_copy[2];
  for (unsigned int i = 0; i < dim; i++)
  {
    error_Lp_deriv_copy[i] = error_Lp_deriv[i];
    error_Lp_deriv[i] = pow(error_Lp_deriv_copy[i], 1.0/param_p);
  }


  //! test
  current_data.energy2 = r2;
  current_data.error_energy1 = abs(current_data.energy1 - u_exact_energy1);
  current_data.error_energy2 = abs(current_data.energy2 - u_exact_energy2);
  //! test


  //! store the (remaining) values in 'current_data'
  current_data.energy = energy;
  current_data.constrained_energy = constrained_energy;

  current_data.error_energy = u_exact_energy - energy;

  current_data.error_Lp = error_Lp;
  current_data.error_Lp_deriv[0] = error_Lp_deriv[0];
  current_data.error_Lp_deriv[1] = error_Lp_deriv[1];

  current_data.error_W1p = error_Lp + error_Lp_deriv[0] + error_Lp_deriv[1];

  current_data.error_L_infty = error_L_infty;
  current_data.error_L_infty_deriv[0] = error_L_infty_deriv[0];
  current_data.error_L_infty_deriv[1] = error_L_infty_deriv[1];

  current_data.error_W1_infty = error_L_infty + error_L_infty_deriv[0] + error_L_infty_deriv[1];

  current_data.rel_error_L_infty = rel_error_L_infty;


  //! compute and save sampled mappings u_approx_sm, abs_error_sm, rel_error_sm and an indexplot

  string filename_u_approx("");
  string filename_abs_error("");
  string filename_rel_error("");
  string filename_indexplot("");
  std::ostringstream ss;


  current_data.u_approx_sm = evaluate(P->basis(), u_epsilon, true, grid_resolution_sm);    //! sampled mapping of approximation u_epsilon

  //! compute absolute error of sampled approximation u_epsilon
  current_data.abs_error_sm = u_exact_sm;
  current_data.abs_error_sm.add(-1., current_data.u_approx_sm);
  current_data.abs_error_sm_max_nrm = maximum_norm(current_data.abs_error_sm.values());             //! compute max norm of absolute error

#if 1

  //! compute relative error of sampled approximation u_epsilon
  Matrix<double> abs_error_matrix(current_data.abs_error_sm.values());
  Matrix<double> rel_error_matrix(u_exact_matrix.row_dimension(), u_exact_matrix.column_dimension() );

  for (Matrix<double>::size_type row = 0; row < u_exact_matrix.row_dimension(); ++row)
  {
    for (Matrix<double>::size_type column = 0; column < u_exact_matrix.column_dimension(); ++column)
    {
      double value = 0.0;
      if ( abs(u_exact_matrix(row, column)) > 1e-14 )
      {
        value = ( abs( abs_error_matrix(row, column) / u_exact_matrix(row, column) ) );
      }
      rel_error_matrix.set_entry(row, column, value);
    }
  }

  Point<2, double> a(0.,0.), b(1., 1.);
  int number_gridpoints( 1UL << grid_resolution_sm);
  Grid<2> grid_sm(a,b, number_gridpoints);

  current_data.rel_error_sm = SampledMapping<2,double>(grid_sm, rel_error_matrix);            //! save relative error as sampled mapping
  current_data.rel_error_sm_max_nrm = maximum_norm(current_data.rel_error_sm.values());         //! compute max norm of relative error


  ss.str("");                                                 //! create filenames ...
  ss << current_iteration;
  filename_u_approx = "./results/plots/u_approx_" + ss.str() + ".m";
  filename_abs_error = "./results/plots/abs_error_" + ss.str() + ".m";
  filename_rel_error = "./results/plots/rel_error_" + ss.str() + ".m";
  filename_indexplot = "./results/plots/indexplot_" + ss.str() + ".m";        //! ... create filenames END

  std::ofstream ofs_u_approx(filename_u_approx.c_str());      //! open file
  current_data.u_approx_sm.matlab_output(ofs_u_approx);                 //! save data (matlab style)
  ofs_u_approx.close();                                       //! close file

  std::ofstream ofs_abs_error(filename_abs_error.c_str());        //! open file
  current_data.abs_error_sm.matlab_output(ofs_abs_error);                   //! save data (matlab style)
  ofs_abs_error.close();                                        //! close file

  std::ofstream ofs_rel_error(filename_rel_error.c_str());    //! open file
  current_data.rel_error_sm.matlab_output(ofs_rel_error);               //! save data (matlab style)
  ofs_rel_error.close();                                        //! close file

  std::ofstream ofs_indexplot(filename_indexplot.c_str());                                             //! open file
  plot_indices_cube(&P->basis(), u_epsilon, j-1, ofs_indexplot, "jet", true, true);      //! save data (matlab style)
  ofs_indexplot.close();                                                                               //! close file


#endif


  //! save ALL computed values in member 'Data'
  Data.push_back(current_data);
}


/*! output of data */
template <class PROBLEM>
void
pPoissonData<PROBLEM>::print_data(const int i, std::ostream& ostr)
{

  ostr << "Iteration: " << i << endl;
  ostr << "-------------" << endl;
  ostr << "   DOF = " << Data[i].DOF  << endl;
  ostr << " Parameters: "  << endl;
  ostr << "   epsilon_CDD1 = " << Data[i].epsilon_CDD1 << endl;
  ostr << "   epsilon_stabilization = " << Data[i].epsilon_stabilization << endl;
  ostr << "   max_lev_wav_CDD1 = " << Data[i].max_lev_wav_CDD1 << endl;
  ostr << "   max_lev_wav_rhs = " << Data[i].max_lev_wav_rhs << endl;
  ostr << "   max_lev_wav_u_exact = " << Data[i].max_lev_wav_u_exact << endl;
  ostr << "   min_res_quad_a = " << Data[i].min_res_quad_a << endl;
  ostr << "   doe_quad_a = " << Data[i].doe_quad_a << endl;
  ostr << "   min_res_quad_f = " << Data[i].min_res_quad_f << endl;
  ostr << "   doe_quad_f = " << Data[i].doe_quad_f << endl;
  ostr << "   gamma_CDD1 = " << Data[i].gamma_CDD1 << endl;
  ostr << "   c1 = " << Data[i].c1 << endl;
  ostr << "   c2 = " << Data[i].c2 << endl;
  ostr << "   Kappa = " << Data[i].kappa << endl;
  ostr << "   F = " << Data[i].F << endl;
  ostr << "   delta_0 = " << Data[i].delta_start << endl;
  ostr << " Errors:" << endl;
  ostr << std::setprecision(10);
  ostr << "   l_2 error of u_epsilon = " << Data[i].error_l2 << endl;
  ostr << "   l_infty error of u_epsilon = " << Data[i].error_l_infty << endl;
  ostr << "   Energy of u_epsilon = " << Data[i].energy << endl;
  ostr << "   Constrained energy of u_epsilon = " << Data[i].constrained_energy << endl;
  ostr << "   Energy of u = " << u_exact_energy << endl;
  ostr << "   Error J(u) - J(u_epsilon) = " << Data[i].error_energy << endl;

  //! test
  ostr << endl << "   TEST: u_exact energy1 = " << u_exact_energy1 << endl;
  ostr << "   TEST: energy1 = " << Data[i].energy1 << endl;
  ostr << "   TEST: error energy1 = " << Data[i].error_energy1 << endl;
  ostr << endl << "   TEST: u_exact energy2 = " << u_exact_energy2 << endl;
  ostr << "   TEST: energy2 = " << Data[i].energy2 << endl;
  ostr << "   TEST: error energy2 = " << Data[i].error_energy2 << endl;
  //!test
  ostr << endl << "   TEST: energy1 of *exact* wavelet expansion = " << energy1_exact_wav << endl;
  ostr << "   TEST: energy2 of *exact* wavelet expansion = " << energy2_exact_wav << endl;
  ostr << "   TEST: energy of *exact* wavelet expansion = " << energy_exact_wav << endl << endl;
  //! test

  ostr << "   L_p-error of u_epsilon = " << Data[i].error_Lp << endl;
  ostr << "   L_p-error of d/dx u_epsilon = " << Data[i].error_Lp_deriv[0] << endl;
  ostr << "   L_p-error of d/dy u_epsilon = " << Data[i].error_Lp_deriv[1] << endl;
  ostr << "   W^1_p-error of u_epsilon = " << Data[i].error_W1p << endl;
  ostr << "   L_infty-error of u_epsilon = " << Data[i].error_L_infty << endl;
  ostr << "   L_infty-error of d/dx u_epsilon = " << Data[i].error_L_infty_deriv[0] << endl;
  ostr << "   L_infty-error of d/dy u_epsilon = " << Data[i].error_L_infty_deriv[1] << endl;
  ostr << "   W^1_infty-error of u_epsilon = " << Data[i].error_W1_infty << endl;
  ostr << "   relative L_infty-error of u_epsilon = " << Data[i].rel_error_L_infty << endl;
#if 0
  ostr << "    max. absolute error of sampled approximation = " << Data[i].abs_error_sm_max_nrm << endl;
  ostr << "    max. relative error of sampled approximation = " << Data[i].rel_error_sm_max_nrm << endl;
#endif
  ostr << endl;
  ostr << std::setprecision(6);

}


/*! save some error plots */
template <class PROBLEM>
void
pPoissonData<PROBLEM>::save_error_plots(std::ostream&  ofs_plot_error_energy,
                                  std::ostream&  ofs_plot_error_l2,
                                  std::ostream&  ofs_plot_error_W1p,
                                  std::ostream&  ofs_plot_error_W1_infty)
{
  ofs_plot_error_energy << "x = [ ";
  ofs_plot_error_l2 << "x = [ ";
  ofs_plot_error_W1p << "x = [ ";
  ofs_plot_error_W1_infty << "x = [ ";

  for (unsigned int i = 0; i < Data.size(); i++)
  {
    ofs_plot_error_energy << i << "  ";
    ofs_plot_error_l2 << i << "  ";
    ofs_plot_error_W1p << i << "  ";
    ofs_plot_error_W1_infty << i << "  ";
  }

  ofs_plot_error_energy << " ]; \n y = [ ";
  ofs_plot_error_l2 << " ]; \n y = [ ";
  ofs_plot_error_W1p << " ]; \n y = [ ";
  ofs_plot_error_W1_infty << " ]; \n y = [ ";

  for (unsigned int i = 0; i < Data.size(); i++)
  {
    ofs_plot_error_energy << abs(Data[i].error_energy) << "  ";
    ofs_plot_error_l2 << abs(Data[i].error_l2) << "  ";
    ofs_plot_error_W1p << abs(Data[i].error_W1p) << "  ";
    ofs_plot_error_W1_infty << abs(Data[i].error_W1_infty) << "  ";
  }

  ofs_plot_error_energy << " ]; \n";
  ofs_plot_error_l2 << " ]; \n";
  ofs_plot_error_W1p << " ]; \n";
  ofs_plot_error_W1_infty << " ]; \n";

}



  template <class PROBLEM>
  double
  pPoissonData<PROBLEM>::dummy_integrate_u(const typename PROBLEM::Index& lambda) const
  {

//  typedef typename PROBLEM::Index Index;

    double r = 0;

    const int DIM = 2;

    const int degree_of_exactness_quadrature_u = 3;
    const int min_res_quadrature_u = 0;

    // first compute supp(psi_lambda)
    typename PROBLEM::WaveletBasis::Support supp;
    support(P->basis(), lambda, supp);

    // setup Gauss points and weights for a composite quadrature formula:
   // const int N_Gauss = 5;
    const int N_Gauss = (degree_of_exactness_quadrature_u + 1) / 2;

    if (min_res_quadrature_u > supp.j)  // adjust resolution of composite quadrature if necessary
    {
      unsigned int diff = (1 << (min_res_quadrature_u - supp.j));
      for (unsigned int i = 0; i < DIM; i++)
      {
        supp.a[i] *= diff;
        supp.b[i] *= diff;
      }
      supp.j = min_res_quadrature_u;
    }

    const double h = ldexp(1.0, -supp.j); // granularity for the quadrature
    FixedArray1D<Array1D<double>,DIM> gauss_points, gauss_weights, v_values;
    for (unsigned int i = 0; i < DIM; i++)
    {
      gauss_points[i].resize(N_Gauss*(supp.b[i]-supp.a[i]));
      gauss_weights[i].resize(N_Gauss*(supp.b[i]-supp.a[i]));
      for (int patch = supp.a[i]; patch < supp.b[i]; patch++)
      {
	    for (int n = 0; n < N_Gauss; n++)
	    {
	      gauss_points[i][(patch-supp.a[i])*N_Gauss+n] = h*(2*patch+1+GaussPoints[N_Gauss-1][n])/2.;
	      gauss_weights[i][(patch-supp.a[i])*N_Gauss+n] = h*GaussWeights[N_Gauss-1][n];
	    }
	  }
    }


    // compute the point values of the integrand (where we use that it is a tensor product)
    for (unsigned int i = 0; i < DIM; i++)
    {
    #if 0
      evaluate(*(P->basis().bases()[i]), 0,
	           typename PROBLEM::WaveletBasis::IntervalBasis::Index(lambda.j(), lambda.e()[i], lambda.k()[i], P->basis().bases()[i]),
	           gauss_points[i],
	           v_values[i]);
    #endif

      evaluate_primbs_22_bc11_v3(lambda.j(), lambda.e()[i], lambda.k()[i], gauss_points[i], v_values[i]);

    }

    // iterate over all points and sum up the integral shares
    int index[DIM]; // current multiindex for the point values
    for (unsigned int i = 0; i < DIM; i++)
    {
      index[i] = 0;
    }

#if 1

    Point<DIM> x;
    while (true)
    {
      for (unsigned int i = 0; i < DIM; i++)
      {
	    x[i] = gauss_points[i][index[i]];
      }
      double share = u_exact->value(x);
      for (unsigned int i = 0; i < DIM; i++)
      {
	    share *= (gauss_weights[i][index[i]] * v_values[i][index[i]]);
      }
      r += share;

      // "++index"
      bool exit = false;
      for (unsigned int i = 0; i < DIM; i++)
      {
	    if (index[i] == N_Gauss*(supp.b[i]-supp.a[i])-1)
	    {
	      index[i] = 0;
	      exit = (i == DIM-1);
	    }
	    else
	    {
	      index[i]++;
	      break;
	    }
      }
      if (exit) break;
    }

#endif

    return r;

  }

}



