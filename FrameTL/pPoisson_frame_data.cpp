// -*- c++ -*-
// +--------------------------------------------------------------------+
// | This file is part of FrameTL - the Frame Template Library      |
// |                                                                    |
// | Copyright (c) 2016-2019                                            |
// | Christoph Hartmann                                                      |
// +--------------------------------------------------------------------+
//
// implementation for pPoisson_frame_data.h

namespace FrameTL
{

/*! constructor */
template <class IBASIS, class PROBLEM>
pPoissonDataFrame<IBASIS, PROBLEM>::pPoissonDataFrame(const AggregatedFrame<IBASIS,2,2>* frame_param, const PROBLEM* prob,
                                                      const Function<2, double>* exact_solution, const Function<2, double>* derivative_exact_solution,
                                                      const int grid_resolution)
  : grid_resolution_sm(grid_resolution), frame(frame_param), P(prob), u_exact(exact_solution), u_exact_deriv(derivative_exact_solution), current_iteration(-1), u_exact_energy(0.0)
{
}

/*! 'precompute_exact_solution'
    computes:   - the exact solution as a sampled mapping,
                - a Matlab plot of the exact solution,
                - the energy of the exact solution
*/
template <class IBASIS, class PROBLEM>
void
pPoissonDataFrame<IBASIS, PROBLEM>::precompute_exact_solution(const int doe_quad, const int res_quad)
{

  //! exact solution as sampled mapping
  u_exact_sm.resize(frame->n_p());

  for (int i = 0; i < frame->n_p(); i++)
  {
    u_exact_sm[i] = SampledMapping<2>(*(frame->atlas()->charts()[i]), grid_resolution_sm);  // all zero
    u_exact_sm[i].add(1.,SampledMapping<2>(u_exact_sm[i], *u_exact));
  }


  //! exact solution as Matlab plot
  std::ofstream u_exact_stream("./results/plots/u_exact.m");
  matlab_output(u_exact_stream, u_exact_sm);
  //gnuplot_output(u_exact_stream, u_exact_sm);
  u_exact_stream.close();


  //! compute energy of exact solution
  Vector<double> deriv_values_u_exact;                      // to store one point evaluation of nabla u_exact
  deriv_values_u_exact.resize(2);
//  double u_exact_in_x;
//  double f_in_x;
  Vector<double> f_values;
  f_values.resize(2);

  int j = res_quad;                                      // resolution of composite quadrature
  double weights;                                        // the current Gauss quadrature weight
  double param_p = P->get_problem()->param_p;


  const unsigned int N_Gauss = (doe_quad+1)/2;
  const double h = ldexp(1.0, -j);                       // granularity for the quadrature (h=2^-j)
  const unsigned int total_patches = (1 << j);

  FixedArray1D<Array1D<double>,2> gauss_points, gauss_weights;

  Point<2> x;

  u_exact_energy = 0.0;
  double integral_f_u = 0.0;
  double share;


    // loop over all non overlapping squares forming the L-shaped domain
    for (unsigned int helpPatch = 0; helpPatch < 3; helpPatch++)
    {
      double a1, a2;
//      double b1, b2;
      switch(helpPatch)
      {
          case 0:
          {
            a1 = 0.;
//            b1 = 1.;
            a2 = -1.;
//            b2 = 0.;
            break;
          }
          case 1:
          {
            a1 = -1.;
//            b1 = 0.;
            a2 = -1.;
//            b2 = 0.;
            break;
          }
          case 2:
          {
            a1 = -1;
//            b1 = 0.;
            a2 = 0.;
//            b2 = 1.;
            break;
          }
      }

      // setup Gauss points and weights for current square
      //##################################################
      for (unsigned int i = 0; i < 2; i++)
      {
	    gauss_points[i].resize(N_Gauss*total_patches);
	    gauss_weights[i].resize(N_Gauss*total_patches);
      }

      for (unsigned int k = 0; k < total_patches ; k++)
      {
	    for (unsigned int n = 0; n < N_Gauss ; n++)
	    {
          double u = a1 + k*h;
          gauss_points[0][k*N_Gauss+n] = (GaussPoints[N_Gauss-1][n]+1)*0.5*h+u;     //! Attention: GaussPoints refer to intervall [-1,1] !!
          gauss_weights[0][k*N_Gauss+n] = h*GaussWeights[N_Gauss-1][n];              //! Attention: GaussWeights refer to intervall [0,1] !!
	    }
      }

      for (unsigned int k = 0; k < total_patches ; k++)
      {
	    for (unsigned int n = 0; n < N_Gauss ; n++)
	    {
          double u = a2 + k*h;
          gauss_points[1][k*N_Gauss+n] = (GaussPoints[N_Gauss-1][n]+1)*0.5*h+u;
          gauss_weights[1][k*N_Gauss+n] = h*GaussWeights[N_Gauss-1][n];
	    }
      }
      //##################################################


      // loop over all Gauss points in current square and sum up the error shares
      //#########################################################################

      for (unsigned int k1 = 0; k1 < total_patches ; k1++)
      {
	    for (unsigned int n1 = 0; n1 < N_Gauss ; n1++)
	    {
	      for (unsigned int k2 = 0; k2 < total_patches ; k2++)
	      {
	        for (unsigned int n2 = 0; n2< N_Gauss ; n2++)
	        {
	          // setup current Gauss point:
              x[0] = gauss_points[0][k1*N_Gauss+n1];
              x[1] = gauss_points[1][k2*N_Gauss+n2];

              // current Gauss weight:
              weights = gauss_weights[0][k1*N_Gauss+n1]*gauss_weights[1][k2*N_Gauss+n2];

              // evaluate first derivative of the exact solution in point x:
              u_exact_deriv->vector_value(x, deriv_values_u_exact);

              // compute |nabla u_epsilon|^2
              share = 0.0;
              for (unsigned int i = 0; i < 2; i++)
              {
                share += deriv_values_u_exact[i]*deriv_values_u_exact[i];
              }

              // compute share for first term of energy (energy1)
              double share_energy = ( pow(share, param_p / 2.0) / param_p );
              u_exact_energy += (weights * share_energy);

              // compute share for second term of energy (energy2)
#if !PPDF_USE_ALTERNATIVE_F         //! if f is smooth

              // evaluate the exact solution in point x:
              u_exact_in_x = u_exact->value(x);

              // evaluate f (given as Function) in x:
              f_in_x = P->get_problem()->get_bvp().f(x);

              // compute share for <f,u>
              integral_f_u += (weights * f_in_x * u_exact_in_x);


#else         //! else if f not smooth

              // evaluate f_ar in x:
              P->get_problem()->get_bvp().f_ar(x, f_values);

              // compute share for <f_ar, \nabla u>
              share = 0.0;
              for (unsigned int i = 0; i < 2; i++)
              {
                share += deriv_values_u_exact[i]*f_values[i];
              }

              integral_f_u += (weights * share);

#endif


	        }
	      }
	    }
      }
      cout << "done patch " << helpPatch << endl;
    }//end for helpPatch

#if 0
  //! Output of parameter mu = sqrt(a(u,u)), where a(u,u) = \int_{Omega} |\nabla u|^p dx
  cout << endl << "PARAMETER mu of u_exact: " << endl;
  cout << "-------------------" << endl;
  cout << std::setprecision(12);
  cout << "Parameter mu = " << sqrt(param_p*u_exact_energy) << endl << endl;
  cout << std::setprecision(6);
#endif

  u_exact_energy1 = u_exact_energy;
  u_exact_energy2 = integral_f_u;
  u_exact_energy -= integral_f_u;
}


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

template <class IBASIS, class PROBLEM>
void
pPoissonDataFrame<IBASIS, PROBLEM>::compute_data(const Array1D<InfiniteVector<double, typename PROBLEM::Index> >& approximations,
                                double param_p, int doe_quad, int min_res_quad, bool compute_constrained_energy,
                                double epsilon_lower, double epsilon_upper, int force_fixed_quad_res)
{
  typedef typename PROBLEM::Index Index;

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

  Matrix<double> evaluation_u_epsilon_grid;  // to store the evaluation of u_epsilon on the hole grid (of gauss points)

  Matrix<double> deriv_x_u_epsilon_grid;    // to store the evaluation of the deriv. (with respect to x) of u_epsilon on the hole grid (of gauss points)
  Matrix<double> deriv_y_u_epsilon_grid;    // to store the evaluation of the deriv. (with respect to y) of u_epsilon on the hole grid (of gauss points)

  FixedArray1D<double,2> deriv_values_u_epsilon;            // to store one point evaluation of nabla u_epsilon
  Vector<double> deriv_values_u_exact;                      // to store one point evaluation of nabla u_exact
  deriv_values_u_exact.resize(2);

  double value_u_epsilon_in_x;          // to store one point evaluation of u_epsilon
  double value_u_exact_in_x;            // to store one point evaluation of u_exact
//  double value_f_in_x;                  // to store one point evaluation of f
  Vector<double> f_values;              // to store one point evaluation of f_ar
  f_values.resize(2);

  double integral_f_u_epsilon = 0.0;

  FixedArray1D<Array1D<double>,2> gauss_points, gauss_weights;


  //! -------------------     compute and save the values of c1, c2, kappa, F, delta_start and DOF     -------------------------


  current_data.c1_patch0 = P->c1_patch0;
  current_data.c2_patch0 = P->c2_patch0;
  current_data.c1_patch1 = P->c1_patch1;
  current_data.c2_patch1 = P->c2_patch1;

  current_data.norm_A = P->norm_A();
  //current_data.kappa = current_data.c2 / current_data.c1;
  current_data.kappa = 99;
  current_data.F = P->F_norm();
  //current_data.delta_start = current_data.F / current_data.c1;
  current_data.delta_start = 99;
  current_data.DOF = approximations[frame->n_p()].size();


  //! -------------------     compute energy, constrained energy, L_p, W^1_p, L_infty, W^1_infty errors     -------------------------


  // set resolution of quadrature
  if (force_fixed_quad_res > 0)   //! fixed resolution for composite quadrature
  {
    j = force_fixed_quad_res;
  }
  else  //! set resolution of quadrature to maximum of: min_res_quad and 'highest wavelet level of u_epsilon' + 1
  {
    typename std::set<Index> supp;
    approximations[frame->n_p()].support(supp);
    typename std::set<Index>::const_iterator it = supp.end();

    if (supp.empty())
    {
      j = frame->j0() + 1;
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
  }

  cout << "'compute_energy': resolution of quadrature = " << j << endl;


  const unsigned int N_Gauss = (doe_quad+1)/2;
  const double h = ldexp(1.0, -j);        // granularity for the quadrature (h=2^-j)
  const unsigned int total_patches = (1 << j);

  Point<2> x, x_patch;

  double share_error_Lp = 0.0;
  double share;


    // loop over all non overlapping squares forming the L-shaped domain
    for (unsigned int helpPatch = 0; helpPatch < 3; helpPatch++)
    {
      double a1, a2;
//      double b1, b2;
      switch(helpPatch)
      {
          case 0:
          {
            a1 = 0.;
//            b1 = 1.;
            a2 = -1.;
//            b2 = 0.;
            break;
          }
          case 1:
          {
            a1 = -1.;
//            b1 = 0.;
            a2 = -1.;
//            b2 = 0.;
            break;
          }
          case 2:
          {
            a1 = -1;
//            b1 = 0.;
            a2 = 0.;
//            b2 = 1.;
            break;
          }
      }


      // setup Gauss points and weights for current square
      //##################################################
      for (unsigned int i = 0; i < 2; i++)
      {
	    gauss_points[i].resize(N_Gauss*total_patches);
	    gauss_weights[i].resize(N_Gauss*total_patches);
      }

      for (unsigned int k = 0; k < total_patches ; k++)
      {
	    for (unsigned int n = 0; n < N_Gauss ; n++)
	    {
          double u = a1 + k*h;
          //double v = a1 + (k+1)*h;
          gauss_points[0][k*N_Gauss+n] = (GaussPoints[N_Gauss-1][n]+1)*0.5*h+u;     //! Attention: GaussPoints refer to intervall [-1,1] !!
          gauss_weights[0][k*N_Gauss+n] = h*GaussWeights[N_Gauss-1][n];              //! Attention: GaussWeights refer to intervall [0,1] !!
	    }
      }

      for (unsigned int k = 0; k < total_patches ; k++)
      {
	    for (unsigned int n = 0; n < N_Gauss ; n++)
	    {
          double u = a2 + k*h;
          //double v = a2 + (k+1)*h;
          gauss_points[1][k*N_Gauss+n] = (GaussPoints[N_Gauss-1][n]+1)*0.5*h+u;
          gauss_weights[1][k*N_Gauss+n] = h*GaussWeights[N_Gauss-1][n];
	    }
      }
      //##################################################


#if 1  //! use this branch if charts *ARE* SIMPLE AFFINE LINEAR!!
      // evaluate u_epsilon on entire grid of gauss points on 'helpPatch'. 'evaluate' makes use of the tensor product structure of the cube basis.
      // -> ONLY FOR SIMPLE AFFINE LINEAR CHARTS!!
      evaluation_u_epsilon_grid = _frame_eval.evaluate(*frame, approximations[frame->n_p()], gauss_points);

      // evaluate first derivatives u_epsilon on entire grid of gauss points on 'helpPatch'. 'evaluate' makes use of the tensor product structure of the cube basis.
      // -> ONLY FOR SIMPLE AFFINE LINEAR CHARTS!!
      _frame_eval.evaluate(*frame, approximations[frame->n_p()], gauss_points, deriv_x_u_epsilon_grid, deriv_y_u_epsilon_grid);
#endif


      // loop over all Gauss points in current square and sum up the error shares
      //#########################################################################

      for (unsigned int k1 = 0; k1 < total_patches ; k1++)
      {
	    for (unsigned int n1 = 0; n1 < N_Gauss ; n1++)
	    {
	      for (unsigned int k2 = 0; k2 < total_patches ; k2++)
	      {
	        for (unsigned int n2 = 0; n2< N_Gauss ; n2++)
	        {
	          // setup current Gauss point:
              x[0] = gauss_points[0][k1*N_Gauss+n1];
              x[1] = gauss_points[1][k2*N_Gauss+n2];

              // current Gauss weight:
              weights = gauss_weights[0][k1*N_Gauss+n1]*gauss_weights[1][k2*N_Gauss+n2];

#if 0  //! use this branch if charts are NOT SIMPLE AFFINE LINEAR!!
              // evaluate u_epsilon in x:
              value_u_epsilon_in_x = 0.0;

              // evaluate the approximate solution given by 'u_epsilon' in point x:
              for(typename InfiniteVector<double,Index>::const_iterator it(approximations[frame->n_p()].begin()),
                  itend(approximations[frame->n_p()].end()); it != itend; ++it)
              {
                if (in_support(*frame, it.index(), x))
                {
                  frame->atlas()->charts()[it.index().p()]->map_point_inv(x,x_patch);
                  double wav_val_x = WaveletTL::evaluate(*(frame->bases()[it.index().p()]->bases()[0]), 0,
                                     typename IBASIS::Index(it.index().j(),
                                                it.index().e()[0],
                                                it.index().k()[0],
                                                frame->bases()[it.index().p()]->bases()[0]),
                                     x_patch[0]);
                  double wav_val_y = WaveletTL::evaluate(*(frame->bases()[it.index().p()]->bases()[1]), 0,
                                     typename IBASIS::Index(it.index().j(),
                                                it.index().e()[1],
                                                it.index().k()[1],
                                                frame->bases()[it.index().p()]->bases()[1]),
                                     x_patch[1]);
                  value_u_epsilon_in_x += (*it)*(wav_val_x*wav_val_y) / frame->atlas()->charts()[it.index().p()]->Gram_factor(x_patch);
                }// end if
              }// end for over coeffs
#endif

#if 1  //! use this branch if charts *ARE* SIMPLE AFFINE LINEAR!!
              // evaluate u_epsilon in x:
              value_u_epsilon_in_x = evaluation_u_epsilon_grid(k1*N_Gauss+n1,k2*N_Gauss+n2);
#endif

              // evaluate the exact solution in point x:
              value_u_exact_in_x = u_exact->value(x);

              // compute share for L_p-error
              share_error_Lp = pow(abs(value_u_epsilon_in_x - value_u_exact_in_x), param_p);  // for L_p -error
              error_Lp += weights*share_error_Lp;

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


#if 0  //! use this branch if charts are NOT SIMPLE AFFINE LINEAR!!
              // evaluate first derivative of u_epsilon in x:
              _frame_eval.evaluate_deriv_tuned( *frame, approximations[frame->n_p()], x, deriv_values_u_epsilon);
#endif

#if 1  //! use this branch if charts *ARE* SIMPLE AFFINE LINEAR!!
              // evaluate first derivative of u_epsilon in x:
              deriv_values_u_epsilon[0] = deriv_x_u_epsilon_grid(k1*N_Gauss+n1,k2*N_Gauss+n2);
              deriv_values_u_epsilon[1] = deriv_y_u_epsilon_grid(k1*N_Gauss+n1,k2*N_Gauss+n2);
#endif


              // evaluate first derivative of the exact solution in point x:
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

              //! compute energy (error)
              // compute first term of energy (energy1)
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

              // compute second term of energy (energy2)
#if !PPDF_USE_ALTERNATIVE_F          //! if f is smooth

              // evaluate f (given as Function) in x:
              value_f_in_x = P->get_problem()->get_bvp().f(x);

              // compute share for <f,u>
              integral_f_u_epsilon += (weights * value_f_in_x * value_u_epsilon_in_x);


#else         //! else if f not smooth

              // evaluate f_ar in x:
              P->get_problem()->get_bvp().f_ar(x, f_values);

              // compute share for <f_ar, \nabla u>
              share = 0.0;
              for (unsigned int i = 0; i < 2; i++)
              {
                share += deriv_values_u_epsilon[i]*f_values[i];
              }

              integral_f_u_epsilon += (weights * share);

#endif

	        }
	      }
	    }
      }
      cout << "done patch " << helpPatch << endl;
    }//end for helpPatch


  double error_Lp_copy = error_Lp;
  error_Lp = pow(error_Lp_copy, 1.0/param_p);


  cout << endl;
  cout << "L_p Error:" << error_Lp << endl;
  cout << "L_infty Error:" << error_L_infty << endl;
  cout << "re. L_infty Error:" << rel_error_L_infty << endl << endl;


  current_data.energy1 = energy;


#if 0
  //! Alternative way to compute <f,u_epsilon>:
  //         <f,u_epsilon> = int_{\Omega} f(x)*u_epsilon(x) dx = int_{\Omega} f(x) * [ sum_{lambda \in u_epsilon} u_epsilon_lambda * Psi_lambda(x) ] dx
  //                       = sum_{lambda \in u_epsilon} u_epsilon_lambda * [ int_{\Omega} f(x) * Psi_lambda(x) dx ]
  //                       = sum_{lambda \in u_epsilon} u_epsilon_lambda * f_lambda

  double r2 = 0.0;
  InfiniteVector<double, typename PROBLEM::Index> dummy;

  // copy P.fcoeffs (an Array1D) into an InfiniteVector (P.fcoeffs contains preconditioned wavelet coefficients of f!)
  for (typename Array1D<std::pair<typename PROBLEM::Index,double> >::const_iterator it((P->get_problem())->fcoeffs.begin()), itend((P->get_problem())->fcoeffs.end()); it != itend; ++it)
  {
    dummy.set_coefficient(it->first, it->second);
  }

  for (typename InfiniteVector<double, typename PROBLEM::Index>::const_iterator it(approximations[frame->n_p()].begin()), itend(approximations[frame->n_p()].end()); it != itend; ++it)
  {
    typename PROBLEM::Index lambda;
    lambda = it.index();
    r2 += ( *it * dummy[lambda] * P->D(lambda) );
  }
  integral_f_u_epsilon = r2;
#endif


  energy -= integral_f_u_epsilon;
  constrained_energy -= integral_f_u_epsilon;


  double error_Lp_deriv_copy[2];
  for (unsigned int i = 0; i < dim; i++)
  {
    error_Lp_deriv_copy[i] = error_Lp_deriv[i];
    error_Lp_deriv[i] = pow(error_Lp_deriv_copy[i], 1.0/param_p);
  }


  //! -------------------------     store the computed values in 'current_data'     -----------------------------


  current_data.energy2 = integral_f_u_epsilon;
  current_data.error_energy1 = abs(current_data.energy1 - u_exact_energy1);
  current_data.error_energy2 = abs(current_data.energy2 - u_exact_energy2);

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


  //! -------------------------     plot approximation, (pointwise) error function, indexplot     ---------------------------


  //! compute and save sampled mappings u_approx_sm, abs_error_sm and an indexplot

  string filename_u_approx("");
  string filename_u_approx_local0("");
  string filename_u_approx_local1("");
  string filename_abs_error("");
  string filename_indexplot("");
  std::ostringstream ss;
  std::ostringstream ss2;

  //! filenames
  ss.str("");                                                                                   // create filenames ...
  ss << current_iteration;
  filename_u_approx = "./results/plots/u_approx_" + ss.str() + ".m";
  filename_u_approx_local0 = "./results/plots/u_approx_" + ss.str() + "_patch0" + ".m";
  filename_u_approx_local1 = "./results/plots/u_approx_" + ss.str() + "_patch1" + ".m";
  filename_abs_error = "./results/plots/abs_error_" + ss.str() + ".m";                          // ... create filenames END


  //! sampled mapping of approximation u_epsilon
  current_data.u_approx_sm = _frame_eval.evaluate(*frame, approximations[frame->n_p()], true, grid_resolution_sm);

  std::ofstream ofs_u_approx(filename_u_approx.c_str());                     // open file
  matlab_output(ofs_u_approx, current_data.u_approx_sm);                     // save data (matlab style)
  ofs_u_approx.close();                                                      // close file


  //! compute absolute error of sampled approximation u_epsilon
  current_data.abs_error_sm.resize(frame->n_p());

  for (int i = 0; i < frame->n_p(); i++)
  {
    current_data.abs_error_sm[i] = SampledMapping<2>(*(frame->atlas()->charts()[i]), grid_resolution_sm);  // all zero
    current_data.abs_error_sm[i].add(1., current_data.u_approx_sm[i]);
    current_data.abs_error_sm[i].add(-1., u_exact_sm[i]);
  }

  std::ofstream ofs_abs_error(filename_abs_error.c_str());                   // open file
  matlab_output(ofs_abs_error, current_data.abs_error_sm);                   // save data (matlab style)
  ofs_abs_error.close();                                                     // close file


  //! sampled mapping of *local* approximations

  Array1D<SampledMapping<2,double> > u_approx_local_sm;

  //! patch 0:
  u_approx_local_sm = _frame_eval.evaluate(*frame, approximations[0], true, grid_resolution_sm);
  std::ofstream ofs_loc0(filename_u_approx_local0.c_str());
  matlab_output(ofs_loc0, u_approx_local_sm);
  ofs_loc0.close();

  //! patch 1:
  u_approx_local_sm = _frame_eval.evaluate(*frame, approximations[1], true, grid_resolution_sm);
  std::ofstream ofs_loc1(filename_u_approx_local1.c_str());
  matlab_output(ofs_loc1, u_approx_local_sm);
  ofs_loc1.close();



#if 1
  //! save indexplots

  Array1D<InfiniteVector<double, CubeIndex<IBASIS,2,MappedCubeBasis<IBASIS,2,2> > > > approximations_cube(frame->n_p());

  //! convert indices to CubeIndices
  for (int i = 0; i < frame->n_p(); i++)
  {
    ss2.str("");
    ss2 << i;

    MappedCubeBasis<IBASIS,2,2>* mapped_basis = frame->bases()[i];
    std::ofstream ofs_indexplot;

    filename_indexplot = "./results/plots/indexplot_" + ss.str() + "_patch" + ss2.str() + ".m";

    ofs_indexplot.open(filename_indexplot.c_str());


    for ( typename InfiniteVector<double, typename PROBLEM::Index >::const_iterator it = approximations[frame->n_p()].begin(),
	   itend = approximations[frame->n_p()].end(); it != itend; ++it)
    {
      if (it.index().p() == i)
      {
        approximations_cube[i].set_coefficient( CubeIndex<IBASIS,2,MappedCubeBasis<IBASIS,2,2> >(it.index().j(),it.index().e(),it.index().k(), mapped_basis),*it);
      }
    }

    plot_indices_cube(frame->bases()[i], approximations_cube[i], frame->jmax(), ofs_indexplot, "jet", false, true);

    ofs_indexplot.close();

  }
#endif


  //! save ALL computed values in member 'Data'
  Data.push_back(current_data);
}


/*! output of data */
template <class IBASIS, class PROBLEM>
void
pPoissonDataFrame<IBASIS, PROBLEM>::print_data(const int i, std::ostream& ostr)
{

  ostr << "--------------------------------------------------------" << endl;
  ostr << "Iteration: " << i <<  endl;
  ostr << "--------------------------------------------------------" << endl;
  ostr << "   DOF = " << Data[i].DOF  << endl << endl;
  ostr << " Parameters: "  << endl;
  ostr << " ---------- "  << endl;
  ostr << "   epsilon_MS = " << Data[i].epsilon_MS << endl;
  ostr << "   epsilon_stabilization = " << Data[i].epsilon_stabilization << endl;
  ostr << "   max_lev_wav_CDD1 = " << Data[i].max_lev_wav_CDD1 << endl;
  ostr << "   max_lev_wav_rhs = " << Data[i].max_lev_wav_rhs << endl;
  ostr << "   max_lev_wav_u_exact = " << Data[i].max_lev_wav_u_exact << endl;
  ostr << "   min_res_quad_a = " << Data[i].min_res_quad_a << endl;
  ostr << "   doe_quad_a = " << Data[i].doe_quad_a << endl;
  ostr << "   min_res_quad_f = " << Data[i].min_res_quad_f << endl;
  ostr << "   doe_quad_f = " << Data[i].doe_quad_f << endl;
  ostr << "   gamma_CDD1 = " << Data[i].gamma_CDD1 << endl;
  ostr << "   c1_patch0 = " << Data[i].c1_patch0 << endl;
  ostr << "   c2_patch0 = " << Data[i].c2_patch0 << endl;
  ostr << "   c1_patch1 = " << Data[i].c1_patch1 << endl;
  ostr << "   c2_patch1 = " << Data[i].c2_patch1 << endl;
  ostr << "   NormA = " << Data[i].norm_A << endl;
  ostr << "   F = " << Data[i].F << endl;
  ostr << "   delta_0 = " << Data[i].delta_start << endl << endl;
  ostr << " Errors:" << endl;
  ostr << " ------ "  << endl;
  ostr << std::setprecision(10);
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
  ostr << "   TEST: error energy2 = " << Data[i].error_energy2 << endl << endl;

  ostr << "   L_p-error of u_epsilon = " << Data[i].error_Lp << endl;
  ostr << "   L_p-error of d/dx u_epsilon = " << Data[i].error_Lp_deriv[0] << endl;
  ostr << "   L_p-error of d/dy u_epsilon = " << Data[i].error_Lp_deriv[1] << endl;
  ostr << "   W^1_p-error of u_epsilon = " << Data[i].error_W1p << endl << endl;
  ostr << "   L_infty-error of u_epsilon = " << Data[i].error_L_infty << endl;
  ostr << "   L_infty-error of d/dx u_epsilon = " << Data[i].error_L_infty_deriv[0] << endl;
  ostr << "   L_infty-error of d/dy u_epsilon = " << Data[i].error_L_infty_deriv[1] << endl;
  ostr << "   W^1_infty-error of u_epsilon = " << Data[i].error_W1_infty << endl << endl;
  ostr << "   relative L_infty-error of u_epsilon = " << Data[i].rel_error_L_infty << endl;

  ostr << endl;
  ostr << std::setprecision(6);

}


/*! save some error plots */
template <class IBASIS, class PROBLEM>
void
pPoissonDataFrame<IBASIS, PROBLEM>::save_error_plots(std::ostream&  ofs_plot_error_energy,
                                  std::ostream&  ofs_plot_error_W1p,
                                  std::ostream&  ofs_plot_error_W1_infty)
{
  ofs_plot_error_energy << "x = [ ";
  ofs_plot_error_W1p << "x = [ ";
  ofs_plot_error_W1_infty << "x = [ ";

  for (unsigned int i = 0; i < Data.size(); i++)
  {
    ofs_plot_error_energy << i << "  ";
    ofs_plot_error_W1p << i << "  ";
    ofs_plot_error_W1_infty << i << "  ";
  }

  ofs_plot_error_energy << " ]; \n y = [ ";
  ofs_plot_error_W1p << " ]; \n y = [ ";
  ofs_plot_error_W1_infty << " ]; \n y = [ ";

  for (unsigned int i = 0; i < Data.size(); i++)
  {
    ofs_plot_error_energy << abs(Data[i].error_energy) << "  ";
    ofs_plot_error_W1p << abs(Data[i].error_W1p) << "  ";
    ofs_plot_error_W1_infty << abs(Data[i].error_W1_infty) << "  ";
  }

  ofs_plot_error_energy << " ]; \n";
  ofs_plot_error_W1p << " ]; \n";
  ofs_plot_error_W1_infty << " ]; \n";

}



  template <class IBASIS, class PROBLEM>
  double
  pPoissonDataFrame<IBASIS, PROBLEM>::dummy_integrate_u(const typename PROBLEM::Index& lambda) const
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



