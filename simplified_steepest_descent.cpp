#ifndef MAX_LOOPS
 #define MAX_LOOPS 10000
#endif
#define COMPRESSION_STRATEGY CDD1
//#define PLOT_U_EPSILON 200


#include <cmath>

#include <adaptive/apply.h>

#include <iostream>
using std::cout;
using std::endl;


#ifdef PLOT_U_EPSILON
#include <fstream>
#include <frame_evaluate.h>
#endif
#ifdef SAVE_ASYMPTOTIC
#include <fstream>
#include <utils/plot_tools.h>
#endif
#ifdef SAVE_LOG
#include <fstream>
#endif


namespace FrameTL
{
#ifdef PLOT_U_EPSILON
  template <class PROBLEM>
  void plot_u_epsilon(const PROBLEM& problem, const InfiniteVector<double, typename PROBLEM::Index>& u_epsilon, const unsigned int iteration)
  {
    cout << "- Evaluating current approximation after " << iteration << " iterations ..." << endl;
    InfiniteVector<double, typename PROBLEM::Index> u(u_epsilon);
    u.scale(&problem, -1); // undo preconditioning
    EvaluateFrame<Basis1D,PROBLEM::space_dimension,PROBLEM::space_dimension> evalObj;
    Array1D<SampledMapping<PROBLEM::space_dimension> > U = evalObj.evaluate(problem.basis(), u, true, 11); // expand in primal basis
    
    cout << "- Plotting current approximation after " << iteration << " iterations ..." << endl;
    std::ostringstream filename;
    filename << "steep_1D_bi_" << BASIS_NAME << "_approx_sol_it" << iteration << ".m";
    std::ofstream ofs_approx(filename.str().c_str());
    matlab_output(ofs_approx,U);
    ofs_approx.close();

    #ifdef PLOT_ERROR
    cout << "- Plotting error after " << iteration << " iterations ..." << endl;
    #if PROBLEM::space_dimension == 1
      #ifdef CONSTANT_RHS
      Polynomial<double> exact_solution(Vector<double>(5, "0 0 16 -32 16"));
      #else
      Biharmonic1D_Solution exact_solution;
      #endif
    #elif PROBLEM::space_dimension == 2
      Point<2> origin;
      origin[0] = 0.0;
      origin[1] = 0.0;
      CornerSingularityBiharmonic exact_solution(origin, 0.5, 1.5);
    #endif
    Array1D<SampledMapping<PROBLEM::space_dimension> > error = evalObj.evaluate_difference(problem.basis(), u, exact_solution, 11);
    std::ostringstream filename_err;
    filename_err << "steep_1D_bi_" << BASIS_NAME << "_error_out_it" << iteration << ".m";
    std::ofstream ofs_error(filename_err.str().c_str());
    matlab_output(ofs_error,error);
    ofs_error.close();
    #endif // PLOT_ERROR
    cout << "- Plotting done." << endl;
  }
#endif // PLOT_U_EPSILON

  template <class PROBLEM>
  void
  simplified_steepest_descent_SOLVE(const PROBLEM& problem, const double epsilon,
                                    InfiniteVector<double, typename PROBLEM::Index>& u_epsilon,
                                    const int jmax)
  {
    unsigned int loop;
    InfiniteVector<double, typename PROBLEM::Index> residual, applied_u, rhs, applied_res;
    double res_norm = 100;
    double alpha;
    double eta = 1.;
    #ifdef SAVE_ASYMPTOTIC
    map<double,double> asymptotic_res;
    map<double,double> asymptotic_time;
    double time;
    #endif

    u_epsilon.clear();

    #ifdef SAVE_LOG
    std::ofstream fs;
    std::ostringstream filename;
    filename << "steep_biharmonic_" << PROBLEM::space_dimension << "D_" << BASIS_NAME << "_jmax" << jmax << ".log";
    fs.open(filename.str().c_str());
    fs << "BASIS = " << BASIS_NAME << endl;
    fs << "iteration\tactive coefficients\tresidual\teta" << endl;
    #endif // SAVE_LOG

    #ifdef SAVE_ASYMPTOTIC
    clock_t time_start = clock();
    clock_t time_now;
    time = 0;
    #endif
                                
    for (loop = 0; loop < MAX_LOOPS && res_norm > epsilon; loop++) {
      // compute residual
      //APPLY(problem, u_epsilon, eta, applied_u, jmax, COMPRESSION_STRATEGY);
      //APPLY_COARSE(problem, u_epsilon, eta, applied_u, 1e-5, jmax, COMPRESSION_STRATEGY);
      InfiniteVector<double, typename PROBLEM::Index> help;
      APPLY(problem, u_epsilon, eta*1e-5, help, jmax, COMPRESSION_STRATEGY);
      cout << help.size() << " actice coefficients after APPLY. Doing COARSE ..." << endl;
      help.COARSE(eta, applied_u);
      //cout << "Au = " << endl << applied_u << endl;
      problem.RHS(eta, rhs);
      //cout << "RHS = " << endl << rhs << endl;
      residual = rhs - applied_u;
      //cout << "residual: " << endl << residual << endl;
      // compute acceleration parameter alpha
      APPLY(problem, residual, eta, applied_res, jmax, COMPRESSION_STRATEGY);
      alpha = residual*residual / (residual*applied_res);
      //cout << "alpha = " << alpha << endl;
      // compute l_2 norm of residual
      res_norm = sqrt(residual*residual);
      // do one iteration step
      u_epsilon += alpha * residual;
      //cout << "u = " << endl << u_epsilon << endl;
      // some output about the current state
      cout << "loop " << loop << ": res_norm = " << res_norm << ", eta = " << eta << ", " << u_epsilon.size() << " active coefficients" << endl;
      #ifdef SAVE_LOG
      fs << loop << "\t" << u_epsilon.size() << "\t" << res_norm << "\t" << eta << endl;
      #endif // SAVE_LOG
      // record some data
      #ifdef SAVE_ASYMPTOTIC
      asymptotic_res[u_epsilon.size()] = res_norm;
      time_now = clock();
      time += (double) (time_now - time_start) / (double)CLOCKS_PER_SEC;
      asymptotic_time[time] = res_norm;
      #endif
      #ifdef PLOT_U_EPSILON
      if (loop % PLOT_U_EPSILON == PLOT_U_EPSILON-1) // plot u_epsilon each PLOT_U_EPSILON steps
        plot_u_epsilon(problem, u_epsilon, loop+1);
      #endif
      // step down accuracy eta
      #ifdef ETA_STEP
        eta *= ETA_STEP;
      #else
        #ifdef BASIS_S
        eta *= 0.9995;
        #else
        //if (PROBLEM::space_dimension == 1)
        //else
        eta *= 0.995;
        #endif
      #endif
    }

    #ifdef SAVE_LOG
    fs.close();
    #endif // SAVE_LOG
    #ifdef SAVE_ASYMPTOTIC
    cout << "Saving asymptotic ..." << endl;
    std::ostringstream filename_asymptotic;
    filename_asymptotic << "steep_biharmonic_" << PROBLEM::space_dimension << "D_" << BASIS_NAME << "_jmax" << jmax << "_asymptotic.m";
    std::ofstream os_asymptotic(filename_asymptotic.str().c_str());
    matlab_output(asymptotic_res, os_asymptotic);
    os_asymptotic.close();
    std::ostringstream filename_asymptotic_time;
    filename_asymptotic_time << "steep_biharmonic_" << PROBLEM::space_dimension << "D_" << BASIS_NAME << "_jmax" << jmax << "_asymptotic_time.m";
    std::ofstream os_asymptotic_time(filename_asymptotic_time.str().c_str());
    matlab_output(asymptotic_time, os_asymptotic_time);
    os_asymptotic_time.close();
    #endif
  }

}
