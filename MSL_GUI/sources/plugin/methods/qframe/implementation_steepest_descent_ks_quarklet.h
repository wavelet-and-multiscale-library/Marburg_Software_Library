/*  -*- c++ -*-

   +-----------------------------------------------------------------------+
   | MSL GUI - A Graphical User Interface for the Marburg Software Library |
   |                                                                       |
   | Copyright (C) 2018 Henning Zickermann                                 |
   | Contact: <zickermann@mathematik.uni-marburg.de>                       |
   +-----------------------------------------------------------------------+

     This file extensively uses code from the Marburg Software Library,
     WaveletTL, which is Copyright (C) 2002-2009 Thorsten Raasch, Manuel Werner.


     This file is part of MSL GUI.

     MSL GUI is free software: you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation, either version 3 of the License, or
     (at your option) any later version.

     MSL GUI is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.

     You should have received a copy of the GNU General Public License
     along with MSL GUI.  If not, see <https://www.gnu.org/licenses/>.
*/


#ifndef IMPLEMENTATION_STEEPEST_DESCENT_KS_QUARKLET_H
#define IMPLEMENTATION_STEEPEST_DESCENT_KS_QUARKLET_H

#include <cmath>
#include <set>
#include <WaveletTL/adaptive/apply.h>

#include "MathTL/algebra/infinite_vector.h"
#include "MathTL/utils/convergence_logger.h"

using std::set;
using MathTL::InfiniteVector;

template <class PROBLEM, typename INDEX>
void steepest_descent_ks_QUARKLET_GuiSolve(const PROBLEM& P,  const double epsilon,
                                           InfiniteVector<double, INDEX>& approximations,
                                           MathTL::AbstractConvergenceLogger& logger,
                                           const int jmax,
                                           const WaveletTL::CompressionStrategy strategy,
                                           const int pmax, const double beta,
                                           const double a, const double b)
{
    unsigned int loops = 0;
    unsigned int niter = 0;

    // #####################################################################################
    // We have to set up the various constants steering the decay of the accuracy.
    // However, most theoretical estimates turn out to be too pessimistic. To obtain
    // a good performance, we have to manually choose them more optimistic.
    // #####################################################################################

    // norm of the pseudo inverse
    double a_inv     = P.norm_Ainv();

    // spectral condition number
    double kappa     = P.norm_A()*a_inv;

    // upper bound for the \ell_2-norm of the exact discrete solution in the range of the
    // stiffness matrix
    double omega_i   = a_inv*P.F_norm();
    cout << "a_inv = " << a_inv << endl;
    cout << "omega_i = " << omega_i << endl;

    //double delta     = 1./(5.*kappa+a_inv);
    double delta = 1.;
    cout << "delta = " << delta << endl;

    //const double A = 1 + delta;
    const double A = 1.;
    //const double C = 1.0 / ((1 - ((kappa*(delta*delta+2.*delta)+a_inv*delta)/((1-delta)*(1-delta))))
    //			    * (((1-delta)*(1-delta))/(a_inv)));
    const double C = 1.0;
    cout << "C = " << C << endl;
    const double B = C * (A*A);
    cout << "B = " << B << endl;

    //double lambda = (kappa-1)/(kappa+1) + P.norm_A()*std::max(3.*A*A*B,C*(1./(1-delta)))*delta;
    //double lambda = ((kappa-1)/(kappa+1)+1.)/2.;
    double lambda = 0.95;
    cout << "lambda = " << lambda << endl;

    const double C3 = B;
    cout << "C3 = " << C3 << endl;

    double mu        = 1.0001; //shall be > 1

    cout << "beta = " << beta << endl;

    //let K be such that beta^K * omega <= epsilon
    unsigned int K   = (int) (log(epsilon/omega_i) / log(beta) + 1);
    //let M be such that lambda^M <= ((1-delta) / (1+delta)) * (beta / ((1+3*mu)*kappa))
    int M            = std::max((int) ((log( ((1-delta)/(1+delta)) * (beta / ((1+3.0*mu)*kappa)) )
                                        / log(lambda)) + 1),1);

    cout << "K = " << K << endl;
    cout << "M = " << M << endl;
    // #####################################################################################
    // End setting up constants.
    // #####################################################################################

    // InfiniteVector's used in the iterative algorithm
    InfiniteVector<double, INDEX> w, tilde_r, help, f, Av;

    double residual_norm = 5.0;
    double dd = 0.5;

    logger.startClock();

    try
    {
        // the adaptive algorithm
        for (unsigned int i = 1; i <= K; i++)
        {
            omega_i *= beta;
            double xi_i = omega_i / ((1+3.0*mu)*C3*M);
            double nu_i = 0.;

            WaveletTL::RES_QUARKLET(P, w, xi_i, delta, omega_i/((1+3.*mu)*a_inv), jmax,
                                    tilde_r, nu_i, niter, strategy, pmax, a, b);

            while ( nu_i > omega_i/((1+3.*mu)*a_inv))
            {
                logger.checkAbortConditions();

                InfiniteVector<double, INDEX> z_i;

                // Instead of using APPLY only, we use a call of apply followed by an
                // immeadiate call of COARSE. This is done to prevent that the number
                // of non-zeros in the iterates grow very quickly.
                WaveletTL::APPLY_QUARKLET_COARSE(P, tilde_r, delta*l2_norm(tilde_r), z_i, jmax, strategy, pmax, a, b, 1e-6);
                //APPLY_COARSE(P, tilde_r, delta*l2_norm(tilde_r), z_i, 0.5, jmax, CDD1);
                double g = z_i*tilde_r;
                if  (g != 0.)
                    dd = (tilde_r*tilde_r)/g;

                w += dd*tilde_r;
                //  	InfiniteVector<double, Index> tmp;
                //  	w.COARSE(1.0/100.0*residual_norm, tmp);
                //  	w = tmp;


                cout << "descent param = " << dd << endl;
                ++loops;
                ++niter;

                WaveletTL::RES_QUARKLET(P, w, xi_i, delta, omega_i/((1+3.*mu)*a_inv), jmax,
                                        tilde_r, nu_i, niter, strategy, pmax, a, b);

                cout << "loop: " << loops << " nu = "
                     << nu_i << " epsilon = " << omega_i/((1+3.*mu)*a_inv) << endl;
                cout << "xi: " << xi_i << endl;

                residual_norm = l2_norm(tilde_r);
                cout << "residual norm = " << residual_norm << endl;
                cout << "active indices: " << w.size() << endl;

            }//end while

            cout << "#######################" << endl;
            cout << "exiting inner loop" << endl;
            cout << "#######################" << endl;
            cout << "number of calls of APPLY = " << niter << endl;
            cout << "Maximal " << K-i << " outer loops to go..." << endl;

            // #####################################################################################
            // Approximate the EXACT residual using a sufficiently small precision
            // and perform output.
            // #####################################################################################

            logger.pauseClock();

            P.RHS(1.0e-6,f);
            WaveletTL::APPLY_QUARKLET(P, w, 1.0e-6, Av, jmax, strategy, pmax, a, b);
            help = f-Av;
            residual_norm = l2_norm(help);

            logger.logConvergenceData(w.size(), residual_norm);

            logger.continueClock();

            if (residual_norm <= epsilon/a_inv)
                break;

            // #####################################################################################
            // Perform the coarsening. This can usually be dropped when COARSE is always applied
            // directly after the APPLY. On the theoretical side, however, it is still mandatory
            // for this algorithm.
            // #####################################################################################
            //       cout << "tolerance for COARSE = " << ((3.*mu*omega_i)/(1+3.*mu)) << endl;
            //       InfiniteVector<double, Index> tmp;
            //       w.COARSE(residual_norm, tmp);
            //       w = tmp;
        }// end for

        // #####################################################################################
        // The algorithm is finished.
        // Collect final approximation
        // #####################################################################################
    }
    catch(...)
    {
        approximations = w;

        throw;
    }

    approximations = w;
}



#endif // IMPLEMENTATION_STEEPEST_DESCENT_KS_QUARKLET_H
