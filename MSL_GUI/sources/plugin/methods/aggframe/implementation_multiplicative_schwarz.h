/*  -*- c++ -*-

   +-----------------------------------------------------------------------+
   | MSL GUI - A Graphical User Interface for the Marburg Software Library |
   |                                                                       |
   | Copyright (C) 2018 Henning Zickermann                                 |
   | Contact: <zickermann@mathematik.uni-marburg.de>                       |
   +-----------------------------------------------------------------------+

     This file extensively uses code from the Marburg Software Library,
     FrameTL, which is Copyright (C) 2002-2010 Thorsten Raasch, Manuel Werner.


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


#ifndef IMPLEMENTATION_MULTIPLICATIVE_SCHWARZ_H
#define IMPLEMENTATION_MULTIPLICATIVE_SCHWARZ_H


// implementation for multiplicative_schwarz_GuiSolve(...)

#include <cmath>
#include <set>
#include <WaveletTL/adaptive/apply.h>
#include <FrameTL/cdd1_local.h>

#include "MathTL/utils/convergence_logger.h"
#include "implementation_additive_schwarz.h"    // for routine "thin_out(...)"

using std::set;
using MathTL::InfiniteVector;



template <class PROBLEM>
void  multiplicative_schwarz_GuiSolve(const PROBLEM& P, const double epsilon,
                                      Array1D<InfiniteVector<double, typename PROBLEM::Index> >& approximations,
                                      MathTL::AbstractConvergenceLogger& logger,
                                      int jmax,
                                      double overlap,
                                      bool sparseVersion,
                                      double mu,    // Estimate for an upper bound of the energy norm of the exact solution
                                      double rho)   // The error reduction rate
{
    typedef typename PROBLEM::Index Index;

    //const double M = 1.0; // we choose the most optimistic case
    //const double M = sqrt(P.norm_A());     //! Christoph
    const double M = 2;             //! Christoph: speedup

    // #####################################################################################
    // Setup of constants.
    // #####################################################################################
    double rho_1 = pow(rho, 1./P.basis().n_p());
    // const double C = 0.5 * rho_1*((1./rho)-1) / (1-rho_1);
    const double C = rho_1 * (1-rho) / (1-rho_1);               // Christoph: Due to Manuel's Diss, p.146 proof of Prop.6.2
    //const double sigma = std::max(1./M, C + 1./M) + 0.1;//10.0
    const double sigma = C/(2*rho*rho) + 1/M + 0.1;//       // Christoph: Due to Manuel's Diss, p.146 (6.1.20)
    cout << "sigma = " << sigma << endl;
    int K = std::max(1,(int)ceil(log(1.0/(2.0 * M * sigma)) / log(rho)));

    int L = (int)ceil( log(epsilon/mu) / log(2.0 * pow(rho,K) * M * sigma) );

    //! Parameters for CDD1
    double c1 = 1.0;   // Manuels's version
    double c2 = 1.0;   //

    cout << endl << "MultSchw parameters:" << endl
         << "  epsilon = " << epsilon << endl
         << "  mu = " << mu << endl
         << "  M = " << M << endl
         << "  rho = " <<  rho << endl
         << "  C = " <<  C << endl
         << "  sigma = " << sigma << endl
         << "  K = " << K << endl
         << "  L = " << L << endl << endl;
    // #####################################################################################
    // End of constant setup.
    // #####################################################################################

    // InfiniteVector's used in the adaptive algorithm
    InfiniteVector<double, Index> f, w, r, tmp, tmp_w;
    InfiniteVector<double, Index> u_k, u_k_sparse,u_k_very_sparse;
    InfiniteVector<double, Index> precond_r_i;

    Array1D<InfiniteVector<double, Index> > xks(P.basis().n_p()); // stores local approximations which are used as starting vectors
    // for the next cycle of local solves

    // initialize u_k and xks with 'guess'
    u_k.clear();

    for (int i = 0; i < P.basis().n_p(); i++)
    {
        xks[i].clear();
        for (typename InfiniteVector<double, Index>::const_iterator it = u_k.begin(), itend = u_k.end(); it != itend; ++it)
        {
            if (it.index().p() == i)
            {
                xks[i].set_coefficient(it.index(),*it);
            }
        }
    }

    // number of patches
    const int m = P.basis().n_p();
    int k = 0;
    double local_eps = 1.0;


    logger.startClock();

    // #####################################################################################
    // The adaptive algorithm.
    // #####################################################################################

    try
    {
        for (int l = 1; l <= L; l++)
        {
            //for (int p = 2; p <= K; p++) {// THAT WAS USED FOR THE 2D PLAIN DD CASE!!!
            for (int p = 1; p <= K; p++)
            {
                for (int i = 0; i < m; i++)
                {
                    logger.checkAbortConditions();

                    k = (l-1)*m*K+(p-1)*m+i+1;
                    cout << "################################" << endl;
                    cout << "number of iteration = " << k << endl;
                    cout << "################################" << endl;

                    // Setup tolerance for the solution of the local problems.
                    // This needs to be manually tuned not to end up with
                    // slow performance.
                    if (m == 4) // --> ring domain
                        local_eps = mu*pow(2.0*pow(rho,K)*M*sigma,l-1)*pow(rho,p)/(m*K)*10;
                    else
                        //local_eps = mu*pow(2.0*pow(rho,K)*M*sigma,l-1)*pow(rho,p)/(m*K)*100;
                        local_eps = mu*pow(2.0*pow(rho,K)*M*sigma,l-1)*pow(rho,p)/(m*K)*50.0;       // Christoph: Test

                    cout << "tolerance for solution of local problem = " << local_eps << endl;

                    precond_r_i.clear();


                    // #####################################################################################
                    // Now follows the solution of the local problems.
                    // #####################################################################################
                    set<Index> Lambda_i;

                    // CDD1 parameter c1,c2
#if 0
                    if(i==0)
                    {
                        c1 = P.c1_patch0;
                        c2 = P.c2_patch0;
                    }
                    else
                    {
                        c1 = P.c1_patch1;
                        c2 = P.c2_patch1;
                    }
#endif


                    if (sparseVersion)
                    {
                        // Preparation
                        thin_out(P.basis(), i, u_k, u_k_sparse, u_k_very_sparse, overlap);
                    }
                    // CDD1
                    cout << "entering CDD solver..." << endl;
                    cout << "jmax = " << jmax << endl;
                    if (sparseVersion)
                    {
                        // Solution of the local problem in case we use the sparse version of the algorithm
                        // as proposed in Stevenson, Werner 2009, where it is proposed to throw away
                        // all degrees of freedom that are contained in the current subdomain before the local solve.
                        FrameTL::CDD1_LOCAL_SOLVE(P, i, local_eps, xks[i], precond_r_i, u_k_very_sparse, c1, c2, jmax, WaveletTL::CDD1);
                    }
                    else
                    {
                        // Solution of the local problem in case we use a plain multiplicative Schwarz adaptive
                        // method, i.e., without removing degrees of freedom in the overlapping region before
                        // the local solve as in the SPARSE-branch.
                        //remove_i<PROBLEM>(i, u_k);
                        if (PROBLEM::space_dimension == 1)
                            FrameTL::CDD1_LOCAL_SOLVE(P, i, 10.0*local_eps, xks[i], precond_r_i, u_k, jmax, WaveletTL::CDD1); // THAT WAS USED FOR THE 1D CASE
                        else
                            FrameTL::CDD1_LOCAL_SOLVE(P, i, local_eps, xks[i], precond_r_i, u_k, jmax, WaveletTL::CDD1); // THAT WAS USED FOR THE 2D CASE

                    }
                    cout << "CDD 1 solve completed, size of output is " << precond_r_i.size() << endl;

                    // #####################################################################################
                    // setup next global iterate
                    // #####################################################################################
                    if (sparseVersion)
                        u_k = precond_r_i + u_k_sparse;
                    else
                        u_k = precond_r_i + u_k;

                    cout << "degrees of freedom: " << u_k.size() << endl;
                    // #####################################################################################

                    // Store the calculated local solution. These are always used as the initial guess in the
                    // next call of CDD1_LOCAL_SOLVE.
                    xks[i].clear();
                    xks[i] = precond_r_i;
                } // end loop over patches
                cout << "############## full cycle of local solves completed ##############"<< endl;

            }// end loop p

            // setup tolerance for coarsening
            double coarse_tol;
            if (m == 4) // --> ring domain
                coarse_tol = (sigma - 1./M)*2.0*pow(rho,K)*mu*pow(2.0*pow(rho,K)*M*sigma,l-1)*0.1;
            else
                //coarse_tol = (sigma - 1./M)*2.0*pow(rho,K)*mu*pow(2.0*pow(rho,K)*M*sigma,l-1);
                coarse_tol = (sigma - 1./M)*2.0*pow(rho,K)*mu*pow(2.0*pow(rho,K)*M*sigma,l-1)*0.5;     //! Christoph: test

            cout << "tolerance for coarsening = " << coarse_tol << endl;
            cout << "norm of u_k = " << l2_norm(u_k) << endl;
            u_k.COARSE(coarse_tol, tmp_w);                //! Christoph: for estimation of rho, comment this out
            u_k = tmp_w;                                  //! Christoph: for estimation of rho, comment this out
            cout << "degrees of freedom after coarsening: " << u_k.size() << endl;


            // #####################################################################################
            //  Approximate global EXACT residual and perform output.
            // #####################################################################################

            logger.pauseClock();

            tmp = u_k;
            tmp.scale(&P,-1);

            //compute global residual
            P.RHS(1.0e-6, f);
            cout << "fsize exact res = " << f.size() << endl;
            tmp_w.clear();
            for (int i = 0; i < P.basis().n_p(); i++) {
                WaveletTL::APPLY(P, i, u_k, 1.0e-6, w, jmax, WaveletTL::CDD1);
                tmp_w += w;
            }
            double residual_norm = l2_norm(f-tmp_w);
            cout << "norm of global residual = " << residual_norm  << endl;

            logger.logConvergenceData(u_k.size(), residual_norm);

            // #####################################################################################
            //  End performing output
            // #####################################################################################

            logger.continueClock();


        }// end loop L


        // #####################################################################################
        // The adaptive algorithm is finished here.
        // #####################################################################################
    }
    catch(...)
    {
        // collect final approximation and its local parts
        approximations[P.basis().n_p()] = u_k;

        for (int i = 0; i < P.basis().n_p(); i++)
        {
            approximations[i].clear();
            for (typename InfiniteVector<double, Index>::const_iterator it = u_k.begin(), itend = u_k.end(); it != itend; ++it)
            {
                if (it.index().p() == i)
                {
                    approximations[i].set_coefficient(it.index(),*it);
                }
            }
        }

        throw;
    }

    // collect final approximation and its local parts
    approximations[P.basis().n_p()] = u_k;

    for (int i = 0; i < P.basis().n_p(); i++)
    {
        approximations[i].clear();
        for (typename InfiniteVector<double, Index>::const_iterator it = u_k.begin(), itend = u_k.end(); it != itend; ++it)
        {
            if (it.index().p() == i)
            {
                approximations[i].set_coefficient(it.index(),*it);
            }
        }
    }


}

//! #####################################################################################
//!  End of MultSchw
//! #####################################################################################


#endif // IMPLEMENTATION_MULTIPLICATIVE_SCHWARZ_H
