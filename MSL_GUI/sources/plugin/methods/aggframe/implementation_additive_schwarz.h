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


#ifndef IMPLEMENTATION_ADDITIVE_SCHWARZ_H
#define IMPLEMENTATION_ADDITIVE_SCHWARZ_H

// implementation for additive_schwarz_GuiSolve(...)


#include <cmath>
#include <set>
#include <WaveletTL/adaptive/apply.h>
#include <FrameTL/cdd1_local.h>
#include <FrameTL/aggregated_frame.h>

#include "MathTL/utils/convergence_logger.h"

using std::set;
using MathTL::InfiniteVector;
using FrameTL::AggregatedFrame;
using WaveletTL::CDD1;


template <class PROBLEM>
void additive_schwarz_GuiSolve(const PROBLEM& P, const double epsilon,
                               Array1D<InfiniteVector<double, typename PROBLEM::Index> >& approximations,
                               MathTL::AbstractConvergenceLogger& logger,
                               int jmax,
                               double overlap,
                               double mu,      // Estimate for an upper bound of the energy norm of the exact solution
                               double rho = sqrt(1 - 1.0/4.0));     // The error reduction rate




// A helper routine for the adaptive algorithm. It works for the case
// of the L-shaped domain (-1,1)^2\setminus[0,1)^2, covered with the two rectangles
// [-OVERLAP,1]\times[-1,0] \cup [-1,0]\times[-1,1].
// We put into u_sparse those coefficients of u that correspond to the patch
// 1-i and that correspond to wavelets not being fully supported in patch i.
// We put into u_very_sparse those coefficients of u that correspond to the patch
// 1-i and that correspond to wavelets which intersect with patch i,
// but which are not fully contained in it.
//
// It also works for the case
// of the rectangular ring-shaped domain as in Sect. 7.2.3 of Manuel's
// PhD thesis, i.e. (-1,2)^2\setminus[0,1]^2, covered with 4 congruent
// rectangles.
// We then put into u_sparse those coefficients of u that correspond to a patch
// different from i and that correspond to wavelets not being fully
// supported in patch i.
// We then put into u_very_sparse those coefficients of u that correspond to a patch
// different from i and that correspond to wavelets which intersect with patch i
// but which are not fully contained in it.
template <class IBASIS, unsigned int DIM>
void thin_out (const AggregatedFrame<IBASIS, DIM>& frame, const int i,
               const InfiniteVector<double, typename AggregatedFrame<IBASIS, DIM>::Index>& u,
               InfiniteVector<double, typename AggregatedFrame<IBASIS, DIM>::Index>& u_sparse,
               InfiniteVector<double, typename AggregatedFrame<IBASIS, DIM>::Index>& u_very_sparse,
               double overlap)
{
    typedef typename AggregatedFrame<IBASIS, DIM>::Support SuppType;
    u_sparse.clear();
    u_very_sparse.clear();
    typename InfiniteVector<double, typename AggregatedFrame<IBASIS, DIM>::Index>::const_iterator it = u.begin();

    if (frame.n_p() == 2)
    {
        if (i==0) {
            for (; it != u.end(); ++it) {
                const SuppType* supp = &(frame.all_patch_supports[it.index().number()]);

                // check whether first line is intersected
                if ((it.index().p() == 1) && (supp->b[0] > -overlap) && (supp->a[1] < 0.) && (0. < supp->b[1]) ) {
                    u_very_sparse.set_coefficient(it.index(), *it);
                }
                // check whether second line is intersected
                if ((it.index().p() == 1) && (supp->a[1] < 0.) && (supp->a[0] < -overlap) && (-overlap < supp->b[0])) {
                    u_very_sparse.set_coefficient(it.index(), *it);
                }
                if ((it.index().p() == 1) && ((supp->a[0] < -overlap) || (supp->b[1] > 0.)) ) {
                    u_sparse.set_coefficient(it.index(), *it);
                }
                // 	if ((it.index().p() == 1) && !((supp->a[0] > -overlap) && (supp->b[1] < 0.)) ) {
                // 	  u_sparse.set_coefficient(it.index(), *it);
                //	}
            }
        }
        else if (i==1) {
            for (; it != u.end(); ++it) {
                const SuppType* supp = &(frame.all_patch_supports[it.index().number()]);
                if ((it.index().p() == 0) && (supp->a[0] < 0.) && (0. < supp->b[0])) {
                    u_very_sparse.set_coefficient(it.index(), *it);
                }
                if ( (it.index().p() == 0) && ( supp->b[0] > 0. ) ) {
                    u_sparse.set_coefficient(it.index(), *it);
                }
            }
        }
    }
    if (frame.n_p() == 3)
    {
        if (i==0) {
            for (; it != u.end(); ++it) {
                const SuppType* supp = &(frame.all_patch_supports[it.index().number()]);
                // check whether first line is intersected
                if ((it.index().p() != 0) && (supp->b[0] > -0.5) && (supp->a[1] < 0.) && (0. < supp->b[1]) ) {
                    u_very_sparse.set_coefficient(it.index(), *it);
                }
                // check whether second line is intersected
                if ((it.index().p() != 0) && (supp->a[1] < 0.) && (supp->a[0] < -0.5) && (-0.5 < supp->b[0])) {
                    u_very_sparse.set_coefficient(it.index(), *it);
                }
                if ((it.index().p() != 0) && ((supp->a[0] < -0.5) || (supp->b[1] > 0.)) ) {
                    u_sparse.set_coefficient(it.index(), *it);
                }
            }
        }
        if (i==1) {
            for (; it != u.end(); ++it) {
                const SuppType* supp = &(frame.all_patch_supports[it.index().number()]);
                // check whether first line is intersected
                if ((it.index().p() != 1) && (supp->b[1] > -0.5) && (supp->a[0] < 0.) && (0. < supp->b[0]) ) {
                    u_very_sparse.set_coefficient(it.index(), *it);
                }
                // check whether second line is intersected
                if ((it.index().p() != 1) && (supp->a[0] < 0.) && (supp->a[1] < -0.5) && (-0.5 < supp->b[1])) {
                    u_very_sparse.set_coefficient(it.index(), *it);
                }
                if ((it.index().p() != 1) && ((supp->a[1] < -0.5) || (supp->b[0] > 0.)) ) {
                    u_sparse.set_coefficient(it.index(), *it);
                }
            }
        }
        if (i==2) {
            for (; it != u.end(); ++it) {
                const SuppType* supp = &(frame.all_patch_supports[it.index().number()]);
                // check whether first line is intersected
                if ((it.index().p() != 2) &&
                        (
                            ((supp->a[0] < 0.) && (0. < supp->b[0]))
                            ||
                            ((supp->a[1] < 0.) && (0. < supp->b[1]))
                            )
                        ) {
                    u_very_sparse.set_coefficient(it.index(), *it);
                }
                if ((it.index().p() != 2) && (((supp->b[1] > 0.) || (supp->b[0] > 0.))) ) {
                    u_sparse.set_coefficient(it.index(), *it);
                }
            }
        }
    }
    if (frame.n_p() == 4) // --> rectangular ring
    {
        switch (i)
        {
        case 0: {
            for (; it != u.end(); ++it) {
                const SuppType* supp = &(frame.all_patch_supports[it.index().number()]);
                if ((it.index().p() != i) && (supp->a[1] < 0.) && (0. < supp->b[1]) ) {
                    u_very_sparse.set_coefficient(it.index(), *it);
                }
                if ((it.index().p() != i) && (supp->b[1] > 0.) ) {
                    u_sparse.set_coefficient(it.index(), *it);
                }
            }
            break;
        }
        case 1: {
            for (; it != u.end(); ++it) {
                const SuppType* supp = &(frame.all_patch_supports[it.index().number()]);
                if ((it.index().p() != i) && (supp->a[0] < 1.0) && (1.0 < supp->b[0]) ) {
                    u_very_sparse.set_coefficient(it.index(), *it);
                }
                if ((it.index().p() != i) && (supp->a[0] < 1.0) ) {
                    u_sparse.set_coefficient(it.index(), *it);
                }
            }
            break;
        }
        case 2: {
            for (; it != u.end(); ++it) {
                const SuppType* supp = &(frame.all_patch_supports[it.index().number()]);
                if ((it.index().p() != i) && (supp->a[1] < 1.0) && (1.0 < supp->b[1]) ) {
                    u_very_sparse.set_coefficient(it.index(), *it);
                }
                if ((it.index().p() != i) && (supp->a[1] < 1.0) ) {
                    u_sparse.set_coefficient(it.index(), *it);
                }
            }
            break;
        }
        case 3: {
            for (; it != u.end(); ++it) {
                const SuppType* supp = &(frame.all_patch_supports[it.index().number()]);
                if ((it.index().p() != i) && (supp->a[0] < 0.0) && (0.0 < supp->b[0]) ) {
                    u_very_sparse.set_coefficient(it.index(), *it);
                }
                if ((it.index().p() != i) && (supp->b[0] > 0.0) ) {
                    u_sparse.set_coefficient(it.index(), *it);
                }
            }
            break;
        }
        }
    }
}



//Overload for DIM = 1:
template <class IBASIS>
void thin_out (const AggregatedFrame<IBASIS, 1>& frame, const int i,
               const InfiniteVector<double, typename AggregatedFrame<IBASIS, 1>::Index>& u,
               InfiniteVector<double, typename AggregatedFrame<IBASIS, 1>::Index>& u_sparse,
               InfiniteVector<double, typename AggregatedFrame<IBASIS, 1>::Index>& u_very_sparse,
               double overlap)
{
    typedef typename AggregatedFrame<IBASIS, 1>::Support SuppType;
    u_sparse.clear();
    u_very_sparse.clear();
    typename InfiniteVector<double, typename AggregatedFrame<IBASIS, 1>::Index>::const_iterator it = u.begin();

    const double patch_width = 0.5 + overlap/2.0;

    if (i==0) {
        Point<1> x(patch_width);
        for (; it != u.end(); ++it) {
            if (it.index().p() == 1) {
                if (in_support(frame, it.index(),x))
                    u_very_sparse.set_coefficient(it.index(), *it);

                const SuppType* supp = &(frame.all_patch_supports[it.index().number()]);
                if (supp->b[0] > patch_width)
                    u_sparse.set_coefficient(it.index(), *it);
            }
        }
    }
    else if (i==1) {
        Point<1> x(1-patch_width);
        for (; it != u.end(); ++it) {
            if (it.index().p() == 0) {
                if (in_support(frame, it.index(),x))
                    u_very_sparse.set_coefficient(it.index(), *it);

                const SuppType* supp = &(frame.all_patch_supports[it.index().number()]);
                if (supp->a[0] < 1-patch_width)
                    u_sparse.set_coefficient(it.index(), *it);
            }
        }
    }
}



// Delete all coefficients of u corresponding to patch i.
template <class PROBLEM>
void remove_i (const int i, InfiniteVector<double, typename PROBLEM::Index>& u)
{
    typename InfiniteVector<double, typename PROBLEM::Index>::const_iterator it = u.begin();
    for (; it != u.end(); ++it) {
        if (it.index().p() == i)
            u.set_coefficient(it.index(), 0.0);
    }
    u.compress();
}



template <class PROBLEM>
void  additive_schwarz_GuiSolve(const PROBLEM& P, const double epsilon,
                                Array1D<InfiniteVector<double, typename PROBLEM::Index> >& approximations,
                                MathTL::AbstractConvergenceLogger& logger,
                                int jmax,
                                double overlap,
                                double mu,
                                double rho)
{
    typedef typename PROBLEM::Index Index;

    //const double M = sqrt(P.norm_A());
    const double M = 1.0; // we choose the most optimistic case

    // number of patches
    const int m = P.basis().n_p();

    // #####################################################################################
    // Setup of constants.
    // #####################################################################################
    const double C = m;

    const double sigma = std::max(1./M, C + 1./M) + 0.1;//10.0
    const int K = std::max(1,(int)ceil(log(1.0/(2.0 * M * sigma)) / log(rho)));

    const int L = (int)ceil(log(epsilon/mu) / log(2.0*pow(rho, K)*M*sigma));
    // #####################################################################################

    cout << "epsilon = " << epsilon << endl
         << "mu = " << mu << endl
         << "M = " << M << endl
         << "rho = " <<  rho << endl
         << "sigma = " << sigma << endl
         << "K = " << K << endl
         << "L = " << L << endl;
    // #####################################################################################
    // End of constant setup.
    // #####################################################################################

    // InfiniteVector's used in the adaptive algorithm
    InfiniteVector<double, Index> f, w, r, tmp, tmp2;
    InfiniteVector<double, Index> u_k, u_k_sparse,u_k_very_sparse;
    InfiniteVector<double, Index> precond_r_i;

    Array1D<InfiniteVector<double, Index> > xks(m); // stores local approximations which are used as starting vectors
    // for the next cycle of local solves
    Array1D<InfiniteVector<double, Index> > uks(m);

    // relaxation parameter
    double alpha;
    if (m == 4) // --> ring domain
    {
        alpha = 0.5;
    }
    else
    {
    // if the relaxation parameter is chosen equal to 1/m, then
    // the right-hand sides for the local auxiliary problems
    // are as sparse as for the multiplicative algorithm, cf. Manuel's
    // PhD thesis page 165/166 and Remark 6.4.
        alpha = 1.0/m;
    }

    int k = 0;
    double local_eps = 1.0;

    logger.startClock();

    // #####################################################################################
    // The adaptive algorithm.
    // #####################################################################################

    try
    {
        for (int l = 1; l < L; l++) {
            for (int p = 1; p <= K; p++) {
                for (int i = 0; i < m; i++) {

                    logger.checkAbortConditions();

                    k = (l-1)*m*K+(p-1)*m+i+1;
                    cout << "################################" << endl;
                    cout << "number of iteration = " << k << endl;
                    cout << "################################" << endl;

                    // Setup tolerance for the solution of the local problems.
                    // This needs to be manually tuned not to end up with
                    // slow performance.
                    if (m == 2)
                        local_eps = 100.0*mu*pow(2.0*pow(rho,K)*M*sigma,l-1)*pow(rho,p)/(alpha*m*K); // Faktor 100 bei Test mit 2 Patches
                    else
                        local_eps = 10.0*mu*pow(2.0*pow(rho,K)*M*sigma,l-1)*pow(rho,p)/(alpha*m*K); // Faktor 10 bei 3 Patches und ring-domain

                    cout << "tolerance for solution of local problem = " << local_eps << endl;
                    precond_r_i.clear();

                    // Preparation
                    thin_out(P.basis(), i, u_k, u_k_sparse, u_k_very_sparse, overlap);

                    tmp = u_k-u_k_sparse;
                    tmp.compress(1.0e-15);
                    tmp = u_k-((1./(m*alpha))*tmp);
                    tmp.compress(1.0e-15);

                    // 	  tmp = u_k;
                    // 	  tmp=u_k-((1./(m*alpha))*(tmp-u_k_sparse));

                    cout << "entering CDD solver..." << endl;

                    // Solution of the local problem in case we use the sparse version of the algorithm
                    // as proposed Manuel's PhD thesis, where it is proposed to throw away
                    // all degrees of freedom that are contained in the current subdomain before the local solve.
                    FrameTL::CDD1_LOCAL_SOLVE(P, i, local_eps, xks[i], precond_r_i, tmp, jmax, CDD1);

                    // Store the calculated local solution. These are always used as the initial guess in the
                    // next call of CDD1_LOCAL_SOLVE.
                    xks[i] = precond_r_i;

                    // setup global(!) intermediate approximation
                    uks[i] = u_k_sparse + (m*alpha*precond_r_i);
                    //xks[i] = precond_r_i;
                } // end loop over the patches

                // compute the average of the global intermediate iterates
                tmp.clear();
                for (int j = 0; j < m; j++) {
                    tmp = tmp + uks[j];
                }
                u_k = (1./m)*tmp;
                cout << "degrees of freedom: " << u_k.size() << endl;
            } // end loop p

            // setup tolerance for coarsening
            double coarse_tol;
            if (m == 4) // --> ring domain
            {
                coarse_tol = (sigma - 1./M)*2.0*pow(rho,K)*mu*pow(2.0*pow(rho,K)*M*sigma,l-1)*0.1;
            }
            else
            {
                coarse_tol = (sigma - 1./M)*2.0*pow(rho,K)*mu*pow(2.0*pow(rho,K)*M*sigma,l-1);
            }

            cout << "tolerance for coarsening = " << coarse_tol << endl;
            cout << "norm of u_k = " << l2_norm(u_k) << endl;
            u_k.COARSE(coarse_tol, tmp2);
            u_k = tmp2;
            cout << "degrees of freedom after coarsening: " << u_k.size() << endl;

            // #####################################################################################
            // Approximate global EXACT residual and perform output.
            // #####################################################################################
            logger.pauseClock();

            tmp = u_k;
            tmp.scale(&P,-1);

            //compute global residual
            tmp2.clear();
            for (int i = 0; i < m; i++) {
                P.RHS(1.0e-8, i, f);
                cout << "fsize exact res = " << f.size() << endl;
                APPLY(P, i, u_k, 1.0e-8, w, jmax, CDD1);
                tmp2 += f-w;
            }
            double residual_norm = l2_norm(tmp2);
            cout << "norm of global residual = " << residual_norm  << endl;

            logger.logConvergenceData(u_k.size(), residual_norm);


            // #####################################################################################
            //  End performing output
            // #####################################################################################

            logger.continueClock();

        }// end loop L
    }
    catch(...)
    {
        // collect final approximation and its local parts
        approximations[P.basis().n_p()] = u_k;
        //approximations[P.basis().n_p()] = uks[0];

        for (int i = 0; i < P.basis().n_p(); i++) {
            approximations[i].clear();
            for (typename InfiniteVector<double, Index>::const_iterator it = u_k.begin(), itend = u_k.end();
                 it != itend; ++it)
                if (it.index().p() == i)
                    approximations[i].set_coefficient(it.index(),*it);
        }

        throw;
    }


    // collect final approximation and its local parts
    approximations[P.basis().n_p()] = u_k;
    //approximations[P.basis().n_p()] = uks[0];

    for (int i = 0; i < P.basis().n_p(); i++) {
        approximations[i].clear();
        for (typename InfiniteVector<double, Index>::const_iterator it = u_k.begin(), itend = u_k.end();
             it != itend; ++it)
            if (it.index().p() == i)
                approximations[i].set_coefficient(it.index(),*it);
    }
}


#endif // IMPLEMENTATION_ADDITIVE_SCHWARZ_H
