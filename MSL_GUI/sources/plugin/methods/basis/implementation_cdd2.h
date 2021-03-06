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


#ifndef IMPLEMENTATION_CDD2_H
#define IMPLEMENTATION_CDD2_H


#include <MathTL/algebra/infinite_vector.h>
#include <WaveletTL/adaptive/compression.h>
#include <cmath>
#include <set>
#include <WaveletTL/adaptive/apply.h>

#include "MathTL/utils/convergence_logger.h"

using MathTL::InfiniteVector;
using std::set;
using std::map;


/*!
    An adaptive solver for the infinite-dimensional problem

      Au = F,

    as developed in [CDD2] and reformulated in [S],[DFR].
    A is assumed to be s.p.d., i.e. the simplifications on [CDD2, p.12] hold.
    Given the problem and a target accuracy \epsilon,
    the algorithm constructs a coefficient vector u_\epsilon, such that

      ||u-u_\epsilon|| <= \epsilon.

    The routine has to be given an estimate of ||u|| <= nu = epsilon_0, which may be
    computed beforehand like nu:=||A^{-1}||*||F||.

    References:
    [CDD2] Cohen/Dahmen/DeVore,
           Adaptive Wavelet Methods II - Beyond the Elliptic Case
    [DFR]  Dahlke/Fornasier/Raasch,
           Adaptive Frame Methods for Elliptic Operator Equations
    [S]    Stevenson,
           Adaptive Solution of Operator Equations using Wavelet Frames
  */


template <class PROBLEM>
void CDD2_GuiSolve(const PROBLEM& P, const double nu, const double epsilon,
                InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& u_epsilon,
                MathTL::AbstractConvergenceLogger& logger,
                const unsigned int maxlevel, WaveletTL::CompressionStrategy strategy)
{
    typedef typename PROBLEM::WaveletBasis::Index Index;

    // compute optimal relaxation parameter omega
    const double omega = 2.0 / (P.norm_A() + 1.0/P.norm_Ainv());
    //const double omega = 0.2;
    cout << "CDD2_SOLVE: omega=" << omega << endl;

    // compute spectral norm rho
    const double cond_A = P.norm_A() * P.norm_Ainv();
    const double rho = (cond_A - 1.0) / (cond_A + 1.0);
    cout << "CDD2_SOLVE: rho=" << rho << endl;

    // desired error reduction factor theta < 1/3
    //     const double theta = 2.0/7.0;
    const double theta = 0.333;
    cout << "CDD2_SOLVE: theta=" << theta << endl;

    // compute minimal K such that 3*rho^K < theta
    const int K = (int) ceil(log10(theta/3.0) / log10(rho));
    cout << "CDD2_SOLVE: K=" << K << endl << endl;

    u_epsilon.clear();

    double epsilon_k = nu, eta;
    unsigned int loops(0);
    InfiniteVector<double,Index> f, v, Av, F;

    P.RHS(1e-6, F);

#if _WAVELETTL_USE_TBASIS == 1
    Array1D<int> jp_guess(0);
#endif

    logger.startClock();

    while (epsilon_k > epsilon)
    {
        logger.checkAbortConditions();

        ++loops;
        cout << "CDD2:: loop: " << loops << endl;
        cout << "CDD2:: u.size() = " << u_epsilon.size() << endl;
        epsilon_k *= 3*pow(rho, K) / theta;
        cout << "CDD2:: epsilon_k=" << epsilon_k << endl;
        //        eta = theta * epsilon_k / (6*omega*K) ;//WAS SOLL HIER DIE 10@PHK
        eta = theta * epsilon_k / (6*omega*K)*10;
        cout << "eta = " << eta << endl;
        P.RHS(eta, f);

        v = u_epsilon;

        for (int j = 1; j <= 1 /*K*/; j++)
        {
#if _WAVELETTL_USE_TBASIS == 1
            WaveletTL::APPLY(P, v, eta, Av, maxlevel, tensor_simple);
            //WaveletTL::APPLY(P, v, eta, jp_guess, Av, maxlevel, tensor_simple);
#else
            WaveletTL::APPLY_COARSE(P, v, eta, Av, 0.5, maxlevel, strategy);
            //          APPLY(P, v, eta, Av, maxlevel, strategy);
#endif
            v += 0.5 *  omega * (f - Av); // the factor 0.5 is needed in case the computed value of normA or normAinv isn't accurate enough
        }

        cout << "CDD2:: v.size() = " << v.size() << endl << endl;

        //      v.COARSE(std::min((1-theta)*epsilon_k,1.0e-6), u_epsilon);
        v.COARSE((1-theta)*epsilon_k, u_epsilon);
        //      v.COARSE(1.0e-6, u_epsilon);

        if (loops%50 == 0) // log convergence data every 50 outer iterations
        {
            logger.pauseClock();
            WaveletTL::APPLY(P, u_epsilon, 1e-6, Av, maxlevel, strategy);
            double residual_norm = l2_norm(F - Av);
            logger.logConvergenceData(u_epsilon.size(), residual_norm);
            cout << "#############################################" << endl;
            cout << "Number of degrees of freedom = " << u_epsilon.size() << endl;
            cout << "current residual error ||f-Av||=" << residual_norm << endl;
            cout << "#############################################" << endl << endl;
            logger.continueClock();
        }
    }

}


#endif // IMPLEMENTATION_CDD2_H
