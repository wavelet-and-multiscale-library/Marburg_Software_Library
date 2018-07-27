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


#ifndef IMPLEMENTATION_RICHARDSON_H
#define IMPLEMENTATION_RICHARDSON_H


#include <cmath>
#include <set>
#include <WaveletTL/adaptive/apply.h>

#include "MathTL/utils/convergence_logger.h"

using std::set;
using MathTL::InfiniteVector;


template <class PROBLEM>
void richardson_GuiSolve(const PROBLEM& P, const double nu, const double epsilon,
                         Array1D<InfiniteVector<double, typename PROBLEM::Index> >& approximations,
                         MathTL::AbstractConvergenceLogger& logger,
                         const int jmax, WaveletTL::CompressionStrategy strategy)
{
//    const double nu = P.norm_Ainv()*P.F_norm();
    typedef typename PROBLEM::Index Index;

    cout << "Rich_SOLVE: nu=" << nu << endl;

    // compute optimal relaxation parameter omega
    //const double omega = 2.0 / (P.norm_A() + 1.0/P.norm_Ainv());

    // ####### 1D #######
    //    const double omega = 0.4;
    const double omega = 0.2;

    // ####### 2D #######
    // bad one
    //    const double omega = 0.05;
    // good one
    //const double omega = 0.25;

    //const double omega = 0.3;


    //const double omega = 2.0/2.47-0.5;
    cout << "Rich_SOLVE: omega=" << omega << endl;

    // compute spectral norm rho
    //const double cond_A = P.norm_A() * P.norm_Ainv();
    //const double rho = (cond_A - 1.0) / (cond_A + 1.0);
    const double rho = 0.8;
    cout << "Rich_SOLVE: rho=" << rho << endl;

    // desired error reduction factor theta < 1/3
    //const double theta = 2.0/7.0;
    const double theta = 0.333;
    cout << "Rich_SOLVE: theta=" << theta << endl;

    // compute minimal K such that 3*rho^K < theta
    const int K = (int) ceil(log(theta/3.0) / log(rho));
    cout << "Rich_SOLVE: K=" << K << endl << endl;

    unsigned int loops = 0;

    double epsilon_k = nu;
    InfiniteVector<double,Index> f, v, Av, F, u_epsilon;

    P.RHS(1e-6, F);

    logger.startClock();

    try
    {
        while (epsilon_k > epsilon)
        {
            logger.checkAbortConditions();

            ++loops;
            epsilon_k *= 3*pow(rho, K) / theta;
            cout << "Rich_SOLVE: loop: " << loops << endl;
            cout << "Rich_SOLVE: u.size() = " << u_epsilon.size() << endl;
            cout << "Rich_SOLVE: epsilon_k=" << epsilon_k << endl;
            double eta = theta * epsilon_k / (6*omega*K);
            cout << "Rich_SOLVE: eta=" << eta << endl << endl;
            P.RHS(eta, f);

            v = u_epsilon;

            for (int j = 1; j <= 1/*K*/; j++)
            {
                WaveletTL::APPLY(P, v, eta, Av, jmax, strategy);
//                WaveletTL::APPLY_COARSE(P, v, eta, Av, 1.0e-6, jmax, strategy);

                v += omega * (f - Av);
            }
            v.COARSE((1-theta)*epsilon_k, u_epsilon);

            if (loops%50 == 0) // log convergence data every 50 outer iterations
            {
                logger.pauseClock();
                WaveletTL::APPLY(P, u_epsilon, 1e-6, Av, jmax, strategy);
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
    catch(...)
    {
        // collect final approximation and its local parts
        approximations[P.basis().n_p()] = u_epsilon;

        for (int i = 0; i < P.basis().n_p(); i++) {
            approximations[i].clear();
            for (typename InfiniteVector<double, Index>::const_iterator it = u_epsilon.begin(), itend = u_epsilon.end();
                 it != itend; ++it)
                if (it.index().p() == i)
                    approximations[i].set_coefficient(it.index(),*it);
        }

        throw;
    }

    // collect final approximation and its local parts
    approximations[P.basis().n_p()] = u_epsilon;

    for (int i = 0; i < P.basis().n_p(); i++) {
        approximations[i].clear();
        for (typename InfiniteVector<double, Index>::const_iterator it = u_epsilon.begin(), itend = u_epsilon.end();
             it != itend; ++it)
            if (it.index().p() == i)
                approximations[i].set_coefficient(it.index(),*it);
    }

}


#endif // IMPLEMENTATION_RICHARDSON_H
