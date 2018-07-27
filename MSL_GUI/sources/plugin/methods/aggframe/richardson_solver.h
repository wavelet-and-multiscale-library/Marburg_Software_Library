/*  -*- c++ -*-

   +-----------------------------------------------------------------------+
   | MSL GUI - A Graphical User Interface for the Marburg Software Library |
   |                                                                       |
   | Copyright (C) 2018 Henning Zickermann                                 |
   | Contact: <zickermann@mathematik.uni-marburg.de>                       |
   +-----------------------------------------------------------------------+

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


#ifndef METHODS_AGGFRAME_RICHARDSON_SOLVER_H
#define METHODS_AGGFRAME_RICHARDSON_SOLVER_H


#include "implementation_richardson.h"
#include "misc/cached_frame_problem_helper.h"
#include "methods/method_base.h"

namespace methods
{
namespace aggframe
{

class RichardsonSolver : public MethodBase<RichardsonSolver, 5, 2>
{
public:
    static QString name()
    {
        return QStringLiteral("Richardson Solver");
    }


    static void addParametersTo(ParameterWidget& widget)
    {

    }


    template <class PROBLEM>
    static void solve(CachedFrameProblemHelper<PROBLEM>& cproblem, const GuiInputData& input,
                      Array1D< InfiniteVector<double, typename PROBLEM::Index> >& coeffs,
                      MathTL::AbstractConvergenceLogger& logger)
    {
        const double nu = cproblem.cachedProblem().norm_Ainv() * cproblem.cachedProblem().F_norm();
        richardson_GuiSolve(cproblem.cachedProblem(), nu, input.epsilon, coeffs, logger, input.jmax, WaveletTL::CDD1);
    }
};

}
}


#endif // METHODS_AGGFRAME_RICHARDSON_SOLVER_H
