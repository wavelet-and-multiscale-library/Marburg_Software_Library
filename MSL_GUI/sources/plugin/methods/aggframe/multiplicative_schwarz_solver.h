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


#ifndef METHODS_AGGFRAME_MULTIPLICATIVE_SCHWARZ_SOLVER_H
#define METHODS_AGGFRAME_MULTIPLICATIVE_SCHWARZ_SOLVER_H


#include "implementation_multiplicative_schwarz.h"
#include "misc/cached_frame_problem_helper.h"
#include "methods/method_base.h"

namespace methods
{
namespace aggframe
{

class MultiplicativeSchwarzSolver : public MethodBase<MultiplicativeSchwarzSolver, 5, 2>
{
public:
    static QString name()
    {
        return QStringLiteral("Multiplicative Schwarz Solver");
    }


    static void addParametersTo(ParameterWidget& widget)
    {
        widget.addDoubleParameter("mu", 1.0, 0.1, 999.9, QChar(0x03BC));
        widget.addDoubleParameter("rho", 0.5, 0.1, 0.9, QChar(0x03C1), 0.1);
        widget.addStringParameter("Version", {"Full", "Sparse"});
    }


    template <class PROBLEM>
    static void solve(CachedFrameProblemHelper<PROBLEM>& cproblem, const GuiInputData& input,
                      Array1D< InfiniteVector<double, typename PROBLEM::Index> >& coeffs,
                      MathTL::AbstractConvergenceLogger& logger)
    {
        bool sparseVersion = false;
        if (stringParameter("Version") == "Sparse")
            sparseVersion = true;

        multiplicative_schwarz_GuiSolve(cproblem.cachedProblemLocal(), input.epsilon, coeffs, logger, input.jmax,
                                        input.overlap,
                                        sparseVersion,
                                        doubleParameter("mu"),
                                        doubleParameter("rho"));
    }
};

}
}


#endif // METHODS_AGGFRAME_MULTIPLICATIVE_SCHWARZ_SOLVER_H
