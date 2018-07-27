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


#ifndef METHODS_QFRAME_STEEPEST_DESCENT_KS_QUARKLET_H
#define METHODS_QFRAME_STEEPEST_DESCENT_KS_QUARKLET_H


#include "implementation_steepest_descent_ks_quarklet.h"

#include "methods/method_base.h"

namespace methods
{
namespace qframe
{

class SteepestDescent_KS_Quarklet : public MethodBase<SteepestDescent_KS_Quarklet, 4, 2>
{
public:
    static QString name()
    {
        return QStringLiteral("Steepest Descent KS Quarklet Solver");
    }


    static void addParametersTo(ParameterWidget& widget)
    {
        widget.addDoubleParameter("beta", 0.85, 0.01, 0.99, QChar(0x03B2), 0.01);
        widget.addDoubleParameter("a", 2.0, 1.0, 99.0);
        widget.addDoubleParameter("b", 2.0, 1.0, 99.0);
    }


    template <class CPROBLEM>
    static void solve(CPROBLEM& cproblem, const GuiInputData& input,
                      MathTL::InfiniteVector<double, typename CPROBLEM::Index>& coeffs,
                      MathTL::AbstractConvergenceLogger& logger)
    {
        MathTL::InfiniteVector<double, int> coeffs_int;

        try
        {
            WaveletTL::CompressionStrategy strategy;
            if (CPROBLEM::space_dimension == 1)
                strategy = WaveletTL::DKR;
            else
                strategy = WaveletTL::tensor_simple;

            steepest_descent_ks_QUARKLET_GuiSolve(cproblem, input.epsilon, coeffs_int, logger, input.jmax,
                                   strategy, input.pmax,
                                   doubleParameter("beta"),
                                   doubleParameter("a"),
                                   doubleParameter("b"));
        }
        catch(...)
        {
            for (typename MathTL::InfiniteVector<double,int>::const_iterator it(coeffs_int.begin()),
                   itend(coeffs_int.end()); it != itend; ++it)
            {
                coeffs.set_coefficient(*(cproblem.frame().get_quarklet(it.index())), *it);
            }

            throw;
        }

        for (typename MathTL::InfiniteVector<double,int>::const_iterator it(coeffs_int.begin()),
               itend(coeffs_int.end()); it != itend; ++it)
        {
            coeffs.set_coefficient(*(cproblem.frame().get_quarklet(it.index())), *it);
        }
    }
};

}
}


#endif // METHODS_QFRAME_STEEPEST_DESCENT_KS_QUARKLET_H
